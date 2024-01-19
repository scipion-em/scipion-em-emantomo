# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import shutil
from enum import Enum
from os.path import exists

from emantomo import Plugin
from emantomo.constants import INIT_MODEL_DIR, INIT_MODEL_NAME, INIT_MODEL_MRC, \
    SYMMETRY_HELP_MSG, REFERENCE_NAME
from emantomo.convert import convertBetweenHdfAndMrc
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from pwem.convert.headers import fixVolume
from pyworkflow.protocol import PointerParam, StringParam, FloatParam, LEVEL_ADVANCED, IntParam, GT, LE
from pyworkflow.utils import Message
from tomo.objects import AverageSubTomogram, SetOfAverageSubTomograms


class OutputsInitModelNew(Enum):
    averages = SetOfAverageSubTomograms
    average = AverageSubTomogram


class EmanProtTomoInitialModelNew(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_sgd_new.py* EMAN2 program.
    It generates an initial model from subtomograms using stochastic gradient descent.
    """

    _label = 'Initial model pppt'
    _possibleOutputs = OutputsInitModelNew

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_SUBTOMOS, PointerParam,
                      pointerClass='EmanSetOfParticles',
                      label='Particles',
                      important=True,
                      help='Select the set of subtomograms to build an initial model')
        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      allowsNull=True,
                      label="Reference volume (opt.)")
        form.addParam('shrink', IntParam,
                      default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Binning factor',
                      help='This option can be used to shrink the input particles by an integer amount '
                           'prior to reconstruction, making them smaller. Default = 1 means no shrinking')

        form.addSection(label='Optimization')
        form.addParam('nIters', IntParam,
                      default=100,
                      label='No. iterations')
        form.addParam('nClasses', IntParam,
                      default=1,
                      label='No. classes')
        form.addParam('symmetry', StringParam,
                      default='c1',
                      label='Symmetry',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('targetRes', FloatParam,
                      default=50,
                      label='Target resolution (Å)')
        form.addParam('batchSize', IntParam,
                      default=12,
                      label='Batch Size',
                      help='SGD batch size. Increasing batch size will use more cores (if you have more than 12), and '
                           'may cause it to converge to the correct answer in fewer iterations, but each iteration '
                           'will not become faster.')
        form.addParam('keptParticles', FloatParam,
                      default=1,
                      label='Fraction of particles to keep',
                      expertLevel=LEVEL_ADVANCED,
                      validators=[GT(0), LE(1)],
                      help='It will actually align more particles and use the number of particles specified by the '
                           'batch size parameter. Default = 1 means that all the particles are kept.')
        form.addParam('learningRate', FloatParam,
                      default=0.2,
                      expertLevel=LEVEL_ADVANCED,
                      label='Learning Rate',
                      help="In the context of stochastic gradient descent (SGD), the learning rate is a hyperparameter "
                           "that determines the step size at each iteration when updating the model's parameters.\n\n"
                           "In other words, the learning rate controls how much the parameters are adjusted in the "
                           "direction of the gradient, which is the direction of steepest descent of the loss "
                           "function. A higher learning rate leads to larger updates and faster convergence, "
                           "but it may also cause the algorithm to overshoot the optimal solution and fail to "
                           "converge. On the other hand, a lower learning rate leads to smaller updates and slower "
                           "convergence, but it may also help the algorithm avoid overshooting and find a more precise "
                           "optimum.\n\n"
                           "Choosing an appropriate learning rate is important for achieving good performance in SGD. "
                           "It typically involves a trade-off between convergence speed and accuracy, and may require "
                           "tuning through trial and error or more advanced optimization techniques such as adaptive "
                           "learning rate methods.")
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.createEmanPrjPostExtractionStep)
        self._insertFunctionStep(self.convertRefVolStep)
        self._insertFunctionStep(self.createInitialModelStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inParticles = self.getAttrib(IN_SUBTOMOS)
        self.inSamplingRate = self.inParticles.getSamplingRate()

    def createInitialModelStep(self):
        # In case of continuing from this step, the previous results dir will be removed to avoid EMAN creating one
        # for each execution (one for each continue)
        initModelDir = self.getInitModelDir()
        if exists(initModelDir):
            shutil.rmtree(initModelDir)
        self.buildEmanSets()
        program = Plugin.getProgram("e2spt_sgd_new.py")
        self.runJob(program, self._genIniModelArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        if self.nClasses.get() == 1:
            averageSubTomogram = self._genOutputAvg()
            # Define outputs and relations
            self._defineOutputs(**{self._possibleOutputs.average.name: averageSubTomogram})
            self._defineSourceRelation(getattr(self, IN_SUBTOMOS), averageSubTomogram)
        else:
            avgSet = SetOfAverageSubTomograms.create(self._getPath(), template='setOfAvgSubTomograms%s.sqlite')
            avgSet.copyInfo(self.inParticles)
            for iVol in range(self.nClasses.get()):
                # Generate the current AverageSubtomogram
                averageSubTomogram = self._genOutputAvg(avgIndex=iVol)
                # Append the current average to the resulting set of averages
                avgSet.append(averageSubTomogram)
            # Define outputs and relations
            self._defineOutputs(**{self._possibleOutputs.averages.name: avgSet})
            self._defineSourceRelation(getattr(self, IN_SUBTOMOS), avgSet)

    # --------------------------- UTILS functions -----------------------------
    def _genIniModelArgs(self):
        args = [f" {self._getLstFile()}"]
        if self.getRefVol():
            args.append(f"--ref {REFERENCE_NAME}.hdf")
        args.extend([f"--shrink {self.shrink.get()}",
                     f"--niter {self.nIters.get()}",
                     f"--ncls {self.nClasses.get()}",
                     f"--sym {self.symmetry.get()}",
                     f"--res {self.targetRes.get():.2f}",
                     f"--batch {self.batchSize.get()}",
                     f"--keep {self.keptParticles.get():.2f}",
                     f"--learnrate {self.learningRate.get():.2f}",
                     f"--parallel thread:{self.numberOfThreads.get()}",
                     "--verbose 9 "])
        return ' '.join(args)

    def getInitialModelHdfFile(self, iVol):
        return self._getExtraPath(INIT_MODEL_DIR, INIT_MODEL_NAME % iVol)

    def getInitialModelMrcFile(self, iVol):
        return self._getExtraPath(INIT_MODEL_MRC % iVol)

    def _genOutputAvg(self, avgIndex=0):
        """Gets the initial volume generated by index from the number of initial volumes requested to be generated,
        converts it from HDF to MRC and generates an AverageSubTomogram object with its metadata filled.

        :param avgIndex: Index of the current initial volume, used when iterating if more than one initial volume
        was requested.
        """
        # Convert the output to MRC
        initModelFile = self.getInitialModelHdfFile(avgIndex)
        outFile = self.getInitialModelMrcFile(avgIndex)
        args = '--apix %.3f' % self.inSamplingRate
        convertBetweenHdfAndMrc(self, initModelFile, outFile, args)
        # Generate the corresponding output
        averageSubTomogram = AverageSubTomogram()
        fixVolume(outFile)  # Fix header to make it interpreted as volume instead of a stack by xmipp
        averageSubTomogram.setFileName(outFile)
        averageSubTomogram.setSamplingRate(self.inSamplingRate)
        return averageSubTomogram

    # -------------------------- INFO functions -------------------------------
