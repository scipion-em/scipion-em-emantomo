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
from enum import Enum
from os.path import join, basename
from emantomo import Plugin
from emantomo.constants import PARTICLES_LST_FILE, SETS_DIR
from emantomo.objects import EmanParticle
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS
from pyworkflow.protocol import PointerParam, StringParam, FloatParam, LEVEL_ADVANCED, IntParam, GT, LE
from pyworkflow.utils import Message
from pyworkflow.utils.path import makePath, createLink
from tomo.objects import AverageSubTomogram


class OutputsInitModelNew(Enum):
    average = AverageSubTomogram


class EmanProtTomoInitialModelNew(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_sgd_new.py* EMAN2 program.
    It generates an initial model from subtomograms using stochastic gradient descent.
    """

    _label = 'New initial model'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._possibleOutputs = OutputsInitModelNew

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_SUBTOMOS, PointerParam,
                      pointerClass='EmanSetOfParticles',
                      label='Particles',
                      important=True,
                      help='Select the set of subtomograms to build an initial model')
        form.addParam('refVol', PointerParam,
                      pointerClass='Volume',
                      allowsNull=True,
                      label="Reference volume (opt.)",
                      help='Specify a 3D volume')
        form.addParam('shrink', IntParam,
                      default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Shrink factor',
                      help='This option can be used to shrink the input particles by an integer amount '
                           'prior to reconstruction, making them smaller. Default = 1 means no shrinking')

        form.addSection(label='Optimization')
        form.addParam('nIters', IntParam,
                      default=100,
                      label='No. iterations')
        form.addParam('symmetry', StringParam,
                      default='c1',
                      label='Symmetry',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('targetRes', FloatParam,
                      default=50,
                      label='Target resolution')
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
                      default=0.1,
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
        self._insertFunctionStep(self.createEmanPrjStep)
        if self.refVol.get():
            self._insertFunctionStep(self.convertRefVolStep)
        self._insertFunctionStep(self.buildEmanSetsStep)
        self._insertFunctionStep(self.createInitialModelStep)

    # --------------------------- STEPS functions -----------------------------
    def createEmanPrjStep(self):
        inSubtomos = getattr(self, IN_SUBTOMOS).get()
        # Create project dir structure
        self.createInitEmanPrjDirs()
        infoDir = self.getInfoDir()
        tsDir = self.getTsDir()
        tomoDir = self.getTomogramsDir()
        stack2dDir = self.getStack2dDir()
        stack3dDir = self.getStack3dDir()
        makePath(self.getSetsDir(), stack2dDir, stack3dDir)
        # Get the unique values of the files to be linked
        dataDict = inSubtomos.getUniqueValues([EmanParticle.INFO_JSON,
                                               EmanParticle.TS_HDF,
                                               EmanParticle.TOMO_HDF,
                                               EmanParticle.STACK_2D_HDF,
                                               EmanParticle.STACK_3D_HDF])
        # Link the files
        for infoJson, tsFile, tomoFile, stack2d, stack3d in zip(dataDict[EmanParticle.INFO_JSON],
                                                                dataDict[EmanParticle.TS_HDF],
                                                                dataDict[EmanParticle.TOMO_HDF],
                                                                dataDict[EmanParticle.STACK_2D_HDF],
                                                                dataDict[EmanParticle.STACK_3D_HDF]):
            createLink(infoJson, join(infoDir, basename(infoJson)))
            createLink(tsFile, join(tsDir, basename(tsFile)))
            createLink(tomoFile, join(tomoDir, basename(tomoFile)))
            createLink(stack2d, join(stack2dDir, basename(stack2d)))
            createLink(stack3d, join(stack3dDir, basename(stack3d)))

    def convertRefVolStep(self):
        inRef = self.refVol.get()
        self.convertOrLink(inRef.getFileName(), 'reference', self._getExtraPath(), inRef.getSamplingRate())

    def buildEmanSetsStep(self):
        program = Plugin.getProgram("e2spt_buildsets.py")
        args = '--allparticles '
        self.runJob(program, args, cwd=self._getExtraPath())

    def createInitialModelStep(self):
        program = Plugin.getProgram("e2spt_sgd_new.py")
        self.runJob(program, self._genIniModelArgs(), cwd=self._getExtraPath())

    # --------------------------- UTILS functions -----------------------------
    def _genIniModelArgs(self):
        args = ' %s ' % join(SETS_DIR, PARTICLES_LST_FILE)
        if self.refVol.get():
            args += '--ref %s ' % self.getRefVolName()
        args += '--shrink %i ' % self.shrink.get()
        args += '--niter %i ' % self.nIters.get()
        args += '--sym %s ' % self.symmetry.get()
        args += '--res  %.2f ' % self.targetRes.get()
        args += '--batch %i ' % self.batchSize.get()
        args += '--keep %.2f ' % self.keptParticles.get()
        args += '--learnrate %.2f ' % self.learningRate.get()
        args += '--parallel thread:%i ' % self.numberOfThreads.get()
        args += '--verbose 9 '
        return args

    def getRefVolName(self):
        return self._getExtraPath('reference.hdf')

    # -------------------------- INFO functions -------------------------------
