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
from emantomo.constants import INIT_MODEL_DIR, INIR_MODEL_NAME_OLD, SYMMETRY_HELP_MSG
from pyworkflow import BETA
from pyworkflow.protocol import params
from pyworkflow.utils.path import makePath, replaceBaseExt
from pwem.protocols import EMProtocol
import emantomo
from emantomo.convert import writeSetOfSubTomograms
from tomo.protocols import ProtTomoBase
from tomo.objects import AverageSubTomogram


class OutputsInitModel(Enum):
    average = AverageSubTomogram


class EmanProtTomoInitialModel(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_sgd.py* EMAN2 program.

    It will take a set of subtomograms (particles) and a subtomogram(reference)
    and build a subtomogram suitable for use as initial models in tomography.
    It also builds a set of subtomograms that contains the original particles
    plus the score, coverage and align matrix per subtomogram .
    """

    _label = 'Initial model'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.volumeFileMrc = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('particles', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Particles", important=True,
                      help='Select the set of subtomograms to build an initial model')

        form.addParam('reference', params.PointerParam,
                      pointerClass='Volume', allowsNull=True,
                      label="Reference volume",
                      help='Specify a 3D volume')

        form.addParam('mask', params.PointerParam,
                      label='Mask',
                      allowsNull=True,
                      pointerClass='VolumeMask',
                      help='Select a 3D Mask to be applied to the initial model')

        form.addSection(label='Optimization')
        form.addParam('symmetry', params.StringParam,
                      default='c1',
                      label='Symmetry',
                      help='Specify the symmetry.\nChoices are: c(n), d(n), '
                           'h(n), tet, oct, icos.\n'
                           'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry\n'
                           'for a detailed description of symmetry in Eman.')
        form.addParam('filterto', params.FloatParam, default=0.02,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Filterto',
                      help='Filter map to frequency after each iteration. Default is 0.02')
        form.addParam('fourier', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Fourier',
                      help='Gradient descent in fourier space')
        form.addParam('batchSize', params.IntParam, default=12,
                      label='Batch Size',
                      help='SGD batch size. Increasing batchsize will use more cores (if you have more than 12), and '
                           'may cause it to converge to the correct answer in fewer iterations, but each iteration '
                           'will not become faster.')
        form.addParam('learningRate', params.FloatParam, default=0.1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Learn Rate',
                      help='Learning Rate. Default is 0.1')
        form.addParam('numberOfIterations', params.IntParam, default=5,
                      label='Number of iterations to perform',
                      help='The total number of refinement iterations to perform.')
        form.addParam('numberOfBatches', params.IntParam, default=10,
                      label='Number of batches',
                      help='Number of batches per iteration')
        form.addParam('shrink', params.IntParam,
                      default=1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Shrink factor',
                      help='Using a box-size >64 is not optimal for making '
                           'initial models. Suggest using this option to '
                           'shrink the input particles by an integer amount '
                           'prior to reconstruction. Default = 1, no shrinking')
        form.addParam('applySim', params.BooleanParam,
                      default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Apply Symmetry',
                      help=SYMMETRY_HELP_MSG)

        # form.addSection(label='Output')
        # form.addParam('returnSubtomos', params.BooleanParam, default=True,
        #               label="Return aligned subtomograms?",
        #               help="Depending on the number of iterations, batch size and number "
        #                    "of batches chosen, it might be possible that some subtomograms "
        #                    "will not contribute to the initial model. By default, algined "
        #                    "subtomograms will be returned, even if the output set might have a "
        #                    "lower number of items.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertImagesStep)
        self._insertFunctionStep(self.createInitialModelStep)
        self._insertFunctionStep(self.converOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    # Get Scipion references to subtomograms and write hdf files for emantomo to process.
    def convertImagesStep(self):
        partSet = self.particles.get()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfSubTomograms(partSet, storePath, alignType=partAlign)

    def createInitialModelStep(self):
        command_params = {
            'symmetry': self.symmetry.get(),
            'filterto': self.filterto.get(),
            'batchSize': self.batchSize.get(),
            'learningRate': self.learningRate.get(),
            'numberOfIterations': self.numberOfIterations.get(),
            'numberOfBatches': self.numberOfBatches.get(),
            'shrink': self.shrink.get(),
            'reference': self.reference.get().getFileName() if self.reference.get() else None,
            'outputPath': self.getOutputPath(),
        }

        args = '%s/*.hdf' % self._getExtraPath("particles")

        if self.mask.get():
            command_params['mask'] = self.mask.get().getFileName()
            args += ' --mask=%(mask)s'

        if command_params['reference']:
            args += ' --reference=%(reference)s'

        args += (' --sym=%(symmetry)s --filterto=%(filterto)f'
                 ' --batchsize=%(batchSize)d --learnrate=%(learningRate)f --niter=%(numberOfIterations)d'
                 ' --nbatch=%(numberOfBatches)d')

        if command_params['shrink'] > 1:
            args += ' --shrink=%(shrink)d'
        if self.fourier.get():
            args += ' --fourier'
        if self.applySim.get():
            args += ' --applysym'

        args += ' --path=%(outputPath)s'

        program = emantomo.Plugin.getProgram("e2spt_sgd.py")
        self.runJob(program, args % command_params)

    def converOutputStep(self):
        # Also fix the sampling rate as it might be set wrong (the value stored in the hdf header may be referred
        # to the original binning, and it will be also in the header of the resulting mrc file
        # 1) Average
        sRate = self.particles.get().getSamplingRate()
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        volumeFileHdf = self._getExtraPath(INIT_MODEL_DIR, INIR_MODEL_NAME_OLD)
        self.volumeFileMrc = self._getExtraPath(replaceBaseExt(volumeFileHdf, "mrc"))
        args = "%s %s --apix=%f" % (volumeFileHdf, self.volumeFileMrc, sRate)
        self.runJob(program, args)

        # # 2) Subtomograms
        # if self.returnSubtomos.get():
        #     particleStacks = glob.glob(self._getExtraPath(PARTICLES_DIR))
        #     for hdfPartStack in particleStacks:
        #         args = ' %s ' % hdfPartStack
        #         args += '%s ' % replaceExt(hdfPartStack, 'mrc')
        #         args += '--apix=%.3f --unstacking' % sRate
        #         self.runJob(program, args)

        # # Remove the hdf files if requested
        # if not self.keepHdfFile.get():
        #     os.remove(hdfFile)

    def createOutputStep(self):
        # outputSubtomos = None
        particles = self.particles.get()
        downSamplingFactor = self.shrink.get()

        # Output 1: Average
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(self.volumeFileMrc)
        averageSubTomogram.setSamplingRate(particles.getSamplingRate() * downSamplingFactor)
        outputDict = {OutputsInitModel.average.name: averageSubTomogram}

        # Output 2: setOfSubTomograms
        # TODO: ¿deberían estas partículas aparecer en el sqlite via setFileName de cada subtomo? Ver la explicación de returnSubtomos
        # if self.returnSubtomos.get():
        #     particleParams = getLastParticlesParams(self.getOutputPath())
        #     outputSubtomos = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
        #     outputSubtomos.setCoordinates3D(particles.getCoordinates3D())
        #     outputSubtomos.copyInfo(particles)
        #     outputSubtomos.setSamplingRate(particles.getSamplingRate() * downSamplingFactor)
        #     updateSetOfSubTomograms(particles, outputSubtomos, particleParams)
        #     outputDict[OutputsInitModel.subtomograms.name] = outputSubtomos

        # Define outputs and relations
        self._defineOutputs(**outputDict)
        self._defineSourceRelation(particles, averageSubTomogram)
        # if self.returnSubtomos.get():
        #     self._defineSourceRelation(particles, outputSubtomos)

    def getOutputPath(self, *args):
        return self._getExtraPath(INIT_MODEL_DIR, *args)

    def _methods(self):
        particles = self.particles.get()
        return [
            "Created an initial model using e2spt_sgd.py (stochastic gradient descent)",
            "A total of %d particles of dimensions %s were used (shrink %d)"
            % (particles.getSize(), particles.getDimensions(), self.shrink.get()),
        ]

    def _summary(self):
        particles = self.particles.get()
        reference = self.reference.get()
        lines = [
            "Particles: %d" % particles.getSize(),
            "Reference file used: %s" % reference.getFileName() if reference else None,
        ]

        return list(filter(bool, lines))
