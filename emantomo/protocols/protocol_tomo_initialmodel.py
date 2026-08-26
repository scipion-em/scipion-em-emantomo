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
    Generates an initial 3D subtomogram model from a collection of aligned or non-aligned subtomograms using stochastic gradient descent strategies implemented in EMAN2. The protocol is intended for cryo-electron tomography workflows where an initial reference structure is required before high-resolution refinement. It produces a biologically meaningful starting map that can later be refined through subtomogram averaging and iterative alignment procedures.

    AI Generated:

    Initial Model Generation (EmanProtTomoInitialModel) - User Manual
        Overview

        The Initial Model Generation protocol creates a low- to medium-resolution
        starting structure from a set of subtomograms. In subtomogram averaging
        workflows, obtaining a reliable initial model is one of the most critical
        steps because the quality of the starting reference strongly influences the
        success of downstream refinement and classification procedures.

        The protocol is designed for situations in which the user has extracted
        subtomograms representing repeated copies of the same biological object
        within tomograms. Typical applications include macromolecular complexes,
        membrane-associated assemblies, viral components, ribosomes, or other
        cellular structures studied directly in situ. The generated model serves
        as an unbiased or weakly biased approximation of the underlying structure
        and provides the foundation for subsequent refinement stages.

        Inputs and Biological Context

        The protocol requires a set of subtomograms as the primary input. These
        particles should ideally represent the same molecular species or conformational
        state. Excessive heterogeneity within the dataset can reduce the quality of
        the generated model and may lead to blurred or biologically ambiguous
        structures.

        An optional reference volume can also be provided. This reference is useful
        when prior structural knowledge exists or when the user wants to guide the
        reconstruction toward a known structural state. In exploratory analyses,
        however, users may prefer to avoid a strong reference in order to minimize
        model bias and allow the reconstruction to emerge more directly from the data.

        Since subtomogram averaging is highly sensitive to voxel size consistency
        and box dimensions, all inputs should share compatible sampling rates and
        dimensions. Mismatched geometry between particles, masks, and references
        may compromise reconstruction quality and biological interpretability.

        Symmetry Considerations

        Symmetry can play an important role in improving the quality of the initial
        model. Biological assemblies with known rotational or point-group symmetry,
        such as viral capsids or oligomeric protein complexes, often benefit from
        symmetry application because averaging equivalent views enhances signal and
        suppresses noise.

        Nevertheless, symmetry should only be applied when biologically justified.
        Incorrect symmetry assumptions may artificially distort the reconstruction
        and hide meaningful asymmetric features. For uncertain systems, beginning
        with no symmetry is often the safest strategy.

        Masking and Structural Focus

        The protocol optionally allows the use of a three-dimensional mask. From a
        biological perspective, masking is highly valuable because it restricts the
        reconstruction focus to the region of interest while suppressing unrelated
        density and surrounding noise.

        This becomes especially important for membrane proteins, flexible complexes,
        or crowded cellular environments where neighboring densities may interfere
        with convergence. A carefully designed mask centered on the stable structural
        core generally improves robustness and produces more interpretable initial
        maps.

        Excessively tight masks, however, may remove biologically relevant regions,
        whereas overly broad masks may fail to suppress background noise effectively.
        Selecting an appropriate mask therefore requires balancing structural focus
        with preservation of meaningful density.

        Stochastic Gradient Descent Reconstruction

        The reconstruction strategy is based on stochastic gradient descent methods,
        which iteratively improve the model using subsets of particles across multiple
        optimization cycles. This approach is computationally efficient and generally
        well suited for large cryo-electron tomography datasets.

        The batch size determines how many subtomograms contribute simultaneously
        during each optimization update. Larger batch sizes may improve stability
        and convergence behavior but also require greater computational resources.
        Smaller batches introduce more stochastic variability, which can sometimes
        help avoid convergence toward suboptimal solutions.

        The learning rate controls how aggressively the reconstruction evolves during
        optimization. Conservative values typically produce more stable convergence,
        whereas excessively large values may destabilize the reconstruction process.
        In most biological workflows, moderate default settings provide a good balance
        between stability and reconstruction speed.

        Iterative Refinement and Filtering

        The number of iterations and batches determines the total extent of optimization.
        Increasing these parameters may improve structural consistency and map quality,
        particularly for noisy datasets, but also increases runtime and computational
        demand.

        Frequency filtering can be applied between iterations to stabilize the early
        stages of reconstruction. This is particularly useful when dealing with highly
        noisy subtomograms or low particle counts. By emphasizing low-resolution
        structural information during early optimization, the protocol promotes the
        emergence of globally consistent shapes before attempting to recover finer
        details.

        Fourier Space Optimization

        The protocol optionally performs optimization in Fourier space. For many
        cryo-electron tomography datasets, Fourier-based optimization improves
        computational efficiency and may enhance convergence stability during early
        model generation stages.

        From a biological perspective, the main objective remains the recovery of
        reliable large-scale structural organization rather than immediate recovery
        of fine high-resolution details. Initial models are expected to represent
        the general architecture of the complex rather than final atomic-quality maps.

        Particle Shrinking and Computational Efficiency

        The protocol allows downsampling of subtomograms before reconstruction. This
        option is often beneficial for very large box sizes because initial model
        generation primarily depends on recovering broad structural features rather
        than fine detail.

        Shrinking substantially reduces computational cost and may improve convergence
        speed. In practical biological workflows, users frequently begin with reduced
        particle sizes during initial model generation and later return to the full
        resolution data during refinement.

        Outputs and Interpretation

        The main output is an averaged subtomogram volume representing the reconstructed
        initial model. This volume can subsequently be used as the reference for
        subtomogram refinement, alignment optimization, or classification workflows.

        Biologically, the resulting structure should be interpreted as an approximate
        consensus representation of the input particles. Flexible regions may appear
        blurred, and heterogeneous populations may produce partially averaged features.
        Visual inspection and iterative refinement are therefore essential for assessing
        model quality and biological plausibility.

        Practical Recommendations

        In most practical tomography workflows, it is advisable to begin with a clean,
        homogeneous subset of particles whenever possible. Strong heterogeneity often
        complicates convergence and reduces interpretability of the resulting model.

        When prior structural information is limited, users commonly start without a
        reference and avoid imposing symmetry. Once a stable initial model is obtained,
        more advanced refinement strategies can progressively improve resolution and
        structural accuracy.

        For noisy datasets, combining moderate filtering, appropriate masking, and
        particle shrinking often yields the most stable reconstructions. Users should
        visually inspect intermediate and final maps to confirm that the emerging
        structure remains biologically meaningful.

        Final Perspective

        Initial model generation is one of the foundational steps in subtomogram
        averaging because it establishes the structural framework used throughout
        later refinement procedures. Reliable initial models improve alignment
        stability, reduce the risk of reconstruction bias, and increase the likelihood
        of recovering biologically accurate structures from noisy tomographic data.

        Successful use of this protocol depends on careful dataset preparation,
        thoughtful symmetry selection, biologically appropriate masking, and realistic
        expectations regarding the low-resolution nature of early-stage reconstructions.
    """

    _label = 'Initial model'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.volumeFileMrc = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('particles', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      strict=True,
                      label="Particles",
                      important=True,
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

    def _validate(self):
        errorMsg = []
        refVol = self.reference.get()
        # Check the dimensions
        inParticles = self.particles.get()
        dimsDict = inParticles._firstDim
        inParticlesDims = (dimsDict[0], dimsDict[1], dimsDict[2])
        # Eman fails if the dimensions are odd numbers
        for iDim in inParticlesDims:
            if iDim % 2 != 0:
                errorMsg.append('The particles dimensions must be an even number')
                break
        if refVol:
            x, y, z = refVol.getDimensions()
            refVolDims = (x, y, z)
            if refVolDims != inParticlesDims:
                errorMsg.append(f'The dimensions of the reference volume {refVolDims} px and the particles '
                                f'{inParticlesDims} px must be the same')
            # Check the sampling rate
            tol = 2e-03
            inParticlesSRate = inParticles.getSamplingRate()
            refVolSRate = refVol.getSamplingRate()
            if abs(inParticlesSRate - refVolSRate) >= tol:
                errorMsg.append(
                    f'The sampling rate of the input particles [{inParticlesSRate} Å/pix] and the reference volume '
                    f'[{refVolSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')
            # Check the mask
            mask = self.mask.get()
            if mask:
                x, y, z = mask.getDimensions()
                maskDims = (x, y, z)
                maskSRate = mask.getSamplingRate()
                if inParticlesDims != maskDims:
                    errorMsg.append(f'The dimensions of the introduced mask {maskDims} px and the particles '
                                    f'{inParticlesDims} px must be the same')
                if abs(inParticlesSRate - maskSRate) >= tol:
                    errorMsg.append(
                        f'The sampling rate of the input particles [{inParticlesSRate} Å/pix] and the mask '
                        f'[{maskSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')

        return errorMsg
