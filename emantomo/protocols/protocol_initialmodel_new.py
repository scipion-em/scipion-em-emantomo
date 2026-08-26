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
from emantomo.protocols.protocol_base import IN_SUBTOMOS, REF_VOL
from emantomo.protocols.protocol_refine_new_base import EmanProtRefineNewBase
from pwem.convert.headers import fixVolume
from pwem.emlib.image import ImageHandler
from pyworkflow.protocol import StringParam, FloatParam, LEVEL_ADVANCED, IntParam, GT, LE
from pyworkflow.utils import Message
from tomo.objects import AverageSubTomogram, SetOfAverageSubTomograms


class OutputsInitModelNew(Enum):
    averages = SetOfAverageSubTomograms
    average = AverageSubTomogram


class EmanProtTomoInitialModelNew(EmanProtRefineNewBase):
    """
    Generates one or more initial 3D subtomogram models using stochastic
    gradient descent optimization within the EMAN2 subtomogram refinement
    framework. The protocol is intended for early stages of cryo-electron
    tomography workflows where no reliable starting structure is available
    or where multiple candidate structural states must be explored before
    refinement and classification.

    AI Generated:

    Initial Model Generation (EmanProtTomoInitialModelNew) — User Manual
        Overview

        The Initial Model Generation protocol reconstructs low-resolution
        three-dimensional reference maps directly from aligned subtomograms.
        Its primary goal is to provide biologically meaningful starting
        structures that can later be refined, classified, or interpreted in
        downstream subtomogram averaging workflows. In cryo-electron tomography,
        generating a stable initial model is often one of the most critical
        steps because the quality of the starting reference strongly influences
        the success of subsequent refinement procedures.

        From a biological perspective, this protocol is especially useful when
        studying complexes with unknown conformations, poorly characterized
        assemblies, or heterogeneous populations for which no high-quality
        external reference exists. It can also be used to generate several
        independent candidate structures in parallel when structural variability
        is expected.

        Inputs and Biological Context

        The protocol requires a set of subtomograms that have already undergone
        preprocessing and approximate alignment. Ideally, the particles should
        correspond to the same biological assembly and should share similar box
        dimensions and voxel sizes. Excessive heterogeneity, severe misalignment,
        or strong contamination may prevent the optimization from converging to
        biologically interpretable structures.

        Optionally, a reference volume may be provided to guide the optimization.
        This can be useful when a low-resolution map, homologous structure, or
        previously reconstructed average already exists. However, biological care
        is required because an overly strong or biased starting reference may
        influence the reconstruction toward a preferred structural solution.

        Stochastic Gradient Descent Optimization

        The protocol uses stochastic gradient descent as the central optimization
        strategy for reconstructing the initial model. Rather than processing the
        entire dataset simultaneously, the method iteratively updates the model
        using subsets of particles. This approach is computationally efficient
        and often improves robustness when handling large subtomogram datasets.

        For biological users, this iterative optimization process allows the
        reconstruction to progressively capture common structural features while
        suppressing random noise. The resulting maps are typically low to medium
        resolution representations that are suitable as starting points for later
        high-resolution refinement.

        Number of Classes and Structural Heterogeneity

        The protocol can generate either a single consensus structure or multiple
        independent classes. Producing several classes is particularly valuable
        when the sample contains conformational variability, compositional
        heterogeneity, or multiple assembly states.

        In practical biological applications, requesting multiple classes may
        reveal distinct structural populations that would otherwise be averaged
        together. However, increasing the number of classes also divides the
        available particle information among several reconstructions, which may
        reduce stability when the dataset is small.

        Symmetry Considerations

        Users may define a symmetry group to constrain the reconstruction. This
        can substantially improve convergence and signal quality when the target
        assembly is known to possess rotational or point-group symmetry.

        From a biological standpoint, symmetry should only be imposed when it is
        strongly supported by prior structural evidence. Applying incorrect
        symmetry may artificially distort flexible regions, suppress asymmetric
        features, or generate misleading structural interpretations.

        Resolution Targets and Model Interpretation

        The target resolution parameter controls the expected level of structural
        detail during optimization. In most initial model generation workflows,
        relatively modest resolutions are sufficient because the objective is to
        recover the global architecture rather than fine molecular detail.

        Biological users should interpret the resulting structures primarily as
        low-resolution organizational maps. Flexible regions, small domains, or
        weakly occupied features may not be accurately represented at this stage.
        The generated models are intended to guide refinement rather than serve
        as final structural interpretations.

        Batch Size and Particle Retention

        The optimization process operates on subsets of particles grouped into
        batches. Larger batches generally improve statistical stability and can
        accelerate convergence when sufficient computational resources are
        available. Smaller batches may introduce additional stochastic behavior
        but can sometimes help avoid poor local optima.

        The protocol also allows retaining only a fraction of particles during
        optimization. Biologically, this may help suppress outliers, damaged
        particles, or rare contaminants that could otherwise bias the initial
        reconstruction. However, overly aggressive particle rejection may discard
        meaningful structural variability.

        Learning Rate and Optimization Stability

        The learning rate controls how strongly each optimization step modifies
        the evolving reconstruction. This parameter strongly influences the
        balance between convergence speed and stability.

        Higher learning rates may accelerate reconstruction but can lead to
        unstable behavior or structural artifacts if updates become too large.
        Lower learning rates generally provide more conservative optimization and
        smoother convergence but may require more iterations. In challenging
        biological datasets, moderate learning rates are often the safest choice.

        Binning and Computational Efficiency

        The protocol optionally supports shrinking or binning the input particles
        prior to reconstruction. Binning reduces the effective box size and
        lowers computational demands, which is especially useful during early
        exploratory stages or when processing very large subtomograms.

        From a biological perspective, binning sacrifices high-frequency detail
        but usually preserves the overall architecture needed for reliable
        initial model estimation. Many workflows begin with strongly binned data
        and later refine the structures using progressively higher resolutions.

        Outputs and Their Interpretation

        Depending on the selected configuration, the protocol produces either a
        single average subtomogram or a set of independently reconstructed
        averages. These outputs represent candidate structural models suitable
        for downstream alignment, refinement, classification, or visualization.

        When multiple classes are generated, each class should be interpreted as
        a potential structural state rather than definitive evidence of
        biological heterogeneity. Visual inspection, refinement consistency, and
        independent biological validation remain essential.

        Practical Recommendations

        In routine cryo-electron tomography workflows, it is generally advisable
        to begin with conservative parameters and a modest target resolution.
        Using a single class is often appropriate for homogeneous datasets,
        whereas heterogeneous samples may benefit from generating several
        candidate structures in parallel.

        Applying known symmetry can substantially improve reconstruction quality,
        but only when biologically justified. If convergence appears unstable,
        increasing the number of iterations, adjusting the learning rate, or
        reducing structural heterogeneity through improved particle selection may
        improve the results.

        Initial models should always be inspected visually before proceeding to
        high-resolution refinement. Features inconsistent with known biology,
        severe fragmentation, or excessive symmetry artifacts may indicate
        problems in preprocessing, alignment, or optimization settings.

        Final Perspective

        Initial model generation is a foundational step in subtomogram averaging
        because it establishes the structural framework for all later refinement
        and interpretation. Careful selection of optimization parameters,
        thoughtful handling of heterogeneity, and biologically realistic symmetry
        assumptions are essential for obtaining reliable starting structures that
        can support meaningful downstream structural analysis.
    """
    _label = 'Initial model pppt'
    _possibleOutputs = OutputsInitModelNew

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        self._addCommonInputParams(form)
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
                     f"--parallel thread:{self.binThreads.get()}",
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
    def _validate(self):
        errorMsg = []
        refVol = self.getAttrib(REF_VOL)
        if refVol:
            # Check the dimensions
            ih = ImageHandler()
            x, y, z, _ = ih.getDimensions(refVol.getFileName())
            refVolDims = (x, y, z)
            inParticles = self.getAttrib(IN_SUBTOMOS)
            inParticlesDims = inParticles.getBoxSize()
            if refVolDims != inParticlesDims:
                errorMsg.append(f'The dimensions of the reference volume {refVolDims} px and the particles '
                                f'{inParticlesDims} px must be the same')
            # Check the sampling rate
            tol = 1e-03
            inParticlesSRate = inParticles.getSamplingRate()
            refVolSRate = refVol.getSamplingRate()
            if abs(inParticlesSRate - refVolSRate) >= tol:
                errorMsg.append(
                    f'The sampling rate of the input particles [{inParticlesSRate} Å/pix] and the reference volume '
                    f'[{refVolSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')
        return errorMsg


