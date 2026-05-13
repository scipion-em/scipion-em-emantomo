# coding=utf-8
# **************************************************************************
# *
# * Authors:     David Herreros  (dherreros@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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
import enum
import glob
import os
import shutil
from os.path import join, basename, abspath
from pyworkflow.protocol import PointerParam, FloatParam, StringParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.utils import Message, makePath, replaceExt, createLink
from .protocol_base import ProtEmantomoBase
from ..constants import SPT_00_DIR, INPUT_PTCLS_LST, THREED_01, SYMMETRY_HELP_MSG, SUBTOMOGRAMS_DIR
from ..convert import writeSetOfSubTomograms, refinement2Json
import emantomo
from tomo.objects import AverageSubTomogram
from ..objects import EmanParticle, EmanSetOfParticles


class OutputsAverageSubtomos(enum.Enum):
    averageSubTomos = AverageSubTomogram


class EmanProtSubTomoAverage(ProtEmantomoBase):
    """
    Computes an averaged subtomogram volume from a selected set of aligned particles in order to enhance
    common structural features and improve the signal-to-noise ratio in cryo-electron tomography studies.

    AI Generated:

    Subtomogram Averaging (EmanProtSubTomoAverage) — User Manual
        Overview

        The Subtomogram Averaging protocol reconstructs an averaged three-dimensional map from a collection
        of subtomograms that already share a common orientation. Its primary objective is to enhance
        reproducible structural information while suppressing random noise, thereby improving the biological
        interpretability of cryo-electron tomography data.

        In cryo-ET workflows, subtomogram averaging is one of the most important approaches for studying
        macromolecular complexes directly within their native cellular environment. Individual subtomograms
        are typically very noisy because of the limited electron dose and the incomplete angular coverage
        inherent to tomography experiments. By combining many particles representing the same biological
        structure, the protocol produces a higher-quality consensus reconstruction that reveals structural
        details not visible in individual particles.

        Biological Context and Typical Applications

        This protocol is commonly used for the structural analysis of ribosomes, membrane proteins, viral
        assemblies, cytoskeletal complexes, and other macromolecular structures imaged in situ. Averaging
        multiple particles enables researchers to study native molecular organization while preserving the
        biological context of the cellular environment.

        In practical applications, the quality of the final average strongly depends on the consistency of
        the input particles. Subtomograms should ideally represent the same structural state and should be
        approximately aligned before averaging. Mixing distinct conformations, compositional states, or
        incorrectly aligned particles can blur structural features and reduce biological interpretability.

        Inputs and General Workflow

        The protocol requires a set of subtomograms as input. These particles are assumed to already contain
        orientation information defining a common spatial reference frame. The protocol combines the particles
        into a consensus volume and generates corresponding half maps that can later be used for resolution
        estimation and validation.

        The workflow is designed to integrate naturally into iterative subtomogram refinement pipelines.
        Users typically perform particle extraction, alignment, classification, and refinement before
        computing a final average suitable for interpretation or visualization.

        Symmetry Considerations

        The protocol allows the application of symmetry during averaging. Symmetry can significantly improve
        the signal-to-noise ratio because equivalent structural information is combined multiple times within
        the reconstruction process. This is particularly beneficial for highly symmetric assemblies such as
        viral capsids, filaments, or oligomeric protein complexes.

        However, symmetry should only be applied when it is biologically justified. Imposing incorrect
        symmetry may artificially distort the reconstruction or obscure biologically relevant asymmetries.
        For many exploratory studies, beginning with no symmetry is often the safest approach until the
        structural organization is well understood.

        Missing Wedge Compensation

        One of the central challenges in cryo-electron tomography is the missing wedge artifact generated by
        incomplete angular sampling during tilt-series acquisition. This anisotropic loss of information can
        introduce directional distortions and reduce the quality of the reconstruction.

        The protocol includes an option to compensate for missing wedge effects during averaging. This
        correction can substantially improve isotropy and reduce reconstruction artifacts, especially in
        datasets acquired with limited tilt ranges. In many biological applications, enabling missing wedge
        compensation improves the interpretability of membrane-associated complexes and elongated structures.

        Nevertheless, the effectiveness of this correction depends on the quality and orientation diversity
        of the input particles. Poor angular coverage or highly biased particle orientations may still limit
        the final reconstruction quality.

        Post-Processing and Validation

        The protocol can optionally skip post-processing operations such as masking, filtering, and Fourier
        shell correlation analysis. This flexibility is useful in workflows where users prefer to perform
        downstream processing separately or apply custom validation strategies.

        In most biological studies, evaluating the resulting half maps is essential for estimating resolution
        and assessing reconstruction reliability. Half-map comparisons help identify overfitting and provide
        an objective measure of structural reproducibility across independent particle subsets.

        The resulting average should always be interpreted in combination with visual inspection and biological
        expectations. High nominal resolution alone does not guarantee biological correctness, particularly in
        heterogeneous or flexible systems.

        Practical Recommendations

        For routine subtomogram averaging projects, users should carefully curate particle sets before
        averaging. Removing misaligned particles, contaminants, or structurally distinct conformations often
        produces larger improvements than simply increasing the number of particles.

        Applying symmetry should be approached conservatively unless supported by strong biological evidence.
        Likewise, missing wedge correction is generally beneficial but should be evaluated alongside the final
        map quality and anisotropy.

        In challenging datasets, iterative cycles of alignment, classification, and averaging frequently yield
        the best results. Producing intermediate averages during refinement can help monitor convergence and
        identify problematic particle populations.

        Outputs and Biological Interpretation

        The protocol generates an averaged subtomogram map together with independent half maps suitable for
        validation and resolution assessment. The final reconstruction represents the consensus structural
        information shared by the input particles.

        Biologically, the resulting map can reveal molecular organization, domain architecture, interaction
        interfaces, and conformational features that are otherwise inaccessible in individual tomograms.
        However, users should remember that averaging emphasizes common structural characteristics while
        suppressing variability. Flexible regions and heterogeneous states may therefore appear blurred or
        absent in the final reconstruction.

        Final Perspective

        Subtomogram averaging is one of the most powerful techniques in cryo-electron tomography for extracting
        high-resolution structural information directly from complex biological environments. Successful results
        depend not only on computational processing but also on careful biological interpretation, thoughtful
        particle selection, and an understanding of the structural heterogeneity present in the sample.
    """

    _label = 'average subtomograms'
    _possibleOutputs = OutputsAverageSubtomos

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.projectPath = None
        self.volumeFileHdf = None
        self.halfEvenFileHdf = None
        self.halfOddFileHdf = None
        self.hdfSubtomosDir = None

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfSubTomogram', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      strict=True,
                      label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('symmetry', StringParam,
                      label='Symmetry',
                      default='c1',
                      allowsNull=False,
                      help=SYMMETRY_HELP_MSG)
        form.addParam('msWedge', FloatParam,
                      default=3,
                      label="Missing wedge threshold",
                      expertLevel=LEVEL_ADVANCED,
                      help="Threshold for identifying missing data in Fourier"
                           " space in terms of standard deviation of each Fourier"
                           " shell. Default 3.0. If set to 0.0, missing wedge correction"
                           " will be skipped")
        form.addParam('skipPostProc', BooleanParam,
                      default=True,
                      label='Skip post process steps (fsc, mask and filters)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('keepHdfFile', BooleanParam,
                      default=False,
                      label='Keep hdf files?',
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to Yes, the generated files will be saved in both HDF and MRC formats. They are '
                           'generated in HDF and then converted into MRC. The HDF files are deleted by default to '
                           'save storage.')
        self._addBinThreads(form)

    # --------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.computeAverageStep)
        self._insertFunctionStep(self.converOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.projectPath = self._getExtraPath(SPT_00_DIR)
        self.hdfSubtomosDir = self._getExtraPath(SUBTOMOGRAMS_DIR)
        makePath(*[self.hdfSubtomosDir, self.projectPath])

    def convertInputStep(self):
        inSubtomos = self.inputSetOfSubTomogram.get()
        if type(inSubtomos) is EmanSetOfParticles:
            # If True, the 3D particles HDF stacks will already exist because they have been extracted with EMAN pppt
            hdfStacks = inSubtomos.getUniqueValues(EmanParticle.STACK_3D_HDF)
            [createLink(abspath(hdfStack), self._getExtraPath(SUBTOMOGRAMS_DIR, basename(hdfStack))) for
             hdfStack in hdfStacks]
        else:
            writeSetOfSubTomograms(inSubtomos, self.hdfSubtomosDir)
        refinement2Json(self, inSubtomos)

        # Generate a virtual stack of particle represented by a .lst file, as expected by EMAN
        program = emantomo.Plugin.getProgram('e2proclst.py')
        particleStacks = [join(SUBTOMOGRAMS_DIR, basename(partStack)) for partStack in glob.glob(join(self.hdfSubtomosDir, '*.hdf'))]
        args = f'{" ".join(particleStacks)} --create {join(SPT_00_DIR, INPUT_PTCLS_LST)}'
        self.runJob(program, args, cwd=self._getExtraPath())

    def computeAverageStep(self):
        args = "--keep 1 --wedgesigma=%f --sym %s --threads %i " % (self.msWedge.get(),
                                                                    self.symmetry.get(),
                                                                    self.binThreads.get())
        if self.skipPostProc.get():
            args += '--skippostp '
        program = emantomo.Plugin.getProgram('e2spt_average.py')
        self.runJob(program, args, cwd=self._getExtraPath())

    def converOutputStep(self):
        # Also fix the sampling rate as it might be set wrong (the value stored in the hdf header may be referred
        # to the original binning, and it will be also in the header of the resulting mrc file
        sRate = self.inputSetOfSubTomogram.get().getSamplingRate()
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        self.volumeFileHdf = self._getExtraPath(SPT_00_DIR, THREED_01 + '.hdf')
        self.halfEvenFileHdf = self._getExtraPath(SPT_00_DIR, THREED_01 + '_even.hdf')
        self.halfOddFileHdf = self._getExtraPath(SPT_00_DIR, THREED_01 + '_odd.hdf')
        filesToConvert = [self.volumeFileHdf, self.halfEvenFileHdf, self.halfOddFileHdf]
        for hdfFile in filesToConvert:
            args = "%s %s --apix %f" % (hdfFile, replaceExt(hdfFile, "mrc"), sRate)
            self.runJob(program, args)
            # Remove the hdf files if requested
            if not self.keepHdfFile.get():
                os.remove(hdfFile)

        # Finally, remove the hdf subtomograms generated in the convert input step, if requested
        if not self.keepHdfFile.get():
            shutil.rmtree(self.hdfSubtomosDir)

    def createOutputStep(self):
        mrcExt = 'mrc'
        inSubtomos = self.inputSetOfSubTomogram.get()
        volume = AverageSubTomogram()
        volume.setFileName(replaceExt(self.volumeFileHdf, mrcExt))
        volume.setHalfMaps([replaceExt(self.halfEvenFileHdf, mrcExt), replaceExt(self.halfOddFileHdf, mrcExt)])
        volume.setSamplingRate(inSubtomos.getSamplingRate())
        volume.fixMRCVolume()

        self._defineOutputs(**{OutputsAverageSubtomos.averageSubTomos.name: volume})
        self._defineSourceRelation(inSubtomos, volume)

    # --------------- INFO functions -------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'averageSubTomos'):
            summary.append("Average not ready yet.")
        else:
            summary.append("Average has been reconstructed.")
        return summary
