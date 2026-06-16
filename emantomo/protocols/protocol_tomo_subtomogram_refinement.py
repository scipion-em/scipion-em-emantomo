# coding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from glob import glob
import re
from os import remove
from os.path import abspath, basename, join
from pwem.objects import SetOfFSCs
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from emantomo.convert import writeSetOfSubTomograms, getLastParticlesParams, updateSetOfSubTomograms, emanFSCsToScipion
import emantomo
import pwem.constants as emcts
from pyworkflow.utils import createLink, Message
from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from .protocol_base import ProtEmantomoBase
from ..constants import SUBTOMOGRAMS_DIR, SPT_00_DIR, SYMMETRY_HELP_MSG


class EmanTomoRefinementOutputs(enum.Enum):
    subtomograms = SetOfSubTomograms
    subtomogramAverage = AverageSubTomogram
    FSCs = SetOfFSCs


class EmanProtTomoRefinement(ProtEmantomoBase):
    """
    Performs iterative subtomogram refinement and averaging using EMAN2.
    The protocol refines a collection of subtomograms against a reference
    volume in order to improve structural consistency, recover higher
    resolution information, and generate biologically meaningful averaged
    reconstructions from heterogeneous tomographic particle datasets.

    AI Generated:

    Subtomogram Refinement (EmanProtTomoRefinement) - User Manual
        Overview

        This protocol performs iterative subtomogram refinement and
        averaging within a cryo-electron tomography workflow using EMAN2.
        Its primary objective is to align and average individual
        subtomograms so that common structural information is reinforced
        while noise and random variability are progressively reduced.

        In cryo-electron tomography, subtomogram averaging is one of the
        most important strategies for improving the resolution of repeated
        macromolecular complexes embedded within tomograms. Because each
        subtomogram typically contains limited signal and substantial
        missing-wedge distortions, iterative refinement becomes essential
        for recovering biologically interpretable structures.

        Biological Purpose

        From a biological perspective, this protocol is designed for the
        structural analysis of macromolecular assemblies observed directly
        in cellular or native environments. Typical applications include
        ribosomes, membrane complexes, viral components, cytoskeletal
        assemblies, and other repeated particles extracted from tomograms.

        The protocol allows refinement of particle orientations and
        positions relative to an initial reference map. As refinement
        progresses, the averaged reconstruction becomes increasingly
        representative of the underlying biological structure while random
        noise and alignment inaccuracies are minimized.

        Inputs and General Workflow

        The protocol requires a set of subtomograms together with an
        initial reference volume. The reference serves as the starting
        structural hypothesis that guides alignment and averaging during
        iterative refinement.

        In most workflows, the subtomograms originate from particle
        extraction procedures performed on reconstructed tomograms. The
        reference may come from an initial model, a low-resolution average,
        a previously refined structure, or an external reconstruction.

        During execution, the protocol repeatedly aligns subtomograms to
        the reference, reconstructs updated averages, and evaluates the
        consistency of the particle population. Through multiple iterations,
        the refinement progressively improves the structural coherence of
        the dataset.

        Iterative Refinement Strategy

        Iterative refinement is central to subtomogram averaging because
        the quality of alignment strongly determines the final achievable
        resolution. Each iteration attempts to improve angular assignment,
        translational positioning, and consistency between particles.

        The number of refinement iterations should be selected according
        to dataset quality and structural complexity. Too few iterations
        may produce under-refined maps, whereas excessive refinement can
        increase the risk of overfitting or amplification of noise.

        In practice, users commonly begin with conservative settings and
        inspect the evolution of the reconstruction quality before extending
        the refinement further.

        Particle Selection and Data Quality

        The protocol includes options for selecting only a fraction of the
        particle population during refinement. This is biologically useful
        because low-quality particles, damaged complexes, or strongly
        heterogeneous conformations can negatively affect the final average.

        Keeping only the most consistent particles often improves the
        interpretability of the reconstruction, especially for datasets
        containing substantial variability or low signal-to-noise ratios.

        From a biological standpoint, particle exclusion should always be
        interpreted carefully. Excessively aggressive filtering may remove
        rare but meaningful structural states, whereas insufficient filtering
        may blur the reconstruction through structural heterogeneity.

        Symmetry and Structural Constraints

        Symmetry definition is one of the most influential parameters in
        subtomogram refinement. When the biological assembly possesses known
        rotational or point-group symmetry, applying the appropriate
        symmetry can substantially improve the effective signal and the
        quality of the final reconstruction.

        However, incorrect symmetry assignment can introduce severe
        structural artifacts and misleading biological interpretations.
        Symmetry should therefore only be imposed when supported by
        experimental evidence or prior structural knowledge.

        The protocol also supports masking strategies that focus refinement
        on biologically relevant regions of the structure. Proper masking
        can improve alignment stability and reduce the influence of solvent
        noise or flexible peripheral domains.

        Gold-Standard Refinement

        The protocol supports gold-standard refinement strategies intended
        to reduce overfitting and provide more reliable resolution
        estimation. In this approach, the dataset is divided into
        independent subsets that are refined separately before comparison.

        Gold-standard refinement is particularly important for publication-
        quality reconstructions because it provides a more rigorous
        assessment of structural reproducibility and resolution.

        Local filtering options may also be used to improve interpretability
        in regions with variable local resolution. Flexible domains or
        poorly sampled regions frequently benefit from adaptive filtering
        approaches.

        Tilt Geometry and Missing Wedge Effects

        Cryo-electron tomography data are intrinsically affected by missing
        wedge artifacts due to limited angular acquisition ranges. The
        protocol allows explicit control over tilt limitations so that
        reconstruction and refinement remain consistent with the geometry
        of the original acquisition.

        From a biological perspective, understanding missing wedge effects
        is essential because anisotropic resolution may influence the
        apparent shape and visibility of structural features. Careful
        interpretation is therefore required, especially for elongated or
        membrane-associated complexes.

        Outputs and Interpretation

        The protocol generates a refined subtomogram average representing
        the consensus structure derived from the aligned particle set. This
        average can be used for visualization, segmentation, interpretation,
        docking, or further high-resolution analysis.

        The workflow also produces updated subtomogram metadata containing
        refined alignment information and particle quality measures. These
        outputs are valuable for downstream classification, particle
        cleaning, or additional refinement cycles.

        Fourier shell correlation curves are additionally generated to
        support resolution estimation and reconstruction quality assessment.
        These measurements help determine the reliability and interpretability
        of structural features observed in the final average.

        Practical Recommendations

        In practical workflows, it is often beneficial to begin refinement
        using a low-resolution reference and relatively conservative angular
        searches. Once alignment stabilizes, refinement parameters can be
        progressively tightened to improve structural detail.

        For heterogeneous datasets, combining careful particle selection
        with biologically meaningful masking usually produces the largest
        improvements in map quality. Conversely, applying overly restrictive
        constraints too early may bias the reconstruction or suppress real
        structural variability.

        Users should visually inspect intermediate averages throughout the
        refinement process to ensure that the reconstruction evolves toward
        biologically plausible conformations.

        Final Perspective

        Subtomogram refinement is one of the central computational steps in
        cryo-electron tomography because it transforms noisy individual
        particles into interpretable structural averages. Successful
        refinement depends not only on computational optimization, but also
        on biologically informed decisions regarding symmetry, masking,
        particle selection, and structural heterogeneity. Careful refinement
        strategy and critical interpretation are therefore essential for
        obtaining reliable in situ structural information.
    """
    _outputClassName = 'SubTomogramRefinement'
    _label = 'subtomogram refinement'
    _possibleOutputs = EmanTomoRefinementOutputs

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.refFileIn = None
        self.alignType = None
        self.refFileOut = None
        self.whole3dStack = None

    # --------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      strict=True,
                      label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('inputRef', params.PointerParam,
                      pointerClass='Volume',
                      allowsNull=False,
                      default=None,
                      label='Input Ref SubTomogram',
                      help='3D reference for initial model generation.'
                           'No reference is used by default.')
        self._addBinThreads(form)

        form.addSection(label='Optimization')
        form.addParam('niter', params.IntParam, default=5,
                      label='Number of iterations',
                      help='The number of iterations to perform.')
        form.addParam('mass', params.FloatParam, default=-1,
                      label='Mass:',
                      help='Mass normalization. Default=-1 Ignores mass')
        form.addParam('pkeep', params.FloatParam, default=0.8,
                      label='Particle keep:',
                      help='Fraction of particles to keep')
        form.addParam('goldstandard', params.IntParam, default=-1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Gold standard:',
                      help='initial resolution for gold standard refinement')
        form.addParam('goldcontinue', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Gold continue',
                      help='continue from an existing gold standard refinement')
        form.addParam('maskFile', params.PointerParam,
                      allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      pointerClass='VolumeMask',
                      label='Mask file',
                      help='Mask file to be applied to initial model')
        form.addParam('setsf', params.PointerParam, allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      pointerClass='VolumeMask', label='Structure factor',
                      help='Select the structure factor')
        form.addParam('sym', params.StringParam, default='c1',
                      label='Symmetry',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('localfilter', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Local filter',
                      help='use tophat local')
        form.addParam('maxtilt', params.FloatParam, default=90.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='maxtilt',
                      help='Explicitly zeroes data beyond specified tilt angle.'
                           'Assumes tilt axis exactly on Y and zero tilt in X-Y'
                           'plane. Default 90 (no limit).')
        form.addParam('useAlign', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Use previous alignments?')
        form.addParam('maxAng', params.FloatParam, default=25,
                      condition='useAlign==True',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum angular change',
                      help='Maximum anglular difference in refine mode (in degrees)')
        form.addParam('extraParams', params.StringParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Extra params",
                      help='Here you can add any extra parameters to run Eman subtomogram refinement. '
                           'Parameters should be written in Eman command line format (--param=val)')

    # --------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.refinementSubtomogram)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.whole3dStack = self._getExtraPath(SUBTOMOGRAMS_DIR, 'whole3dstack.hdf')
        if self.useAlign.get():
            self.alignType = emcts.ALIGN_3D
        else:
            self.alignType = emcts.ALIGN_NONE

    def convertInputStep(self):
        inSubtomos = self.inputSetOfSubTomogram.get()
        storePath = self._getExtraPath(SUBTOMOGRAMS_DIR)
        pwutils.makePath(storePath)
        writeSetOfSubTomograms(inSubtomos, storePath, alignType=self.alignType)

        # Fix the sampling rate as it might be set wrong
        sampling_rate = inSubtomos.getSamplingRate()
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        for file in sorted(glob(join(storePath, '*.hdf'))):  # It's critical that the HDF stack keep the sqlite order
            # to later generate the output correctly in terms of transformation matrix assignment

            # Bypass the problematic lst generation with this classical approach by making one HDF file containing all
            # the particles from all the tomograms
            args = "--apix %f %s %s --append" % (sampling_rate, file, self.whole3dStack)
            self.runJob(program, args)
            remove(file)

        # The same with the reference volume
        self.refFileIn = self.inputRef.get().getFileName()
        self.refFileOut = self._getExtraPath('refVol.hdf')
        if self.refFileIn.endswith('.mrc'):
            args = "--apix %f %s %s" % (sampling_rate, self.refFileIn, self.refFileOut)
            self.runJob(program, args)
        else:
            createLink(self.refFileIn, self.refFileOut)

    def refinementSubtomogram(self):
        """Run the Subtomogram refinement. """
        args = ' %s' % self.whole3dStack
        args += ' --verbose=9'
        args += (' --reference=%s ' % abspath(self.refFileOut))
        args += (' --mass=%f' % self.mass.get())
        args += ' --goldstandard=%d ' % self.goldstandard.get()
        args += ' --pkeep=%f ' % self.pkeep.get()
        args += ' --sym=%s ' % self.sym.get()
        args += ' --maxtilt=%s ' % self.maxtilt
        args += ' --path=%s ' % abspath(self.getOutputPath())
        args += ' --niter=%d' % self.niter.get()
        if self.goldcontinue:
            args += ' --goldcontinue '
        if self.maskFile.get():
            args += ' --mask=%s' % abspath(self.maskFile.get().getFileName())
        if self.localfilter:
            args += ' --localfilter '
        if self.alignType != emcts.ALIGN_NONE:
            args += ' --refine --maxang=%f' % self.maxAng.get()
        if self.extraParams.get():
            args += ' ' + self.extraParams.get()
        args += ' --threads=%d' % self.binThreads.get()

        program = emantomo.Plugin.getProgram('e2spt_refine.py')
        self.runJob(program, args)

    def convertOutputStep(self):
        self.hdfToMrc("threed_\d+.hdf", self.getAverageFn())
        self.hdfToMrc("threed_even_unmasked.hdf", self.getEvenFn())
        self.hdfToMrc("threed_odd_unmasked.hdf", self.getOddFn())

    def createOutputStep(self):
        inputSetOfSubTomograms = self.inputSetOfSubTomogram.get()
        # Output 1: AverageSubTomogram
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(self.getAverageFn())
        averageSubTomogram.setHalfMaps([self.getEvenFn(), self.getOddFn()])
        averageSubTomogram.setSamplingRate(inputSetOfSubTomograms.getSamplingRate())
        # Output 2: setOfSubTomograms
        particleParams = getLastParticlesParams(self.getOutputPath())
        outputSetOfSubTomograms = self._createSet(SetOfSubTomograms, 'subtomograms%s.sqlite', "particles")
        outputSetOfSubTomograms.copyInfo(inputSetOfSubTomograms)
        outputSetOfSubTomograms.setCoordinates3D(inputSetOfSubTomograms.getCoordinates3D())
        updateSetOfSubTomograms(inputSetOfSubTomograms, outputSetOfSubTomograms, particleParams)
        # Output 3: FSCs
        fscs= self._createSet(SetOfFSCs, 'fsc%s.sqlite', "")
        fscMasked = self.getLastFromOutputPath('fsc_masked_\d+.txt')
        fscUnmasked = self.getLastFromOutputPath('fsc_unmasked_\d+.txt')
        fscTight = self.getLastFromOutputPath('fsc_maskedtight_\d+.txt')

        emanFSCsToScipion(fscs, fscMasked, fscUnmasked, fscTight)
        outputs = {EmanTomoRefinementOutputs.subtomogramAverage.name: averageSubTomogram,
                   EmanTomoRefinementOutputs.subtomograms.name: outputSetOfSubTomograms,
                   EmanTomoRefinementOutputs.FSCs.name: fscs}

        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputSetOfSubTomogram, averageSubTomogram)
        self._defineSourceRelation(self.inputSetOfSubTomogram, outputSetOfSubTomograms)

    # --------------- INFO functions -------------------------
    def _summary(self):
        summary = ["Subtomograms source: %s" % (self.inputSetOfSubTomogram.get().getFileName())]

        if self.inputRef.get() is not None:
            summary.append("Referenced Tomograms source: %s" % (self.inputRef.get().getFileName()))

        if self.getOutputsSize() >= 1:
            summary.append("Subtomogram refinement Completed")
        else:
            summary.append("Subtomogram refinement not ready yet.")

        return summary

    def _methods(self):
        inputSetOfSubTomgrams = self.inputSetOfSubTomogram.get()
        return [
            "Applied refinement using e2spt_refine (stochastic gradient descent)",
            "A total of %d particles of dimensions %s were used"
            % (inputSetOfSubTomgrams.getSize(), inputSetOfSubTomgrams.getDimensions()),
        ]

    # --------------------------- UTILS functions ------------------------------
    def hdfToMrc(self, pattern, mrcName):
        # Fix the sampling rate as it might be set wrong
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        lastImage = self.getLastFromOutputPath(pattern)
        args = "--apix %f %s %s" % (self.inputSetOfSubTomogram.get().getSamplingRate(),
                                    lastImage, mrcName)
        self.runJob(program, args)

    def getAverageFn(self):
        return self._getExtraPath("Average_refined.mrc")

    def getEvenFn(self):
        return self._getExtraPath("even.mrc")

    def getOddFn(self):
        return self._getExtraPath("odd.mrc")

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob(self.getOutputPath("*"))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]

    def getOutputPath(self, *args):
        return join(self._getExtraPath(SPT_00_DIR, *args))

