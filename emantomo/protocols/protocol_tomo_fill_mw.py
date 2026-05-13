# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from os.path import basename, join, splitext
from random import randint

from emantomo.constants import TOMOGRAMS_DIR
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TOMOS
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils import makePath, createLink, replaceExt
from pyworkflow.utils.properties import Message
from tomo.objects import Tomogram, SetOfTomograms
import emantomo


class fillMWOutputs(Enum):
    tomograms = SetOfTomograms


class MWDataObj:

    def __init__(self, tomo=None, inLinkedScipionFile=None, inEmanFile=None, outEmanFile=None, outScipionFile=None):
        self.tomo = tomo
        self.inLinkedScipionFile = inLinkedScipionFile
        self.inEmanFile = inEmanFile
        self.outEmanFile = outEmanFile
        self.outScipionFile = outScipionFile


class EmanProtTomoFillMW(ProtEmantomoBase):
    """
    Restores missing wedge information in tomographic reconstructions using a
    deep learning strategy based on EMAN2 tools. The protocol improves the
    interpretability of tomograms affected by incomplete angular sampling by
    generating more isotropic structural information across different viewing
    directions.

    AI Generated:

    Tomogram Missing Wedge Filling (EmanProtTomoFillMW) — User Manual
        Overview

        The Tomogram Missing Wedge Filling protocol is designed to reduce the
        visual and structural artifacts produced by the missing wedge effect in
        cryo-electron tomography datasets. In tomographic reconstruction, the
        limited tilt range imposed during acquisition leaves regions of Fourier
        space unsampled, creating anisotropic resolution and distortions that
        can complicate visualization, segmentation, interpretation, and
        downstream analysis.

        This protocol applies a deep learning based strategy that attempts to
        infer biologically meaningful information within the missing wedge
        region. Rather than performing a classical reconstruction refinement,
        the method learns structural patterns from existing tomographic content
        and transfers this information into poorly sampled orientations. The
        resulting tomograms often display improved continuity of membranes,
        filaments, macromolecular assemblies, and cellular densities.

        Biological Purpose and Typical Applications

        In biological cryo-electron tomography workflows, the missing wedge can
        strongly distort elongated or directional structures. Membranes may
        appear artificially stretched, filaments can become fragmented, and
        structural interpretation along the Z direction may become unreliable.
        This protocol is intended to alleviate these limitations and provide
        tomograms that are easier to inspect visually and more suitable for
        segmentation, annotation, particle picking, or qualitative biological
        interpretation.

        Typical applications include cellular tomography, membrane remodeling
        studies, cytoskeletal organization, organelle visualization, and
        exploratory subtomogram analysis. The protocol is particularly useful
        when downstream tasks depend heavily on visual continuity and isotropic
        appearance of structural features.

        General Workflow

        The protocol takes a set of tomograms as input. One tomogram from the
        dataset is automatically selected for neural network training, and the
        learned representation is subsequently applied to the remaining
        tomograms in the set. This approach assumes that the tomograms share
        sufficiently similar imaging characteristics and biological content.

        During training, the method extracts many small volumetric regions from
        the selected tomogram. These examples are used to learn correlations
        between well-sampled and poorly sampled directions. Once training is
        completed, the model predicts improved information for the missing wedge
        regions of all tomograms in the dataset.

        The protocol preserves the original tomographic geometry and sampling
        while generating corrected output volumes suitable for visualization and
        further processing.

        Training Volume Selection

        Because only one tomogram is used for training, choosing a homogeneous
        dataset is biologically important. The best results are usually obtained
        when all tomograms originate from similar acquisition conditions,
        contain related biological structures, and share comparable contrast and
        noise properties.

        If the dataset contains highly heterogeneous specimens, the learned
        representation may not generalize well across all tomograms. In such
        cases, splitting the dataset into biologically consistent subsets before
        processing may improve the quality of the results.

        Box Size and Structural Context

        The box size parameter determines the size of the volumetric regions
        used during neural network training. Smaller box sizes focus on local
        texture and fine structural details, whereas larger box sizes capture
        broader contextual information and extended biological features.

        For compact macromolecular environments, moderate box sizes are often
        sufficient. Larger structures such as membranes, vesicles, or filament
        networks may benefit from larger contextual windows. Excessively large
        boxes, however, may increase computational cost and reduce training
        efficiency.

        Number of Samples and Learning Stability

        The number of extracted samples influences how extensively the protocol
        explores structural variability within the tomogram. Increasing the
        number of samples generally improves robustness and representation of
        biological diversity, particularly in crowded cellular environments.

        However, larger sample counts also increase processing time and GPU
        memory usage. In many routine workflows, intermediate values provide a
        practical balance between reconstruction quality and computational cost.

        Learning Rate Considerations

        The learning rate controls how rapidly the neural network adapts during
        optimization. Conservative learning rates typically produce more stable
        and biologically coherent results, especially in noisy tomographic
        datasets. Larger values may accelerate convergence but can also produce
        unstable or overly artificial reconstructions.

        For most biological applications, the default configuration is a
        reasonable starting point. Advanced users may adjust the learning rate
        when working with particularly noisy data or highly specialized imaging
        conditions.

        GPU Usage and Computational Requirements

        The protocol relies on GPU acceleration for neural network training and
        prediction. Processing speed depends strongly on GPU memory and hardware
        capabilities. Since training and inference can be computationally
        intensive, users should ensure that sufficient GPU resources are
        available before execution.

        The protocol is intended primarily for exploratory enhancement and
        interpretation rather than direct quantitative reconstruction accuracy.
        Biological conclusions should therefore always be validated against the
        original tomograms whenever possible.

        Outputs and Interpretation

        The protocol produces a new set of tomograms with reduced missing wedge
        artifacts. These corrected tomograms retain the metadata and structural
        organization of the original dataset while presenting more isotropic
        visual information.

        The outputs are often easier to interpret visually and may improve
        downstream segmentation or annotation tasks. Nevertheless, users should
        remember that the restored information is computationally inferred
        rather than experimentally measured. Features introduced by the neural
        network should therefore be interpreted cautiously, especially in
        regions with very limited original information.

        Practical Recommendations

        In routine biological practice, it is advisable to begin with default
        parameters and visually compare the corrected tomograms against the
        original data. The protocol is particularly valuable when the missing
        wedge strongly interferes with visualization of membranes, filaments,
        or directional assemblies.

        For homogeneous datasets acquired under stable conditions, the protocol
        can substantially improve interpretability. In contrast, highly diverse
        datasets may require separate processing groups to avoid inconsistent
        predictions.

        Final Perspective

        Missing wedge correction represents an important complementary tool in
        cryo-electron tomography analysis. Although it does not replace the
        original experimental information, it can significantly enhance the
        usability and interpretability of tomographic datasets. Careful
        biological validation, comparison with raw reconstructions, and
        thoughtful interpretation remain essential for reliable scientific
        conclusions.
    """
    _label = 'tomo fill missing wedge'
    _possibleOutputs = fillMWOutputs
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomoFiles = []
        self.sRate = None
        self.emanMdList = []

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addLine('IMPORTANT: This EMAN protocol uses GPU Id 0 and does not allow to select another.')
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Tomograms for training",
                      important=True,
                      help='One of the tomograms from the set will be used for to train and the result will be '
                           'applied to all the tomograms from the introduced set.')
        form.addParam('boxSize', IntParam,
                      default=64,
                      label='Box size of the training volumes.')
        form.addParam('nSamples', IntParam,
                      default=2000,
                      label='Number of samples to extract.')
        form.addParam('learnRate', FloatParam,
                      default=2e-4,
                      label='Learning rate',
                      expertLevel=LEVEL_ADVANCED)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.fillMWStep)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTomos = self.getAttrib(IN_TOMOS)
        inTomoList = [tomo.clone() for tomo in inTomos]
        for tomo in inTomoList:
            inEmanFile = self._getFilePathForEman(tomo.getFileName())
            outEmanFile = splitext(inEmanFile)[0] + '_mw.hdf'
            self.emanMdList.append(
                MWDataObj(
                    tomo=tomo,
                    inLinkedScipionFile=self._getExtraPath(inEmanFile),
                    inEmanFile=inEmanFile,
                    outEmanFile=outEmanFile,
                    outScipionFile=self._getExtraPath(replaceExt(outEmanFile, 'mrc'))
                )
            )
            self.tomoFiles.append(inEmanFile)
        self.sRate = inTomos.getSamplingRate()
        makePath(self.getTomogramsDir())

    def convertInputStep(self):
        [createLink(mdObj.tomo.getFileName(), mdObj.inLinkedScipionFile) for mdObj in self.emanMdList]

    def fillMWStep(self):
        args = '--train %s ' % self._getTomoFile4Training()
        args += '--applyto %s ' % ','.join(self.tomoFiles)
        args += '--boxsz %i ' % self.boxSize.get()
        args += '--nsample %i ' % self.nSamples.get()
        args += '--learnrate %s ' % self.learnRate.get()
        program = emantomo.Plugin.getProgram("e2tomo_mwfill.py")
        self.runJob(program, args, cwd=self._getExtraPath())

    def convertOutputStep(self):
        [self.convertBetweenHdfAndMrc(mdObj.outEmanFile, replaceExt(mdObj.outEmanFile, 'mrc'), extraArgs=f'--apix {self.sRate:.3f}')
         for mdObj in self.emanMdList]

    def createOutputStep(self):
        inTomosPointer = self.getAttrib(IN_TOMOS, getPointer=True)
        outTomos = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        outTomos.copyInfo(inTomosPointer.get())
        for mdObj in self.emanMdList:
            outTomo = Tomogram()
            outTomo.copyInfo(mdObj.tomo)
            outTomo.setFileName(mdObj.outScipionFile)
            outTomos.append(outTomo)

        self._defineOutputs(**{self._possibleOutputs.tomograms.name: outTomos})
        self._defineSourceRelation(inTomosPointer, outTomos)

    # --------------------------- UTIL functions -----------------------------------
    @staticmethod
    def _getFilePathForEman(iFile):
        return join(TOMOGRAMS_DIR, basename(iFile))

    def _getTomoFile4Training(self):
        """Pick a random tomograms from the set if it has more than one element"""
        nTomos = len(self.tomoFiles)
        return self.tomoFiles[0] if nTomos == 1 else self.tomoFiles[randint(0, nTomos - 1)]
