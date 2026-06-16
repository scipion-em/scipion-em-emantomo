# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *             Scipion Team (scipion@cnb.csic.es) [1]
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
from os.path import join, basename
from emantomo import Plugin
from emantomo.constants import TOMOGRAMS_DIR
from emantomo.protocols.protocol_base import IN_TOMOS, ProtEmantomoBase
from pwem.convert.headers import setMRCSamplingRate
from pwem.objects import Transform
from pyworkflow.protocol.params import PointerParam, IntParam, GT
from pyworkflow.utils import runJob, makePath, createLink
from pyworkflow.utils.properties import Message
from tomo.objects import SetOfTomograms, Tomogram


class clipTomogramsOuts(Enum):
    tomograms = SetOfTomograms


class EmanProtTomoClip(ProtEmantomoBase):
    """
    Adjusts tomogram dimensions and spatial positioning by clipping or padding volumetric data while preserving
    the biological region of interest within a consistent coordinate system.

    AI Generated:

    Tomogram Clipping and Re-Centering (EmanProtTomoClip) — User Manual
        Overview

        The Tomogram Clipping protocol modifies the dimensions and center position of tomograms in order to
        standardize datasets for downstream cryo-electron tomography analysis. Its main objective is to generate
        volumes with consistent sizes and spatial alignment while preserving the biologically relevant structures
        contained within the tomograms.

        In practical cryo-ET workflows, tomograms frequently originate from different acquisition sessions,
        reconstruction pipelines, or preprocessing strategies. As a result, their dimensions and coordinate
        systems may vary considerably. This protocol provides a controlled way to crop excessive regions,
        reduce unnecessary empty space, or expand volumes through padding so that subsequent processing steps
        can operate on standardized inputs.

        Biological Motivation

        Standardizing tomogram dimensions is particularly important when combining datasets for subtomogram
        averaging, particle picking, template matching, neural network segmentation, or comparative structural
        analysis. Large differences in volume dimensions often complicate automated processing pipelines and
        may introduce inconsistencies during alignment or classification.

        Re-centering is equally important from a biological perspective. Cellular structures, macromolecular
        assemblies, or regions of interest are not always positioned near the geometric center of the tomogram.
        By defining a new center position, the protocol allows the user to focus the resulting volume on the
        biologically relevant area while excluding unnecessary regions dominated by noise or empty solvent.

        Inputs and General Workflow

        The protocol requires a set of tomograms as input. Each tomogram is processed independently, allowing
        datasets with heterogeneous dimensions to be standardized within the same workflow. The protocol
        preserves the original voxel sampling rate so that downstream quantitative interpretation remains
        biologically meaningful.

        Users may optionally define new dimensions along the X, Y, and Z axes. When dimensions are reduced,
        the tomogram is clipped around the selected center region. When dimensions are increased, the protocol
        pads the tomogram to create a larger output volume. If dimensions are left undefined, the original
        size along those axes is preserved.

        The protocol also allows users to define a new center coordinate. This capability is especially useful
        when the structure of interest is displaced from the original tomogram center or when preparing focused
        subregions for detailed analysis. If no new center is specified, the original center of the tomogram
        is retained.

        Choosing Appropriate Dimensions

        Selecting appropriate output dimensions is biologically important because excessive clipping may remove
        meaningful structural information. In cellular tomography, relevant contextual features such as membranes,
        neighboring complexes, or cytoskeletal elements may extend beyond the immediate target region. Users
        should therefore ensure that the selected dimensions preserve the biological context required for the
        intended analysis.

        Conversely, unnecessarily large tomograms increase storage requirements and computational cost. Reducing
        empty regions often improves efficiency in downstream subtomogram extraction, alignment, and averaging.
        A balanced strategy is typically recommended, preserving sufficient contextual information while avoiding
        excessive empty space.

        Re-Centering and Coordinate Interpretation

        When a new center is introduced, the protocol updates the spatial origin of the resulting tomogram so
        that coordinate consistency is maintained throughout the workflow. This behavior is particularly relevant
        for particle coordinates, segmentation masks, or annotations that depend on accurate spatial referencing.

        From a biological perspective, careful interpretation of the new center is important. Re-centering does
        not modify the biological structure itself, but it changes the spatial frame in which the tomogram is
        represented. This is especially relevant when comparing tomograms across experiments or integrating
        information with external coordinate systems.

        Practical Recommendations

        In most cryo-ET workflows, it is advisable to inspect tomograms visually before selecting clipping
        parameters. Regions containing fiducials, reconstruction artifacts, or large empty solvent regions are
        often suitable candidates for exclusion. At the same time, users should avoid aggressive cropping that
        may truncate flexible domains, membrane boundaries, or neighboring structural features.

        For subtomogram averaging projects, generating tomograms with standardized dimensions can simplify
        extraction and improve reproducibility across datasets. For focused biological studies, re-centering
        around the region of interest often facilitates visualization and downstream interpretation.

        Final Perspective

        Tomogram clipping and re-centering are not merely geometric manipulations but important preprocessing
        operations that influence the biological interpretability and computational efficiency of cryo-electron
        tomography workflows. Careful selection of dimensions and spatial centers helps preserve meaningful
        structural information while producing cleaner and more standardized datasets for downstream analysis.
    """
    IN_TOMOGRAMS_DIR = 'inTomograms'
    _possibleOutputs = clipTomogramsOuts
    _label = 'clip tomograms'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.outSet = None
        self.program = None
        self.samplingRate = None
        self.tomoSizeAndCenterDict = {}

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input tomograms",
                      important=True)

        line = form.addLine('New center coordinates [px] (opt.)',
                            help='For each empty coordinate, the correspònding original center coordinate '
                                 'will be used.')
        line.addParam('xc', IntParam,
                      label="cx",
                      allowsNull=True)
        line.addParam('yc', IntParam,
                      label="cy",
                      allowsNull=True)
        line.addParam('zc', IntParam,
                      label="cz",
                      allowsNull=True)

        line = form.addLine('New dimensions [px] (opt.)',
                            help='For each empty dimension, the original corresponding dimension will be used.')
        line.addParam('xDim', IntParam,
                      label="dx",
                      allowsNull=True)
        line.addParam('yDim', IntParam,
                      label="dy",
                      allowsNull=True)
        line.addParam('zDim', IntParam,
                      label="dz",
                      allowsNull=True)

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        tomoList = self._initilize()
        for tomo in tomoList:
            self._insertFunctionStep(self.convertInputStep, tomo)
            self._insertFunctionStep(self.clipTomogramsStep, tomo)
            self._insertFunctionStep(self.createOutputStep, tomo)
        self._insertFunctionStep(self.closingStep)

    def _initilize(self):
        makePath(self.getInTomosDir(), self.getTomogramsDir())
        inTomos = self.getObjByName(IN_TOMOS)
        self.program = Plugin.getProgram("e2proc3d.py")
        self.samplingRate = inTomos.getSamplingRate()
        newDimsAndCenter = [self.xDim.get(), self.yDim.get(), self.zDim.get(),
                            self.xc.get(), self.yc.get(), self.zc.get()]
        tomoList = [tomo.clone() for tomo in inTomos]
        # We'll do tomo by tomo instead of considering all the tomograms from the set to be of the same size
        for tomo in tomoList:
            self.tomoSizeAndCenterDict[tomo.getTsId()] = self._getSizeAndCenter(tomo, newDimsAndCenter)
        return tomoList

    def convertInputStep(self, tomo):
        inFName = tomo.getFileName()
        createLink(inFName, self._genLinkedFileName(inFName))

    def clipTomogramsStep(self, tomo):
        inFile = tomo.getFileName()
        sizeAndCenter = self.tomoSizeAndCenterDict[tomo.getTsId()]
        args = f'{self._getInFileName(inFile)} {self._getOutFileName(inFile)} ' \
               f'--clip {",".join([str(int(val)) for val in sizeAndCenter])}'  # EMAN expects these values as integers
        runJob(None, self.program, args, env=Plugin.getEnviron(), cwd=self._getExtraPath())

    def createOutputStep(self, tomo):
        sr = self.samplingRate
        # Create the output set if it does not exist yet
        if not self.outSet:
            inTomos = self.getObjByName(IN_TOMOS)
            self.outSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            self.outSet.copyInfo(inTomos)
            self.outSet.setSamplingRate(sr)
        outTomo = Tomogram()
        outTomo.copyInfo(tomo)
        outFName = self._getExtraPath(self._getOutFileName(tomo.getFileName()))
        setMRCSamplingRate(outFName, sr)
        outTomo.setFileName(outFName)
        # Update the origin if necessary
        if self._shiftCenter():
            sizeAndCenter = self.tomoSizeAndCenterDict[tomo.getTsId()]
            sx, sy, sz = sizeAndCenter[3::]
            origin = Transform()
            origin.setShifts(-sx * sr,
                             -sy * sr,
                             -sz * sr)
            outTomo.setOrigin(origin)
        self.outSet.append(outTomo)

    def closingStep(self):
        self._defineOutputs(**{self._possibleOutputs.tomograms.name: self.outSet})
        self._defineSourceRelation(getattr(self, IN_TOMOS), self.outSet)

    # --------------------------- INFO functions -----------------------------
    def _methods(self):
        methodsMsgs = ["Tomogram clipping using e2proc3d.py"]
        return methodsMsgs

    def _validate(self):
        errorMsg = []
        inVals = [self.xDim.get(), self.yDim.get(), self.zDim.get(), self.xc.get(), self.yc.get(), self.zc.get()]
        if all(x is None for x in inVals):
            errorMsg.append('At least one of the new center coordinates or the new dimensions must be filled.')
        else:
            for val in inVals:
                if val is not None:
                    if val < 1:
                        errorMsg.append('The introduced values must be greater than 1.')
                        break
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _shiftCenter(self):
        return True if self.xc.get() or self.yc.get() or self.zc.get() else False

    @staticmethod
    def _getSizeAndCenter(tomo, newDimsAndCenter):
        # sr = self.samplingRate
        origXDim, origYDim, origZDim = tomo.getDimensions()  # In pixels
        origVals = [origXDim, origYDim, origZDim,
                    origXDim / 2, origYDim / 2, origZDim / 2]  # In pixels
        valSize = len(origVals)
        finalVals = [None] * valSize  # List pre-allocating
        for i, newVal in enumerate(newDimsAndCenter):
            finalVals[i] = newVal if newVal else origVals[i]
        return finalVals

    def getInTomosDir(self):
        return self._getExtraPath(self.IN_TOMOGRAMS_DIR)

    def _genLinkedFileName(self, inFile):
        return join(self.getInTomosDir(), basename(inFile))

    @staticmethod
    def _getOutFileName(inFile):
        return join(TOMOGRAMS_DIR, basename(inFile))

    def _getInFileName(self, inFile):
        return join(self.IN_TOMOGRAMS_DIR, basename(inFile))
