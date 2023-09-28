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
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils import runJob, makePath, createLink
from pyworkflow.utils.properties import Message
from tomo.objects import SetOfTomograms, Tomogram


class clipTomogramsOuts(Enum):
    tomograms = SetOfTomograms


class EmanProtTomoClip(ProtEmantomoBase):
    """Make the output have this size by padding/clipping using e2proc3d.py"""
    IN_TOMOGRAMS_DIR = 'inTomograms'
    _possibleOutputs = clipTomogramsOuts
    _label = 'clip tomograms'
    _devStatus = BETA

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

        group = form.addGroup('New size')
        group.addParam('xDim', IntParam,
                       label="New X dimension (pix.)",
                       allowsNull=True,
                       help='If empty, the original dimension will be used.')
        group.addParam('yDim', IntParam,
                       label="New Y dimension (pix.)",
                       allowsNull=True,
                       help='If empty, the original dimension will be used.')
        group.addParam('zDim', IntParam,
                       label="New Z dimension (pix.)",
                       allowsNull=True,
                       help='If empty, the original dimension will be used.')

        group = form.addGroup('New center')
        group.addParam('xc', IntParam,
                       label="New center x coordinate (pix.)",
                       allowsNull=True,
                       help='If empty, the original center current coordinate will be used.')
        group.addParam('yc', IntParam,
                       label="New center y coordinate (pix.)",
                       allowsNull=True,
                       help='If empty, the original center current coordinate will be used.')
        group.addParam('zc', IntParam,
                       label="New center z coordinate (pix.)",
                       allowsNull=True,
                       help='If empty, the original center current coordinate will be used.')

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
        if not self.outSet:
            inTomos = self.getObjByName(IN_TOMOS)
            self.outSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            self.outSet.copyInfo(inTomos)
            self.outSet.setSamplingRate(sr)
        outTomo = Tomogram()
        outTomo.setFileName(self._getExtraPath(self._getOutFileName(tomo.getFileName())))
        outTomo.setSamplingRate(sr)
        ix, iy, iz = tomo.getShiftsFromOrigin()
        sx, sy, sz = self.tomoSizeAndCenterDict[tomo.getTsId()][3::]  # The first 3 elements are the new size,
        # while the 3 last are the new center
        # TODO: check if an invert Y axis options should be included for the de-centering, depending on the ref system were the use measured the new center
        outTomo.setShiftsInOrigin(x=(ix + sx) * sr,
                                  y=(iy + sy) * sr,
                                  z=(iz + sz) * sr)
        self.outSet.append(outTomo)

    def closingStep(self):
        self._defineOutputs(**{self._possibleOutputs.tomograms.name: self.outSet})
        self._defineSourceRelation(getattr(self, IN_TOMOS), self.outSet)

    # --------------------------- INFO functions -----------------------------
    def _methods(self):
        methodsMsgs = ["Tomogram clipping using e2proc3d.py"]
        return methodsMsgs

    # def _summary(self):
    #     summary = []
    #     if self.getOutputsSize() < 1:
    #         summary.append(Message.TEXT_NO_OUTPUT_CO)
    #     else:
    #         summary.append("A total of %d Tomograms were resized to dimensions (%d,%d,%d)"
    #                        % (self.resizedTomograms.getSize(), self.xDim.get(),
    #                           self.yDim.get(), self.zDim.get()))
    #     return summary

    # --------------------------- UTILS functions -----------------------------

    def _getSizeAndCenter(self, tomo, newDimsAndCenter):
        sr = self.samplingRate
        origXDim, origYDim, origZDim = tomo.getDimensions()  # In pixels
        origVals = [origXDim, origYDim, origZDim,
                    origXDim / 2, origYDim / 2, origZDim / 2]   # In pixels
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
