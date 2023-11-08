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
from pwem.objects import Transform
from pyworkflow.protocol.params import PointerParam, IntParam, GT
from pyworkflow.utils import runJob, makePath, createLink
from pyworkflow.utils.properties import Message
from tomo.objects import SetOfTomograms, Tomogram


class clipTomogramsOuts(Enum):
    tomograms = SetOfTomograms


class EmanProtTomoClip(ProtEmantomoBase):
    """Make the output have this size by padding/clipping and re-centered using e2proc3d.py. Both
    new center and new dimensions are referred to the original position and dimensions, respectively."""
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

        line = form.addLine('New center coordinates [pix.] (opt.)',
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

        line = form.addLine('New dimensions [pix.] (opt.)',
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
        outTomo.setFileName(self._getExtraPath(self._getOutFileName(tomo.getFileName())))
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
