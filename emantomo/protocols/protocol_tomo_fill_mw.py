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
    This protocol is a wrapper of *e2tomo_mwfill*.

    This program can be used to fill in the missing wedge of a SetOfTomograms with useful information
    based on a neural network.

    Fill the missing wedge of Tomograms with somewhat meaningful information using a deep learning based tool.
    The idea is similar to a "style transform" that makes the features in the x-z 2D slice views similar to
    the x-y slice views.'
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
