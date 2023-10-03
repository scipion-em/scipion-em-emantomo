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
import os
from enum import Enum
from os.path import basename, join

from emantomo.constants import TOMOGRAMS_DIR
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TOMOS
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils import makePath, createLink, removeBaseExt
from pyworkflow.utils.properties import Message
from pwem.protocols import EMProtocol
from tomo.objects import Tomogram, SetOfTomograms
import emantomo


class fillMWOutputs(Enum):
    tomograms = SetOfTomograms


class EmanProtTomoFillMW(ProtEmantomoBase):
    """
    This protocol is a wrapper of *e2tomo_mwfill*.

    This program can be used to fill in the missing wedge of a SetOfTomograms with useful information
    based on a neural network
    """
    _label = 'tomo fill missing wedge'
    _possibleOutputs = fillMWOutputs
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.trainTomoFiles = None
        self.tomos2DenoiseFiles = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Tomograms for training",
                      important=True)
        form.addParam('tomos2denoise', PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Tomograms to be denoised (opt.)",
                      allowsNull=True,
                      help='If empty, they will be the same as the tomograms used for training.')
        form.addParam('boxSize', IntParam,
                      default=64,
                      label='Box size of the training volumes')
        form.addParam('nSamples', IntParam,
                      default=2000,
                      label='Number of samples to extract')
        form.addParam('learnRate', FloatParam,
                      default=2e-4,
                      label='Learning rate',
                      expertLevel=LEVEL_ADVANCED)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.fillMWStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        makePath(self.getTomogramsDir())

    def convertInputStep(self):
        tomos2Denoise = self.tomos2denoise.get()
        trainTomos = self.getAttrib(IN_TOMOS)
        sRate = trainTomos.getSamplingRate()
        trainTomoFiles = [tomo.getFileName() for tomo in trainTomos]
        tomo2DenoiseFiles = [tomo.getFileName() for tomo in tomos2Denoise] if tomos2Denoise else trainTomoFiles[:]
        files2link = set(trainTomoFiles + tomo2DenoiseFiles)
        [self.convertOrLink(iFile, removeBaseExt(iFile), TOMOGRAMS_DIR, sRate) for iFile in files2link]
        # [createLink(iFile, self._getExtraPath(TOMOGRAMS_DIR, basename(iFile))) for iFile in files2link]
        self.trainTomoFiles = [self._getFilePathForEman(iFile) for iFile in trainTomoFiles]
        self.tomos2DenoiseFiles = [self._getFilePathForEman(iFile) for iFile in tomo2DenoiseFiles] if not (
            tomos2Denoise) else self.trainTomoFiles[:]

    def fillMWStep(self):
        args = '--train %s ' % ','.join(self.trainTomoFiles)
        args += '--applyto %s ' % ','.join(self.tomos2DenoiseFiles)
        args += '--boxsz %i ' % self.boxSize.get()
        args += '--nsample %i ' % self.nSamples.get()
        args += '--learnrate %s ' % self.learnRate.get()
        program = emantomo.Plugin.getProgram("e2tomo_mwfill.py")
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        pass
        # tilt_series = self.tiltSeries.get()
        #
        # sampling_rate = tilt_series.getSamplingRate()
        # if self.outsize.get() == 0:
        #     sampling_rate *= 4
        # elif self.outsize.get() == 1:
        #     sampling_rate *= 2
        #
        # # Output 1: Main tomograms
        # tomograms_paths = self._getOutputTomograms()
        # tomograms = self._createSet(SetOfTomograms, 'tomograms%s.sqlite', "")
        # tomograms.copyInfo(tilt_series)
        # tomograms.setSamplingRate(sampling_rate)
        #
        # for tomogram_path in tomograms_paths:
        #     self._log.info('Main tomogram: ' + tomogram_path)
        #     tomogram = Tomogram()
        #     tomogram.setFileName(tomogram_path)
        #     tomogram.copyInfo(tilt_series)
        #     tomogram.setSamplingRate(sampling_rate)
        #     tomograms.append(tomogram)
        #
        # self._defineOutputs(tomograms=tomograms)
        # self._defineSourceRelation(self.tiltSeries, tomograms)

    # --------------------------- INFO functions -----------------------------
    def _methods(self):
        return [
            'Fill the missing wedge of Tomograms with somewhat meaningful information using a '
            'deep learning based tool. ',
            'The idea is similar to a "style transform" that makes the features in the x-z '
            '2D slice views similar to the x-y slice views'
        ]

    def getInfo(self, output):
        msg = '\t A total of *%d* were process to fill their missing wedge' % output.getSize()
        return msg

    def _summary(self):
        summary = []
        if self.getOutputsSize() < 1:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        else:
            for key, output in self.iterOutputAttributes():
                msg = self.getInfo(output)
                summary.append("%s: \n %s" % (self.getObjectTag(output), msg))
        return summary

    @staticmethod
    def _getFilePathForEman(iFile):
        return join(TOMOGRAMS_DIR, removeBaseExt(iFile) + '.hdf')
        # return join(TOMOGRAMS_DIR, basename(iFile))
