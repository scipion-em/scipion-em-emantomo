# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
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
import glob

from pyworkflow import BETA
from pyworkflow.protocol import params
from pyworkflow.utils.properties import Message

from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms

import emantomo


class EmanProtTomoFillMW(ProtTomoBase, EMProtocol):
    """
    This protocol is a wrapper of *e2tomo_mwfill*.

    This program can be used to fill in the missing wedge of a SetOfTomograms with useful information
    based on a neural network
    """
    _label = 'tomo fill mw'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input Tomograms", important=True,
                      help='Select the SetOfTomograms to be processed')

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createCommandStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def createCommandStep(self):
        tomograms = self.inputTomograms.get()
        files = [os.path.abspath(tomo) for tomo in tomograms.getFiles()]
        files_join = ",".join(files)
        args = "--train %s --applyto %s" % (files[0], files_join)
        program = emantomo.Plugin.getProgram("e2tomo_mwfill.py")
        self._log.info('Launching: ' + program + ' ' + args)
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
