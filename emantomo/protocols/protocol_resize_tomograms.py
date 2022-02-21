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


from pwem.protocols import EMProtocol

from pyworkflow import BETA
import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils.properties import Message

from tomo.protocols import ProtTomoBase

import emantomo


class EmanProtTomoResize(ProtTomoBase, EMProtocol):
    """ Tomogram resize using e2proc3d """
    _label = 'resize tomos'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTomograms', PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input tomograms", important=True,
                      help='Select the tomograms to be resized')
        form.addParam('xDim', IntParam, important=True, label="New X dimension",
                      help='Select the new X dimension for the Tomogram')
        form.addParam('yDim', IntParam, important=True, label="New Y dimension",
                      help='Select the new Y dimension for the Tomogram')
        form.addParam('zDim', IntParam, important=True, label="New Z dimension",
                      help='Select the new Z dimension for the Tomogram')

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.resizeTomograms)
        self._insertFunctionStep(self.createOutputStep)

    def resizeTomograms(self):
        program = emantomo.Plugin.getProgram("e2proc3d.py")
        xDim = self.xDim.get()
        yDim = self.yDim.get()
        zDim = self.zDim.get()
        for tomo in self.inputTomograms.get().iterItems():
            tomo_file = tomo.getFileName()
            out_file = self._getExtraPath(pwutils.removeBaseExt(tomo_file) + ".mrc")
            args = "%s %s --process normalize --clip %d,%d,%d" \
                   % (tomo_file, out_file, xDim, yDim, zDim)
            pwutils.runJob(None, program, args, env=emantomo.Plugin.getEnviron())

    def createOutputStep(self):
        input_tomos = self.inputTomograms.get()
        sr = input_tomos.getSamplingRate()
        out_tomos = self._createSetOfTomograms()
        out_tomos.copyInfo(input_tomos)
        out_tomos.setSamplingRate(sr)

        for tomo in input_tomos.iterItems():
            out_tomo = tomo.clone()
            tomo_file = self._getExtraPath(pwutils.removeBaseExt(tomo.getFileName()) + ".mrc")
            out_tomo.setFileName(tomo_file)
            out_tomo.setSamplingRate(sr)
            out_tomos.append(out_tomo)

        self._defineOutputs(resizedTomograms=out_tomos)
        self._defineSourceRelation(self.inputTomograms, out_tomos)

    # --------------------------- INFO functions -----------------------------
    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("Tomogram resizing by e2proc3d.")
        return methodsMsgs

    def _summary(self):
        summary = []
        if self.getOutputsSize() < 1:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        else:
            summary.append("A total of %d Tomograms were resized to dimensions (%d,%d,%d)"
                           % (self.resizedTomograms.getSize(), self.xDim.get(),
                              self.yDim.get(), self.zDim.get()))
        return summary