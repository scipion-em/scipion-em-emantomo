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
from os.path import join, exists
from typing import List

from emantomo import Plugin
from emantomo.constants import TS_DIR
from emantomo.convert import ts2Json_
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TS
from emantomo.utils import genJsonFileName
from pyworkflow.protocol import FloatParam, IntParam, BooleanParam, PointerParam
from pyworkflow.utils import Message
from tomo.objects import SetOfCTFTomoSeries


class EstimateCtfOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class EmanProtEstimateCTFBase(ProtEmantomoBase):
    """
    Protocol for CTF estimation from tilt series using e2spt_tomoctf.py
    """
    _label = 'ctf estimation'
    _possibleOutputs = EstimateCtfOutputs
    program = Plugin.getProgram('e2spt_tomoctf.py')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series",
                      important=True)
        self._defineCTFParams(form)
        self._addBinThreads(form)

    @staticmethod
    def _defineCTFParams(form):
        lineDefocus = form.addLine('Defocus range (μm)',
                                   help="Search range of defocus (start, end, step). Note that "
                                        "they must be introduced in microns.")

        lineDefocus.addParam('minDefocus', FloatParam, default=2.0, label='min')
        lineDefocus.addParam('maxDefocus', FloatParam, default=7.0, label='max')
        lineDefocus.addParam('stepDefocus', FloatParam, default=0.02, label='Step')
        form.addParam('doPhaseShiftSearch', BooleanParam,
                      label='Do phase shift search?',
                      default=False)
        linePhaseShift = form.addLine('Phase shift range (deg.)',
                                      condition='doPhaseShiftSearch',
                                      help="Search range of the phase shift (start, end, step). To avoid"
                                           "the phase shift search use min 0.0, max 1.0, and step 1.0.")
        linePhaseShift.addParam('minPhaseShift', FloatParam, default=0, label='min', condition='doPhaseShiftSearch')
        linePhaseShift.addParam('maxPhaseShift', FloatParam, default=1, label='max', condition='doPhaseShiftSearch')
        linePhaseShift.addParam('stepPhaseShift', FloatParam, default=1, label='Step', condition='doPhaseShiftSearch')

        form.addParam('tilesize', IntParam,
                      label='Size of tile to calculate FFT',
                      default=256)
        form.addParam('nref', IntParam,
                      label='Number of references',
                      default=15,
                      help='Using N tilt images near the center tilt to estimate the range of defocus for all images.')
        form.addParam('stepx', IntParam,
                      label='Step in X direction',
                      default=20,
                      help='Number of tiles to generate on x-axis (different defocus)')
        form.addParam('stepy', IntParam,
                      label='Step in Y direction',
                      default=40,
                      help='Number of tiles to generate on y-axis (same defocus)')

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inTsSet = self.getAttrib(IN_TS)
        self.createInitEmanPrjDirs()

    def writeData2JsonFileStep(self, tsId: str):
        ts = self.getCurrentTs(tsId)
        jsonFile = genJsonFileName(self.getInfoDir(), tsId)
        mode = 'a' if exists(jsonFile) else 'w'
        ts2Json_(ts, jsonFile, mode=mode)

    def estimateCtfStep(self, tsId: str):
        args = ' '.join(self._genCtfEstimationArgs(tsId))
        self.runJob(self.program, args, cwd=self._getExtraPath())

    # --------------------------- UTILS functions -----------------------------
    def _genCtfEstimationArgs(self, tsId: str) -> List[str]:
        ts = self.getCurrentTs(tsId)
        acq = ts.getAcquisition()
        args = ['%s ' % join(TS_DIR, f'{tsId}.hdf'),
                f'--dfrange {self.minDefocus.get()},{self.maxDefocus.get()},{self.stepDefocus.get()}',
                f'--psrange {self.minPhaseShift.get()},{self.maxPhaseShift.get()},{self.stepPhaseShift.get()}',
                '--tilesize %i ' % self.tilesize.get(),
                '--voltage %i ' % acq.getVoltage(),
                '--cs %.2f ' % acq.getSphericalAberration(),
                '--nref %i ' % self.nref.get(),
                '--stepx %i ' % self.stepx.get(),
                '--stepy %i ' % self.stepy.get(),
                '--threads %i ' % self.binThreads.get(),
                '--verbose 9 ']
        return args
