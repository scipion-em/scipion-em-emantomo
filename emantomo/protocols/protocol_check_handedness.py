# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import logging
import typing
from os.path import join
import numpy as np
from emantomo.constants import TS_DIR
from emantomo.objects import EmanMetaData
from emantomo.protocols.protocol_base import IN_TS
from emantomo.protocols.protocol_estimate_ctf_base import EmanProtEstimateCTFBase
from emantomo.utils import genJsonFileName
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message, cyanStr

logger = logging.getLogger(__name__)


class EmanProtEstimateHandedness(EmanProtEstimateCTFBase):
    """
    Protocol for checking handedness by CTF using e2spt_tomoctf.py
    """
    _label = 'check handedness'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS, PointerParam,
                      pointerClass='TiltSeries',
                      label="Tilt Series",
                      important=True)
        self._defineCTFParams(form)
        self._addBinThreads(form)
        form.addParallelSection(threads=1, mpi=0)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        tsId = self.getAttrib(IN_TS).getTsId()
        self._insertFunctionStep(self.convertTsStep, tsId, needsGPU=False)
        self._insertFunctionStep(self.writeData2JsonFileStep, tsId, needsGPU=False)
        self._insertFunctionStep(self.estimateCtfStep, tsId, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        ts = self.getAttrib(IN_TS)
        tsId = ts.getTsId()
        self.createInitEmanPrjDirs()
        # Get the required acquisition data
        acq = ts.getAcquisition()
        self.sphAb = acq.getSphericalAberration()
        self.voltage = acq.getVoltage()
        # Do this to take advantage of the super-class methods, prepared to work with a dict of EmanMetadata objects
        self.mdObjDict[tsId] = EmanMetaData(tsId=tsId,
                                            ts=ts,
                                            tsHdfName=join(TS_DIR, f'{tsId}.hdf'),
                                            jsonFile=genJsonFileName(self.getInfoDir(), tsId))

    def estimateCtfStep(self, tsId: str):
        logger.info(cyanStr(f'===> tsId = {tsId}: estimating the handedness...'))
        super().estimateCtfStep(tsId)

    def createOutputStep(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genCtfEstimationArgs(self, mdObj: EmanMetaData):
        args = super()._genCtfEstimationArgs(mdObj)
        args += '--checkhand --writetmp'
        return args

    # --------------------------- INFO functions ------------------------------
    def _validate(self) -> typing.List[str]:
        msg = []
        firstTi = self.getAttrib(IN_TS).getFirstItem()
        transform = firstTi.getTransform()
        errorMsg = 'The tilt-series introduced must contain alignment data.'
        if transform:
            trMatrix = transform.getMatrix()
            if trMatrix is not None:
                idMatrix = np.eye(3)
                if np.array_equal(trMatrix, idMatrix):
                    msg.append(errorMsg)
            else:
                msg.append(errorMsg)
        else:
            msg.append(errorMsg)
        return msg
