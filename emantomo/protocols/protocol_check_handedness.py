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
import re
import subprocess
import typing
from enum import Enum

import numpy as np
from emantomo.protocols.protocol_base import IN_TS
from emantomo.protocols.protocol_estimate_ctf_base import EmanProtEstimateCTFBase
from pyworkflow import BETA
from pyworkflow.object import Boolean, String
from pyworkflow.protocol import StringParam
from pyworkflow.utils import cyanStr, greenStr
from tomo.objects import TiltSeries

logger = logging.getLogger(__name__)

class EstimatedHandednessOutputs(Enum):
    EmanHandedness = Boolean


class EmanProtEstimateHandedness(EmanProtEstimateCTFBase):
    """
    Protocol for checking handedness by CTF using e2spt_tomoctf.py
    """
    _label = 'check handedness'
    _possibleOutputs = EstimatedHandednessOutputs
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.parsedMsg = String()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        super()._defineParams(form)
        form.addParam('chosenTsId', StringParam,
                      label='Tilt-series id',
                      help='Tilt-series identifier of the tilt-series that will be used to '
                           'estimate the handedness. The wizard provided by this parameter displays a '
                           'list with the tilt-series ids of the tilt-series contained in the set of '
                           'tilt-series introduced.')

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        tsId = self.chosenTsId.get()
        self._insertFunctionStep(self.convertTsStep, tsId, needsGPU=False)
        self._insertFunctionStep(self.writeData2JsonFileStep, tsId, needsGPU=False)
        self._insertFunctionStep(self.estimateCtfStep, tsId, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def estimateCtfStep(self, tsId: str):
        logger.info(cyanStr(f'===> tsId = {tsId}: estimating the handedness...'))
        args = ' '.join(self._genCtfEstimationArgs(tsId))
        cmd = f'{self.program} {args}'
        logger.info(greenStr(cmd))
        res = subprocess.run(cmd, cwd=self._getExtraPath(), capture_output=True, text=True, shell=True)
        stderrMsg = res.stderr
        if stderrMsg:
            raise Exception(stderrMsg)
        parsedMsg, isFlipped = self._parseStdOut(res.stdout)
        self.parsedMsg.set(parsedMsg)
        isFlipped = Boolean(isFlipped)
        self._store(self.parsedMsg)
        self._defineOutputs(**{self._possibleOutputs.EmanHandedness.name: isFlipped})
        self._defineSourceRelation(getattr(self, IN_TS), isFlipped)

    # --------------------------- UTILS functions -----------------------------
    def _genCtfEstimationArgs(self, tsId: str) -> typing.List[str]:
        args = super()._genCtfEstimationArgs(tsId)
        args.append('--checkhand')
        args.append('--writetmp')
        return args

    @staticmethod
    def _parseStdOut(stdoutMsg: str) -> typing.Tuple[str, bool]:
        """Parses the stdout of the EMAN handedness estimation program and returns a substring with
        the relevant information.

        Example of favorable case:

        [...]
        Average score: Current hand - 0.452, flipped hand - 0.379
        Defocus std: Current hand - 0.753, flipped hand - 0.744
        Current hand is better than the flipped hand in 82.5% tilt images
        tiltseries/TS_03.hdf: Result: the handedness (--tltax=85.3) seems to be correct.
        Rerun CTF estimation without the checkhand option to finish the process.

        Example of non-favorable case:
        Average score: Current hand - 0.382, flipped hand - 0.445
        Defocus std: Current hand - 0.888, flipped hand - 0.530
        Current hand is better than the flipped hand in 15.0% tilt images
        tiltseries/TS_03.hdf: Result: the handedness seems to be flipped.
        Consider rerun the tomogram reconstruction with --tltax=-274.7
        then rerun the CTF estimation.
        [...]
        """
        pattern = r"(current hand is.*?rerun)"
        matchList = re.findall(pattern, stdoutMsg.lower(), re.DOTALL)
        parsedMsg = matchList[0].replace(' consider rerun', '').replace(' rerun', '')
        isFlipped = False if 'correct' in parsedMsg else True
        return parsedMsg, isFlipped

    # --------------------------- INFO functions ------------------------------
    def _summary(self) -> typing.List[str]:
        msg = []
        if self.isFinished():
            msg.append(f'*IsFlipped: {getattr(self, self._possibleOutputs.EmanHandedness.name).get()}*\n'
                       f'{self.parsedMsg.get()}\n\n'
                       f'NOTE: The result relies on the fitting of CTF rings in the micrographs, and can be '
                       f'unreliable if there is not enough signal. Only believe the result when more than the '
                       f'80% of the tilt angles have the same handedness.')
        return msg

    def _validate(self) -> typing.List[str]:
        msg = []
        ts = self.getAttrib(IN_TS).getItem(TiltSeries.TS_ID_FIELD, self.chosenTsId.get())
        transform = ts.getFirstItem().getTransform()
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
