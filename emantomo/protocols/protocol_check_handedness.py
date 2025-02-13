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

from os.path import exists, join
from emantomo import Plugin
from emantomo.constants import TS_DIR
from emantomo.objects import EmanMetaData
from emantomo.protocols.protocol_base import IN_TS
from emantomo.protocols.protocol_estimate_ctf_base import EmanProtEstimateCTFBase
from emantomo.utils import getPresentTsIdsInSet, genJsonFileName
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message
from tomo.objects import SetOfCTFTomoSeries, CTFTomoSeries, CTFTomo
from emantomo.convert import loadJson, ts2Json


class EmanProtEstimateCTF(EmanProtEstimateCTFBase):
    """
    Protocol for checking handedness by CTF using e2spt_tomoctf.py
    """
    _label = 'Check handedness by CTF'

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
        form.addParallelSection(threads=1, mpi=0)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertTsStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
            self._insertFunctionStep(self.estimateCtfStep, mdObj)
        self._insertFunctionStep(self.createOutputStep, mdObjDict)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTsSet = self.getAttrib(IN_TS)
        self.createInitEmanPrjDirs()
        # Manage the TS
        presentTsIds = set(getPresentTsIdsInSet(inTsSet))
        tsIdsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inTsSet if ts.getTsId() in presentTsIds}
        # Get the required acquisition data
        self.sphAb = inTsSet.getAcquisition().getSphericalAberration()
        self.voltage = inTsSet.getAcquisition().getVoltage()
        # Split all the data into subunits (EmanMetaData objects) referred to the same tsId
        mdObjDict = {}
        for tomoId, ts in tsIdsDict.items():
            mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
                                             ts=ts,
                                             tsHdfName=join(TS_DIR, f'{tomoId}.hdf'),
                                             jsonFile=genJsonFileName(self.getInfoDir(), tomoId))
        return mdObjDict

    @staticmethod
    def writeData2JsonFileStep(mdObj):
        mode = 'a' if exists(mdObj.jsonFile) else 'w'
        ts2Json(mdObj, mode=mode)

    def estimateCtfStep(self, mdObj):
        program = Plugin.getProgram('e2spt_tomoctf.py')
        self.runJob(program, self._genCtfEstimationArgs(mdObj), cwd=self._getExtraPath())

    def createOutputStep(self, mdObjDict):
        inTsSet = self.getAttrib(IN_TS)
        outCtfSet = SetOfCTFTomoSeries.create(self._getPath(), template='CTFmodels%s.sqlite')
        outCtfSet.setSetOfTiltSeries(inTsSet)

        for tsId, mdObj in mdObjDict.items():
            ts = mdObj.ts
            newCTFTomoSeries = CTFTomoSeries()
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            newCTFTomoSeries.setObjId(ts.getObjId())
            newCTFTomoSeries.setTsId(tsId)
            outCtfSet.append(newCTFTomoSeries)
            jsonDict = loadJson(mdObj.jsonFile)
            defocus = jsonDict['defocus']
            phase_shift = jsonDict['phase']
            for idx, tiltImage in enumerate(ts.iterItems()):
                defocusU = defocusV = 10000.0 * defocus[idx]
                newCTFTomo = CTFTomo()
                newCTFTomo.setIndex(idx + 1)
                newCTFTomo.setAcquisitionOrder(tiltImage.getAcquisitionOrder())
                if phase_shift[idx] != 0:
                    newCTFTomo.setPhaseShift(phase_shift[idx])
                newCTFTomo.setDefocusU(defocusU)
                newCTFTomo.setDefocusV(defocusV)
                newCTFTomo.setDefocusAngle(0)
                newCTFTomoSeries.append(newCTFTomo)

        self._defineOutputs(**{self._possibleOutputs.CTFs.name: outCtfSet})
        self._defineSourceRelation(getattr(self, IN_TS), outCtfSet)

    # --------------------------- UTILS functions -----------------------------
    def _genCtfEstimationArgs(self, mdObj):
        args = '%s ' % mdObj.tsHdfName
        args += '--dfrange %s ' % ','.join([str(self.minDefocus.get()),
                                            str(self.maxDefocus.get()),
                                            str(self.stepDefocus.get())])
        args += '--psrange %s ' % ','.join([str(self.minPhaseShift.get()),
                                            str(self.maxPhaseShift.get()),
                                            str(self.stepPhaseShift.get())])
        args += '--tilesize %i ' % self.tilesize.get()
        args += '--voltage %i ' % self.voltage
        args += '--cs %.2f ' % self.sphAb
        args += '--nref %i ' % self.nref.get()
        args += '--stepx %i ' % self.stepx.get()
        args += '--stepy %i ' % self.stepy.get()
        args += '--threads %i ' % self.numberOfThreads.get()
        args += '--verbose 9 '
        return args


