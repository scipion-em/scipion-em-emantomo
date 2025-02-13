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
from emantomo.protocols.protocol_base import IN_TS
from emantomo.protocols.protocol_estimate_ctf_base import EmanProtEstimateCTFBase
from pyworkflow.object import Set
from pyworkflow.protocol import PointerParam, STEPS_PARALLEL
from pyworkflow.utils import Message
from tomo.objects import SetOfCTFTomoSeries, CTFTomoSeries, CTFTomo
from emantomo.convert import loadJson


class EstimateCtfOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class EmanProtEstimateCTF(EmanProtEstimateCTFBase):
    """
    Protocol for CTF estimation from tilt series using e2spt_tomoctf.py
    """
    _label = 'ctf estimation'
    _possibleOutputs = EstimateCtfOutputs
    stepsExecutionMode = STEPS_PARALLEL

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
        closeSetStepDeps = []
        self._initialize()
        for tsId in self.mdObjDict.keys():
            cId = self._insertFunctionStep(self.convertTsStep, tsId,
                                           prerequisites=[],
                                           needsGPU=False)
            wJsonId = self._insertFunctionStep(self.writeData2JsonFileStep, tsId,
                                               prerequisites=cId,
                                               needsGPU=False)
            ctfId = self._insertFunctionStep(self.estimateCtfStep, tsId,
                                             prerequisites=wJsonId,
                                             needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=ctfId,
                                              needsGPU=False)
            closeSetStepDeps.append(cOutId)
        self._insertFunctionStep(self._closeOutputSet,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, tsId: str):
        with self._lock:
            mdObj = self.mdObjDict[tsId]
            outCtfSet = self.getOutputCtfTomoSet()
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

            newCTFTomoSeries.write()
            outCtfSet.update(newCTFTomoSeries)
            outCtfSet.write()
            self._store(outCtfSet)

    # --------------------------- UTILS functions -----------------------------
    def getOutputCtfTomoSet(self) -> SetOfCTFTomoSeries:
        outCtfSet = getattr(self, self._possibleOutputs.CTFs.name, None)
        if outCtfSet:
            outCtfSet.enableAppend()
        else:
            inTsSetPointer = self.getAttrib(IN_TS, getPointer=True)
            outCtfSet = SetOfCTFTomoSeries.create(self._getPath(),
                                                  template='ctfTomoSeries%s.sqlite')
            outCtfSet.setSetOfTiltSeries(inTsSetPointer)
            outCtfSet.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{self._possibleOutputs.CTFs.name: outCtfSet})
            self._defineSourceRelation(inTsSetPointer, outCtfSet)

        return outCtfSet
