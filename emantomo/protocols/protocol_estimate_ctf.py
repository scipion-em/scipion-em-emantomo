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
from os.path import exists, join
from emantomo import Plugin
from emantomo.constants import TS_DIR
from emantomo.objects import EmanMetaData
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TS
from emantomo.utils import getPresentTsIdsInSet, genJsonFileName
from pyworkflow.protocol import PointerParam, FloatParam, IntParam, BooleanParam
from tomo.objects import SetOfCTFTomoSeries, CTFTomoSeries, CTFTomo
from emantomo.convert import loadJson, ts2Json


class EstimateCtfOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class EmanProtEstimateCTF(ProtEmantomoBase):
    """
    Protocol for CTF estimation from tilt series using e2spt_tomoctf.py
    """
    _label = 'ctf estimation'
    _possibleOutputs = EstimateCtfOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series",
                      important=True)
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
        form.addParallelSection(threads=4, mpi=0)

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


