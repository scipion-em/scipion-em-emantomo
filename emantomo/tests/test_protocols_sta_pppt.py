# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
from emantomo.protocols import EmanProtTsAlignTomoRec
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.objects import TomoAcquisition
from tomo.protocols import ProtImportTs, ProtImportTsBase
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestEmanTsAlignAndTomoRec(TestBaseCentralizedLayer):

    ds = None
    importedTs = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.importedTs = cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=DataSetRe4STATuto.exclusionWords.value,
                                       anglesFrom=ProtImportTsBase.ANGLES_FROM_TLT,
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=DataSetRe4STATuto.unbinnedPixSize.value,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        cls.assertIsNotNone(tsImported, "There was a problem importing the tilt series")
        return tsImported

    def test_ts_align_tomo_rec(self):
        print(magentaStr("\n==> Aligning the tilt series and reconstructing the tomograms:"))
        protTsAlignTomoRec = self.newProtocol(EmanProtTsAlignTomoRec,
                                              inputTS=self.importedTs,
                                              boxSizeTrk=100,
                                              tltStep=3,
                                              clipz=250,
                                              numberOfThreads=8)
        self.launchProtocol(protTsAlignTomoRec)
        tsAligned = getattr(protTsAlignTomoRec, protTsAlignTomoRec._possibleOutputs.tiltSeries.name, None)
        tomosRec = getattr(protTsAlignTomoRec, protTsAlignTomoRec._possibleOutputs.tomograms.name, None)
        self.assertIsNotNone(tsAligned, "There was a problem aligning the tilt series")
        self.assertIsNotNone(tomosRec, "There was a problem reconstructing the tomograms")
        # Check the tilt series
        testAcq = TomoAcquisition(voltage=DataSetRe4STATuto.voltage.value,
                                  sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                  amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                  magnification=DataSetRe4STATuto.magnification.value,
                                  tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value,
                                  doseInitial=DataSetRe4STATuto.initialDose.value,
                                  dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                  accumDose=DataSetRe4STATuto.accumDose.value
                                  )
        self.checkTiltSeries(tsAligned)
        # Check the tomograms
