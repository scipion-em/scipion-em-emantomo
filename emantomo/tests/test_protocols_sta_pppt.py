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
from cistem.protocols import CistemProtTsImportCtf

from emantomo.protocols import EmanProtTsAlignTomoRec, EmanProtEstimateCTF, EmanProtTSExtraction
from imod.protocols import ProtImodImportTransformationMatrix
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME
from pwem import ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols import ProtImportCoordinates3DFromStar
from tomo.objects import TomoAcquisition
from tomo.protocols import ProtImportTs, ProtImportTomograms
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestEmanBasePPPT(TestBaseCentralizedLayer):
    ds = None
    importedTs = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)

    @classmethod
    def _runImportTs(cls, exclusionWords=DataSetRe4STATuto.exclusionWordsTs03.value):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=exclusionWords,
                                       anglesFrom=2,  # From tlt file
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


class TestEmanTsAlignAndTomoRec(TestEmanBasePPPT):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.importedTs = cls._runImportTs()

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
        self.checkTiltSeries(tsAligned,
                             expectedSetSize=1,
                             expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value,
                             expectedDimensions=[3710, 3838, 40],
                             testAcqObj=testAcq,
                             alignment=ALIGN_2D,
                             hasAlignment=True,
                             anglesCount=40,
                             )
        # Check the tomograms
        self.checkTomograms(tomosRec,
                            expectedSetSize=1,
                            expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value * 4,  # Bin 4
                            expectedDimensions=[1050, 1000, 250]
                            )


class TestEmanEstimateCtf(TestEmanBasePPPT):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.importedTs = cls._runImportTs()

    def test_estimate_ctf(self):
        print(magentaStr("\n==> Estimating the CTF:"))
        protEstimateCtf = self.newProtocol(EmanProtEstimateCTF,
                                           inputTS=self.importedTs,
                                           minDefocus=2,
                                           maxDefocus=4,
                                           stepDefocus=0.02,
                                           tilesize=512,
                                           numberOfThreads=8)
        self.launchProtocol(protEstimateCtf)
        outCtfs = getattr(protEstimateCtf, protEstimateCtf._possibleOutputs.CTFs.name, None)
        self.assertIsNotNone(outCtfs, "There was a problem estimating the CTFs")
        # Check the CTFs
        self.checkCTFs(outCtfs, expectedSetSize=1)


class TestEmanExtractParticlesFromTs(TestEmanBasePPPT):
    tsWithAlignment = None
    importedCtfs = None
    importedTomos = None
    importedCoords = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs(exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value)
        cls.tsWithAlignment = cls._runImportTrMatrix()
        cls.importedCtfs = cls._runImportCtfs()
        cls.importedTomos = cls._runImportTomograms()
        cls.importedCoords = cls._runImportCoordinatesFromStar()

    @classmethod
    def _runImportTrMatrix(cls):
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=DataSetRe4STATuto.transformPattern.value,
                                             inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        cls.assertIsNotNone(outTsSet, "There was a problem importing the transformation matrices")
        return outTsSet

    @classmethod
    def _runImportCtfs(cls):
        print(magentaStr("\n==> Importing the TS' CTFs with Cistem:"))
        protImportCtfs = cls.newProtocol(CistemProtTsImportCtf,
                                         filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                         filesPattern=DataSetRe4STATuto.ctfPattern.value,
                                         inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportCtfs)
        outCtfSet = getattr(protImportCtfs, CistemProtTsImportCtf._possibleOutputs.outputCTFs.name, None)
        cls.assertIsNotNone(outCtfSet, "There was a problem importing the estimated CTFs")
        return outCtfSet

    @classmethod
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=DataSetRe4STATuto.tomosPattern.value,
                                          samplingRate=DataSetRe4STATuto.tomosSRate.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        cls.assertIsNotNone(outTomos, "There was a problem importing the tomograms")
        return outTomos

    @classmethod
    def _runImportCoordinatesFromStar(cls):
        print(magentaStr("\n==> Importing the coordinates with Relion:"))
        protImportCoords = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                           starFile=cls.ds.getFile(DataSetRe4STATuto.coordsStar.value),
                                           inTomos=cls.importedTomos,
                                           samplingRate=DataSetRe4STATuto.unbinnedPixSize.value)
        cls.launchProtocol(protImportCoords)
        outCoords = getattr(protImportCoords, ProtImportCoordinates3DFromStar._possibleOutputs.coordinates.name, None)
        cls.assertIsNotNone(outCoords, "There was a problem importing the coordinates with Relion")
        return outCoords

    @classmethod
    def _runExtractParticlesFromTs(cls):
        print(magentaStr("\n==> Extracting the particles from the TS:"))
        protExtractFromTs = cls.newProtocol(EmanProtTSExtraction,
                                            inputSubtomos=cls.importedCoords,
                                            inputCTF=cls.importedCtfs,
                                            inputTS=cls.tsWithAlignment,
                                            boxSize=120,
                                            shrink=4,
                                            numberOfThreads=8)
        cls.launchProtocol(protExtractFromTs)
        # outCoords = getattr(protExtractFromTs, ProtImportCoordinates3DFromStar._possibleOutputs.coordinates.name, None)
        # cls.assertIsNotNone(outCoords, "There was a problem importing the coordinates with Relion")
        # return outCoords

    def test_extract_particles_from_ts(self):
        self._runExtractParticlesFromTs()
