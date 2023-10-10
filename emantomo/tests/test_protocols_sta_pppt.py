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
from os.path import exists, basename

from cistem.protocols import CistemProtTsImportCtf

from emantomo.constants import TOMOBOX
from emantomo.protocols import EmanProtTsAlignTomoRec, EmanProtEstimateCTF, EmanProtTSExtraction, \
    EmanProtTomoInitialModelNew
from emantomo.utils import getPresentPrecedents
from imod.protocols import ProtImodImportTransformationMatrix
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME
from pwem import ALIGN_2D
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols import ProtImportCoordinates3DFromStar
from tomo.constants import TR_EMAN
from tomo.objects import TomoAcquisition, Coordinate3D
from tomo.protocols import ProtImportTs, ProtImportTomograms
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestEmanBasePPPT(TestBaseCentralizedLayer):
    ds = None
    importedTs = None
    particlesUnbinnedBoxSize = 120
    particlesExtractedBoxSize = 32

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


class TestEmanPPPTRefineCycle(TestEmanBasePPPT):
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
        return outCtfSet

    @classmethod
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=DataSetRe4STATuto.tomosPattern.value,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    @classmethod
    def _runImportCoordinatesFromStar(cls):
        print(magentaStr("\n==> Importing the coordinates with Relion:"))
        protImportCoords = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                           starFile=cls.ds.getFile(DataSetRe4STATuto.coordsStarSubset.value),
                                           inTomos=cls.importedTomos,
                                           samplingRate=DataSetRe4STATuto.unbinnedPixSize.value)
        cls.launchProtocol(protImportCoords)
        outCoords = getattr(protImportCoords, ProtImportCoordinates3DFromStar._possibleOutputs.coordinates.name, None)
        return outCoords

    @classmethod
    def _runExtractParticlesFromTs(cls):
        print(magentaStr("\n==> Extracting the particles from the TS:"))
        protExtractFromTs = cls.newProtocol(EmanProtTSExtraction,
                                            inputSubtomos=cls.importedCoords,
                                            inputCTF=cls.importedCtfs,
                                            inputTS=cls.tsWithAlignment,
                                            boxSize=cls.particlesUnbinnedBoxSize,
                                            shrink=4,
                                            numberOfThreads=8)
        cls.launchProtocol(protExtractFromTs)
        outParticles = getattr(protExtractFromTs, protExtractFromTs._possibleOutputs.subtomograms.name, None)
        return outParticles

    @classmethod
    def _runNewInitialModel(cls, inEmanParticles):
        print(magentaStr("\n==> Generating the initial model:"))
        protNewInitModel = cls.newProtocol(EmanProtTomoInitialModelNew,
                                           inputSubtomos=inEmanParticles,
                                           nIters=10,
                                           nClasses=1,
                                           symmetry='c6',
                                           numberOfThreads=8)
        cls.launchProtocol(protNewInitModel)
        outAvgs = getattr(protNewInitModel, EmanProtTomoInitialModelNew._possibleOutputs.averages.name, None)
        return outAvgs[1]  # There's only one element in the set of averages as the no. classes requested is 1

    def checkEmanParticles(self, inCoords, outParticles, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                           convention=TR_EMAN, orientedParticles=False, hasEman2dAliFile=False, hasEman3dAliFile=False):
        """Checks exhaustively the Eman particles generated after having carried out a particle extraction from TS
        with EMAN.

        :param inCoords: SetOf3DCoordinates introduced for the subtomo extraction.
        :param outParticles: the resulting SetOfSubTomograms.
        :param expectedSetSize: expected set site to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py
        :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be an
        eye matrix (False) or not (True).
        :param hasEman2dAliFile: used to indicate if the 2d particles have been aligned by EMAN or not.
        :param hasEman3dAliFile: the same as hasEmanAli2dFile, but for 3d particles.
        """

        def _checkEmanFile(inFile, expectedBaseName):
            self.assertTrue(exists(inFile), msg=f'{inFile} does not exist')
            self.assertEqual(basename(inFile), expectedBaseName)

        # Check the subtomograms attributes as EmanSetOfParticles inherits from SetOfSubTomograms
        self.checkExtractedSubtomos(inCoords, outParticles,
                                    expectedSetSize=expectedSetSize,
                                    expectedSRate=expectedSRate,
                                    expectedBoxSize=expectedBoxSize,
                                    convention=convention,
                                    orientedParticles=orientedParticles)
        # Check the specific attributes of EmanSetOfParticles and EmanParticle objects
        ali2dFile = outParticles.getAli2dLstFile()
        self.assertTrue(exists(ali2dFile)) if hasEman2dAliFile else self.assertIsNone(ali2dFile)
        ali3dFile = outParticles.getAli3dLstFile()
        self.assertTrue(exists(ali3dFile)) if hasEman3dAliFile else self.assertIsNone(ali3dFile)
        presentTsIds = inCoords.getUniqueValues(Coordinate3D.TOMO_ID_ATTR)
        for tomo in getPresentPrecedents(inCoords, presentTsIds):
            tomoId = tomo.getTsId()
            expectedInfoJson = f'{tomoId}_info.json'
            expectedTsBaseName = f'{tomoId}.hdf'  # Converted into HDF
            expectedHdfStacksBaseName = f'{tomoId}__{TOMOBOX}_bin4.hdf'
            for particle in outParticles.iterSubtomos(volume=tomo):
                infoJson = particle.getInfoJson()
                _checkEmanFile(infoJson, expectedInfoJson)
                tsHdfFile = particle.getTsHdf()
                _checkEmanFile(tsHdfFile, expectedTsBaseName)
                stack2dFile = particle.getStack2dHdf()
                _checkEmanFile(stack2dFile, expectedHdfStacksBaseName)
                stack3dFile = particle.getStack3dHdf()
                _checkEmanFile(stack3dFile, expectedHdfStacksBaseName)

    def test_pppt_refine_cycle(self):
        # [1] EXTRACT PARTICLES FROM THE TS
        outParticles = self._runExtractParticlesFromTs()
        # Check the results
        expectedSetSize = DataSetRe4STATuto.nCoordsFromTs03.value + DataSetRe4STATuto.nCoordsFromTs54.value
        self.checkEmanParticles(self.importedCoords, outParticles,
                                expectedSetSize=expectedSetSize,
                                expectedSRate=DataSetRe4STATuto.sRateBin4.value,
                                expectedBoxSize=self.particlesExtractedBoxSize,
                                convention=TR_EMAN,
                                orientedParticles=True,  # The tutorial dataset provides them oriented
                                hasEman2dAliFile=False,
                                hasEman3dAliFile=False,  # False as EMAN has not aligned them. The orientations
                                # provided will be used as initial condition for the refinement
                                )
        # [2] NEW INITIAL VOLUME
        initModel = self._runNewInitialModel(outParticles)
        self.checkAverage(initModel,
                          expectedSRate=DataSetRe4STATuto.sRateBin4.value,
                          expectedBoxSize=self.particlesExtractedBoxSize,
                          hasHalves=False)  # EMAN's e2spt_sgd_new.py does not generate the halves
