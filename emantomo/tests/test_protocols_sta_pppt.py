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
from emantomo.constants import TOMOBOX
from emantomo.protocols import EmanProtTsAlignTomoRec, EmanProtEstimateCTF, EmanProtTSExtraction, \
    EmanProtTomoInitialModelNew, EmanProtTomoRefinementNew, EmanProtMultiRefinementNew
from emantomo.utils import getPresentPrecedents
from imod.protocols import ProtImodImportTransformationMatrix
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME
from pwem import ALIGN_2D
from pwem.protocols import ProtImportVolumes
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols import ProtImportCoordinates3DFromStar
from tomo.constants import TR_EMAN
from tomo.objects import TomoAcquisition, Coordinate3D, SetOfCoordinates3D
from tomo.protocols import ProtImportTs, ProtImportTomograms, ProtImportTsCTF
from tomo.protocols.protocol_import_ctf import ImportChoice
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestEmanBasePPPT(TestBaseCentralizedLayer):
    ds = None
    importedTs = None
    particlesUnbinnedBoxSize = 256
    particlesExtractedBoxSize = 64

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


class TestBaseRefineCyclePPPT(TestEmanBasePPPT):
    expectedSetSize = DataSetRe4STATuto.nCoordsFromTs03.value + DataSetRe4STATuto.nCoordsFromTs54.value
    tsWithAlignment = None
    importedCtfs = None
    importedTomos = None
    importedCoords = None
    extractedParticles = None

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
        protImportCtfs = cls.newProtocol(ProtImportTsCTF,
                                         importFrom=ImportChoice.CTFFIND.value,
                                         filesPath=cls.ds.getFile(DataSetRe4STATuto.cistemFilesPath.value),
                                         filesPattern='*.txt',
                                         inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportCtfs)
        outCtfSet = getattr(protImportCtfs, protImportCtfs._possibleOutputs.CTFs.name, None)
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

    def checkEmanParticles(self, inSet, outParticles, hasEman2dAliFile=False, hasEman3dAliFile=False):
        """Check the specific attributes of EmanSetOfParticles and EmanParticle objects.

        :param inSet: SetOf3DCoordinates or SetOfSubTomograms introduced.
        :param outParticles: the resulting SetOfSubTomograms.
        :param hasEman2dAliFile: used to indicate if the 2d particles have been aligned by EMAN or not.
        :param hasEman3dAliFile: the same as hasEmanAli2dFile, but for 3d particles.
        """

        def _checkEmanFile(inFile, expectedBaseName):
            self.assertTrue(exists(inFile), msg=f'{inFile} does not exist')
            self.assertEqual(basename(inFile), expectedBaseName)

        inCoords = inSet if type(inSet) is SetOfCoordinates3D else inSet.getCoordinates3D()
        ali2dFile = outParticles.getAli2dLstFile()
        if hasEman2dAliFile:
            self.assertTrue(exists(ali2dFile), msg=f'{ali2dFile} does not exist')
        else:
            self.assertIsNone(ali2dFile)
        ali3dFile = outParticles.getAli3dLstFile()
        if hasEman3dAliFile:
            self.assertTrue(exists(ali3dFile), msg=f'{ali3dFile} does not exist')
        else:
            self.assertIsNone(ali3dFile)
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


class TestEmanExtractParticlesFromTS(TestBaseRefineCyclePPPT):

    def test_extract_particle_from_ts(self):
        extractedParticles = self._runExtractParticlesFromTs()
        # Check the results
        expectedSetSize = DataSetRe4STATuto.nCoordsFromTs03.value + DataSetRe4STATuto.nCoordsFromTs54.value
        expectedSRate = DataSetRe4STATuto.sRateBin4.value
        expectedBozxSize = self.particlesExtractedBoxSize
        # Check the subtomograms attributes as EmanSetOfParticles inherits from SetOfSubTomograms
        self.checkExtractedSubtomos(self.importedCoords, extractedParticles,
                                    expectedSetSize=expectedSetSize,
                                    expectedSRate=expectedSRate,
                                    expectedBoxSize=expectedBozxSize,
                                    convention=TR_EMAN,
                                    orientedParticles=True)  # The tutorial dataset provides them oriented
        self.checkEmanParticles(self.importedCoords, extractedParticles,
                                hasEman2dAliFile=False,
                                hasEman3dAliFile=False,  # False as EMAN has not aligned them. The orientations
                                # provided will be used as initial condition for the refinement
                                )


class TestEmanInitialVolumeNew(TestBaseRefineCyclePPPT):

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.extractedParticles = cls._runExtractParticlesFromTs()

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
        outAvgs = getattr(protNewInitModel, EmanProtTomoInitialModelNew._possibleOutputs.average.name, None)
        return outAvgs

    def test_new_initial_volume(self):
        initModel = self._runNewInitialModel(self.extractedParticles)
        # Check the results
        expectedSRate = DataSetRe4STATuto.sRateBin4.value
        expectedBozxSize = self.particlesExtractedBoxSize
        self.checkAverage(initModel,
                          expectedSRate=expectedSRate,
                          expectedBoxSize=expectedBozxSize,
                          hasHalves=False)  # EMAN's e2spt_sgd_new.py does not generate the halves


class TestEmanRefinementAndClassif(TestBaseRefineCyclePPPT):
    importedRef = None

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.extractedParticles = cls._runExtractParticlesFromTs()
        cls.importedRef = cls._runImportReference()

    @classmethod
    def _runImportReference(cls):
        print(magentaStr("\n==> Importing the reference volume:"))
        protImportRef = cls.newProtocol(ProtImportVolumes,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.initVolByEman.name),
                                        samplingRate=DataSetRe4STATuto.sRateBin4.value)
        cls.launchProtocol(protImportRef)
        return getattr(protImportRef, ProtImportVolumes._possibleOutputs.outputVolume.name, None)

    @classmethod
    def _runSubtomoRefineNew(cls, inEmanParticles, refVol):
        print(magentaStr("\n==> Refining the particles:"))
        protNewRefine = cls.newProtocol(EmanProtTomoRefinementNew,
                                        inputSubtomos=inEmanParticles,
                                        refVol=refVol,
                                        iters='p2,t',  # Only 3 to avoid the test taking too long
                                        symmetry='c6',
                                        numberOfThreads=8)
        cls.launchProtocol(protNewRefine)
        outAvg = getattr(protNewRefine, EmanProtTomoRefinementNew._possibleOutputs.average.name, None)
        outParticles = getattr(protNewRefine, EmanProtTomoRefinementNew._possibleOutputs.subtomograms.name, None)
        outFscs = getattr(protNewRefine, EmanProtTomoRefinementNew._possibleOutputs.FSCs.name, None)
        return outAvg, outParticles, outFscs

    @classmethod
    def _runMultiRefClassif(cls, inEmanParticles, doAlignment=None, labelMsg='', nClasses=-1, ref=None):
        print(magentaStr(f"\n==> Running the multi ref classification {labelMsg}:"))
        inputDict = {
            'inputSubtomos': inEmanParticles,
            'maxRes': 30,
            'symmetry': 'c6',
            'numberOfThreads': 8
        }
        if nClasses > 0:
            inputDict['nClasses'] = nClasses
        if ref:
            inputDict['refVol'] = ref
        if doAlignment:
            inputDict['doAlignment'] = True
        protMultiCl = cls.newProtocol(EmanProtMultiRefinementNew, **inputDict)
        cls.launchProtocol(protMultiCl)
        protMultiCl.setObjLabel(f'MultiRefCl {labelMsg}')
        outParticles = getattr(protMultiCl, EmanProtMultiRefinementNew._possibleOutputs.subtomograms.name, None)
        outClasses = getattr(protMultiCl, EmanProtMultiRefinementNew._possibleOutputs.classes.name, None)
        return outParticles, outClasses

    def checkAliEmanParticles(self, inParticles, outParticles):
        expectedSRate = DataSetRe4STATuto.sRateBin4.value
        expectedBozxSize = self.particlesExtractedBoxSize
        self.checkRefinedSubtomograms(inParticles, outParticles,
                                      expectedSetSize=self.expectedSetSize,
                                      expectedBoxSize=expectedBozxSize,
                                      expectedSRate=expectedSRate,
                                      convention=TR_EMAN,
                                      orientedParticles=True
                                      )
        self.checkEmanParticles(inParticles, outParticles,
                                hasEman2dAliFile=True,
                                hasEman3dAliFile=True)

    def test_refine_and_classify(self):
        expectedSRate = DataSetRe4STATuto.sRateBin4.value
        expectedBozxSize = self.particlesExtractedBoxSize
        # [1] REFINEMENT ###############################################################################################
        refinedAvg, refinedParticles, outFscs = self._runSubtomoRefineNew(self.extractedParticles, self.importedRef)
        # Check the average
        self.checkAverage(refinedAvg,
                          expectedSRate=expectedSRate,
                          expectedBoxSize=expectedBozxSize,
                          hasHalves=True)
        # Check the particles
        self.checkAliEmanParticles(self.extractedParticles, refinedParticles)
        # Check the FSCs
        self.assertSetSize(outFscs, size=3)

        # [2A] CLASSIFY WITHOUT ALIGNMENT, N CLASSES 1, NO REF #########################################################
        nClasses = 1
        outParticles, outClasses = self._runMultiRefClassif(refinedParticles,
                                                            doAlignment=False,
                                                            nClasses=nClasses,
                                                            labelMsg='no ali, 1 cl, no ref')
        # Check the particles
        self.checkAliEmanParticles(refinedParticles, outParticles)
        # Check the classes
        self.checkClasses(outClasses, expectedSetSize=nClasses, expectedSRate=expectedSRate)

        # [2B] CLASSIFY WITHOUT ALIGNMENT, N CLASSES 3, NO REF #########################################################
        nClasses = 3
        outParticles, outClasses = self._runMultiRefClassif(refinedParticles,
                                                            doAlignment=False,
                                                            nClasses=nClasses,
                                                            labelMsg='no ali, 3 cl, no ref')
        # Check the particles
        self.checkAliEmanParticles(refinedParticles, outParticles)
        # Check the classes
        self.checkClasses(outClasses, expectedSetSize=nClasses, expectedSRate=expectedSRate)

        # [2C] CLASSIFY WITHOUT ALIGNMENT WITH 1 REF ###################################################################
        nClasses = 1
        outParticles, outClasses = self._runMultiRefClassif(refinedParticles,
                                                            doAlignment=False,
                                                            ref=refinedAvg,
                                                            labelMsg='no ali, with 1 ref')
        # Check the particles
        self.checkAliEmanParticles(refinedParticles, outParticles)
        # Check the classes
        self.checkClasses(outClasses, expectedSetSize=nClasses, expectedSRate=expectedSRate)

        # [2D] CLASSIFY WITH ALIGNMENT AND REF #########################################################################
        nClasses = 1
        outParticles, outClasses = self._runMultiRefClassif(refinedParticles,
                                                            doAlignment=True,
                                                            ref=refinedAvg,
                                                            labelMsg='Ali with ref')
        # Check the particles
        self.checkAliEmanParticles(refinedParticles, outParticles)
        # Check the classes
        self.checkClasses(outClasses, expectedSetSize=nClasses, expectedSRate=expectedSRate)

        # [2E] CLASSIFY WITH ALIGNMENT, N CLASSES 1 ####################################################################
        nClasses = 1
        outParticles, outClasses = self._runMultiRefClassif(refinedParticles,
                                                            doAlignment=True,
                                                            nClasses=nClasses,
                                                            labelMsg='Ali, 1 cl, no ref')
        # Check the particles
        self.checkAliEmanParticles(refinedParticles, outParticles)
        # Check the classes
        self.checkClasses(outClasses, expectedSetSize=nClasses, expectedSRate=expectedSRate)
