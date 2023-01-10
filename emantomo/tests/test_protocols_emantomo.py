# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import os

from pyworkflow.tests import (BaseTest, setupTestProject, DataSet)
from pwem.protocols import (ProtImportMicrographs, ProtImportParticles, ProtImportVolumes,
                            ProtImportAverages)
import pyworkflow.utils as pwutils
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from ..constants import EMAN_COVERAGE, EMAN_SCORE
from ..protocols import *
import tomo.protocols
from tomo.constants import TR_EMAN

import emantomo
from ..protocols.protocol_tomo_extraction_from_tomo import OutputExtraction
from ..protocols.protocol_tomo_subtomogram_refinement import EmanTomoRefinementOutputs


class TestEmanBase(BaseTest):

    origSRate = 5

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')

    @classmethod
    def setData(cls, projectData='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(projectData)
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification, sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: passe the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        if samplingRate is not None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=0, filesPath=pattern,
                                             samplingRate=samplingRate, magnification=magnification,
                                             voltage=voltage, sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=1, filesPath=pattern,
                                             scannedPixelSize=scannedPixelSize,
                                             voltage=voltage, magnification=magnification,
                                             sphericalAberration=sphericalAberration)

        cls.proj.launchProtocol(cls.protImport, wait=True)
        cls.assertIsNotNone(cls.protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")

        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern, samplingRate=samplingRate,
                                         checkStack=checkStack)
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return cls.protImport

    @classmethod
    def runImportParticlesSqlite(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         importFrom=4,
                                         sqliteFile=pattern, samplingRate=samplingRate)
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes,
                                     filesPath=pattern, samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def runImportAverages(cls, pattern, samplingRate):
        """ Run an Import averages protocol. """
        cls.protImportAvg = cls.newProtocol(ProtImportAverages,
                                            filesPath=pattern,
                                            samplingRate=samplingRate,
                                            checkStack=True)
        cls.launchProtocol(cls.protImportAvg)
        return cls.protImportAvg


class TestEmanTomoBase(TestEmanBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def setData(cls, projectData='tomo-em'):
        from tomo.tests import DataSet
        cls.dataset = DataSet.getDataSet(projectData)
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
        cls.coords3D_Large = cls.dataset.getFile('overview_wbp_large.txt')
        cls.inputSetOfSubTomogram = cls.dataset.getFile('subtomo')
        cls.smallTomogram = cls.dataset.getFile('coremask_normcorona.mrc')


class TestEmanTomoInitialModel(TestEmanTomoBase):
    """This class check if the protocol to extract particles
    in Relion works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")

        # if emantomo.Plugin.getActiveVersion(versions=[emantomo.V2_39]):
        boxSize = 128
        # else:
        #     boxSize = 32
        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=boxSize,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                              inputTomograms=protImportTomogram.outputTomograms,
                                              inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                              tomoSource=0,
                                              doInvert=False,
                                              doNormalize=False,
                                              boxSize=boxSize)

        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(getattr(protTomoExtraction, OutputExtraction.subtomograms.name),
                             "There was a problem with SetOfSubtomogram output")
        return protTomoExtraction

    def _performFinalValidation(self, protInitialModel):
        averageSubTomogram = protInitialModel.subtomogramAverage
        self.assertEqual(os.path.basename(averageSubTomogram.getFirstItem().getFileName()), "output.hdf")
        self.assertEqual(averageSubTomogram.getFirstItem().getSamplingRate(), 20.0)
        self.assertEqual(averageSubTomogram.getSamplingRate(), 20.0)

        setOfSubTomograms = protInitialModel.outputParticles
        self.assertEqual(setOfSubTomograms.getSize(), 5)
        self.assertEqual(setOfSubTomograms.getCoordinates3D().getObjValue().getSize(), 5)
        self.assertEqual(setOfSubTomograms.getSamplingRate(), 20.0)

        for subTomogram in setOfSubTomograms:
            self.assertEqual(subTomogram.getSamplingRate(), 20)
            self.assertTrue(hasattr(subTomogram, "coverage"))
            self.assertTrue(hasattr(subTomogram, "score"))
            matrix = subTomogram.getTransform(convention=TR_EMAN).getMatrix()
            self.assertEqual(matrix.shape, (4, 4))

        summary = protInitialModel.summary()
        self.assertIsNotNone(summary)
        self.assertEqual(len(summary), 2)
        self.assertEqual(summary[0], "Particles: 5")
        self.assertTrue(summary[1].startswith("Reference file used:"))

    def _runTomoSubtomogramInitialModelWithSubtomo(self):
        protTomoExtraction = self._runPreviousProtocols()

        particles = getattr(protTomoExtraction,OutputExtraction.subtomograms.name)

        protInitialModel = self.newProtocol(EmanProtTomoInitialModel,
                                            particles=particles,
                                            reference=protTomoExtraction,
                                            symmetry="c1",
                                            gaussFilter=-1.5,
                                            filterto=0.03,
                                            fourier=False,
                                            batchSize=20,
                                            learningRate=1,
                                            numberOfIterations=2,
                                            numberOfBatches=1,
                                            shrink=4,
                                            applySim=False)
        protInitialModel.reference.setExtended("%s.1" % OutputExtraction.subtomograms.name)

        self.launchProtocol(protInitialModel)

        self.assertIsNotNone(protInitialModel.subtomogramAverage,
                             "There was a problem with subTomograms output")
        self.assertIsNotNone(protInitialModel.outputParticles,
                             "There was a problem with particles output")

        return protInitialModel

    def _runTomoSubtomogramInitialModelWithVolume(self):
        protTomoExtraction = self._runPreviousProtocols()

        particles = getattr(protTomoExtraction, OutputExtraction.subtomograms.name)

        self.dataset = DataSet.getDataSet('eman')
        self.vol = self.dataset.getFile('volume')
        self.protImportVol = self.runImportVolumes(self.vol, 3.5)

        self.assertIsNotNone(self.protImportVol.outputVolume,
                             "There was a problem with SetOfSubtomogram output")

        protInitialModel = self.newProtocol(EmanProtTomoInitialModel,
                                            particles=particles,
                                            reference=self.protImportVol.outputVolume,
                                            symmetry="c1",
                                            gaussFilter=-1.5,
                                            filterto=0.03,
                                            fourier=False,
                                            batchSize=20,
                                            learningRate=1,
                                            numberOfIterations=2,
                                            numberOfBatches=1,
                                            shrink=4,
                                            applySim=False)

        self.launchProtocol(protInitialModel)

        self.assertIsNotNone(protInitialModel.subtomogramAverage,
                             "There was a problem with subTomograms output")
        self.assertIsNotNone(protInitialModel.outputParticles,
                             "There was a problem with particles output")

        return protInitialModel

    def test_initialModelOutputWithSubtomo(self):
        protInitialModel = self._runTomoSubtomogramInitialModelWithSubtomo()

        self._performFinalValidation(protInitialModel)

        return protInitialModel

    def test_initialModelOutputWithVolume(self):
        protInitialModel = self._runTomoSubtomogramInitialModelWithVolume()

        self._performFinalValidation(protInitialModel)

        return protInitialModel


class TestEmanTomoTempMatch(TestEmanTomoBase):
    """This class check if the program Template Matching
    from Eman works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runTomoTempMatch(self):
        protImportTomogramBig = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                                 filesPath=self.tomogram,
                                                 samplingRate=5)

        protImportTomogramSmall = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                                   filesPath=self.smallTomogram,
                                                   samplingRate=5)

        self.launchProtocol(protImportTomogramBig)

        self.launchProtocol(protImportTomogramSmall)

        self.assertSetSize(protImportTomogramBig.outputTomograms, size=1,
                           msg="There was a problem with tomogram output")

        self.assertSetSize(protImportTomogramSmall.outputTomograms, size=1, msg="There was a problem with tomogram "
                                                                                "output")

        self.dataset = DataSet.getDataSet('eman')
        self.vol = self.dataset.getFile('volume')
        self.protImportVol = self.runImportVolumes(self.vol, 3.5)

        self.assertIsNotNone(self.protImportVol.outputVolume,
                             "There was a problem with SetOfSubtomogram output")

        protTomoTempMatchBig = self.newProtocol(EmanProtTomoTempMatch,
                                                inputSet=protImportTomogramBig.outputTomograms,
                                                ref=self.protImportVol.outputVolume,
                                                boxSize=128,
                                                sym="c1")

        protTomoTempMatchSmall = self.newProtocol(EmanProtTomoTempMatch,
                                                  inputSet=protImportTomogramSmall.outputTomograms,
                                                  ref=self.protImportVol.outputVolume,
                                                  boxSize=128,
                                                  sym="d1")

        self.launchProtocol(protTomoTempMatchBig)
        self.launchProtocol(protTomoTempMatchSmall)
        return protTomoTempMatchBig, protTomoTempMatchSmall

    def test_TempMatch(self):
        # if EmanProtTomoTempMatch.isDisabled():
        #     print("Test Cancelled. Template Matching is not supported in Eman 2.31")
        #     return
        protTomoTempMatch = self._runTomoTempMatch()

        outputCoordsBig = protTomoTempMatch[0].output3DCoordinates
        if emantomo.Plugin.isVersion(emantomo.constants.V_CB):
            self.assertAlmostEqual(outputCoordsBig.getSize(), 154, delta=1)
        else:
            self.assertAlmostEqual(outputCoordsBig.getSize(), 500, delta=1)
        self.assertEqual(outputCoordsBig.getBoxSize(), 128)
        self.assertEqual(outputCoordsBig.getSamplingRate(), 5)

        outputCoordsSmall = protTomoTempMatch[1].output3DCoordinates
        if emantomo.Plugin.isVersion(emantomo.constants.V_CB):
            self.assertAlmostEqual(outputCoordsBig.getSize(), 154, delta=1)
        else:
            self.assertAlmostEqual(outputCoordsSmall.getSize(), 145, delta=1)
        self.assertEqual(outputCoordsSmall.getBoxSize(), 128)
        self.assertEqual(outputCoordsSmall.getSamplingRate(), 5)

        return protTomoTempMatch

class TestEmanTomoClassifySubtomos(TestEmanTomoBase):
    """This class check if the protocol classify subtomos works properly.
    """

    @classmethod
    def setUpClass(cls):
        from tomo.tests import DataSet
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.inputSetOfSubTomogram = cls.dataset.getFile('subtomo')
        cls.coords3D_Large = cls.dataset.getFile('overview_wbp_large.txt')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')

    def _runTomoClassifySubtomos(self):
        from tomo.protocols import ProtImportCoordinates3D, ProtImportTomograms
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        coords = pwutils.removeBaseExt(self.coords3D)
        coords = protImportTomogram._getExtraPath(coords + '.txt')
        pwutils.createAbsLink(self.coords3D_Large, coords)
        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto='eman',
                                                   filesPath=coords,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")

        doInvert = False
        doNormalize = False
        boxSize = 32
        protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                              inputTomograms=protImportTomogram.outputTomograms,
                                              inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                              downsampleType=0,
                                              doInvert=doInvert,
                                              doNormalize=doNormalize,
                                              boxSize=boxSize)

        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(getattr(protTomoExtraction, OutputExtraction.subtomograms.name),
                             "There was a problem with SetOfSubtomogram output")

        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)
        self.launchProtocol(protImportTomogram)

        protClassifySubtomos = self.newProtocol(EmanProtTomoClassifySubtomos,
                                                inputSetOfSubTomogram=getattr(protTomoExtraction, OutputExtraction.subtomograms.name))

        self.launchProtocol(protClassifySubtomos)
        self.assertEqual(protClassifySubtomos.outputVolumes.getSize(), 2,
                         "Wrong number of volumes registered")
        self.assertEqual(protClassifySubtomos.outputClasses.getSize(), 2,
                         "Wrong number of classes registered")
        self.assertEqual(protClassifySubtomos.outputClasses[1].getSize(), 12,
                         "Wrong number of particles in first class")
        self.assertEqual(protClassifySubtomos.outputClasses[2].getSize(), 3,
                         "Wrong number of particles in second class")
        return protClassifySubtomos

    def test_classify_subtomos(self):
        protClassifySubtomos = self._runTomoClassifySubtomos()
        return protClassifySubtomos


class TestEmanTomoReconstruction(TestEmanTomoBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _getProtImportTs(self):
        return self.newProtocol(
            tomo.protocols.ProtImportTs,
            filesPath=TestEmanTomoReconstruction.dataset.getPath(),
            filesPattern='tomo{TS}.hdf',
            minAngle=-55,
            maxAngle=60,
            stepAngle=2,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=1.35,
            doseInitial=0,
            dosePerFrame=0.3,
            tiltAxisAngle=87.1
        )

    def _getProtAlignTs(self, protImportTs):
        return self.newProtocol(
            EmanProtAlignTs,
            tiltSeries=protImportTs.outputTiltSeries,
            tiltStep=2.0,
            niter='1,1,1,1',
            bxsz=64)

    def _getProtReconstruct(self, protAlignTs):
        return self.newProtocol(
            EmanProtTomoReconstruction,
            tiltSeries=protAlignTs.alignedTiltSeries)

    def _runPreviousProtocols(self):
        protImportTs = self._getProtImportTs()
        self.launchProtocol(protImportTs)
        self.assertIsNotNone(protImportTs.outputTiltSeries, "Output tilt series not found")

        protAlignTs = self._getProtAlignTs(protImportTs)
        self.launchProtocol(protAlignTs)
        self.assertIsNotNone(protImportTs.outputTiltSeries, "Output aligned tilt series not found")

        protReconstruct = self._getProtReconstruct(protAlignTs)
        self.launchProtocol(protReconstruct)

        return protReconstruct

    def _validateOutput(self, protReconstruct):
        tomograms = list(protReconstruct.tomograms)
        # 1 tomogram per input file
        self.assertEqual(len(tomograms), 1)
        for tomogram in tomograms:
            self.assertEqual(tomogram.getSamplingRate(), 5.4)
            self.assertEqual(tomogram.getDimensions(), (1120, 1120, 256))

    def test_protocol(self):
        protTomoExtraction = self._runPreviousProtocols()
        self._validateOutput(protTomoExtraction)
        self.assertTrue(protTomoExtraction.summary())
        self.assertTrue(protTomoExtraction.methods())
