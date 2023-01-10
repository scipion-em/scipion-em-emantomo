# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
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
import numpy as np
from emantomo.protocols import EmanProtTomoExtraction
from emantomo.protocols.protocol_tomo_extraction_from_tomo import SAME_AS_PICKING, OutputExtraction, OTHER
from imod.protocols import ProtImodTomoNormalization
from imod.protocols.protocol_base import OUTPUT_TOMOGRAMS_NAME
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtImportTomograms, ProtImportCoordinates3DFromScipion
from tomo.tests import EMD_10439, DataSetEmd10439


class TestEmantomoBase(BaseTest):
    nParticles = 39
    boxSize = 44
    binnedBoxSize = 22
    origSRate = 13.68
    binnedSRate = 27.36
    ds = DataSet.getDataSet(EMD_10439)

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def runImportTomograms(cls):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.origSRate, )

        cls.launchProtocol(protImportTomogram)
        tomoImported = protImportTomogram.outputTomograms
        cls.assertIsNotNone(tomoImported, "There was a problem with tomogram output")
        return tomoImported

    @classmethod
    def runImport3dCoords(cls, tomoImported):
        # Import coordinates
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        protImportCoordinates3d = cls.newProtocol(ProtImportCoordinates3DFromScipion,
                                                  sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
                                                  importTomograms=tomoImported,
                                                  boxSize=cls.boxSize)

        cls.launchProtocol(protImportCoordinates3d)
        coordsImported = protImportCoordinates3d.outputCoordinates
        cls.assertIsNotNone(coordsImported, "There was a problem with the 3D coordinates output")
        return coordsImported

    @classmethod
    def runBinTomograms(cls, tomoImported, binning=2):
        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram binning:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=tomoImported,
                                                binning=binning)

        cls.launchProtocol(protTomoNormalization)
        tomosBinned = getattr(protTomoNormalization, OUTPUT_TOMOGRAMS_NAME, None)
        cls.assertIsNotNone(tomosBinned, 'No tomograms were genetated in tomo normalization.')
        return tomosBinned

    @classmethod
    def runExtractSubtomograms(cls, coordsImported, tomoSource=SAME_AS_PICKING, tomograms=None, boxSize=None):
        # Extract subtomograms
        print(magentaStr("\n==> Extracting the subtomograms:"))
        protLabel = 'Extraction - same as picking'
        argsDict = {'inputCoordinates': coordsImported,
                    'tomoSource': tomoSource,
                    'boxSize': boxSize,
                    'doInvert': True}
        if tomoSource != SAME_AS_PICKING:
            argsDict['tomoSource'] = OTHER
            argsDict['inputTomograms'] = tomograms
            protLabel = 'Extraction - another tomo source'

        protTomoExtraction = cls.newProtocol(EmanProtTomoExtraction, **argsDict)
        protTomoExtraction.setObjLabel(protLabel)
        cls.launchProtocol(protTomoExtraction)
        subtomosExtracted = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        cls.assertIsNotNone(subtomosExtracted, "There was a problem with the subtomograms extraction")
        return subtomosExtracted

    def check3dTransformMatrix(self, outMatrix):
        transfMatrixShape = (4, 4)
        self.assertIsNotNone(outMatrix)
        if type(outMatrix) != np.ndarray:
            outMatrix = np.array(outMatrix)
        self.assertIsNotNone(outMatrix)
        self.assertTrue(outMatrix.shape, transfMatrixShape)
        self.assertFalse(np.array_equal(outMatrix, np.eye(4)))

    def checkShiftsScaling(self, inTransform, outTransform, scaleFactor):
        # Check if the shifts have been scaled properly
        sx, sy, sz = inTransform.getShifts()
        osx, osy, osz = outTransform.getShifts()
        for inShift, outShift in zip([sx, sy, sz], [osx, osy, osz]):
            self.assertEqual(outShift, scaleFactor * inShift)

    @staticmethod
    def getMinAndMaxCoordValuesFromSet(inSet):
        if type(inSet) == SetOfSubTomograms:
            inSet = inSet.getCoordinates3D()

        dataDict = inSet.aggregate(['MAX'], '_tomoId', ['_x', '_y', '_z'])
        xcoords, ycoords, zcoords = zip(*[(d['_x'], d['_y'], d['_z']) for d in dataDict])
        return np.array([min(xcoords),
                         max(xcoords),
                         min(ycoords),
                         max(ycoords),
                         min(zcoords),
                         max(zcoords)])
