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
from xmipp3.constants import MASK3D_CYLINDER
from xmipp3.protocols import XmippProtCreateMask3D
from xmipp3.protocols.protocol_preprocess.protocol_create_mask3d import SOURCE_GEOMETRY
from emantomo.protocols import EmanProtTomoExtraction, EmanProtSubTomoAverage
from emantomo.protocols.protocol_average_subtomos import OutputsAverageSubtomos
from emantomo.protocols.protocol_extraction_from_tomo import SAME_AS_PICKING, OutputExtraction, OTHER
from imod.protocols import ProtImodTomoNormalization
from imod.protocols.protocol_base import OUTPUT_TOMOGRAMS_NAME
from pyworkflow.tests import DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTomograms, ProtImportCoordinates3DFromScipion
from tomo.protocols.protocol_import_coordinates_from_scipion import outputObjs
from tomo.tests import EMD_10439, DataSetEmd10439
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestEmantomoStaClassicBase(TestBaseCentralizedLayer):
    ds = None
    coordsImported = None
    tomosBinned = None
    tomoImported = None
    subtomosExtracted = None
    nParticles = 39
    boxSize = 44
    binnedBoxSize = 22
    origSRate = 13.68
    binnedSRate = 27.36

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(EMD_10439)

    @classmethod
    def runPreviousProtocols(cls):
        cls.tomoImported = cls.runImportTomograms()  # Import tomograms
        cls.tomosBinned = cls.runBinTomograms(cls.tomoImported)  # Bin the tomogram to make it smaller
        cls.coordsImported = cls.runImport3dCoords(
            cls.tomosBinned)  # Import the coordinates from the binned tomogram
        # Extract subtomograms
        cls.subtomosExtracted = cls.runExtractSubtomograms(cls.coordsImported,
                                                           tomoSource=SAME_AS_PICKING,
                                                           boxSize=cls.binnedBoxSize)

    # --------------- RUNS functions -----------------------
    @classmethod
    def runCreate3dMask(cls):
        print(magentaStr("\n==> Generating the 3D mask:"))
        protCreateParticleMask = cls.newProtocol(XmippProtCreateMask3D,
                                                 source=SOURCE_GEOMETRY,
                                                 samplingRate=cls.binnedSRate,
                                                 size=cls.binnedBoxSize,
                                                 geo=MASK3D_CYLINDER,
                                                 radius=6,
                                                 shiftCenter=True,
                                                 centerZ=3,
                                                 height=15,
                                                 doSmooth=True)
        cls.launchProtocol(protCreateParticleMask)
        genMask = getattr(protCreateParticleMask, 'outputMask', None)
        cls.assertIsNotNone(genMask, 'the 3D mask was not generated')
        return genMask

    @classmethod
    def runImportTomograms(cls):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.origSRate, )

        cls.launchProtocol(protImportTomogram)
        tomoImported = protImportTomogram.Tomograms
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
        coordsImported = getattr(protImportCoordinates3d, outputObjs.coordinates.name, None)
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
    def runExtractSubtomograms(cls, coordsImported, tomoSource=SAME_AS_PICKING, tomograms=None, boxSize=None, label=None):
        # Extract subtomograms
        print(magentaStr("\n==> Extracting the subtomograms:"))
        protLabel = 'Extraction - same as picking' if label is None else label
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

    @classmethod
    def runAverageSubtomograms(cls):
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(EmanProtSubTomoAverage, inputSetOfSubTomogram=cls.subtomosExtracted)
        cls.launchProtocol(protAvgSubtomo)
        avgSubtomo = getattr(protAvgSubtomo, OutputsAverageSubtomos.averageSubTomos.name, None)
        cls.assertIsNotNone(avgSubtomo, "There was a problem calculating the average subtomogram")
        return avgSubtomo
