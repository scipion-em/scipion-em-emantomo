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

import numpy as np

from imod.protocols import ProtImodTomoNormalization
from imod.protocols.protocol_base import OUTPUT_TOMOGRAMS_NAME
from pyworkflow.tests import setupTestProject, DataSet, BaseTest
from pyworkflow.utils import magentaStr
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from . import EMD_10439, DataSetEmd10439
from ..constants import EMAN_COVERAGE, EMAN_SCORE
import tomo.protocols
from tomo.constants import TR_EMAN
from ..protocols.protocol_average_subtomos import OutputsAverageSubtomos, EmanProtSubTomoAverage
from ..protocols.protocol_tomo_extraction_from_tomo import OutputExtraction, SAME_AS_PICKING, EmanProtTomoExtraction
from ..protocols.protocol_tomo_subtomogram_refinement import EmanTomoRefinementOutputs, EmanProtTomoRefinement


class TestEmanTomoSubtomogramRefinement(BaseTest):
    """This class check if the protocol Subtomogram refinement works properly.
    """

    nParticles = 39
    boxSize = 32
    origSRate = 13.68
    ds = DataSet.getDataSet(EMD_10439)

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.subtomosExtracted, cls.avgSubtomo = cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(tomo.protocols.ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.origSRate, )

        cls.launchProtocol(protImportTomogram)
        tomoImported = protImportTomogram.outputTomograms
        cls.assertIsNotNone(tomoImported, "There was a problem with tomogram output")

        # Import coordinates
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        protImportCoordinates3d = cls.newProtocol(tomo.protocols.ProtImportCoordinates3DFromScipion,
                                                  sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
                                                  importTomograms=tomoImported,
                                                  boxSize=cls.boxSize)

        cls.launchProtocol(protImportCoordinates3d)
        coordsImported = protImportCoordinates3d.outputCoordinates
        cls.assertIsNotNone(coordsImported, "There was a problem with the 3D coordinates output")

        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram normalization:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=protImportTomogram.outputTomograms,
                                                binning=2)

        cls.launchProtocol(protTomoNormalization)
        tomosBinned = getattr(protTomoNormalization, OUTPUT_TOMOGRAMS_NAME, None)
        cls.assertIsNotNone(tomosBinned, 'No tomograms were genetated in tomo normalization.')

        # Extract subtomograms
        print(magentaStr("\n==> Extracting the subtomograms:"))
        protTomoExtraction = cls.newProtocol(EmanProtTomoExtraction,
                                             inputTomograms=tomosBinned,
                                             inputCoordinates=coordsImported,
                                             tomoSource=SAME_AS_PICKING,
                                             boxSize=cls.boxSize)

        cls.launchProtocol(protTomoExtraction)
        subtomosExtracted = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        cls.assertIsNotNone(subtomosExtracted, "There was a problem with the subtomograms extraction")

        # Average the subtomograms
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(EmanProtSubTomoAverage, inputSetOfSubTomogram=subtomosExtracted)
        cls.launchProtocol(protAvgSubtomo)
        avgSubtomo = getattr(protAvgSubtomo, OutputsAverageSubtomos.averageSubTomos.name, None)
        cls.assertIsNotNone(avgSubtomo, "There was a problem calculating the average subtomogram")

        return subtomosExtracted, avgSubtomo

    def _performFinalValidation(self, protTomoSubtomogramRefinement):
        outputSetOfSubTomograms = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.subtomograms.name)
        averageSubTomogram = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name)

        self.assertTrue(os.path.exists(averageSubTomogram.getFileName()))

        # Check average is a mrc file
        self.assertTrue(averageSubTomogram.getFileName().endswith(".mrc"))

        self.assertEqual(averageSubTomogram.getSamplingRate(), self.origSRate)

        self.assertEqual(outputSetOfSubTomograms.getDimensions(), (self.boxSize, self.boxSize, self.boxSize))
        self.assertEqual(outputSetOfSubTomograms.getSize(), self.nParticles)
        self.assertEqual(outputSetOfSubTomograms.getCoordinates3D().getSize(), self.nParticles)

        for subTomogram in outputSetOfSubTomograms:
            self.assertEqual(subTomogram.getSamplingRate(), self.origSRate)
            self.assertTrue(hasattr(subTomogram, EMAN_COVERAGE))
            self.assertTrue(hasattr(subTomogram, EMAN_SCORE))
            matrix = subTomogram.getTransform(convention=TR_EMAN).getMatrix()
            self.assertEqual(matrix.shape, (4, 4))
            self.assertFalse((matrix == np.eye(4)).all())

        # FSCs
        fscs = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.FSCs.name)
        self.assertSetSize(fscs, 3, msg="FSCs not registered properly")

    def _runTomoSubtomogramRefinement(self, niter=2, mass=500.0, threads=1, pkeep=1, goldstandard=-1,
                                      goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0, refVol=None):

        print(magentaStr("\n==> Refining the subtomograms:"))
        inputDict = {'inputSetOfSubTomogram': self.subtomosExtracted,
                     'niter': niter,
                     'mass': mass,
                     'threads': threads,
                     'pkeep': pkeep,
                     'goldstandard': goldstandard,
                     'goldcontinue': goldcontinue,
                     'sym': sym,
                     'localfilter': localfilter,
                     'maxtilt': maxtilt,
                     'numberOfThreads': 10}

        if refVol:
            inputDict['inputRef'] = refVol

        protTomoRefinement = self.newProtocol(EmanProtTomoRefinement, **inputDict)
        self.launchProtocol(protTomoRefinement)

        average = getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name)
        self.assertIsNotNone(average, "There was a problem with average output")
        self.assertTrue(os.path.exists(average.getFileName()), "Average %s does not exists" % average.getFileName())
        self.assertTrue(average.hasHalfMaps(), "Halves not registered")

        half1, half2 = average.getHalfMaps().split(',')
        self.assertTrue(os.path.exists(half1), msg="Average 1st half %s does not exists" % half1)
        self.assertTrue(os.path.exists(half2), msg="Average 2nd half %s does not exists" % half2)

        self.assertSetSize(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomograms.name), None,
                           "There was a problem with particles output")

        return protTomoRefinement

    def test_defaultSubTomogramRefinementWithSubTomo(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinement()
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement

    def test_defaultSubTomogramRefinementWithVolume(self):
        # Generate the reference volume averaging the extracted particles
        # Refine the subtomograms
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinement(refVol=self.avgSubtomo)
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement
