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
from os.path import exists
import numpy as np
from xmipp3.constants import MASK3D_CYLINDER
from xmipp3.protocols import XmippProtCreateMask3D
from xmipp3.protocols.protocol_preprocess.protocol_create_mask3d import SOURCE_GEOMETRY

from imod.protocols import ProtImodTomoNormalization
from imod.protocols.protocol_base import OUTPUT_TOMOGRAMS_NAME
from pyworkflow.tests import setupTestProject, DataSet, BaseTest
from pyworkflow.utils import magentaStr
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from tomo.tests import EMD_10439, DataSetEmd10439
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
    boxSize = 44
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

    def testSubTomoRefinement(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinement()
        self._performFinalValidation(protTomoSubtomogramRefinement)

    def testSubTomoRefinementWithMask(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinement(mask=self._runCreate3dMask())
        self._performFinalValidation(protTomoSubtomogramRefinement)

    def _runTomoSubtomogramRefinement(self, mask=None):

        print(magentaStr("\n==> Refining the subtomograms:"))
        inputDict = {'inputSetOfSubTomogram': self.subtomosExtracted,
                     'inputRef': self.avgSubtomo,
                     'pkeep': 1,
                     'niter': 2,
                     'numberOfThreads': 10}

        objLabel = 'Subtomo refinement'
        if mask:
            inputDict['maskFile'] = mask
            objLabel += ' with mask'

        protTomoRefinement = self.newProtocol(EmanProtTomoRefinement, **inputDict)
        protTomoRefinement.setObjLabel(objLabel)
        self.launchProtocol(protTomoRefinement)

        average = getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name)
        self.assertIsNotNone(average, "There was a problem with average output")
        self.assertTrue(exists(average.getFileName()), "Average %s does not exists" % average.getFileName())
        self.assertTrue(average.hasHalfMaps(), "Halves not registered")

        half1, half2 = average.getHalfMaps().split(',')
        self.assertTrue(exists(half1), msg="Average 1st half %s does not exists" % half1)
        self.assertTrue(exists(half2), msg="Average 2nd half %s does not exists" % half2)

        self.assertSetSize(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomograms.name), None,
                           "There was a problem with particles output")

        return protTomoRefinement

    def _performFinalValidation(self, protRefineSubtomos):
        refinedSubtomos = getattr(protRefineSubtomos, EmanTomoRefinementOutputs.subtomograms.name, None)
        avgSubtomo = getattr(protRefineSubtomos, EmanTomoRefinementOutputs.subtomogramAverage.name, None)
        testBoxSize = (self.boxSize, self.boxSize, self.boxSize)

        # Generation checking
        self.assertIsNotNone(refinedSubtomos)
        self.assertIsNotNone(avgSubtomo)

        # Avg checking
        self.assertTrue(exists(avgSubtomo.getFileName()))
        self.assertTrue(avgSubtomo.getFileName().endswith(".mrc"))
        self.assertEqual(avgSubtomo.getSamplingRate(), self.origSRate)
        self.assertEqual(avgSubtomo.getDimensions(), testBoxSize)

        # Subtomos checking
        self.assertEqual(refinedSubtomos.getDimensions(), testBoxSize)
        self.assertSetSize(refinedSubtomos, self.nParticles, msg='Expected size of the generated subtomograms is '
                                                                 'different than expected: %i != %i' %
                                                                 (refinedSubtomos.getSize(), self.nParticles))
        self.assertEqual(refinedSubtomos.getCoordinates3D().getSize(), self.nParticles)
        for subTomogram in refinedSubtomos:
            self.assertEqual(subTomogram.getSamplingRate(), self.origSRate)
            self.assertTrue(hasattr(subTomogram, EMAN_COVERAGE))
            self.assertTrue(hasattr(subTomogram, EMAN_SCORE))
            matrix = subTomogram.getTransform(convention=TR_EMAN).getMatrix()
            self.assertEqual(matrix.shape, (4, 4))
            self.assertFalse((matrix == np.eye(4)).all())  # Then it contains the angles and shitfs

        # FSCs checking
        fscs = getattr(protRefineSubtomos, EmanTomoRefinementOutputs.FSCs.name)
        self.assertSetSize(fscs, 3, msg="FSCs not registered properly")

    @classmethod
    def _runCreate3dMask(cls):
        print(magentaStr("\n==> Generating the 3D mask:"))
        protCreateParticleMask = cls.newProtocol(XmippProtCreateMask3D,
                                                 source=SOURCE_GEOMETRY,
                                                 samplingRate=cls.origSRate,
                                                 size=cls.boxSize,
                                                 geo=MASK3D_CYLINDER,
                                                 radius=12,
                                                 shiftCenter=True,
                                                 centerZ=6,
                                                 height=30,
                                                 doSmooth=True)
        cls.launchProtocol(protCreateParticleMask)
        genMask = getattr(protCreateParticleMask, 'outputMask', None)
        cls.assertIsNotNone(genMask, 'the 3D mask was not generated')
        return genMask
