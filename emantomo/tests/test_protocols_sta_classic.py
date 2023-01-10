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
from os.path import exists
import numpy as np
from xmipp3.constants import MASK3D_CYLINDER
from xmipp3.protocols import XmippProtCreateMask3D
from xmipp3.protocols.protocol_preprocess.protocol_create_mask3d import SOURCE_GEOMETRY
from pyworkflow.utils import magentaStr
from .test_eman_base import TestEmantomoBase
from ..constants import EMAN_COVERAGE, EMAN_SCORE
from tomo.constants import TR_EMAN
from ..protocols import EmanProtTomoInitialModel
from ..protocols.protocol_average_subtomos import OutputsAverageSubtomos, EmanProtSubTomoAverage
from ..protocols.protocol_tomo_extraction_from_tomo import SAME_AS_PICKING
from ..protocols.protocol_tomo_initialmodel import OutputsInitModel
from ..protocols.protocol_tomo_subtomogram_refinement import EmanTomoRefinementOutputs, EmanProtTomoRefinement


class TestEmanTomoSubtomogramRefinement(TestEmantomoBase):
    """This class check if the protocol Subtomogram refinement works properly.
    """

    avgSubtomo = None
    subtomosExtracted = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.subtomosExtracted, cls.avgSubtomo = cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        tomoImported = super().runImportTomograms()  # Import tomograms
        tomosBinned = super().runBinTomograms(tomoImported)  # Bin the tomogram to make it smaller
        coordsImported = super().runImport3dCoords(tomosBinned)  # Import the coordinates from the binned tomogram
        # Extract subtomograms
        subtomosExtracted = super().runExtractSubtomograms(coordsImported,
                                                           tomoSource=SAME_AS_PICKING,
                                                           boxSize=super().binnedBoxSize)
        # Average the subtomograms
        avgSubtomo = cls.runAverageSubtomograms(subtomosExtracted)
        return subtomosExtracted, avgSubtomo

    @classmethod
    def runAverageSubtomograms(cls, subtomosExtracted):
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(EmanProtSubTomoAverage, inputSetOfSubTomogram=subtomosExtracted)
        cls.launchProtocol(protAvgSubtomo)
        avgSubtomo = getattr(protAvgSubtomo, OutputsAverageSubtomos.averageSubTomos.name, None)
        cls.assertIsNotNone(avgSubtomo, "There was a problem calculating the average subtomogram")
        return avgSubtomo

    @classmethod
    def runGenInitialModel(cls):
        print(magentaStr("\n==> Generating the initial model:"))
        protInitModel = cls.newProtocol(EmanProtTomoInitialModel,
                                        particles=cls.subtomosExtracted,
                                        reference=cls.avgSubtomo,
                                        numberOfIterations=2)
        cls.launchProtocol(protInitModel)
        initModel = getattr(protInitModel, OutputsInitModel.average.name, None)
        cls.assertIsNotNone(initModel, "There was a problem calculating the initial model")
        return initModel

    @classmethod
    def runCreate3dMask(cls):
        print(magentaStr("\n==> Generating the 3D mask:"))
        protCreateParticleMask = cls.newProtocol(XmippProtCreateMask3D,
                                                 source=SOURCE_GEOMETRY,
                                                 samplingRate=super().binnedSRate,
                                                 size=super().binnedBoxSize,
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
    def runTomoSubtomogramRefinement(cls, mask=None):

        print(magentaStr("\n==> Refining the subtomograms:"))
        inputDict = {'inputSetOfSubTomogram': cls.subtomosExtracted,
                     'inputRef': cls.avgSubtomo,
                     'pkeep': 1,
                     'niter': 2,
                     'numberOfThreads': 10}

        objLabel = 'Subtomo refinement'
        if mask:
            inputDict['maskFile'] = mask
            objLabel += ' with mask'

        protTomoRefinement = cls.newProtocol(EmanProtTomoRefinement, **inputDict)
        protTomoRefinement.setObjLabel(objLabel)
        cls.launchProtocol(protTomoRefinement)
        return protTomoRefinement

    def test_averageSubtomograms(self):
        self.checkAverage(self.avgSubtomo, boxSize=super().binnedBoxSize)

    def test_genInitialModel(self):
        initModel = self.runGenInitialModel()
        self.checkAverage(initModel, boxSize=super().binnedBoxSize, halvesExpected=False)

    def test_subtomoRefinement(self):
        protTomoSubtomogramRefinement = self.runTomoSubtomogramRefinement()
        self.checkRefinementResults(protTomoSubtomogramRefinement)

    def test_subTomoRefinementWithMask(self):
        protTomoSubtomogramRefinement = self.runTomoSubtomogramRefinement(mask=self.runCreate3dMask())
        self.checkRefinementResults(protTomoSubtomogramRefinement)

    def checkAverage(self, avg, boxSize=None, halvesExpected=True):
        testBoxSize = (boxSize, boxSize, boxSize)
        self.assertTrue(exists(avg.getFileName()), "Average %s does not exists" % avg.getFileName())
        self.assertTrue(avg.getFileName().endswith(".mrc"))
        # The imported coordinates correspond to a binned 2 tomogram
        self.assertEqual(avg.getSamplingRate(), super().binnedSRate)
        self.assertEqual(avg.getDimensions(), testBoxSize)
        # Check the halves
        if halvesExpected:
            self.assertTrue(avg.hasHalfMaps(), "Halves not registered")
            half1, half2 = avg.getHalfMaps().split(',')
            self.assertTrue(exists(half1), msg="Average 1st half %s does not exists" % half1)
            self.assertTrue(exists(half2), msg="Average 2nd half %s does not exists" % half2)

    def checkRefinementResults(self, protRefineSubtomos):
        refinedSubtomos = getattr(protRefineSubtomos, EmanTomoRefinementOutputs.subtomograms.name, None)
        avgSubtomo = getattr(protRefineSubtomos, EmanTomoRefinementOutputs.subtomogramAverage.name, None)
        testDims = (super().binnedBoxSize, super().binnedBoxSize, super().binnedBoxSize)

        # Generation checking
        self.assertIsNotNone(refinedSubtomos)
        self.assertIsNotNone(avgSubtomo)

        # Avg checking
        self.checkAverage(avgSubtomo, boxSize=super().binnedBoxSize)

        # Subtomos checking
        self.assertEqual(refinedSubtomos.getDimensions(), testDims)
        self.assertSetSize(refinedSubtomos, self.nParticles, msg='Expected size of the generated subtomograms is '
                                                                 'different than expected: %i != %i' %
                                                                 (refinedSubtomos.getSize(), self.nParticles))
        self.assertEqual(refinedSubtomos.getCoordinates3D().getSize(), self.nParticles)
        for subTomogram in refinedSubtomos:
            self.assertEqual(subTomogram.getSamplingRate(), super().binnedSRate)
            self.assertTrue(hasattr(subTomogram, EMAN_COVERAGE))
            self.assertTrue(hasattr(subTomogram, EMAN_SCORE))
            matrix = subTomogram.getTransform(convention=TR_EMAN).getMatrix()
            self.assertEqual(matrix.shape, (4, 4))
            self.assertFalse(np.array_equal(matrix, np.eye(4)))  # Then it contains the angles and shitfs

        # FSCs checking
        fscs = getattr(protRefineSubtomos, EmanTomoRefinementOutputs.FSCs.name)
        self.assertSetSize(fscs, 3, msg="FSCs not registered properly")

