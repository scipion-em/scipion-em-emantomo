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
# from os.path import exists
from pyworkflow.utils import magentaStr
from .test_eman_sta_classic_base import TestEmantomoStaClassicBase
# from ..constants import SPTCLS_00_DIR
from tomo.constants import TR_EMAN
from ..protocols import EmanProtTomoInitialModel, EmanProtTemplateMatching
# from emantomo.protocols.deprecated_20230914.protocol_pca_kmeans_classify_subtomos import pcaOutputObjects, EmanProtPcaKMeansClassifySubtomos
from ..protocols.protocol_tomo_initialmodel import OutputsInitModel
from ..protocols.protocol_tomo_subtomogram_refinement import EmanTomoRefinementOutputs, EmanProtTomoRefinement


class TestEmanTomoAverageSubtomogramsStaClassic(TestEmantomoStaClassicBase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        try:
            super().runPreviousProtocols()
        except Exception as e:
            raise Exception('Some of the previous protocols failed --> \n%s' % e)

    def test_averageSubtomograms(self):
        avgSubtomo = super().runAverageSubtomograms()
        super().checkAverage(avgSubtomo,
                             expectedSRate=self.binnedSRate,
                             expectedBoxSize=self.binnedBoxSize)


class TestEmanTomoInitialModelStaClassic(TestEmantomoStaClassicBase):
    avgSubtomo = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        try:
            super().runPreviousProtocols()
            cls.avgSubtomo = super().runAverageSubtomograms()
        except Exception as e:
            raise Exception('Some of the previous protocols failed --> \n%s' % e)

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

    def test_genInitialModel(self):
        initModel = self.runGenInitialModel()
        super().checkAverage(initModel,
                             expectedSRate=self.binnedSRate,
                             expectedBoxSize=self.binnedBoxSize,
                             hasHalves=False)


class TestEmanTomoSubtomogramRefinementStaClassic(TestEmantomoStaClassicBase):
    mask = None
    avgSubtomo = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        try:
            cls.mask = super().runCreate3dMask()  # Generate a mask
            super().runPreviousProtocols()
            cls.avgSubtomo = super().runAverageSubtomograms()
        except Exception as e:
            raise Exception('Some of the previous protocols failed --> \n%s' % e)

    @classmethod
    def runSubtomogramRefinement(cls, inputRef, mask=None):
        print(magentaStr("\n==> Refining the subtomograms:"))
        inputDict = {'inputSetOfSubTomogram': cls.subtomosExtracted,
                     'inputRef': inputRef,
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

    def test_subtomoRefinement(self):
        protSubtomoRefinement = self.runSubtomogramRefinement(self.avgSubtomo)
        refinedSubtomos = getattr(protSubtomoRefinement, EmanTomoRefinementOutputs.subtomograms.name, None)
        self._checkRefinmentResults(refinedSubtomos)

    def test_subTomoRefinementWithMask(self):
        protSubtomoRefinement = self.runSubtomogramRefinement(self.avgSubtomo, mask=self.mask)
        refinedSubtomos = getattr(protSubtomoRefinement, EmanTomoRefinementOutputs.subtomograms.name, None)
        self._checkRefinmentResults(refinedSubtomos)

    def _checkRefinmentResults(self, outSubtomos):
        self.checkRefinedSubtomograms(self.subtomosExtracted, outSubtomos,
                                      expectedSetSize=self.nParticles,
                                      expectedBoxSize=self.binnedBoxSize,
                                      expectedSRate=self.binnedSRate,
                                      convention=TR_EMAN,
                                      orientedParticles=True)  # The coords imported were picked with PySeg


class TestEmanTemplateMatchingStaClassic(TestEmantomoStaClassicBase):
    avgSubtomo = None
    outNoParticles = 100

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        try:
            super().runPreviousProtocols()
            cls.avgSubtomo = super().runAverageSubtomograms()
        except Exception as e:
            raise Exception('Some of the previous protocols failed --> \n%s' % e)

    @classmethod
    def runTemplateMatching(cls):
        print(magentaStr("\n==> Running the template matching:"))
        protTempMatch = cls.newProtocol(EmanProtTemplateMatching,
                                        inputTomograms=cls.tomosBinned,
                                        refVol=cls.avgSubtomo,
                                        nptcl=cls.outNoParticles,
                                        delta=45,
                                        dthr=10)
        cls.launchProtocol(protTempMatch)
        return getattr(protTempMatch, protTempMatch._possibleOutputs.coordinates.name, None)

    def test_template_matching(self):
        # The imported coordinates were referred to the tomo used for the picking, which was at bin 2
        extractedCoords = self.runTemplateMatching()
        self.checkCoordinates(extractedCoords,
                              expectedSetSize=self.outNoParticles,
                              expectedSRate=self.binnedSRate,
                              expectedBoxSize=self.binnedBoxSize,
                              orientedParticles=False)  # The template matching generates non-oriented coords

