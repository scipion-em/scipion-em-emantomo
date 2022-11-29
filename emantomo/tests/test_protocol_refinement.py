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

from emantomo.tests.test_protocols_emantomo import TestEmanTomoBase
from imod.protocols import ProtImodTomoNormalization
from pyworkflow.tests import setupTestProject, DataSet
import pyworkflow.utils as pwutils
from pyworkflow.utils import magentaStr
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from ..constants import EMAN_COVERAGE, EMAN_SCORE
from ..protocols import *
import tomo.protocols
from tomo.constants import TR_EMAN
from ..protocols.protocol_average_subtomos import OutputsAverageSubtomos
from ..protocols.protocol_tomo_extraction_from_tomo import OutputExtraction, OTHER
from ..protocols.protocol_tomo_subtomogram_refinement import EmanTomoRefinementOutputs


class TestEmanTomoSubtomogramRefinement(TestEmanTomoBase):
    """This class check if the protocol Subtomogram refinement works properly.
    """

    boxSize = 12

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()
        cls.protTomoExtraction = cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(tomo.protocols.ProtImportTomograms,
                                             filesPath='/ExtraDisk/extraDataSets/overview_wbp.mrc',
                                             samplingRate=super().origSRate)

        cls.launchProtocol(protImportTomogram)
        cls.assertIsNotNone(protImportTomogram.outputTomograms, "There was a problem with tomogram output")

        # Import coordinates
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        coords = pwutils.removeBaseExt(super().coords3D)
        coords = protImportTomogram._getExtraPath(coords + '.txt')
        pwutils.createAbsLink(super().coords3D_Large, coords)
        protImportCoordinates3d = cls.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                  auto=IMPORT_FROM_EMAN,
                                                  filesPath=coords,
                                                  importTomograms=protImportTomogram.outputTomograms,
                                                  filesPattern='',
                                                  boxSize=cls.boxSize,
                                                  samplingRate=super().origSRate)

        cls.launchProtocol(protImportCoordinates3d)
        cls.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                            "There was a problem with the 3D coordinates output")

        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram normalization:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=protImportTomogram.outputTomograms,
                                                binning=4)

        cls.launchProtocol(protTomoNormalization)
        outputTomosBinned = getattr(protTomoNormalization, 'Tomograms', None)
        cls.assertIsNotNone(outputTomosBinned, 'No tomograms were genetated in tomo normalization.')

        # Extract subtomograms
        print(magentaStr("\n==> Extracting the subtomograms:"))
        doInvert = False
        doNormalize = False
        protTomoExtraction = cls.newProtocol(EmanProtTomoExtraction,
                                             inputTomograms=outputTomosBinned,
                                             inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                             tomoSource=OTHER,
                                             doInvert=doInvert,
                                             doNormalize=doNormalize,
                                             boxSize=cls.boxSize)

        cls.launchProtocol(protTomoExtraction)
        cls.assertIsNotNone(getattr(protTomoExtraction, OutputExtraction.subtomograms.name),
                            "There was a problem with the subtomograms output")

        return protTomoExtraction

    def _performFinalValidation(self, protTomoSubtomogramRefinement):
        outputSetOfSubTomograms = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.subtomograms.name)
        averageSubTomogram = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name)

        self.assertTrue(os.path.exists(averageSubTomogram.getFileName()))

        # Check average is a mrc file
        self.assertTrue(averageSubTomogram.getFileName().endswith(".mrc"))

        self.assertEqual(averageSubTomogram.getSamplingRate(), self.origSRate)

        self.assertEqual(outputSetOfSubTomograms.getDimensions(), (self.boxSize, self.boxSize, self.boxSize))
        self.assertEqual(outputSetOfSubTomograms.getSize(), 15)
        self.assertEqual(outputSetOfSubTomograms.getCoordinates3D().getSize(), 15)

        for subTomogram in outputSetOfSubTomograms:
            self.assertEqual(subTomogram.getSamplingRate(), self.origSRate)
            self.assertTrue(hasattr(subTomogram, EMAN_COVERAGE))
            self.assertTrue(hasattr(subTomogram, EMAN_SCORE))
            matrix = subTomogram.getTransform(convention=TR_EMAN).getMatrix()
            self.assertEqual(matrix.shape, (4, 4))

        # FSCs
        fscs = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.FSCs.name)
        self.assertSetSize(fscs, 3, msg="FSCs not registered properly")

    # def _runTomoSubtomogramRefinementWithSubtomo(self, niter=2, mass=500.0, threads=1, pkeep=1, goldstandard=-1,
    #                                              goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0):
    #     protTomoExtraction = self._runPreviousProtocols()
    #     protTomoRefinement = self.newProtocol(EmanProtTomoRefinement,
    #                                           inputSetOfSubTomogram=getattr(protTomoExtraction,
    #                                                                         OutputExtraction.subtomograms.name),
    #                                           inputRef=protTomoExtraction,
    #                                           niter=niter,
    #                                           mass=mass,
    #                                           threads=threads,
    #                                           pkeep=pkeep,
    #                                           goldstandard=goldstandard,
    #                                           goldcontinue=goldcontinue,
    #                                           sym=sym,
    #                                           localfilter=localfilter,
    #                                           maxtilt=maxtilt)
    #     protTomoRefinement.inputRef.setExtended("%s.1" % OutputExtraction.subtomograms.name)
    #
    #     self.launchProtocol(protTomoRefinement)
    #
    #     self.assertIsNotNone(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name),
    #                          "There was a problem with subTomograms output")
    #     self.assertSetSize(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomograms.name), None,
    #                        "There was a problem with particles output")
    #
    #     return protTomoRefinement

    def _runTomoSubtomogramRefinement(self, niter=2, mass=500.0, threads=1, pkeep=1, goldstandard=-1,
                                      goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0, refVol=None):

        print(magentaStr("\n==> Refining the subtomograms:"))
        particles = getattr(self.protTomoExtraction, OutputExtraction.subtomograms.name)
        inputDict = {'inputSetOfSubTomogram': particles,
                     'niter': niter,
                     'mass': mass,
                     'threads': threads,
                     'pkeep': pkeep,
                     'goldstandard': goldstandard,
                     'goldcontinue': goldcontinue,
                     'sym': sym,
                     'localfilter': localfilter,
                     'maxtilt': maxtilt}

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
        refVol = self._runAverageSubtomos()
        # Refine the subtomograms
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinement(refVol=refVol)
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement

    def _runAverageSubtomos(self):
        extractedSubtomos = getattr(self.protTomoExtraction, OutputExtraction.subtomograms.name, None)
        protSubtomoAvg = self.newProtocol(EmanProtSubTomoAverage, inputSetOfSubTomogram=extractedSubtomos)
        self.launchProtocol(protSubtomoAvg)
        outputAvg = getattr(protSubtomoAvg, OutputsAverageSubtomos.averageSubTomos.name, None)
        self.assertIsNotNone(outputAvg, 'There was a problem when calculating the avg subtomogram')
        return outputAvg
