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
from pyworkflow.tests import  setupTestProject, DataSet
import pyworkflow.utils as pwutils
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from ..constants import EMAN_COVERAGE, EMAN_SCORE
from ..protocols import *
import tomo.protocols
from tomo.constants import TR_EMAN

from ..protocols.protocol_tomo_extraction_from_tomo import OutputExtraction
from ..protocols.protocol_tomo_subtomogram_refinement import EmanTomoRefinementOutputs



class TestEmanTomoSubtomogramRefinement(TestEmanTomoBase):
    """This class check if the protocol Subtomogram refinement works properly.
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

        coords = pwutils.removeBaseExt(self.coords3D)
        coords = protImportTomogram._getExtraPath(coords + '.txt')
        pwutils.createAbsLink(self.coords3D_Large, coords)
        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_EMAN,
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
                                              tomoSource=0,
                                              doInvert=doInvert,
                                              doNormalize=doNormalize,
                                              boxSize=boxSize)

        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(getattr(protTomoExtraction, OutputExtraction.subtomograms.name),
                             "There was a problem with SetOfSubtomogram output")

        return protTomoExtraction

    def _performFinalValidation(self, protTomoSubtomogramRefinement):
        outputSetOfSubTomograms = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.subtomograms.name)
        averageSubTomogram = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name)

        self.assertTrue(os.path.exists(averageSubTomogram.getFileName()))

        # Check average is a mrc file
        self.assertTrue(averageSubTomogram.getFileName().endswith(".mrc") )

        self.assertEqual(averageSubTomogram.getSamplingRate(), 5.0)

        self.assertEqual(outputSetOfSubTomograms.getDimensions(), (32, 32, 32))
        self.assertEqual(outputSetOfSubTomograms.getSize(), 15)
        self.assertEqual(outputSetOfSubTomograms.getCoordinates3D().getSize(), 15)

        for subTomogram in outputSetOfSubTomograms:
            self.assertEqual(subTomogram.getSamplingRate(), 5)
            self.assertTrue(hasattr(subTomogram, EMAN_COVERAGE))
            self.assertTrue(hasattr(subTomogram, EMAN_SCORE))
            matrix = subTomogram.getTransform(convention=TR_EMAN).getMatrix()
            self.assertEqual(matrix.shape, (4, 4))

        # FSCs
        fscs = getattr(protTomoSubtomogramRefinement, EmanTomoRefinementOutputs.FSCs.name)
        self.assertSetSize(fscs, 3, msg="FSCs not registered properly")

    def _runTomoSubtomogramRefinementWithSubtomo(self, niter=2, mass=500.0, threads=1, pkeep=1, goldstandard=-1,
                                                 goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0):
        protTomoExtraction = self._runPreviousProtocols()
        protTomoRefinement = self.newProtocol(EmanProtTomoRefinement,
                                              inputSetOfSubTomogram=getattr(protTomoExtraction, OutputExtraction.subtomograms.name),
                                              inputRef=protTomoExtraction,
                                              niter=niter,
                                              mass=mass,
                                              threads=threads,
                                              pkeep=pkeep,
                                              goldstandard=goldstandard,
                                              goldcontinue=goldcontinue,
                                              sym=sym,
                                              localfilter=localfilter,
                                              maxtilt=maxtilt)
        protTomoRefinement.inputRef.setExtended("%s.1" % OutputExtraction.subtomograms.name)

        self.launchProtocol(protTomoRefinement)

        self.assertIsNotNone(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name),
                             "There was a problem with subTomograms output")
        self.assertSetSize(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomograms.name), None,
                             "There was a problem with particles output")

        return protTomoRefinement

    def _runTomoSubtomogramRefinementWithVolume(self, niter=2, mass=500.0, threads=1, pkeep=1, goldstandard=-1,
                                                goldcontinue=False, sym="c1", localfilter=False, maxtilt=90.0):
        protTomoExtraction = self._runPreviousProtocols()

        particles = getattr(protTomoExtraction, OutputExtraction.subtomograms.name)

        self.dataset = DataSet.getDataSet('eman')
        self.vol = self.dataset.getFile('volume')
        self.protImportVol = self.runImportVolumes(self.vol, 3.5)

        self.assertIsNotNone(self.protImportVol.outputVolume,
                             "There was a problem with SetOfSubtomogram output")

        protTomoRefinement = self.newProtocol(EmanProtTomoRefinement,
                                              inputSetOfSubTomogram=particles,
                                              inputRef=self.protImportVol.outputVolume,
                                              niter=niter,
                                              mass=mass,
                                              threads=threads,
                                              pkeep=pkeep,
                                              goldstandard=goldstandard,
                                              goldcontinue=goldcontinue,
                                              sym=sym,
                                              localfilter=localfilter,
                                              maxtilt=maxtilt)

        self.launchProtocol(protTomoRefinement)

        average = getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomogramAverage.name)
        self.assertIsNotNone(average,  "There was a problem with average output")
        self.assertTrue(os.path.exists(average.getFileName()), "Average %s does not exists" % average.getFileName())

        self.assertTrue(average.hasHalfMaps(), "Halves not registered")

        half1, half2 = average.getHalfMaps().split(',')
        self.assertTrue(os.path.exists(half1), msg="Average 1st half %s does not exists" % half1)
        self.assertTrue(os.path.exists(half2), msg="Average 2nd half %s does not exists" % half2)

        self.assertSetSize(getattr(protTomoRefinement, EmanTomoRefinementOutputs.subtomograms.name), None,
                             "There was a problem with particles output")

        return protTomoRefinement

    def test_defaultSubTomogramRefinementWithSubTomo(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinementWithSubtomo()
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement

    def test_defaultSubTomogramRefinementWithVolume(self):
        protTomoSubtomogramRefinement = self._runTomoSubtomogramRefinementWithVolume()
        self._performFinalValidation(protTomoSubtomogramRefinement)

        return protTomoSubtomogramRefinement
