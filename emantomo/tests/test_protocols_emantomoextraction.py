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
from emantomo.tests.test_protocols_emantomo import TestEmanTomoBase
from pyworkflow.tests import (setupTestProject)
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from ..protocols import *
import tomo.protocols

from ..protocols.protocol_tomo_extraction_from_tomo import OutputExtraction


class TestEmanTomoExtraction(TestEmanTomoBase):
    """This class check if the protocol to extract subtomograms
    in Eman works properly.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestEmanTomoBase.setData()

    def _runTomoExtraction(self, tomoSource=0, doInvert=False, doNormalize=False, boxSize=32, downFactor=1):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertSetSize(protImportTomogram.outputTomograms, 1,
                           "There was a problem with tomogram output")
        self.assertSetSize(protImportCoordinates3d.outputCoordinates, 5,
                           "There was a problem with coordinates 3d output")

        if tomoSource == 1:
            protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                                  inputTomograms=protImportTomogram.outputTomograms,
                                                  inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                                  tomoSource=tomoSource,
                                                  doInvert=doInvert,
                                                  doNormalize=doNormalize,
                                                  boxSize=boxSize,
                                                  downFactor=downFactor)
        else:
            protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                                  inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                                  tomoSource=tomoSource,
                                                  doInvert=doInvert,
                                                  doNormalize=doNormalize,
                                                  boxSize=boxSize,
                                                  downFactor=downFactor)
        self.launchProtocol(protTomoExtraction)
        self.assertSetSize(getattr(protTomoExtraction, OutputExtraction.subtomograms.name), 5,
                           "There was a problem with SetOfSubtomogram output")
        return protTomoExtraction

    def test_extractParticlesSameAsPicking(self):
        protTomoExtraction = self._runTomoExtraction()
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def test_extractParticlesOther(self):
        protTomoExtraction = self._runTomoExtraction(tomoSource=1)
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def test_extractParticlesWithDoInvert(self):
        protTomoExtraction = self._runTomoExtraction(doInvert=True)
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def test_extractParticlesWithDoNormalize(self):
        protTomoExtraction = self._runTomoExtraction(doNormalize=True)
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def test_extractParticlesModifiedDownFactor(self):
        protTomoExtraction = self._runTomoExtraction(downFactor=2)
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def test_extractParticlesModifiedBoxSize(self):
        protTomoExtraction = self._runTomoExtraction(boxSize=64)
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def test_extractParticlesWithAllOptions(self):
        protTomoExtraction = self._runTomoExtraction(boxSize=64, downFactor=2, doNormalize=True, doInvert=True)
        output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
        self.assessOutput(output)

    def assessOutput(self, outputSet, size=5):

        self.assertSetSize(outputSet, size)
        self.assertTrue(outputSet.hasCoordinates3D())

