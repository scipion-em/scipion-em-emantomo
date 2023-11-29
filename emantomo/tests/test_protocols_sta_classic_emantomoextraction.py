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
from pyworkflow.utils import magentaStr
from tomo.constants import TR_EMAN
from tomo.protocols import ProtTomoExtractCoords
from tomo.protocols.protocol_extract_coordinates import Output3dCoordExtraction
from .test_eman_sta_classic_base import TestEmantomoStaClassicBase
from ..protocols.protocol_extraction_from_tomo import OTHER


class TestEmanTomoExtractionStaClassic(TestEmantomoStaClassicBase):
    """This class check if the protocol to extract subtomograms
    in Eman works properly.
    """

    tomosImported = None
    tomosBinned = None
    coordsImported = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tomosImported, cls.coordsImported, cls.tomosBinned = cls.runPreviousProtocols()
        cls.coordExtremeValsBin2 = super().getMinAndMaxCoordValuesFromSet(cls.coordsImported)
        cls.subtomosSameAsPicking = super().runExtractSubtomograms(cls.coordsImported, boxSize=super().boxSize / 2)
        cls.subtomosAnotherTomo = super().runExtractSubtomograms(cls.coordsImported,
                                                                 tomoSource=OTHER,
                                                                 tomograms=cls.tomosImported,
                                                                 boxSize=super().boxSize)

        cls.subtomosFromSubtomos = super().runExtractSubtomograms(cls.subtomosSameAsPicking,
                                                                  boxSize=super().boxSize / 2,
                                                                  label="Extraction from subtomograms")

    @classmethod
    def runPreviousProtocols(cls):
        tomoImported = super().runImportTomograms()  # Import tomograms
        tomosBinned = super().runBinTomograms(tomoImported)  # Bin the tomogram to make it smaller
        coordsImported = super().runImport3dCoords(tomosBinned)  # Import coordinates

        return tomoImported, coordsImported, tomosBinned

    @classmethod
    def runExtract3dCoords(cls, inputSubTomos=None, inputTomos=None, boxSize=None):
        print(magentaStr("\n==> Extracting the 3D coordinates:"))
        protExtract3dCoords = cls.newProtocol(ProtTomoExtractCoords,
                                              inputSubTomos=inputSubTomos,
                                              inputTomos=inputTomos,
                                              boxSize=boxSize)

        cls.launchProtocol(protExtract3dCoords)
        coordsExtracted = getattr(protExtract3dCoords, Output3dCoordExtraction.coordinates3d.name, None)
        cls.assertIsNotNone(coordsExtracted, "There was a problem with the 3d coordinates extraction")
        return coordsExtracted

    def test_extractParticlesSameAsPicking(self):
        # The imported 3d coordinates were picked from the binned tomogram
        self.checkExtractedSubtomos(self.coordsImported, self.subtomosSameAsPicking,
                                    expectedSetSize=self.nParticles,
                                    expectedSRate=self.binnedSRate,
                                    expectedBoxSize=self.binnedBoxSize,
                                    convention=TR_EMAN,
                                    orientedParticles=True)  # The coords imported were picked with PySeg

    def test_extractParticlesSameAsPickingSubtomos(self):
        # The imported 3d coordinates were picked from the binned tomogram
        self.checkExtractedSubtomos(self.coordsImported, self.subtomosFromSubtomos,
                                    expectedSetSize=self.nParticles,
                                    expectedSRate=self.binnedSRate,
                                    expectedBoxSize=self.binnedBoxSize,
                                    convention=TR_EMAN,
                                    orientedParticles=True)  # The coords imported were picked with PySeg

    def test_extractParticlesOtherTomoSource(self):
        self.checkExtractedSubtomos(self.coordsImported, self.subtomosAnotherTomo,
                                    expectedSetSize=self.nParticles,
                                    expectedSRate=self.origSRate,
                                    expectedBoxSize=self.boxSize,
                                    convention=TR_EMAN,
                                    orientedParticles=True)  # The coords imported were picked with PySeg

    # __________________________________________________________________________________________________________________
    # NOTE:
    # Although the coordinates extraction is not a part of the plugin emantomo, a part of its functionality
    # will be tested here as it has a very direct relation with the particle extraction. In the coordinates
    # extraction a new set of coordinates is generated, so in this case not only the shifts but the coordinates
    # are expected to be scaled properly according to the sampling rate of the introduced tomograms.
    # __________________________________________________________________________________________________________________

    def test_extract3dCoordsToBiggerTomo(self):
        """Subtomos extracted from the same tomo used for the picking, which was at bin 2. Coordinates will
        be extracted to the original size (unbinned)."""
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosSameAsPicking,
                                                  inputTomos=self.tomosImported,
                                                  boxSize=self.boxSize)

        self.checkExtracted3dCoordinates(self.coordsImported, extractedCoords,
                                         expectedSetSize=self.nParticles,
                                         expectedSRate=self.origSRate,
                                         expectedBoxSize=self.boxSize,
                                         convention=TR_EMAN,
                                         orientedParticles=True)  # The coords imported were picked with PySeg

    def test_extract3dCoordsToSmallerTomo(self):
        """Subtomos extracted from the another tomo source, which was unbinned. Coordinates will
        be extracted to bin 2."""
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosAnotherTomo,
                                                  inputTomos=self.tomosBinned,
                                                  boxSize=self.binnedBoxSize)

        self.checkExtracted3dCoordinates(self.coordsImported, extractedCoords,
                                         expectedSetSize=self.nParticles,
                                         expectedSRate=self.binnedSRate,
                                         expectedBoxSize=self.binnedBoxSize,
                                         convention=TR_EMAN,
                                         orientedParticles=True)  # The coords imported were picked with PySeg

    def test_extract3dCoordsToTheSameTomo(self):
        """Subtomos extracted from the same tomo used for the picking, which was at bin 2. Coordinates will
               be extracted to the same tomogram."""
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosSameAsPicking,
                                                  inputTomos=self.tomosBinned,
                                                  boxSize=self.binnedBoxSize)

        self.checkExtracted3dCoordinates(self.coordsImported, extractedCoords,
                                         expectedSetSize=self.nParticles,
                                         expectedSRate=self.binnedSRate,
                                         expectedBoxSize=self.binnedBoxSize,
                                         convention=TR_EMAN,
                                         orientedParticles=True)  # The coords imported were picked with PySeg

