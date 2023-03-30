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
from pyworkflow.utils import magentaStr
from tomo.constants import TR_EMAN
from tomo.protocols import ProtTomoExtractCoords
from tomo.protocols.protocol_extract_coordinates import Output3dCoordExtraction
from .test_eman_sta_classic_base import TestEmantomoStaClassicBase
from ..protocols.protocol_tomo_extraction_from_tomo import OTHER


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

        cls.subtomosFromSubtomos = super().runExtractSubtomograms(cls.subtomosSameAsPicking, boxSize=super().boxSize / 2, label="Extraction from subtomograms")
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

    def checkExtractedSubtomos(self, subtomograms, setSize=-1, sRate=-1, boxSize=-1, tomo4extraction=None):
        scaleFactor = subtomograms.getSamplingRate() / super().origSRate
        # Check the critical properties of the set
        self.assertSetSize(subtomograms, setSize)
        self.assertEqual(subtomograms.getSamplingRate(), sRate)
        self.assertEqual(subtomograms.getDimensions(), (boxSize, boxSize, boxSize))
        self.assertTrue(subtomograms.hasCoordinates3D())
        # Check the subtomograms that compose the set
        for subtomo in subtomograms:
            subtomoTr = subtomo.getTransform(convention=TR_EMAN)
            subtomoMatrix = subtomoTr.getMatrix()
            coordinate = subtomo.getCoordinate3D()
            coordTr = coordinate._eulerMatrix
            coordMatrix = coordinate.getMatrix(convention=TR_EMAN)
            self.assertTrue(exists(subtomo.getFileName()))
            self.assertEqual(subtomo.getSamplingRate(), sRate)
            # The shifts in the subtomograms transformation matrix should have been scaled properly
            super().checkShiftsScaling(coordTr, subtomoTr, scaleFactor)
            # Imported coordinates were picked using PySeg, so they must have an orientation
            super().check3dTransformMatrix(subtomoMatrix)
            # Also, at this point the transformation matrix should be the same as the coordinate matrix as the angles
            # have not been refined yet
            super().check3dTransformMatrix(coordMatrix)
            self.assertTrue(np.array_equal(subtomoMatrix, coordMatrix))
            # Check the tomoId
            self.assertEqual(coordinate.getTomoId(), tomo4extraction.getTsId())

        # Check that the coordinates remain the same (the scaling is only applied to the shifts of the
        # transformation matrix, while the coordinates are only scaled in the coordinates extraction protocol
        # from the plugin scipion-em-tomo
        currentCoordsExtremes = super().getMinAndMaxCoordValuesFromSet(subtomograms)
        unbinnedCoordsExtremes = self.coordExtremeValsBin2
        self.assertTrue(np.array_equal(currentCoordsExtremes, unbinnedCoordsExtremes))

    def checkExtracted3dCoordinates(self, extractedCoords, coordScaleFactor=1.0, shiftsScaleFactor=1.0,
                                    tomoId=None, tomoSRate=-1):
        binned2CoordsExtremes = self.coordExtremeValsBin2
        currentCoordsExtremes = super().getMinAndMaxCoordValuesFromSet(extractedCoords)
        # Check the coordinate extremes
        self.assertTrue(np.array_equal(currentCoordsExtremes, coordScaleFactor * binned2CoordsExtremes))
        # Check the sampling rate
        self.assertEqual(extractedCoords.getSamplingRate(), tomoSRate)
        # Check the set size
        self.assertSetSize(extractedCoords, super().nParticles)
        # Check the box size
        self.assertEqual(extractedCoords.getBoxSize(), super().boxSize)
        # Other checks per coordinate
        for inSubtomos, outCoord in zip(self.subtomosSameAsPicking, extractedCoords):
            # Check the transformation matrices and shifts
            subtomoTr = inSubtomos.getTransform(convention=TR_EMAN)
            coordTr = outCoord._eulerMatrix
            super().check3dTransformMatrix(outCoord.getMatrix(convention=TR_EMAN))
            super().checkShiftsScaling(subtomoTr, coordTr, shiftsScaleFactor)
            # Check the tomoId
            self.assertEqual(outCoord.getTomoId(), tomoId)

    def test_extractParticlesSameAsPicking(self):
        # The imported 3d coordinates were picked from the binned tomogram
        self.checkExtractedSubtomos(self.subtomosSameAsPicking,
                                    setSize=super().nParticles,
                                    sRate=super().origSRate * 2,
                                    boxSize=super().boxSize / 2,
                                    tomo4extraction=self.tomosImported.getFirstItem())

    def test_extractParticlesSameAsPickingSubtomos(self):
        # The imported 3d coordinates were picked from the binned tomogram
        self.checkExtractedSubtomos(self.subtomosFromSubtomos,
                                    setSize=super().nParticles,
                                    sRate=super().origSRate * 2,
                                    boxSize=super().boxSize / 2,
                                    tomo4extraction=self.tomosImported.getFirstItem())

    def test_extractParticlesOtherTomoSource(self):
        self.checkExtractedSubtomos(self.subtomosAnotherTomo,
                                    setSize=super().nParticles,
                                    sRate=super().origSRate,
                                    boxSize=super().boxSize,
                                    tomo4extraction=self.tomosBinned.getFirstItem())

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
                                                  boxSize=super().boxSize)

        self.checkExtracted3dCoordinates(extractedCoords,
                                         coordScaleFactor=2,
                                         shiftsScaleFactor=2,
                                         tomoId=self.tomosImported.getFirstItem().getTsId(),
                                         tomoSRate=super().origSRate)

    def test_extract3dCoordsToSmallerTomo(self):
        """Subtomos extracted from the another tomo source, which was unbinned. Coordinates will
        be extracted to bin 2."""
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosAnotherTomo,
                                                  inputTomos=self.tomosBinned,
                                                  boxSize=super().boxSize)

        self.checkExtracted3dCoordinates(extractedCoords,
                                         coordScaleFactor=1,  # They'll remain the same
                                         shiftsScaleFactor=0.5,
                                         tomoId=self.tomosImported.getFirstItem().getTsId(),
                                         tomoSRate=super().origSRate * 2)

    def test_extract3dCoordsToTheSameTomo(self):
        """Subtomos extracted from the same tomo used for the picking, which was at bin 2. Coordinates will
               be extracted to the same tomogram."""
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosSameAsPicking,
                                                  inputTomos=self.tomosBinned,
                                                  boxSize=super().boxSize)

        self.checkExtracted3dCoordinates(extractedCoords,
                                         coordScaleFactor=1,
                                         shiftsScaleFactor=1,
                                         tomoId=self.tomosImported.getFirstItem().getTsId(),
                                         tomoSRate=super().origSRate * 2)








