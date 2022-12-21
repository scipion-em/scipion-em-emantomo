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

from tomo.constants import TR_EMAN
from .test_eman_base import TestEmantomoBase
from ..protocols.protocol_tomo_extraction_from_tomo import OTHER


class TestEmanTomoExtraction(TestEmantomoBase):
    """This class check if the protocol to extract subtomograms
    in Eman works properly.
    """

    tomosImported = None
    tomosBinned = None
    coordsImported = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tomosImported, cls.coordsImported, cls.tomosBinned = cls._runPreviousProtocols()
        cls.coordExtremeValsNoBin = super().getMinAndMaxCoordValuesFromSet(cls.coordsImported,
                                                                           cls.tomosImported.getFirstItem())

    @classmethod
    def _runPreviousProtocols(cls):
        tomoImported = super().runImportTomograms()  # Import tomograms
        coordsImported = super().runImport3dCoords(tomoImported)  # Import coordinates
        tomosBinned = super().runBinTomograms(tomoImported)  # Bin the tomogram to make it smaller
        return tomoImported, coordsImported, tomosBinned

    def test_extractParticlesSameAsPicking(self):
        subtomosExtracted = super().runExtractSubtomograms(self.coordsImported, boxSize=super().boxSize)
        self.checkExtractedSubtomos(subtomosExtracted,
                                    setSize=super().nParticles,
                                    sRate=super().origSRate,
                                    boxSize=super().boxSize,
                                    tomo4extraction=self.tomosImported.getFirstItem())

    def test_extractParticlesOtherTomoSource(self):
        subtomosExtracted = super().runExtractSubtomograms(self.coordsImported,
                                                           tomoSource=OTHER,
                                                           tomograms=self.tomosBinned,
                                                           boxSize=super().boxSize)

        self.checkExtractedSubtomos(subtomosExtracted,
                                    setSize=super().nParticles,
                                    sRate=super().origSRate * 2,
                                    boxSize=super().boxSize,
                                    tomo4extraction=self.tomosBinned.getFirstItem())

    # def test_extractParticlesWithDoInvert(self):
    #     protTomoExtraction = self._runTomoExtraction(doInvert=True)
    #     output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
    #     self.checkExtractedSubtomos(output)
    #
    # def test_extractParticlesWithDoNormalize(self):
    #     protTomoExtraction = self._runTomoExtraction(doNormalize=True)
    #     output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
    #     self.checkExtractedSubtomos(output)
    #
    # def test_extractParticlesModifiedDownFactor(self):
    #     protTomoExtraction = self._runTomoExtraction(downFactor=2)
    #     output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
    #     self.checkExtractedSubtomos(output)
    #
    # def test_extractParticlesModifiedBoxSize(self):
    #     protTomoExtraction = self._runTomoExtraction(boxSize=64)
    #     output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
    #     self.checkExtractedSubtomos(output)
    #
    # def test_extractParticlesWithAllOptions(self):
    #     protTomoExtraction = self._runTomoExtraction(boxSize=64, downFactor=2, doNormalize=True, doInvert=True)
    #     output = getattr(protTomoExtraction, OutputExtraction.subtomograms.name, None)
    #     self.checkExtractedSubtomos(output)

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
        currentCoordsExtremes = super().getMinAndMaxCoordValuesFromSet(subtomograms, tomo4extraction)
        unbinnedCoordsExtremes = self.coordExtremeValsNoBin
        self.assertTrue(np.array_equal(currentCoordsExtremes, unbinnedCoordsExtremes))

