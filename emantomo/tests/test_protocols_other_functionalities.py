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
from .test_eman_sta_classic_base import TestEmantomoStaClassicBase
from ..protocols import EmanProtTomoClip


class TestEmanTomoClipTomograms(TestEmantomoStaClassicBase):
    tomosImported = None
    newCenter = [645, 783, 251]
    newDimensions = [389, 335, 240]
    biggerXdim = 1200

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tomosImported = super().runImportTomograms()
        x, y, z = cls.tomosImported.getDimensions()
        cls.origDimensions = [x, y, z]

    @classmethod
    def _runClipTomograms(cls, cx=None, cy=None, cz=None, dx=None, dy=None, dz=None, objectLabel=None):
        protClipTomograms = cls.newProtocol(EmanProtTomoClip,
                                            inputTomograms=cls.tomosImported,
                                            xc=cx,
                                            yc=cy,
                                            zc=cz,
                                            xDim=dx,
                                            yDim=dy,
                                            zDim=dz)
        if objectLabel:
            protClipTomograms.setObjLabel(objectLabel)
        cls.launchProtocol(protClipTomograms)
        outTomos = getattr(protClipTomograms, protClipTomograms._possibleOutputs.tomograms.name, None)
        cls.assertIsNotNone(outTomos, 'The tomograms were not generated.')
        return outTomos

    def test_clip_tomograms_new_dims_01(self):
        print(magentaStr("\n==> Clipping the tomograms with smaller new dimensions:"))
        outTomos = self._runClipTomograms(dx=self.newDimensions[0],
                                          dy=self.newDimensions[1],
                                          objectLabel='smaller xDim and yDim')
        self.checkTomograms(outTomos,
                            expectedSetSize=1,
                            expectedSRate=self.origSRate,  # This protocol does not change the sampling rate
                            expectedDimensions=[self.newDimensions[0], self.newDimensions[1], self.origDimensions[2]],
                            expectedOriginShifts=[-self.origDimensions[0] / 2 * self.origSRate,
                                                  -self.origDimensions[1] / 2 * self.origSRate,
                                                  -self.origDimensions[2] / 2 * self.origSRate])

    def test_clip_tomograms_new_dims_02(self):
        # TODO--> review whis one, I think the expected ortigin should be updated
        print(magentaStr("\n==> Clipping the tomograms with bigger new dimensions:"))
        outTomos = self._runClipTomograms(dx=self.biggerXdim,
                                          objectLabel='Bigger xDim')
        self.checkTomograms(outTomos,
                            expectedSetSize=1,
                            expectedSRate=self.origSRate,  # This protocol does not change the sampling rate
                            expectedDimensions=[self.biggerXdim, self.origDimensions[1], self.origDimensions[2]],
                            expectedOriginShifts=[-self.origDimensions[0] / 2 * self.origSRate,
                                                  -self.origDimensions[1] / 2 * self.origSRate,
                                                  -self.origDimensions[2] / 2 * self.origSRate])

    def test_tomograms_recenter(self):
        print(magentaStr("\n==> Re-centering the tomograms:"))
        outTomos = self._runClipTomograms(cx=self.newCenter[0],
                                          cy=self.newCenter[1],
                                          cz=self.newCenter[2],
                                          objectLabel='Re-centered only')
        self.checkTomograms(outTomos,
                            expectedSetSize=1,
                            expectedSRate=self.origSRate,  # This protocol does not change the sampling rate
                            expectedDimensions=[self.origDimensions[0], self.origDimensions[1], self.origDimensions[2]],
                            expectedOriginShifts=[-self.newCenter[0] * self.origSRate,
                                                  -self.newCenter[1] * self.origSRate,
                                                  -self.newCenter[2] * self.origSRate])

    def test_clip_tomograms_and_recenter(self):
        print(magentaStr("\n==> Clipping and re-centering the tomograms:"))
        outTomos = self._runClipTomograms(cx=self.newCenter[0],
                                          cy=self.newCenter[1],
                                          cz=self.newCenter[2],
                                          dx=self.newDimensions[0],
                                          dy=self.newDimensions[1],
                                          dz=self.newDimensions[2],
                                          objectLabel='Clipped and Re-centered')
        self.checkTomograms(outTomos,
                            expectedSetSize=1,
                            expectedSRate=self.origSRate,  # This protocol does not change the sampling rate
                            expectedDimensions=[self.newDimensions[0], self.newDimensions[1], self.newDimensions[2]],
                            expectedOriginShifts=[-self.newCenter[0] * self.origSRate,
                                                  -self.newCenter[1] * self.origSRate,
                                                  -self.newCenter[2] * self.origSRate])

