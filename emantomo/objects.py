# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
import h5py
from pyworkflow.object import String
from tomo.objects import SetOfSubTomograms, SubTomogram


class EmanMetaData:

    def __init__(self, tsId=None, inTomo=None, tomoHdfName=None, ts=None, tsHdfName=None, ctf=None, coords=None,
                 particles=None, jsonFile=None, processingInd=None):

        self.tsId = tsId
        self.inTomo = inTomo
        self.tomoHdfName = tomoHdfName
        self.ts = ts
        self.tsHdfName = tsHdfName
        self.ctf = ctf
        self.coords = coords
        self.particles = particles
        self.jsonFile = jsonFile
        self.processingInd = processingInd


class EmanHdf5Handler:

    def __init__(self, fileName):
        self._imgObjList = None
        self.fileName = fileName

    @property
    def fileName(self):
        return self._fileName

    @fileName.setter
    def fileName(self, fileName):
        try:
            f = h5py.File(fileName, 'r')
            self._imgObjList = f['MDF']['images']
            self._fileName = fileName
        except FileNotFoundError as e:
            raise e

    def getProjsFrom2dStack(self):
        """Generate a list of elements in which each element is a list containing the following data:
         [tiltId, particleId, xCoord, yCoord]. In a 2d stack, each image corresponds to a particle cropped
         in a tilt image"""
        projsList = []
        for particleId in range(len(self._imgObjList)):
            img = self._imgObjList[str(particleId)]
            projsList.append([img.attrs['EMAN.tilt_id'][0], particleId] +
                             img.attrs['EMAN.ptcl_source_coord'][:-1].tolist())
        return projsList


class EmanParticle(SubTomogram):

    INFO_JSON = '_infoJson'
    TS_HDF = '_tsHdf'
    TOMO_HDF = '_tomoHdf'
    STACK_2D_HDF = '_stack2dHdf'
    STACK_3D_HDF = '_stack3dHdf'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._infoJson = None
        self._tsHdf = None
        self._tomoHdf = None
        self._stack2dHdf = None
        self._stack3dHdf = None

    def setInfoJson(self, val):
        self._infoJson = String(val)

    def setTsHdf(self, val):
        self._tsHdf = String(val)

    def setTomoHdf(self, val):
        self._tomoHdf = String(val)

    def setStack2dHdf(self, val):
        self._stack2dHdf = String(val)

    def setStack3dHdf(self, val):
        self._stack3dHdf = String(val)

    def getInfoJson(self):
        return self._infoJson.get()

    def getTsHdf(self):
        return self._tsHdf.get()

    def getTomoHdf(self):
        return self._tomoHdf.get()

    def getStack2dHdf(self):
        return self._stack2dHdf.get()

    def getStack3dHdf(self):
        return self._stack3dHdf.get()


class EmanSetOfParticles(SetOfSubTomograms):

    ITEM_TYPE = EmanParticle
    ALI_2D = '_ali2d'
    ALI_3D = '_ali3d'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._ali2d = String()
        self._ali3d = String()

    def setAli2d(self, val):
        self._ali2d.set(val)

    def setAli3d(self, val):
        self._ali3d.set(val)

    def getAli2d(self):
        return self._ali2d.get()

    def getAli3d(self):
        return self._ali3d.get()

    def copyInfo(self, other):
        super().copyInfo(other)
        moreAttrs2Copy = [self.ALI_2D, self.ALI_3D]
        for attr in moreAttrs2Copy:
            if hasattr(other, attr):
                self.copyAttributes(other, attr)
