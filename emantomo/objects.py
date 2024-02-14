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
from pyworkflow.object import String, Float, Integer
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

    def getProjsFrom2dStack(self, shrink=1):
        """Generate a list of elements in which each element is a list containing the following data:
         [tiltId, xCoord, yCoord]. In a 2d stack, each image corresponds to a particle cropped
         in a tilt image"""
        projsList = []
        for particleId in range(len(self._imgObjList)):
            img = self._imgObjList[str(particleId)]
            projsList.append([img.attrs['EMAN.tilt_id'][0]] +
                             (shrink * img.attrs['EMAN.ptcl_source_coord'][:-1]).tolist())
        return projsList

    def getSamplingRate(self):
        """Reads the sampling rate from the HDF header"""
        return self._imgObjList['0'].attrs['EMAN.apix_x']


class EmanParticle(SubTomogram):
    INFO_JSON = '_infoJson'
    TS_HDF = '_tsHdf'
    TOMO_HDF = '_tomoHdf'
    STACK_2D_HDF = '_stack2dHdf'
    STACK_3D_HDF = '_stack3dHdf'
    EMAN_SCORE = '_emanScore'
    ABS_INDEX = '_absIndex'  # Absolute index within the whole particles from all the tomograms (required for the 2d/3d

    # particle matching and sub-setting)

    def __init__(self, infoJson=None, tsHdf=None, tomoHdf=None, stack2dHdf=None, stack3dHdf=None,
                 emanScore=0, absIndex=None, **kwargs):
        super().__init__(**kwargs)
        self._infoJson = String(infoJson)
        self._tsHdf = String(tsHdf)
        self._tomoHdf = String(tomoHdf)
        self._stack2dHdf = String(stack2dHdf)
        self._stack3dHdf = String(stack3dHdf)
        self._emanScore = Float(emanScore)
        self._absIndex = Integer(absIndex)

    def setInfoJson(self, val):
        self._infoJson.set(val)

    def setTsHdf(self, val):
        self._tsHdf.set(val)

    def setTomoHdf(self, val):
        self._tomoHdf.set(val)

    def setStack2dHdf(self, val):
        self._stack2dHdf.set(val)

    def setStack3dHdf(self, val):
        self._stack3dHdf.set(val)

    def setEmanScore(self, val):
        self._emanScore.set(val)

    def setAbsIndex(self, val):
        self._absIndex.set(val)

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

    def getEmanScore(self):
        return self._emanScore.get()

    def getAbsIndex(self):
        return self._absIndex.get()


class EmanSetOfParticles(SetOfSubTomograms):
    ITEM_TYPE = EmanParticle
    ALI_2D = '_emanAli2dFile'
    ALI_3D = '_emanAli3dFile'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._emanAli2dFile = String()
        self._emanAli3dFile = String()

    def setAli2dLstFile(self, val):
        self._emanAli2dFile.set(val)

    def setAli3dLstFile(self, val):
        self._emanAli3dFile.set(val)

    def getAli2dLstFile(self):
        return self._emanAli2dFile.get()

    def getAli3dLstFile(self):
        return self._emanAli3dFile.get()

    def copyInfo(self, other):
        super().copyInfo(other)
        moreAttrs2Copy = [self.ALI_2D, self.ALI_3D]
        for attr in moreAttrs2Copy:
            if hasattr(other, attr):
                self.copyAttributes(other, attr)
