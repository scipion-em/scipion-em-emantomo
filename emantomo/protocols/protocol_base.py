# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]             
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import glob
from os.path import join, abspath, basename
from emantomo import Plugin
from emantomo.constants import INFO_DIR, TOMOGRAMS_DIR, TS_DIR, SETS_DIR, PARTICLES_DIR, PARTICLES_3D_DIR, \
    REFERENCE_NAME, TOMOBOX
from emantomo.objects import EmanParticle
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer
from pyworkflow.utils import makePath, createLink
from tomo.protocols import ProtTomoBase
from tomo.utils import getNonInterpolatedTsFromRelations

IN_TS = 'inputTS'
IN_COORDS = 'inputCoordinates'
IN_CTF = 'inputCTF'
IN_SUBTOMOS = 'inputSubtomos'
IN_BOXSIZE = 'boxSize'
IN_TOMOS = 'inputTomograms'
REF_VOL = 'refVol'


class ProtEmantomoBase(EMProtocol, ProtTomoBase):
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inSamplingRate = -1
        self.scaleFactor = 1
        self.voltage = 300
        self.sphAb = 2.7

    # --------------------------- STEPS functions -----------------------------
    def convertTsStep(self, mdObj):
        # The converted TS must be unbinned, because EMAN will read the sampling rate from its header. This is why the
        # the TS associated to the CTF is the one considered first. Later, when generating the json, the TS alignment
        # parameters are read from the introduced TS and the shifts are scaled to at the unbinned scale
        if mdObj.ctf:
            inTsFName = mdObj.ctf.getTiltSeries().getFirstItem().getFileName()
        else:
            inTsFName = mdObj.ts.getFirstItem().getFileName()
        dirName = TS_DIR
        sRate = mdObj.ts.getSamplingRate()
        outFile = self.convertOrLink(inTsFName, mdObj.tsId, dirName, sRate)
        # Store the tsHdfName in the current mdObj
        mdObj.tsHdfName = outFile

    def convertTomoStep(self, mdObj):
        inTomoFName = mdObj.inTomo.getFileName()
        dirName = TOMOGRAMS_DIR
        sRate = mdObj.inTomo.getSamplingRate()
        outFile = self.convertOrLink(inTomoFName, mdObj.tsId, dirName, sRate)
        # Store the tomoHdfName in the current mdObj
        mdObj.tomoHdfName = outFile

    def createEmanPrjPostExtractionStep(self):
        inSubtomos = getattr(self, IN_SUBTOMOS).get()
        self.inSamplingRate = inSubtomos.getSamplingRate()
        # Create project dir structure
        self.createInitEmanPrjDirs()
        infoDir = self.getInfoDir()
        tsDir = self.getTsDir()
        tomoDir = self.getTomogramsDir()
        stack2dDir = self.getStack2dDir()
        stack3dDir = self.getStack3dDir()
        makePath(self.getSetsDir(), stack2dDir, stack3dDir)
        # Get the unique values of the files to be linked
        dataDict = inSubtomos.getUniqueValues([EmanParticle.INFO_JSON,
                                               EmanParticle.TS_HDF,
                                               EmanParticle.TOMO_HDF,
                                               EmanParticle.STACK_2D_HDF,
                                               EmanParticle.STACK_3D_HDF])
        # Link the files
        for infoJson, tsFile, tomoFile, stack2d, stack3d in zip(dataDict[EmanParticle.INFO_JSON],
                                                                dataDict[EmanParticle.TS_HDF],
                                                                dataDict[EmanParticle.TOMO_HDF],
                                                                dataDict[EmanParticle.STACK_2D_HDF],
                                                                dataDict[EmanParticle.STACK_3D_HDF]):
            createLink(infoJson, join(infoDir, basename(infoJson)))
            createLink(tsFile, join(tsDir, basename(tsFile)))
            createLink(tomoFile, join(tomoDir, basename(tomoFile)))
            createLink(stack2d, join(stack2dDir, basename(stack2d)))
            createLink(stack3d, join(stack3dDir, basename(stack3d)))

    def convertRefVolStep(self):
        inRef = self.getRefVol()
        if inRef:
            self.convertOrLink(inRef.getFileName(), REFERENCE_NAME, '', inRef.getSamplingRate())

    def buildEmanSetsStep(self):
        program = Plugin.getProgram("e2spt_buildsets.py")
        args = '--allparticles '
        self.runJob(program, args, cwd=self._getExtraPath())

    # --------------------------- UTILS functions ----------------------------------
    def getObjByName(self, name):
        """Return an object, from a protocol, named 'name' instead of a pointer."""
        obj = getattr(self, name, None)
        if obj and type(obj) == Pointer:
            return obj.get()
        else:
            return obj

    def getTs(self):
        """If the user provides a set of tilt series, use them. If not (expected behaviour), get the non-interpolated
        tilt series from the introduced coordinates."""
        tsSet = getattr(self, IN_TS).get()
        return tsSet if tsSet else getNonInterpolatedTsFromRelations(getattr(self, IN_COORDS), self)

    def getBoxSize(self):
        """If the user provides a box size value. If not (expected behaviour), read it from the set of subtomograms
        introduced."""
        boxSizeFromForm = getattr(self, IN_BOXSIZE).get()
        return boxSizeFromForm if boxSizeFromForm else getattr(self, IN_SUBTOMOS).get().getCoordinates3D().getBoxSize()

    def getRefVol(self):
        refVol = getattr(self, REF_VOL, None)
        return refVol.get() if refVol else None

    def createInitEmanPrjDirs(self):
        """Create in the current protocol's extra path the initial directory structure of an EMAN tomo project: info,
        tiltseries and tomograms dirs"""
        makePath(self.getInfoDir(), self.getTomogramsDir(), self.getTsDir())

    def getInfoDir(self):
        return self._getExtraPath(INFO_DIR)

    def getTomogramsDir(self):
        return self._getExtraPath(TOMOGRAMS_DIR)

    def getTsDir(self):
        return self._getExtraPath(TS_DIR)

    def getStack2dDir(self):
        return self._getExtraPath(PARTICLES_DIR)

    def getStack3dDir(self):
        return self._getExtraPath(PARTICLES_3D_DIR)

    def getSetsDir(self):
        return self._getExtraPath(SETS_DIR)

    def convertOrLink(self, inFile, tsId, outDir, sRate):
        """Fill the simulated EMAN project directories with the expected data at this point of the pipeline.
        Also convert the precedent tomograms into HDF files if they are not. The converted filename will be the tsId,
        avoiding added suffixes this way"""
        hdf = '.hdf'
        outFile = join(outDir, tsId + hdf)
        program = Plugin.getProgram('e2proc3d.py')
        if inFile.endswith(hdf):
            createLink(inFile, self._getExtraPath(outFile))
        else:
            args = '%s %s --apix %.2f ' % (abspath(inFile), outFile, sRate)
            self.runJob(program, args, cwd=self._getExtraPath())
        return outFile

    def convertBetweenHdfAndMrc(self, inFile, outFile, extraArgs=''):
        program = Plugin.getProgram("e2proc3d.py")
        args = '%s %s ' % (inFile, outFile)
        self.runJob(program, args + extraArgs, cwd=self._getExtraPath())

    def _getLstFile(self):
        lstFile = glob.glob(join(self.getSetsDir(), TOMOBOX + '*.lst'))[0]
        return join(SETS_DIR, basename(lstFile))
