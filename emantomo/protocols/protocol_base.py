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
import re
from os.path import join, abspath, basename, realpath, relpath
from emantomo import Plugin
from emantomo.constants import INFO_DIR, TOMOGRAMS_DIR, TS_DIR, SETS_DIR, PARTICLES_DIR, PARTICLES_3D_DIR, \
    REFERENCE_NAME, TOMOBOX, SPT_00_DIR, THREED, ALI3D_BASENAME, ALI2D_BASENAME, FSC_MASKED_BNAME, FSC_UNMASKED_BNAME, \
    FSC_MASKED_TIGHT_BNAME, LST_LINE, PARTICLE_IND, PARTICLE_FILE, PART3D_ID
from emantomo.convert import emanFSCsToScipion
from emantomo.convert.lstAlignConvert import EmanLstReader
from emantomo.objects import EmanParticle, EmanSetOfParticles
from pwem.objects import SetOfFSCs
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
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
        self.inParticles = None
        self.inSamplingRate = -1
        self.scaleFactor = 1
        self.voltage = 300
        self.sphAb = 2.7

    # --------------------------- STEPS functions -----------------------------
    def convertTsStep(self, mdObj):
        # The converted TS must be unbinned, because EMAN will read the sampling rate from its header. This is why
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
        # Prevent sub-setting: only the 3d HDF stacks are linked when initializing the corresponding protocols,
        # covering partially the sub-setting functionality (as the stacks from which no particles are present in
        # the subset won't appear in the EMAN project simulated in the current protocol). But, some particles can
        # have been removed from the present stacks. Thus, a lst file is generated for all the stacks removing then
        # the lines corresponding to non-present particle indices
        present3dStacks = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '*.hdf'))
        presentParticlesDict = {}
        for i, stack3d in enumerate(present3dStacks):
            stack3dRealPath = relpath(realpath(stack3d))
            pathKey = join(PARTICLES_3D_DIR, basename(stack3d))
            indList = []
            for particle in self.inParticles.iterItems(where=f'_stack3dHdf="{stack3dRealPath}"'):
                indList.append(str(particle.getIndex()))
            presentParticlesDict[pathKey] = indList

        buildSetProgram = Plugin.getProgram("e2spt_buildsets.py")
        self.runJob(buildSetProgram, '--allparticles', cwd=self._getExtraPath())
        particlesFile = join(self.getSetsDir(), f'{TOMOBOX}.lst')
        with open(particlesFile, 'r') as inLst:
            lines = inLst.readlines()

        with open(particlesFile, 'w') as outLst:
            for line in lines:
                if line[0] == '#':
                    outLst.write(line)
                else:
                    lineContentsList = line.strip().split('\t')
                    particleInd = lineContentsList[0]
                    filePath = lineContentsList[1]
                    if particleInd in presentParticlesDict[filePath]:
                        outLst.write(line)

        # Do the same with the 2d and 3d align files if necessary
        align3dFile = getattr(self.inParticles, EmanSetOfParticles.ALI_3D, String()).get()
        align2dFile = getattr(self.inParticles, EmanSetOfParticles.ALI_2D, String()).get()
        if align3dFile:
            new3dAlignFile = self.getNewAliFile()
            dataDictList = EmanLstReader.read(align3dFile)
            with open(new3dAlignFile, 'w') as out3dAlign:
                for dataDict in dataDictList:
                    line = dataDict[LST_LINE]
                    if line[0] == '#':
                        out3dAlign.write(line)
                    else:
                        particleInd = dataDict[PARTICLE_IND]
                        filePath = dataDict[PARTICLE_FILE]
                        if particleInd in presentParticlesDict[filePath]:
                            out3dAlign.write(line)

        if align2dFile:
            new2dAlignFile = self.getNewAliFile(is3d=False)
            dataDictList = EmanLstReader.read(align2dFile)
            with open(new2dAlignFile, 'w') as out2dAlign:
                for dataDict in dataDictList:
                    line = dataDict[LST_LINE]
                    if line[0] == '#':
                        out2dAlign.write(line)
                    else:
                        particle3dInd = dataDict[PART3D_ID]
                        filePath = dataDict[PARTICLE_FILE].replace(PARTICLES_DIR, PARTICLES_3D_DIR)
                        if particle3dInd in presentParticlesDict[filePath]:
                            out2dAlign.write(line)

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

    def getRefineDir(self):
        return self._getExtraPath(SPT_00_DIR)

    def getNewAliFile(self, is3d=True):
        fName = 'iniAlign3d.lst' if is3d else 'iniAlign2d.lst'
        return join(self.getRefineDir(), fName)

    def getAttrib(self, attribName, getPointer=False):
        attribPointer = getattr(self, attribName)
        return attribPointer if getPointer else attribPointer.get()

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

    def getRefinedAverageFn(self, iterNumber, ext='hdf', half=None):
        # Example, for 9 iterations, the resulting file would be called threed_09.hdf
        pattern = self._getExtraPath(SPT_00_DIR, THREED + f'_{iterNumber:02d}')
        return pattern + f'_{half}.{ext}' if half else pattern + f'.{ext}'

    def getRefineEvenFn(self, iterNumber):
        # Example, for 9 iterations, the resulting file would be called threed_09_even.hdf
        return self.getRefinedAverageFn(iterNumber, half='even')

    def getRefineOddFn(self, iterNumber):
        # Example, for 9 iterations, the resulting file would be called threed_09_odd.hdf
        return self.getRefinedAverageFn(iterNumber, half='odd')

    def getAli3dFile(self, iterNum):
        return self._getExtraPath(SPT_00_DIR, f'{ALI3D_BASENAME}{iterNum:02d}.lst')

    def getAli2dFile(self, iterNum):
        return self._getExtraPath(SPT_00_DIR, f'{ALI2D_BASENAME}{iterNum:02d}.lst')

    def genFscs(self, iterNum):
        sptPath = self._getExtraPath(SPT_00_DIR)
        fscs = SetOfFSCs.create(self._getPath(), template='fscs%s.sqlite')
        fscMasked = join(sptPath, f'{FSC_MASKED_BNAME}{iterNum:02d}.txt')
        fscUnmasked = join(sptPath, f'{FSC_UNMASKED_BNAME}{iterNum:02d}.txt')
        fscTight = join(sptPath, f'{FSC_MASKED_TIGHT_BNAME}{iterNum:02d}.txt')
        emanFSCsToScipion(fscs, fscMasked, fscUnmasked, fscTight)
        return fscs

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob.glob(join(self.getRefineDir(), '*'))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]
