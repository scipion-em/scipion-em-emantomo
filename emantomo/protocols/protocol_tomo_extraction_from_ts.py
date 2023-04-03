# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es)
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
from enum import Enum
from os.path import abspath, join, basename
from emantomo import Plugin
from emantomo.constants import INFO_DIR, TOMO_ID, GROUP_ID, TS_ID, PARTICLES_3D_DIR, PARTICLES_DIR, TOMOGRAMS_DIR
from emantomo.utils import getPresentPrecedents, getBoxSize, getFromPresentObjects, genEmanGrouping
from emantomo.objects import EmanMetaData
from pwem.objects import Transform
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, FloatParam, BooleanParam, LEVEL_ADVANCED, GE, LE, GT, IntParam
from pyworkflow.utils import makePath, Message, replaceBaseExt, createLink
from emantomo.convert import coords2Json, ts2Json, ctfTomo2Json
from tomo.constants import TR_EMAN
from tomo.objects import SetOfSubTomograms, SubTomogram
from tomo.protocols import ProtTomoBase
from tomo.utils import getNonInterpolatedTsFromRelations

SAME_AS_PICKING = 0
OTHER = 1


class outputObjects(Enum):
    subtomograms = SetOfSubTomograms


class EmanProtTSExtraction(EMProtocol, ProtTomoBase):
    """Extract 2D subtilt particles from the tilt series, and reconstruct 3D subvolumes."""

    _label = 'particles extraction from TS'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.groupIds = None
        self.emanDict = None
        self.sphAb = None
        self.voltage = None
        self.presentTsIds = None
        self.scaleFactor = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputCoordinates', PointerParam,
                      label="Coordinates",
                      important=True,
                      pointerClass='SetOfCoordinates3D',
                      help='The corresponding tomograms data will be accessed from the provided coordinates.')
        form.addParam('inputCTF', PointerParam,
                      label="CTF tomo series",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Estimated CTF for the tilt series associates to the tomograms used to pick the input '
                           'coordinates. The corresponding tilt series data will be also accessed through them.')
        form.addParam('inputTS', PointerParam,
                      help="Tilt series with alignment (non interpolated) used in the tomograms reconstruction. "
                           "To be deprecated!!",
                      pointerClass='SetOfTiltSeries',
                      label="Input tilt series",
                      important=True,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True)
        form.addParam('boxSize', FloatParam,
                      label='Box size unbinned (pix.)',
                      help='The subtomograms are extracted as a cubic box of this size. The wizard selects same '
                           'box size as picking')
        form.addParam('maxTilt', IntParam,
                      default=100,
                      label='Max tilt',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('tltKeep', FloatParam,
                      default=1.0,
                      label='Tilt fraction to keep',
                      expertLevel=LEVEL_ADVANCED,
                      validators=[GT(0), LE(1)],
                      help='Keep a fraction of tilt images with good score determined from tomogram reconstruction')
        form.addParam('rmThr', FloatParam,
                      default=-1,
                      label='Contrast threshold for 2D particle removal',
                      expertLevel=LEVEL_ADVANCED,
                      help='Remove 2d particles with high contrast object beyond N sigma at 100Å. Note that '
                           'this may result in generating fewer particles than selected. Default is -1 '
                           '(include all particles). 0.5 might be a good choice for removing gold beads.')
        form.addParam('paddingFactor', FloatParam,
                      label='Padding factor',
                      default=2,
                      allowsNull=False,
                      validators=[GE(0)],
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to 0, no padding will be considered. If your particles are deeply buried in other '
                           'densities, using a bigger padtwod may help, but doing so may significantly increase the '
                           'memory usage and slow down the process.')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertInputStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
            self._insertFunctionStep(self.extractParticlesStep, mdObj)
            self._insertFunctionStep(self.convertOutputStep, mdObj)
        self._insertFunctionStep(self.createOutputStep, mdObjDict)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        infoDir = self._getExtraPath(INFO_DIR)
        makePath(infoDir, self._getExtraPath(TOMOGRAMS_DIR))

        # Get the group ids and the emanDict to have the correspondence between the previous classes and
        # how EMAN will refer them
        coords = self.inputCoordinates.get()
        uniqueTomoValsDict = getFromPresentObjects(coords, [TOMO_ID, GROUP_ID])
        self.groupIds = uniqueTomoValsDict[GROUP_ID]
        self.emanDict = genEmanGrouping(self.groupIds)

        # Manage the TS and CTF tomo Series
        inCtfSet = self.inputCTF.get()
        inTsSet = self.inputTS.get()
        uniqueCtfValsDict = getFromPresentObjects(inCtfSet, [TS_ID])
        tsIdsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inTsSet if uniqueCtfValsDict[TS_ID]}
        presentTsIds = list(tsIdsDict.keys())
        ctfIdsDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inCtfSet if ctf.getTsId() in presentTsIds}

        # Get the required acquisition data
        self.sphAb = inTsSet.getAcquisition().getSphericalAberration()
        self.voltage = inTsSet.getAcquisition().getVoltage()

        # Split all the data into subunits referred to the same tsId
        mdObjDict = {}
        self.presentTsIds = []
        self.scaleFactor = coords.getSamplingRate() / inTsSet.getSamplingRate()
        for tomo in getPresentPrecedents(coords, uniqueTomoValsDict[TOMO_ID]):
            tomoId = tomo.getTsId()
            self.presentTsIds.append(tomoId)
            mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
                                             inTomo=tomo,
                                             ts=tsIdsDict[tomoId],
                                             ctf=ctfIdsDict[tomoId],
                                             coords=[coord.clone() for coord in coords.iterCoordinates(volume=tomo)],
                                             jsonFile=join(infoDir, '%s_info.json' % tomoId)
                                             # jsonFile=genTomoJsonFileName(tomo.getFileName(), infoDir)
                                             )
        return mdObjDict

    def convertInputStep(self, mdObj):
        """Convert the precedent tomograms into HDF files if they are not"""
        program = Plugin.getProgram('e2proc3d.py')
        inTomoFName = mdObj.inTomo.getFileName()
        if inTomoFName.endswith('.hdf'):
            outFile = join(TOMOGRAMS_DIR, basename(inTomoFName))
            createLink(inTomoFName, outFile)
        else:
            inFile = abspath(inTomoFName)
            outFile = join(TOMOGRAMS_DIR, replaceBaseExt(inTomoFName, 'hdf'))
            args = '%s %s --apix %i ' % (inFile, outFile, mdObj.inTomo.getSamplingRate())
            self.runJob(program, args, cwd=self._getExtraPath())
        # Store the tomoHdfName in the current mdObj
        mdObj.tomoHdfName = outFile

    def writeData2JsonFileStep(self, mdObj):
        coords2Json(mdObj, self.emanDict, self.groupIds, getBoxSize(self))
        ts2Json(mdObj, mode='a')
        # tltParams2Json([mdObj.jsonFile], mdObj.ts, mode="w")
        # If CTF tomo series are introduced, the defocus data is read and added to the corresponding json file
        if mdObj.ctf:
            ctfTomo2Json(mdObj, self.sphAb, self.voltage, mode='a')

    def extractParticlesStep(self, mdObj):
        # TODO: update the command with the new parameters and considerations
        program = Plugin.getProgram('e2spt_extract.py')
        self.runJob(program, self._genExtractArgs(mdObj), cwd=self._getExtraPath())

    def convertOutputStep(self, mdObj):
        # stacks2d = glob.glob(self._getExtraPath(PARTICLES_DIR, '*.hdf'))
        stacks3d = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.hdf' % mdObj.inTomo.getTsId()))
        # self.unstackParticles(stacks2d)
        self.unstackParticles(stacks3d)
        # cleanPattern(hdfFile)

    def createOutputStep(self, mdObjDict):
        sRate = self.inputCoordinates.get().getSamplingRate()
        subtomoSet = SetOfSubTomograms.create(self._getPath(), template='subtomograms%s.sqlite')
        subtomoSet.setSamplingRate(sRate)
        subtomoSet.setCoordinates3D(self.inputCoordinates.get())
        for tsId in self.presentTsIds:
            mdObj = mdObjDict[tsId]
            coords = mdObj.coords
            tomoFile = mdObj.inTomo.getFileName()
            subtomoFiles = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.mrc' % tsId))
            tsSubStack = glob.glob(self._getExtraPath(PARTICLES_DIR, '%s*.hdf' % tsId))[0]
            for coord, subtomoFile in zip(coords, subtomoFiles):
                subtomogram = SubTomogram()
                transform = Transform()
                subtomogram.setFileName(subtomoFile)
                subtomogram.setLocation(subtomoFile)
                subtomogram.setSamplingRate(sRate)
                subtomogram.setCoordinate3D(coord)
                M = coord.getMatrix()
                shift_x = M[0, 3]
                shift_y = M[1, 3]
                shift_z = M[2, 3]
                transform.setMatrix(M)
                transform.setShifts(self.scaleFactor * shift_x,
                                    self.scaleFactor * shift_y,
                                    self.scaleFactor * shift_z)
                subtomogram.setTransform(transform, convention=TR_EMAN)
                subtomogram.setVolName(tomoFile)
                subtomogram._tsSubStack = String(tsSubStack)
                subtomoSet.append(subtomogram)

        self._defineOutputs(**{outputObjects.subtomograms.name: subtomoSet})
        self._defineSourceRelation(self.inputCoordinates.get(), subtomoSet)
        self._defineSourceRelation(self.inputCTF.get(), subtomoSet)
        if self.inputTS.get():
            self._defineSourceRelation(self.inputTS.get(), subtomoSet)

        # # Output pseudosubtomograms --> set of volumes for visualization purposes
        # outputSet = SetOfEmanPseudoSubtomograms.create(self._getPath(), EMAN_PSUBTOMOS_SQLITE)
        # outputSet.setSamplingRate(sRate)
        # for mdObj in mdObjList:
        #     tsId = mdObj.tsId
        #     stackFile, stackLen = emanPrj.getParticlesStackInfo(tsId)
        #     for i in range(1, stackLen + 1):
        #         emanPSubtomo = EmanPSubtomogram(fileName=stackFile,
        #                                         index=i,
        #                                         tsId=tsId,
        #                                         samplingRate=sRate)
        #         outputSet.append(emanPSubtomo)
        #
        # self._defineOutputs(**{outputObjects.emanPSubtomograms.name: outputSet,
        #                        outputObjects.emanProject.name: emanPrj})
        # self._defineSourceRelation(self.inputCoordinates.get(), outputSet)
        # self._defineSourceRelation(self.inputCTF.get(), outputSet)
        # self._defineSourceRelation(self.inputCoordinates.get(), emanPrj)
        # self._defineSourceRelation(self.inputCTF.get(), emanPrj)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        valMsg = []
        if not self.inputTS.get():
            try:
                getNonInterpolatedTsFromRelations(self.inputCoordinates.get(), self)
            except:
                valMsg.append('Unable to go via relations from the introduced coordinates to the '
                              'corresponding non-interpolated tilt series. Please introduce them using the '
                              'advanced parameter "Tilt series with alignment..."')
        return valMsg

    # --------------------------- UTILS functions ----------------------------------
    def _genExtractArgs(self, mdObj):
        args = '%s ' % mdObj.tomoHdfName
        args += '--rmbeadthr=%.2f ' % self.rmThr.get()
        args += '--shrink=%.2f ' % self.downFactor.get()
        args += '--boxsz_unbin=%i ' % getBoxSize(self)
        args += '--tltkeep=%.2f ' % self.tltKeep.get()
        args += '--padtwod=%.2f ' % self.paddingFactor.get()
        args += '--threads=%i ' % self.numberOfThreads.get()
        args += '--append '
        if self.doSkipCtfCorrection.get():
            args += '--noctf '
        if self.skip3dRec.get():
            args += '--skip3d '
        return args

    def unstackParticles(self, stackList, outExt='mrc'):
        """Unstacks and coverts a list of stack files into separate particles"""
        program = Plugin.getProgram('e2proc3d.py')
        outDir = PARTICLES_3D_DIR if PARTICLES_3D_DIR in stackList[0] else PARTICLES_DIR
        for hdfFile in stackList:
            args = ' --unstacking'
            args += ' %s' % abspath(hdfFile)
            args += ' %s' % join(outDir, replaceBaseExt(hdfFile, outExt))
            self.runJob(program, args, cwd=self._getExtraPath())