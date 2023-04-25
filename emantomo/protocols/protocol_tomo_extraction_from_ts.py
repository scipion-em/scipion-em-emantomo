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
from os.path import abspath, join, exists
from emantomo import Plugin
from emantomo.constants import TOMO_ID, GROUP_ID, PARTICLES_3D_DIR, PARTICLES_DIR
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_COORDS, IN_CTF, IN_TS, IN_BOXSIZE
from emantomo.utils import getFromPresentObjects, genEmanGrouping
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, FloatParam, LEVEL_ADVANCED, GE, LE, GT, IntParam
from pyworkflow.utils import Message, replaceBaseExt
from emantomo.convert import coords2Json, ts2Json, ctfTomo2Json
from tomo.constants import TR_EMAN
from tomo.objects import SetOfSubTomograms, SubTomogram
from tomo.protocols import ProtTomoBase
from tomo.utils import getNonInterpolatedTsFromRelations

SAME_AS_PICKING = 0
OTHER = 1


class outputObjects(Enum):
    subtomograms = SetOfSubTomograms


class EmanProtTSExtraction(ProtEmantomoBase):
    """Extract 2D subtilt particles from the tilt series, and reconstruct 3D subvolumes."""

    _label = 'particles extraction from TS'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.groupIds = None
        self.emanDict = None
        self.sphAb = None
        self.voltage = None
        self.scaleFactor = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_COORDS, PointerParam,
                      label="Coordinates",
                      important=True,
                      pointerClass='SetOfCoordinates3D',
                      help='The corresponding tomograms data will be accessed from the provided coordinates.')
        form.addParam(IN_CTF, PointerParam,
                      label="CTF tomo series",
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      help='Estimated CTF for the tilt series associates to the tomograms used to pick the input '
                           'coordinates. The corresponding tilt series data will be also accessed through them.')
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt series with alignment, non-interpolated',
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction. '
                           'To be deprecated!!')
        form.addParam(IN_BOXSIZE, FloatParam,
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
            self._insertFunctionStep(self.convertTsStep, mdObj)
            self._insertFunctionStep(self.convertTomoStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
            self._insertFunctionStep(self.extractParticlesStep, mdObj)
            self._insertFunctionStep(self.convertOutputStep, mdObj)
        self._insertFunctionStep(self.createOutputStep, mdObjDict)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        coords = getattr(self, IN_COORDS).get()
        inCtfSet = getattr(self, IN_CTF).get()
        inTsSet = self.getTs()
        # Get the group ids and the emanDict to have the correspondence between the previous classes and
        # how EMAN will refer them
        uniqueTomoValsDict = getFromPresentObjects(coords, [TOMO_ID, GROUP_ID])
        self.groupIds = uniqueTomoValsDict[GROUP_ID]
        self.emanDict = genEmanGrouping(self.groupIds)
        # Calculate the scale factor
        self.scaleFactor = coords.getSamplingRate() / inTsSet.getSamplingRate()
        # Generate the md object data units as a dict
        return self.genMdObjDict(inTsSet, inCtfSet, coords=coords)

    def writeData2JsonFileStep(self, mdObj):
        mode = 'a' if exists(mdObj.jsonFile) else 'w'
        coords2Json(mdObj, self.emanDict, self.groupIds, self.getBoxSize(), mode=mode)
        ts2Json(mdObj, mode='a')
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
        for tsId, mdObj in mdObjDict.items():
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

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        valMsg = []
        if not self.inputTS.get():
            try:
                getNonInterpolatedTsFromRelations(self.inputCoordinates.get(), self)
            except:
                valMsg.append('Unable to go through the relations from the introduced coordinates to the '
                              'corresponding non-interpolated tilt series. Please introduce them using the '
                              'advanced parameter "Tilt series with alignment, non-interpolated."')
        return valMsg

    # --------------------------- UTILS functions ----------------------------------
    def _genExtractArgs(self, mdObj):
        args = '%s ' % mdObj.tomoHdfName
        args += '--boxsz_unbin=%i ' % self.getBoxSize()
        args += '--maxtilt=%i ' % self.maxTilt.get()
        args += '--tltkeep=%.2f ' % self.tltKeep.get()
        args += '--padtwod=%.2f ' % self.paddingFactor.get()
        args += '--rmbeadthr=%.2f ' % self.rmThr.get()
        args += '--threads=%i ' % self.numberOfThreads.get()
        args += '--newlabel=%s ' % mdObj.tsId
        args += '--append '
        args += '--verbose=9 '
        # if self.doSkipCtfCorrection.get():
        #     args += '--noctf '
        # if self.skip3dRec.get():
        #     args += '--skip3d '
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


