# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
import logging
from enum import Enum
from os.path import basename, join, exists
import numpy as np
from emantomo import Plugin
from emantomo.constants import SYMMETRY_HELP_MSG, REFERENCE_NAME, TOMOGRAMS_DIR
from emantomo.convert import loadJson, readCoordinate3D
from emantomo.protocols.protocol_base import ProtEmantomoBase, REF_VOL, IN_TOMOS
from pwem.emlib.image.image_readers import MRCImageReader
from pyworkflow.object import Set, String
from pyworkflow.protocol import PointerParam, StringParam, FloatParam, LEVEL_ADVANCED, IntParam, BooleanParam, GE, LE, \
    STEPS_PARALLEL
from pyworkflow.utils import Message, replaceExt, redStr, cyanStr
from tomo.objects import SetOfCoordinates3D, Tomogram

logger = logging.getLogger(__name__)


class OutputsTemplateMatch(Enum):
    coordinates = SetOfCoordinates3D


class EmanProtTemplateMatching(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_tempmatch.py* EMAN2 program.
    It is a reference-based picking (template matching).
    """

    _label = 'Template matching picking'
    _possibleOutputs = OutputsTemplateMatch
    stepsExecutionMode = STEPS_PARALLEL
    program = Plugin.getProgram("e2spt_tempmatch.py")

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inTomosDict = None
        self.badTsIds = String('')
        self.zeroCoordsTsIds = String('')

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      important=True,
                      help='Specify tomograms containing reference-like particles to be exctracted. These '
                           'should be "dark contrast".')
        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Reference volume",
                      help="This should be 'light contrast'.")
        self._addBinThreads(form)

        form.addSection(label='options')
        form.addParam('nptcl', IntParam,
                      default=500,
                      label='Maximum no. particles picked among the tomograms',
                      help='If a higher number of particles is detected, the program will take the best N, '
                           'being N the value of the current parameter.')
        form.addParam('delta', FloatParam,
                      default=30,
                      label='Anglular sampling to rotate the reference (deg.',
                      help='The lower value, the higher number of orientations that will be checked.')
        form.addParam('dthr', FloatParam,
                      default=-1,
                      label='Distance threshold (Å)',
                      help='Particles closer than the value introduced will be removed. If default =-1, '
                           'it will be considered as half of the box size of the reference.')
        form.addParam('vthr', FloatParam,
                      default=5,
                      label='Template matching threshold (n sigma)',
                      validators=[GE(1), LE(10)],
                      help='Particles with score lower than the introduced value will be removed. Admitted values '
                           'are [1, 10], where 1 means very insensitive picking (pick a lot of particles with a high '
                           'probability of having a a high number of false positives) and 10 means very sensitive '
                           'picking (strict picking, with less particles picked and very low presence of false '
                           'positives, but with high probability of some or even a lot of particles not to be picked.')
        form.addParam('minvol', IntParam,
                      default=-1,
                      label='Minimum peak volume',
                      expertLevel=LEVEL_ADVANCED,
                      help='If default=-1, this filter is not applied.')
        form.addParam('maxvol', IntParam,
                      default=-1,
                      label='Maximum peak volume',
                      expertLevel=LEVEL_ADVANCED,
                      help='If default=-1, this filter is not applied.')
        form.addParam('symmetry', StringParam,
                      default='c1',
                      label='Symmetry of the reference',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('rmedge', BooleanParam,
                      default=True,
                      label='Remove particles on the edge?')
        form.addParam('rmgold', BooleanParam,
                      default=True,
                      label='Remove particles near gold fiducials?')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        cInputId = self._insertFunctionStep(self.convertRefVolStep,
                                            prerequisites=[],
                                            needsGPU=False)
        for tsId in self.inTomosDict.keys():
            prjId = self._insertFunctionStep(self.prepareEmanPrj, tsId,
                                             prerequisites=[cInputId],
                                             needsGPU=False)
            tmId = self._insertFunctionStep(self.templateMatchingStep, tsId,
                                            prerequisites=prjId,
                                            needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=tmId,
                                              needsGPU=False)
            closeSetStepDeps.append(cOutId)
        self._insertFunctionStep(self.closeOutputSet,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inTomosDict = {tomo.getTsId(): tomo.clone() for tomo in getattr(self, IN_TOMOS).get()}

    def prepareEmanPrj(self, tsId: str):
        logger.info(cyanStr(f'tsId = {tsId}: preparing the EMAN project...'))
        super().createInitEmanPrjDirs()
        tomo = self.inTomosDict[tsId]
        sRate = tomo.getSamplingRate()
        inTomoName = tomo.getFileName()
        # Required to ensure that the sampling rate is correct in the header, as EMAN reads and compares it
        # with the sampling rate of the reference volume to scale the data
        self.convertOrLink(inTomoName, tsId, TOMOGRAMS_DIR, sRate)

    def templateMatchingStep(self, tsId: str):
        try:
            logger.info(cyanStr(f'tsId = {tsId}: running the template matching...'))
            self.runJob(self.program, self._genTempMatchArgs(tsId), cwd=self._getExtraPath())
        except Exception as e:
            tomo = self.inTomosDict[tsId]
            logger.error(redStr(f'--------->Template matching failed for:'
                         f'\n\ttsId = {tsId}, tomoName = {tomo.getFileName()}'), exc_info=e)
            self.badTsIds.set(self.badTsIds.get() + f' {tsId}')

    def createOutputStep(self, tsId: str):
        with self._lock:
            logger.info(cyanStr(f'tsId = {tsId}: registering the picked coordinates...'))
            tomo = self.inTomosDict[tsId]
            outCoords = self.createOutputSet()
            tomoJsonFile = join(self.getInfoDir(), f'{tomo.getTsId()}_info.json')
            if exists(tomoJsonFile):
                jsonBoxDict = loadJson(tomoJsonFile)
                boxes = jsonBoxDict.get("boxes_3d", None)
                if boxes:
                    for box in boxes:
                        newCoord = readCoordinate3D(box, tomo)
                        outCoords.append(newCoord)
                    outCoords.write()
                    self._store(outCoords)
                else:
                    logger.error(redStr(f'tsId = {tsId} --> No coordinates picked ({tomoJsonFile})'))
                    self.zeroCoordsTsIds.set(self.zeroCoordsTsIds.get() + f' {tsId}')
            else:
                logger.error(redStr(f'tsId = {tsId} --> Json file not found ({tomoJsonFile})'))

    def closeOutputSet(self):
        self._closeOutputSet()
        # Throw an exception if no coordinates were registered
        if len(getattr(self, self._possibleOutputs.coordinates.name, '')) == 0:
            raise ValueError('ERROR!!! No coordinates were registered.')

    # --------------------------- UTILS functions -----------------------------
    def _genTempMatchArgs(self, tsId: str) -> str:
        convertedOrLinkedTomo = join(self.getTomogramsDir(), tsId + '.hdf')
        args = [f" {self._getTomoName(convertedOrLinkedTomo)}",
                f"--ref {REFERENCE_NAME}.hdf ",
                f"--nptcl {self.nptcl.get()}",
                f"--dthr {self.dthr.get():.2f}",
                f"--vthr {self.vthr.get():.2f}",
                f"--minvol {self.minvol.get()}",
                f"--maxvol {self.maxvol.get()}",
                f"--delta {self.delta.get():.2f}",
                f"--sym {self.symmetry.get()}",
                f"--threads {self.binThreads.get()}"]
        if self.rmedge.get():
            args.append("--rmedge")
        if self.rmgold.get():
            args.append('--rmgold')
        return ' '.join(args)

    def createOutputSet(self) -> SetOfCoordinates3D:
        outCoords = getattr(self, self._possibleOutputs.coordinates.name, None)
        if outCoords:
            outCoords.enableAppend()
        else:
            outCoords = SetOfCoordinates3D.create(self._getPath(), template="coordinates%s")
            outCoords.setPrecedents(self.inTomosDict)
            outCoords.setSamplingRate(self.inTomosDict.getSamplingRate())
            outCoords.setBoxSize(self.getRefVol().getDim()[0])
            outCoords.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{self._possibleOutputs.coordinates.name: outCoords})
            self._defineSourceRelation(getattr(self, IN_TOMOS), outCoords)

        return outCoords

    @staticmethod
    def _getTomoName(tomoFile: str) -> str:
        return join(TOMOGRAMS_DIR, replaceExt(basename(tomoFile), 'hdf'))

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        tol = 1e-03
        inTomos = self.getAttrib(IN_TOMOS)
        reference = self.getAttrib(REF_VOL)
        tomosSRate = inTomos.getSamplingRate()
        refSRate = reference.getSamplingRate()
        refDims = MRCImageReader.getDimensions(reference.getFileName())
        tomoXDims, tomoYDims, _, _ = zip(*[MRCImageReader.getDimensions(tomo.getFileName())
                                                   for tomo in inTomos])
        tomoXDims = np.array(tomoXDims)
        tomoYDims = np.array(tomoYDims)
        if abs(tomosSRate - refSRate) >= tol:
            errorMsg.append(f'The sampling rate of the input tomograms [{tomosSRate} Å/pix] and the reference volume '
                            f'[{refSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')
        if refDims[0] < 32:
            errorMsg.append(f'Box sizes below 32 px are not allowed.')
        if np.any(tomoXDims < 900) or np.any(tomoYDims < 900):
            errorMsg.append('The lowest value admitted for X, Y dimensions of the introduced tomograms is 900 px. '
                            'This is because the lowest tomogram reconstruction allowed by native EMAN is around 1k '
                            'px.')
        return errorMsg

    def _summary(self):
        msg = []
        badTsIds = self.badTsIds.get()
        zeroPickedTsIds = self.zeroCoordsTsIds.get()
        if badTsIds:
            msg.append(f'Some tomograms were not processed.\n'
                       f'*{badTsIds}*'
                       f'Check the log for more info\n')
        if zeroPickedTsIds:
            msg.append('Zero coordinates were picked in tomograms:\n'
                       f'{zeroPickedTsIds}')
        return msg
