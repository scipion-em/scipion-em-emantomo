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
import logging
from enum import Enum
from os.path import join, exists, basename, abspath
from typing import Dict, Union, Set, Any, List

import numpy as np
from emantomo import Plugin
from emantomo.constants import GROUP_ID, PARTICLES_3D_DIR, PARTICLES_DIR, TOMOBOX, TOMOGRAMS_DIR, TS_DIR
from emantomo.objects import EmanHdf5Handler, EmanSetOfParticles, EmanParticle, EmanMetaData
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_CTF, IN_TS, IN_BOXSIZE, IN_SUBTOMOS
from emantomo.utils import getFromPresentObjects, genEmanGrouping, getPresentTsIdsInSet, genJsonFileName, \
    getPresentPrecedents
from pwem.convert.headers import fixVolume
from pwem.objects import Transform
from pyworkflow.protocol import PointerParam, FloatParam, LEVEL_ADVANCED, GE, LE, GT, IntParam, BooleanParam, \
    STEPS_PARALLEL
from pyworkflow.utils import Message, replaceExt, removeExt, createLink, cyanStr, yellowStr
from emantomo.convert import coords2Json, ctfTomo2Json, ts2Json
from tomo.constants import TR_EMAN, SCIPION
from tomo.objects import SetOfLandmarkModels, LandmarkModel, Coordinate3D, SetOfSubTomograms, SetOfCoordinates3D, \
    SetOfMeshes, SetOfCTFTomoSeries, SetOfTiltSeries, Tomogram, CTFTomoSeries

logger = logging.getLogger(__name__)


class outputObjects(Enum):
    subtomograms = EmanSetOfParticles
    projected2DCoordinates = SetOfLandmarkModels


class EmanProtTSExtraction(ProtEmantomoBase):
    """Extract 2D subtilt particles from the tilt series, and reconstruct 3D subvolumes."""

    _label = 'Extraction from TS'
    _possibleOutputs = outputObjects
    stepsExecutionMode = STEPS_PARALLEL
    program = Plugin.getProgram('e2spt_extract.py')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # self.mdObjDict = {}
        self.coords = None
        self.inCtfSet = None
        self.inTsSet = None
        self.commonTsIdList = None
        self.procIndDict = {}
        self.groupIds = None
        self.emanDict = None
        self.sphAb = None
        self.scaleFactor = None
        self.projectionsDict = {}
        self.currentBoxSize = None
        self.currentSRate = None
        self.fiducialSize = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_SUBTOMOS, PointerParam,
                      label="Coordinates or 3D particles",
                      important=True,
                      pointerClass='SetOfCoordinates3D, SetOfSubTomograms',
                      help='The corresponding tomograms data will be accessed from the provided coordinates or the '
                           'coordinates associated to the 3D particles.')
        form.addParam(IN_CTF, PointerParam,
                      label="CTF tomo series",
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      help='Estimated CTF for the tilt series associates to the tomograms used to pick the input '
                           'coordinates. The corresponding tilt series data will be also accessed through them.')
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt series',
                      # expertLevel=LEVEL_ADVANCED,
                      important=True,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction.')
        form.addParam('doFlipZInTomo', BooleanParam,
                      default=True,
                      important=True,
                      label='Flip Z axis in tomogram?',
                      help='If the reconstruction was carried out with EMAN, it would be set to No.')
        form.addParam(IN_BOXSIZE, IntParam,
                      allowsNull=False,
                      important=True,
                      label='Box size unbinned (px)',
                      help='The subtomograms are extracted as a cubic box of this size.')
        form.addParam('shrink', FloatParam,
                      default=1,
                      allowsNull=False,
                      important=True,
                      label='Particles binning factor',
                      help='For example, if the unbinned box size is 160 pix and the particles binning factor '
                           'introduced is 4, the 2D tilt particles will be cropped on the tilt series with a box of '
                           '160 x 160 pix, and then shrink to 160 / 4 = 40 pix. Thus, both resulting 2D and 3D sets '
                           'of particles will be of size 40 pix.')
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
        form.addParam('minDist', FloatParam,
                      label='Minimum distance between particles (Å)',
                      default=10,
                      validators=[GE(0)],
                      expertLevel=LEVEL_ADVANCED)
        self._addBinThreads(form)
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        pIds = []
        self._initialize()
        for tsId in self.commonTsIdList:
            convertInId = self._insertFunctionStep(self.convertTsStep, tsId, prerequisites=[])
            writeJsonId = self._insertFunctionStep(self.writeData2JsonFileStep, tsId, prerequisites=convertInId)
            extractId = self._insertFunctionStep(self.extractParticlesStep, tsId, prerequisites=writeJsonId)
            convertOutId = self._insertFunctionStep(self.convertOutputStep, tsId, prerequisites=extractId)
            createOutputId = self._insertFunctionStep(self.createOutputStep, tsId, prerequisites=convertOutId)
            pIds.append(createOutputId)
        self._insertFunctionStep(self._closeOutputSet, prerequisites=pIds)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        logger.info(cyanStr('Initializing...'))
        self.createInitEmanPrjDirs()
        inParticles = self.getAttrib(IN_SUBTOMOS)
        self.inCtfSet = self.getAttrib(IN_CTF)
        self.inTsSet = self.getAttrib(IN_TS)
        if type(inParticles) in [SetOfCoordinates3D, SetOfMeshes]:  # Extraction from coords
            self.coords = inParticles
        else:  # Extraction from particles
            self.coords = inParticles.getCoordinates3D()
        binFactor = self.shrink.get()
        self.currentBoxSize = self.getAttrib(IN_BOXSIZE) / binFactor
        self.currentSRate = self.inTsSet.getSamplingRate() * binFactor

        # Get the group ids and the emanDict to have the correspondence between the previous classes and
        # how EMAN will refer them
        uniqueTomoValsDict = getFromPresentObjects(self.coords, [Coordinate3D.TOMO_ID_ATTR, GROUP_ID])
        self.groupIds = uniqueTomoValsDict[GROUP_ID]
        self.emanDict = genEmanGrouping(self.groupIds)
        # Calculate the scale factor for coords and shift scaling: according to EMAN's behavior, it is the rate between
        # the sampling rate of the input coordinates and the TS unbinned sampling rate multiplied by the binning factor
        # introduced
        self.scaleFactor = inParticles.getSamplingRate() / (self.inTsSet.getSamplingRate() * binFactor)
        # The fiducial size is the diameter in nm
        self.fiducialSize = 0.1 * self.getAttrib(IN_BOXSIZE) * self.getAttrib(IN_TS).getSamplingRate() / 2
        # Generate the list of the common tsIds:
        # Sorting the items before the iteration based on the tsId (key) will ensure the alphabetical order in
        # all the data generated, preventing data mismatching between the alignment files generated by EMAN in
        # posterior refinements and its corresponding version in Scipion format
        self.commonTsIdList = sorted(self._genCommonTsIdList(self.inTsSet, self.inCtfSet, self.coords))
        processingInd = 0
        for tsId in self.commonTsIdList:
            self.procIndDict[tsId] = processingInd
            tomo = self.getCurrentTomo(tsId, doLock=False)
            processingInd += len([coord for coord in self.coords.iterCoordinates(volume=tomo)])

    def convertTsStep(self, tsId: str):
        # Create project dir structure
        tomogramsDir = self.getTomogramsDir()
        # TS must be in HDF (in EMAN's native code it is searched in hardcoded HDF format based on the tomogram
        # basename), so it will be converted later. However, the tomograms can be linked in MRC without any problem
        with self._lock:
            currentTomo = self.getCurrentTomo(tsId, doLock=False)
            currentTs = self.getCurrentTs(tsId, doLock=False)
        tomoFName = currentTomo.getFileName()
        createLink(tomoFName, join(tomogramsDir, basename(tomoFName)))
        # Convert the TS into HDF
        # The converted TS must be unbinned, because EMAN will read the sampling rate from its header. This is why
        # the TS associated to the CTF is the one considered first. Later, when generating the json, the TS alignment
        # parameters are read from the introduced TS and the shifts are scaled to at the unbinned scale
        logger.info(cyanStr(f'===> tsId = {tsId}: converting the tilt-series into HDF...'))
        inTsFName = currentTs.getFirstItem().getFileName()
        sRate = currentTs.getSamplingRate()
        self.convertOrLink(inTsFName, tsId, TS_DIR, sRate)

    def writeData2JsonFileStep(self, tsId: str):
        logger.info(cyanStr(f'===> tsId = {tsId}: writing the json files...'))
        with self._lock:
            currentTs = self.getCurrentTs(tsId, doLock=False)
            currentCtf = self.getCurrentCtf(tsId, doLock=False)
            currentTomo = self.getCurrentTomo(tsId, doLock=False)
        acq = currentTs.getAcquisition()
        jsonFile = genJsonFileName(self.getInfoDir(), tsId)
        coords = [coord.clone() for coord in self.coords.iterCoordinates(volume=currentTomo)]
        mode = 'a' if exists(jsonFile) else 'w'

        coords2Json(coords, currentTomo, currentTs.getSamplingRate(), jsonFile, self.emanDict,
                    self.groupIds, self.getBoxSize(), doFlipZ=self.doFlipZInTomo.get(),
                    mode=mode)
        ts2Json(currentTs, jsonFile, currentCtf, mode='a')
        ctfTomo2Json(currentCtf, jsonFile, acq.getSphericalAberration(), acq.getVoltage(), mode='a')

    def extractParticlesStep(self, tsId: str):
        logger.info(cyanStr(f'===> tsId = {tsId}: extracting the particles...'))
        self.runJob(self.program, self._genExtractArgs(tsId), cwd=self._getExtraPath())

    def convertOutputStep(self, tsId: str):
        logger.info(f'===> {tsId}: converting and unstacking the 3D particles into MRC...')
        fName = self._getEmanFName(tsId)
        # Unstack the 3d particles HDF stack into individual MRC files
        stack3d = join(PARTICLES_3D_DIR, fName)
        self.unstackParticles(stack3d)
        fixVolume(glob.glob(join(self.getStack3dDir(), f"{removeExt(fName)}*.mrc")))
        # # Convert the 2d particles HDF stack into MRC for visualization purposes (compatibility with viewers)
        # stack2d = join(PARTICLES_DIR, fName)
        # inFile = join(PARTICLES_DIR, basename(stack2d))
        # outFile = replaceExt(inFile, 'mrc')
        # self.convertBetweenHdfAndMrc(inFile, outFile)

    def createOutputStep(self, tsId: str):
        logger.info(f'===> {tsId}: creating the corresponding outputs...')
        with self._lock:
            currentTs = self.getCurrentTs(tsId, doLock=False)
            currentTomo = self.getCurrentTomo(tsId, doLock=False)

        jsonFile = genJsonFileName(self.getInfoDir(), tsId)
        tsHdfName = join(TS_DIR, f'{tsId}.hdf')
        subtomoSet = self.getOutSetOfParticles()
        # Generate the fiducial model (for data visualization purpose)
        fiducialModelGaps = self.getOutSetIfFiducials()

        # Get the projections
        tsProjections = self.getProjections(tsId)
        landmarkModelGaps = self.createLandmarkModelGaps(currentTs)

        # Because the order is changed when the particles are processed with EMAN, if they have transformation matrix,
        # there will be a data mismatch. To avoid it, The EMAN coordinates will be read from the 3d HDF stack, then the
        # corresponding coordinate will be located among the input set, so the transformation is correctly matched, and
        # finally, the corresponding output unstacked MRC file will be located and added to the output particle.
        emanCoords = self.getEmanCoordsFromHdfStack(tsId)

        _, _, nImgs = currentTs.getDimensions()
        subtomoFiles = sorted(glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.mrc' % tsId)))
        stack2d = glob.glob(self._getExtraPath(PARTICLES_DIR, '%s*.hdf' % tsId))[0]
        stack3d = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.hdf' % tsId))[0]
        coords = [coord.clone() for coord in self.coords.iterCoordinates(volume=currentTomo)]

        nTotalParticles = self.procIndDict[tsId]
        particleCounter = 0
        for emanCoord, subtomoFile in zip(emanCoords, subtomoFiles):
            self.fillLandmarkModel(landmarkModelGaps, tsProjections, particleCounter, nImgs)
            subtomogram = EmanParticle()
            transform = Transform()
            subtomogram.setFileName(subtomoFile)
            subtomogram.setLocation(particleCounter, subtomoFile)  # The index is stored to be used as
            # the position of the particle in the 3d HDF stack (to handle subsets - LST file gen.)
            subtomogram.setSamplingRate(self.currentSRate)
            subtomogram.setVolName(tsId)
            scipionCoord, inCoords = self.getEmanMatchingCoord(emanCoord, coords)
            subtomogram.setCoordinate3D(scipionCoord)
            M = scipionCoord.getMatrix()
            shift_x = M[0, 3]
            shift_y = M[1, 3]
            shift_z = M[2, 3]
            transform.setMatrix(M)
            transform.setShifts(self.scaleFactor * shift_x,
                                self.scaleFactor * shift_y,
                                self.scaleFactor * shift_z)
            subtomogram.setTransform(transform, convention=TR_EMAN)
            # Fill EmanParticle own attributes
            subtomogram.setInfoJson(jsonFile)
            subtomogram.setTsHdf(self._getExtraPath(tsHdfName))
            subtomogram.setStack2dHdf(stack2d)
            subtomogram.setStack3dHdf(stack3d)
            subtomogram.setAbsIndex(nTotalParticles)
            subtomoSet.append(subtomogram)
            particleCounter += 1
            nTotalParticles += 1
        subtomoSet.write()
        self._store(subtomoSet)
        fiducialModelGaps.append(landmarkModelGaps)
        fiducialModelGaps.write()
        self._store(fiducialModelGaps)

    # --------------------------- INFO functions -----------------------------------
    def _warnings(self):
        warnMsg = []
        if not (self.getAttrib(IN_TS).hasAlignment() and not self.getAttrib(IN_TS).interpolated()):
            warnMsg.append('The introduced tilt series do not have an alignment transformation associated.')
        return warnMsg

    # --------------------------- UTILS functions ----------------------------------
    def _genCommonTsIdList(self,
                     inTsSet: SetOfTiltSeries,
                     inCtfSet: SetOfCTFTomoSeries,
                     coords: Union[SetOfCoordinates3D, SetOfSubTomograms]) -> Union[Set[str], None]:
        self.createInitEmanPrjDirs()
        if type(coords) is SetOfSubTomograms:
            coords = coords.getCoordinates3D()
        # Considering the possibility of subsets, let's find the tsIds present in all the sets ob objects introduced,
        # which means the intersection of the tsId lists
        presentCtfTsIds = set(inCtfSet.getTSIds())
        presentTsSetTsIds = set(inTsSet.getTSIds())
        presentTomoSetTsIds = set(coords.getTSIds())
        presentTsIds = presentCtfTsIds & presentTsSetTsIds & presentTomoSetTsIds
        logger.info(cyanStr(f"TsIds present the introduced sets of tilt-series, CTFs and coordinates are: "
                            f"{presentTsIds}" ))
        nonPresentTsIds = self.getNonCommonTsIds(presentTsSetTsIds, presentCtfTsIds, presentTomoSetTsIds)
        if nonPresentTsIds:
            logger.warning(yellowStr(f'Some tsIds are not present in all the sets introduced (tilt-series, '
                                     f'CTFs and coordinates): {nonPresentTsIds}'))
        return presentTsIds

    # def genMdObjDict(self,
    #                  inTsSet: SetOfTiltSeries,
    #                  inCtfSet: SetOfCTFTomoSeries,
    #                  coords: Union[SetOfCoordinates3D, SetOfSubTomograms]) -> Dict[Any, EmanMetaData]:
    #     self.createInitEmanPrjDirs()
    #     if type(coords) is SetOfSubTomograms:
    #         coords = coords.getCoordinates3D()
    #     # Considering the possibility of subsets, let's find the tsIds present in all the sets ob objects introduced,
    #     # which means the intersection of the tsId lists
    #     presentCtfTsIds = set(inCtfSet.getTSIds())
    #     presentTsSetTsIds = set(inTsSet.getTSIds())
    #     presentTomoSetTsIds = set(coords.getTSIds())
    #     presentTsIds = presentCtfTsIds & presentTsSetTsIds & presentTomoSetTsIds
    #     logger.info(cyanStr(f"TsIds present the introduced sets of tilt-series, CTFs and coordinates are: "
    #                         f"{presentTsIds}" ))
    #     nonPresentTsIds = self.getNonCommonTsIds(presentTsSetTsIds, presentCtfTsIds, presentTomoSetTsIds)
    #     if nonPresentTsIds:
    #         logger.warning(yellowStr(f'Some tsIds are not present in all the sets introduced (tilt-series, '
    #                                  f'CTFs and coordinates): {nonPresentTsIds}'))
    #     # The tomograms are obtained as the coordinates precedents. Operating this way, the code considers the case of
    #     # subsets of coordinates
    #     presentTomograms = getPresentPrecedents(coords, presentTsIds)
    #
    #     # Manage the TS, CTF tomo Series and tomograms
    #     tsIdsDict = {ts.getTsId(): ts.clone() for ts in inTsSet if ts.getTsId() in presentTsIds}
    #     ctfIdsDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inCtfSet if ctf.getTsId() in presentTsIds}
    #     tomoIdsDict = {tomo.getTsId(): tomo for tomo in presentTomograms}
    #
    #     # Get the required acquisition data
    #     self.sphAb = inTsSet.getAcquisition().getSphericalAberration()
    #     self.voltage = inTsSet.getAcquisition().getVoltage()
    #
    #     # Split all the data into subunits (EmanMetaData objects) referred to the same tsId
    #     mdObjDict = {}
    #     procInd = 0
    #     for tomoId, ts in sorted(tsIdsDict.items()):
    #         # Sorting the items before the iteration based on the tsId (key) will ensure the alphabetical order in
    #         # all the data generated, preventing data mismatching between the alignment files generated by EMAN in
    #         # posterior refinements and its corresponding version in Scipion format
    #         tomo = tomoIdsDict[tomoId]
    #         hdfFileBaseName = f'{tomoId}.hdf'
    #         iCoords = [coord.clone() for coord in coords.iterCoordinates(volume=tomo)] if coords else None
    #         iParticles = [particle.clone() for particle in self.inParticles.iterSubtomos(volume=tomo)] if \
    #             self.inParticles else None
    #         mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
    #                                          inTomo=tomo,
    #                                          tomoHdfName=join(TOMOGRAMS_DIR, hdfFileBaseName),
    #                                          ts=ts,
    #                                          tsHdfName=join(TS_DIR, hdfFileBaseName),
    #                                          ctf=ctfIdsDict[tomoId],
    #                                          coords=iCoords,
    #                                          particles=iParticles,
    #                                          jsonFile=genJsonFileName(self.getInfoDir(), tomoId),
    #                                          processingInd=procInd)
    #         procInd += len(iCoords)
    #     if not mdObjDict:
    #         raise Exception("There isn't any common tilt series among the coordinates, the CTF tomo series and "
    #                         "the tilt series chosen.")
    #     return dict(sorted(mdObjDict.items()))

    def _genExtractArgs(self, tsId: str):
        # align2dFile = self.getNewAliFile(is3d=False)
        currentTomo = self.getCurrentTomo(tsId, doLock=True)
        tomoFiles = join(TOMOGRAMS_DIR, basename(currentTomo.getFileName()))
        align3dFile = self.getNewAliFile()
        args = [
            # The particles can be read using the input tomograms or directly the particles contained in the ali3d file
            # in clase it exists
            f'--jsonali={align3dFile}' if exists(self._getExtraPath(align3dFile)) else f'{tomoFiles}',
            f'--boxsz_unbin={self.getBoxSize()}',
            f'--shrink={self.shrink.get():.2f}',
            f'--maxtilt={self.maxTilt.get()}',
            f'--tltkeep={self.tltKeep.get():.2f}',
            f'--padtwod={self.paddingFactor.get():.2f}',
            f'--rmbeadthr={self.rmThr.get():.2f}',
            f'--mindist={self.minDist.get():.2f}',
            f'--threads={self.binThreads.get()}',
            '--append',
            '--verbose=9'
        ]
        return ' '.join(args)

    def unstackParticles(self, stackFile, outExt='mrc'):
        """Unstacks and coverts a list of stack files into separate particles"""
        program = Plugin.getProgram('e2proc3d.py')
        args = [
            '--unstacking',
            f'{stackFile}',
            f'{replaceExt(stackFile, outExt)}'
        ]
        self.runJob(program, ' '.join(args), cwd=self._getExtraPath())

    def _getEmanFName(self, tsId):
        # if self.isReExtraction:
        #     partFiles = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, f'{tsId}*_reextract.hdf'))
        #     return basename(partFiles[0])
        # else:
        shrink = self.shrink.get()
        pattern = f'__{TOMOBOX}'
        if shrink:
            if shrink > 1:
                pattern += f'_bin{int(self.shrink.get())}'
        return tsId + pattern + '.hdf'

    def getProjections(self, tsId: str):
        """load the corresponding HDF 2d stack file and read from the header the required data to create later a
        landmark model. To do that, the h5py module is used.
        """
        stack2dHdf = abspath(self._getExtraPath(PARTICLES_DIR, self._getEmanFName(tsId)))
        eh = EmanHdf5Handler(stack2dHdf)
        projList = eh.getProjsFrom2dStack(shrink=self.shrink.get())
        # Add the required tsId as the first element of each sublist
        list(map(lambda sublist: sublist.insert(0, tsId), projList))  # More efficient than comprehension for huge
        # lists of lists, as expected
        return projList

    def createLandmarkModelGaps(self, ts):
        tsId = ts.getTsId()
        landmarkModelGapsFilePath = self._getExtraPath(str(tsId) + "_gaps.sfid")
        landmarkModelGaps = LandmarkModel(tsId=tsId,
                                          tiltSeriesPointer=ts,
                                          fileName=landmarkModelGapsFilePath,
                                          modelName=None,
                                          size=self.fiducialSize,
                                          applyTSTransformation=False)
        landmarkModelGaps.setTiltSeries(ts)
        return landmarkModelGaps

    @staticmethod
    def fillLandmarkModel(landmarkModelGaps, tsProjections, particle3dInd, nImgs):
        ind = particle3dInd * nImgs
        for i in range(ind, ind + nImgs):
            try:
                # Found cases in which some 2d tilt images seem to have been excluded, e.g., TS with 58 images, 819
                # coordinates, and the stack 2d is of size 48251, which is not divisible by 59 (48251 / 59 = 817.8135).
                # But in the end it's not critical as this 2d landmark models are only used for visualization and
                # result checking purposes
                projection = tsProjections[i]
                tiltId = projection[1] + 1
                partId = particle3dInd + 1
                xCoord = round(projection[2])
                yCoord = round(projection[3])
                landmarkModelGaps.addLandmark(xCoord, yCoord, tiltId, partId, 0, 0)
            except Exception as e:
                continue

    def getOutSetOfParticles(self):
        outParticles = getattr(self, self._possibleOutputs.subtomograms.name, None)
        if outParticles:
            outParticles.enableAppend()
        else:
            inTsPointer = getattr(self, IN_TS)
            inParticlesPointer = getattr(self, IN_SUBTOMOS)
            outParticles = EmanSetOfParticles.create(self._getPath(), template='emanParticles%s.sqlite')
            if isinstance(inParticlesPointer.get(), SetOfSubTomograms):
                outParticles.copyInfo(inParticlesPointer.get())
            else:
                # The input are coordinates
                outParticles.setCoordinates3D(inParticlesPointer)
            outParticles.setSamplingRate(self.currentSRate)
            # Define outputs and relations
            self._defineOutputs(**{self._possibleOutputs.subtomograms.name: outParticles})
            self._defineSourceRelation(getattr(self, IN_SUBTOMOS), outParticles)
            self._defineSourceRelation(getattr(self, IN_CTF), outParticles)
            self._defineSourceRelation(inTsPointer, outParticles)

        return outParticles

    def getOutSetIfFiducials(self):
        outFiducials = getattr(self, self._possibleOutputs.projected2DCoordinates.name, None)
        if outFiducials:
            outFiducials.enableAppend()
        else:
            inTsPointer = getattr(self, IN_TS)
            inTs = inTsPointer.get()
            outFiducials = SetOfLandmarkModels.create(self._getPath(), suffix='Gaps')
            outFiducials.copyInfo(inTs)
            outFiducials.setSetOfTiltSeries(inTsPointer)
            # Define outputs and relations
            self._defineOutputs(**{self._possibleOutputs.projected2DCoordinates.name: outFiducials})
            self._defineSourceRelation(inTsPointer, outFiducials)
        return outFiducials

    def getEmanCoordsFromHdfStack(self, tsId):
        fName = self._getEmanFName(tsId)
        # Unstack the 3d particles HDF stack into individual MRC files
        stack3d = abspath(join(self.getStack3dDir(), fName))
        eh = EmanHdf5Handler(stack3d)
        return eh.get3dCoordsFrom3dStack(invertZ=self.doFlipZInTomo.get())

    def getEmanMatchingCoord(self, emanCoord, scipionCoords):
        def euclideanDist(coord1, coord2):
            return np.linalg.norm(np.array(coord1) * self.scaleFactor - np.array(coord2))

        distances = [euclideanDist(scipionCoord.getPosition(SCIPION), emanCoord) for scipionCoord in scipionCoords]
        minDistInd = np.argmin(distances)
        logger.info(f'MIN DIST = {distances[minDistInd]}')
        matchingCoord = scipionCoords[minDistInd]
        scaledCoords = np.array(matchingCoord.getPosition(SCIPION)) * self.scaleFactor
        matchingCoord.setPosition(scaledCoords[0], scaledCoords[1], scaledCoords[2], SCIPION)
        # Remove the matching coord from the list, so in an iterative project, the list becomes smaller with each
        # iteration, being more efficient
        reducedList = scipionCoords[:minDistInd] + scipionCoords[minDistInd+1:]
        return matchingCoord, reducedList

    @staticmethod
    def getNonCommonTsIds(tsIds: Set[str],
                          ctfTsIds: Set[str],
                          tomoTsIds: Set[str]) -> Set[str]:
        # Union
        union_todos = tsIds | ctfTsIds | tomoTsIds
        # Intersection
        interseccion_todos = tsIds & ctfTsIds & tomoTsIds
        # Elements present in one or two sets, but not in all the three sets
        return union_todos - interseccion_todos

    def getCurrentTomo(self, tsId: str, doLock: bool = True):
        tomograms = self.coords.getPrecedents()
        if doLock:
            with self._lock:
                return tomograms.getItem(Tomogram.TS_ID_FIELD, tsId)
        else:
            return tomograms.getItem(Tomogram.TS_ID_FIELD, tsId)

    def getCurrentCtf(self, tsId: str, doLock: bool = True):
        if doLock:
            with self._lock:
               return self.inCtfSet.getItem(CTFTomoSeries.TS_ID_FIELD, tsId)
        else:
            return self.inCtfSet.getItem(CTFTomoSeries.TS_ID_FIELD, tsId)

