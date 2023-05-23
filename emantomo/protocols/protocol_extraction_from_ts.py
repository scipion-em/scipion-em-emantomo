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
from os.path import join, exists, basename, abspath
from emantomo import Plugin
from emantomo.constants import GROUP_ID, PARTICLES_3D_DIR, PARTICLES_DIR, TOMOBOX
from emantomo.objects import EmanHdf5Handler, EmanSetOfParticles, EmanParticle, EmanMetaData
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_COORDS, IN_CTF, IN_TS, IN_BOXSIZE
from emantomo.utils import getFromPresentObjects, genEmanGrouping, getPresentTsIdsInSet, genJsonFileName
from pwem.objects import Transform
from pyworkflow.protocol import PointerParam, FloatParam, LEVEL_ADVANCED, GE, LE, GT, IntParam, BooleanParam
from pyworkflow.utils import Message, replaceExt
from emantomo.convert import coords2Json, ts2Json, ctfTomo2Json
from tomo.constants import TR_EMAN
from tomo.objects import SetOfLandmarkModels, LandmarkModel, Coordinate3D


class outputObjects(Enum):
    subtomograms = EmanSetOfParticles
    projected2DCoordinates = SetOfLandmarkModels


class EmanProtTSExtraction(ProtEmantomoBase):
    """Extract 2D subtilt particles from the tilt series, and reconstruct 3D subvolumes."""

    _label = 'particle extraction from TS'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.groupIds = None
        self.emanDict = None
        self.sphAb = None
        self.voltage = None
        self.scaleFactor = None
        self.projectionsDict = {}

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
                      # expertLevel=LEVEL_ADVANCED,
                      important=True,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction.')
        form.addParam('doFlipZInTomo', BooleanParam,
                      default=True,
                      important=True,
                      label='Flip Z axis in tomogram?',
                      help='If the reconstruction was carried out with EMAN, it would be set to No.')
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
            self._insertFunctionStep(super().convertTsStep, mdObj)
            self._insertFunctionStep(super().convertTomoStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
            self._insertFunctionStep(self.extractParticlesStep, mdObj)
            self._insertFunctionStep(self.convertOutputStep, mdObj)
            self._insertFunctionStep(self.getProjectionsStep, mdObj)
        self._insertFunctionStep(self.createOutputStep, mdObjDict)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        coords = getattr(self, IN_COORDS).get()
        inCtfSet = getattr(self, IN_CTF).get()
        inTsSet = self.getTs()
        # Get the group ids and the emanDict to have the correspondence between the previous classes and
        # how EMAN will refer them
        uniqueTomoValsDict = getFromPresentObjects(coords, [Coordinate3D.TOMO_ID_ATTR, GROUP_ID])
        self.groupIds = uniqueTomoValsDict[GROUP_ID]
        self.emanDict = genEmanGrouping(self.groupIds)
        # Calculate the scale factor
        self.scaleFactor = coords.getSamplingRate() / inTsSet.getSamplingRate()
        # Generate the md object data units as a dict
        return self.genMdObjDict(inTsSet, inCtfSet, coords=coords)

    def writeData2JsonFileStep(self, mdObj):
        mode = 'a' if exists(mdObj.jsonFile) else 'w'
        coords2Json(mdObj, self.emanDict, self.groupIds, self.getBoxSize(), doFlipZ=self.doFlipZInTomo.get(), mode=mode)
        ts2Json(mdObj, mode='a')
        ctfTomo2Json(mdObj, self.sphAb, self.voltage, mode='a')

    def extractParticlesStep(self, mdObj):
        program = Plugin.getProgram('e2spt_extract.py')
        self.runJob(program, self._genExtractArgs(mdObj), cwd=self._getExtraPath())

    def convertOutputStep(self, mdObj):
        fName = self._getEmanFName(mdObj.tsId)
        # Unstack the 3d particles HDF stack into individual MRC files
        stack3d = join(PARTICLES_3D_DIR, fName)
        self.unstackParticles(stack3d)
        # Convert the 2d particles HDF stack into MRC for visualization purposes (compatibility with viewers)
        stack2d = join(PARTICLES_DIR, fName)
        inFile = join(PARTICLES_DIR, basename(stack2d))
        outFile = replaceExt(inFile, 'mrc')
        self.convertBetweenHdfAndMrc(inFile, outFile)
        # cleanPattern(hdfFile)

    def getProjectionsStep(self, mdObj):
        """load the corresponding HDF 2d stack file and read from the header the required data to create later a
        landmark model. To do that, the h5py module is used.
        """
        tsId = mdObj.tsId
        stack2dHdf = abspath(self._getExtraPath(PARTICLES_DIR, self._getEmanFName(mdObj.tsId)))
        eh = EmanHdf5Handler(stack2dHdf)
        projList = eh.getProjsFrom2dStack()
        # Add the required tsId as the first element of each sublist
        list(map(lambda sublist: sublist.insert(0, tsId), projList))  # More efficient than comprehension for huge
        # lists of lists, as expected
        self.projectionsDict[tsId] = projList

    def createOutputStep(self, mdObjDict):
        inCoordsPointer = getattr(self, IN_COORDS)
        inCoords = inCoordsPointer.get()
        sRate = inCoords.getSamplingRate()
        inTsPointer = getattr(self, IN_TS)
        inTs = inTsPointer.get()
        subtomoSet = EmanSetOfParticles.create(self._getPath(), template='emanParticles%s.sqlite')
        subtomoSet.setSamplingRate(sRate)
        subtomoSet.setCoordinates3D(inCoords)
        for tsId, mdObj in mdObjDict.items():
            coords = mdObj.coords
            tomoFile = mdObj.inTomo.getFileName()
            subtomoFiles = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.mrc' % tsId))
            stack2d = glob.glob(self._getExtraPath(PARTICLES_DIR, '%s*.hdf' % tsId))[0]
            stack3d = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.hdf' % tsId))[0]

            for coord, subtomoFile in zip(coords, subtomoFiles):
                subtomogram = EmanParticle()
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
                # Fill EmanParticle own attributes
                subtomogram.setInfoJson(mdObj.jsonFile)
                subtomogram.setTsHdf(mdObj.tsHdfName)
                subtomogram.setTomoHdf(mdObj.tomoHdfName)
                subtomogram.setStack2dHdf(stack2d)
                subtomogram.setStack3dHdf(stack3d)
                subtomoSet.append(subtomogram)

        # Generate the fiducial model (for data visualization purpose)
        fiducialModelGaps = SetOfLandmarkModels.create(self._getPath(), suffix='Gaps')
        fiducialModelGaps.copyInfo(inTs)
        fiducialModelGaps.setSetOfTiltSeries(inTsPointer)

        fiducialSize = round((inCoords.getBoxSize() / 2 * sRate) / 10)  # Box size is too large, a tenth of the radius
        for ts in inTs:
            tsId = ts.getTsId()
            landmarkModelGapsFilePath = self._getExtraPath(str(tsId) + "_gaps.sfid")
            landmarkModelGaps = LandmarkModel(tsId=tsId,
                                              tiltSeriesPointer=ts,
                                              fileName=landmarkModelGapsFilePath,
                                              modelName=None,
                                              size=fiducialSize,
                                              applyTSTransformation=False)
            landmarkModelGaps.setTiltSeries(ts)
            tsProjections = self.projectionsDict[tsId]
            for projection in tsProjections:
                tiltId = projection[1] + 1
                partId = projection[2] + 1
                xCoord = round(projection[3])
                yCoord = round(projection[4])
                landmarkModelGaps.addLandmark(xCoord, yCoord, tiltId, partId, 0, 0)
            fiducialModelGaps.append(landmarkModelGaps)

        # Define outputs and relations
        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: subtomoSet,
                               self._possibleOutputs.projected2DCoordinates.name: fiducialModelGaps})
        self._defineSourceRelation(inCoordsPointer, subtomoSet)
        self._defineSourceRelation(getattr(self, IN_CTF), subtomoSet)
        self._defineSourceRelation(inTsPointer, subtomoSet)
        self._defineSourceRelation(inTsPointer, fiducialModelGaps)

    # --------------------------- INFO functions --------------------------------
    # def _validate(self):
    #     valMsg = []
    #     if not self.inputTS.get():
    #         try:
    #             getNonInterpolatedTsFromRelations(self.inputCoordinates.get(), self)
    #         except:
    #             valMsg.append('Unable to go through the relations from the introduced coordinates to the '
    #                           'corresponding non-interpolated tilt series. Please introduce them using the '
    #                           'advanced parameter "Tilt series with alignment, non-interpolated."')
    #     return valMsg

    # --------------------------- UTILS functions ----------------------------------
    def genMdObjDict(self, inTsSet, inCtfSet, tomograms=None, coords=None):
        self.createInitEmanPrjDirs()

        # Considering the possibility of subsets, let's find the tsIds present in all the sets ob objects introduced,
        # which means the intersection of the tsId lists
        tomograms = tomograms if tomograms and not coords else coords.getPrecedents()
        presentCtfTsIds = set(getPresentTsIdsInSet(inCtfSet))
        presentTsSetTsIds = set(getPresentTsIdsInSet(inTsSet))
        presentTomoSetTsIds = set(getPresentTsIdsInSet(tomograms))
        presentTsIds = presentCtfTsIds & presentTsSetTsIds & presentTomoSetTsIds

        # Manage the TS, CTF tomo Series and tomograms
        tsIdsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inTsSet if ts.getTsId() in presentTsIds}
        ctfIdsDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inCtfSet if ctf.getTsId() in presentTsIds}
        tomoIdsDict = {tomo.getTsId(): tomo.clone() for tomo in tomograms if tomo.getTsId() in presentTsIds}

        # Get the required acquisition data
        self.sphAb = inTsSet.getAcquisition().getSphericalAberration()
        self.voltage = inTsSet.getAcquisition().getVoltage()

        # Split all the data into subunits (EmanMetaData objects) referred to the same tsId
        mdObjDict = {}
        for tomoId, ts in tsIdsDict.items():
            tomo = tomoIdsDict[tomoId]
            mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
                                             inTomo=tomo,
                                             ts=ts,
                                             ctf=ctfIdsDict[tomoId],
                                             coords=[coord.clone() for coord in
                                                     coords.iterCoordinates(volume=tomo)] if coords else None,
                                             jsonFile=genJsonFileName(self.getInfoDir(), tomoId))
        return mdObjDict

    def _genExtractArgs(self, mdObj):
        args = '%s ' % mdObj.tomoHdfName
        args += '--boxsz_unbin=%i ' % self.getBoxSize()
        args += '--maxtilt=%i ' % self.maxTilt.get()
        args += '--tltkeep=%.2f ' % self.tltKeep.get()
        args += '--padtwod=%.2f ' % self.paddingFactor.get()
        args += '--rmbeadthr=%.2f ' % self.rmThr.get()
        args += '--threads=%i ' % self.numberOfThreads.get()
        # args += '--newlabel=%s ' % mdObj.tsId
        args += '--append '
        args += '--verbose=9 '
        # if self.doSkipCtfCorrection.get():
        #     args += '--noctf '
        # if self.skip3dRec.get():
        #     args += '--skip3d '
        return args

    def unstackParticles(self, stackFile, outExt='mrc'):
        """Unstacks and coverts a list of stack files into separate particles"""
        program = Plugin.getProgram('e2proc3d.py')
        args = ' --unstacking'
        args += ' %s' % stackFile
        args += ' %s' % replaceExt(stackFile, outExt)
        self.runJob(program, args, cwd=self._getExtraPath())

    @staticmethod
    def _getEmanFName(tsId):
        return tsId + '__%s.hdf' % TOMOBOX
