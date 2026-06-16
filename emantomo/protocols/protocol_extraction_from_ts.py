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
from pyworkflow.utils import Message, replaceExt, removeExt, createLink
from emantomo.convert import coords2Json, ts2Json, ctfTomo2Json
from tomo.constants import TR_EMAN, SCIPION
from tomo.objects import SetOfLandmarkModels, LandmarkModel, Coordinate3D, SetOfSubTomograms, SetOfCoordinates3D, \
    SetOfMeshes


logger = logging.getLogger(__name__)


class outputObjects(Enum):
    subtomograms = EmanSetOfParticles
    projected2DCoordinates = SetOfLandmarkModels


class EmanProtTSExtraction(ProtEmantomoBase):
    """
    Extracts 2D tilt-series particle projections and reconstructs 3D subtomograms from cryo-electron tomography data.

    AI Generated:

    Extraction from Tilt Series (EmanProtTSExtraction) — User Manual
        Overview

        The Extraction from Tilt Series protocol is designed to generate
        subtomograms from cryo-electron tomography experiments by extracting
        particle projections directly from aligned tilt series and reconstructing
        them into individual 3D subvolumes. The protocol integrates tomography
        geometry, contrast transfer information, and particle coordinates in
        order to produce particle-centered reconstructions suitable for
        downstream subtomogram averaging, classification, or structural analysis.

        In practical cryo-ET workflows, this protocol is commonly used after
        tomogram reconstruction and particle picking. Users typically begin with
        either a set of 3D particle coordinates or an existing set of
        subtomograms associated with tomographic volumes. The protocol then
        revisits the original tilt-series data to extract the corresponding 2D
        particle images and reconstruct cleaner or standardized 3D particles
        using EMAN tomography tools.

        Biological Purpose and Typical Applications

        Subtomogram extraction is an essential step in structural analysis of
        macromolecular complexes in situ. By reconstructing individual particles
        from the original tilt images, the protocol preserves experimental
        geometry and allows subsequent averaging workflows to improve signal and
        reveal structural details.

        This approach is particularly valuable when studying heterogeneous
        assemblies, membrane-associated complexes, viral particles, ribosomes,
        cytoskeletal filaments, or macromolecular machinery embedded in crowded
        cellular environments. Because the extraction is tied directly to the
        original tilt-series geometry, the resulting particles maintain spatial
        consistency with the tomographic acquisition.

        Inputs and Experimental Consistency

        The protocol requires three principal inputs: particle coordinates or
        existing subtomograms, aligned tilt series, and corresponding CTF tilt
        series information. These datasets must refer to the same tomographic
        acquisitions. Correct correspondence between coordinates, tilt series,
        and CTF estimations is essential for successful reconstruction.

        The tilt series should already contain alignment information and should
        normally represent non-interpolated aligned data. Using improperly
        aligned tilt series may introduce geometric inconsistencies that reduce
        reconstruction quality or generate inaccurate particle orientations.

        Coordinates may originate from manual picking, template matching,
        segmentation workflows, or previous subtomogram processing steps. When
        subtomograms are provided as input, the protocol reuses their associated
        coordinates to maintain spatial continuity.

        Particle Size and Binning

        One of the most important parameters is the extraction box size, which
        determines the physical region reconstructed around each particle. The
        selected box should fully contain the biological structure of interest
        together with enough surrounding contextual density to avoid edge
        artifacts.

        The protocol also allows particle binning during extraction. Binning
        reduces image dimensions and computational cost while increasing signal-
        to-noise ratio. This is especially useful during exploratory analyses or
        when processing large datasets. Smaller binned particles are faster to
        reconstruct and easier to align, although excessive binning may remove
        high-resolution information.

        In biological practice, moderate binning is commonly used during initial
        subtomogram averaging iterations, whereas final refinements often employ
        lower binning or fully unbinned particles.

        Tilt Selection and Data Quality

        The protocol includes options to restrict the angular range of the tilt
        images used during reconstruction. Excluding extremely tilted images may
        improve reconstruction stability because very high tilts usually contain
        lower signal quality, stronger distortions, and increased effective
        sample thickness.

        Users may also retain only a fraction of the best tilt images according
        to reconstruction quality criteria. This selective approach can improve
        particle consistency, particularly in datasets affected by motion,
        contamination, charging, or variable acquisition quality.

        From a biological perspective, careful tilt selection can substantially
        improve the interpretability of flexible or low-contrast complexes.

        Removal of Fiducials and High-Contrast Features

        Cryo-electron tomography datasets frequently contain fiducial gold beads
        or other highly contrasted contaminants used during tilt-series
        alignment. The protocol provides a contrast-based filtering mechanism to
        exclude problematic 2D particle projections that overlap with such
        features.

        This option is especially useful when particles are located close to
        fiducials or dense cellular structures. Removing contaminated
        projections can significantly improve the quality of the reconstructed
        subtomograms and reduce artifacts during averaging.

        However, overly aggressive filtering may discard valid particle data and
        reduce particle completeness. Biological users should therefore inspect
        reconstruction quality carefully when enabling this feature.

        Padding and Particle Isolation

        The extraction workflow optionally supports additional padding around
        particle projections. Padding may help when particles are deeply buried
        within crowded environments or neighboring densities interfere with
        reconstruction quality.

        Increasing padding improves contextual isolation of the particle but also
        increases memory consumption and computational cost. In practice, modest
        padding values are often sufficient for most cellular cryo-ET datasets.

        Coordinate Consistency and Geometric Accuracy

        During processing, the protocol preserves the relationship between
        particle coordinates, tomograms, and reconstructed subtomograms. This is
        particularly important for workflows involving particle transformations,
        iterative refinements, or coordinate-based biological interpretation.

        The protocol also manages coordinate scaling associated with binning and
        tilt-series sampling differences. This ensures that reconstructed
        particles remain geometrically consistent with the original tomographic
        coordinate system.

        In addition, optional inversion of the tomogram Z axis is available to
        accommodate differences between reconstruction conventions used by
        various tomography software packages.

        Landmark and Projection Information

        Besides producing reconstructed subtomograms, the protocol generates
        landmark-style projection information derived from the extracted 2D tilt
        particles. These landmarks are mainly intended for visualization,
        validation, and inspection purposes.

        Such projection models can help users verify that particle extraction
        correctly follows the experimental geometry across the tilt series and
        can reveal problematic particles, missing projections, or alignment
        inconsistencies.

        Outputs and Downstream Analysis

        The primary output is a set of reconstructed subtomograms suitable for
        downstream cryo-ET analysis. These particles can subsequently be used
        for subtomogram alignment, classification, averaging, or structural
        interpretation within larger tomography workflows.

        Each output particle preserves its relationship to the originating tilt
        series and tomogram, facilitating traceability and reproducibility in
        large-scale studies.

        Additional visualization-oriented outputs containing projection landmark
        information are also generated to support quality assessment and data
        inspection.

        Practical Recommendations

        For routine cryo-ET processing, users should begin with conservative
        binning and moderate box sizes while validating reconstruction quality
        visually. If particles are located in crowded cellular environments,
        enabling padding and carefully tuning fiducial filtering often improves
        downstream averaging.

        When processing heterogeneous or flexible assemblies, maintaining
        accurate coordinate consistency is especially important because small
        geometric mismatches can strongly affect alignment quality in later
        subtomogram averaging steps.

        Biological interpretation should always consider the missing wedge,
        variable tilt-image quality, and the possibility of reconstruction
        artifacts introduced by incomplete angular coverage.

        Final Perspective

        The Extraction from Tilt Series protocol provides a complete bridge
        between tomographic particle coordinates and reconstructed subtomograms.
        By combining experimental geometry, tilt-series alignment, CTF
        information, and particle extraction into a unified workflow, the
        protocol enables robust generation of subtomograms suitable for modern
        in situ structural biology analyses.

        In cryo-electron tomography workflows, accurate extraction is a critical
        prerequisite for meaningful subtomogram averaging and reliable biological
        interpretation. Careful selection of extraction parameters, tilt-image
        quality, and particle geometry therefore has a direct impact on the
        final structural results.
    """

    _label = 'Extraction from TS'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.groupIds = None
        self.emanDict = None
        self.sphAb = None
        self.voltage = None
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
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            createPrjId = self._insertFunctionStep(self.createExtractionEmanPrjStep, mdObj)
            convertInId = self._insertFunctionStep(self.convertTsStep, mdObj, prerequisites=createPrjId)
            writeJsonId = self._insertFunctionStep(self.writeData2JsonFileStep, mdObj, prerequisites=convertInId)
            extractId = self._insertFunctionStep(self.extractParticlesStep, mdObj, prerequisites=writeJsonId)
            convertOutId = self._insertFunctionStep(self.convertOutputStep, mdObj, prerequisites=extractId)
            createOutputId = self._insertFunctionStep(self.createOutputStep, mdObj, prerequisites=convertOutId)
            pIds.extend([createPrjId, convertInId, writeJsonId, extractId, convertOutId, createOutputId])
        self._insertFunctionStep(self._closeOutputSet, prerequisites=pIds)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        logger.info('Initializing...')
        self.createInitEmanPrjDirs()
        inParticles = getattr(self, IN_SUBTOMOS).get()
        inCtfSet = getattr(self, IN_CTF).get()
        inTsSet = self.getAttrib(IN_TS)
        binFactor = self.shrink.get()
        typeInParticles = type(inParticles)
        if typeInParticles in [SetOfCoordinates3D, SetOfMeshes]:  # Extraction from coords
            coords = inParticles
        else:  # Extraction from particles
            coords = inParticles.getCoordinates3D()
        self.currentBoxSize = self.getAttrib(IN_BOXSIZE) / binFactor
        self.currentSRate = inTsSet.getSamplingRate() * binFactor

        # Get the group ids and the emanDict to have the correspondence between the previous classes and
        # how EMAN will refer them
        uniqueTomoValsDict = getFromPresentObjects(coords, [Coordinate3D.TOMO_ID_ATTR, GROUP_ID])
        self.groupIds = uniqueTomoValsDict[GROUP_ID]
        self.emanDict = genEmanGrouping(self.groupIds)
        # Calculate the scale factor for coords and shift scaling: according to EMAN's behavior, it is the rate between
        # the sampling rate of the input coordinates and the TS unbinned sampling rate multiplied by the binning factor
        # introduced
        self.scaleFactor = inParticles.getSamplingRate() / (inTsSet.getSamplingRate() * self.shrink.get())
        # The fiducial size is the diameter in nm
        self.fiducialSize = 0.1 * self.getAttrib(IN_BOXSIZE) * self.getAttrib(IN_TS).getSamplingRate() / 2
        # Generate the md object data units as a dict
        return self.genMdObjDict(inTsSet, inCtfSet, coords)

    def createExtractionEmanPrjStep(self, mdObj):
        logger.info(f'Creating the project for tsIf = {mdObj.tsId}...')
        # Create project dir structure
        tomogramsDir = self.getTomogramsDir()
        # TS must be in HDF (in EMAN's native code it is searched in hardcoded HDF format based on the tomogram
        # basename), so it will be converted later. However, the tomograms can be linked in MRC without any problem
        tomoFName = mdObj.inTomo.getFileName()
        createLink(tomoFName, join(tomogramsDir, basename(tomoFName)))

    def writeData2JsonFileStep(self, mdObj):
        logger.info(f'Writing the json files for TS {mdObj.tsId}...')
        mode = 'a' if exists(mdObj.jsonFile) else 'w'
        coords2Json(mdObj, self.emanDict, self.groupIds, self.getBoxSize(), doFlipZ=self.doFlipZInTomo.get(),
                    mode=mode)
        ts2Json(mdObj, mode='a')
        ctfTomo2Json(mdObj, self.sphAb, self.voltage, mode='a')

    def extractParticlesStep(self, mdObj):
        logger.info(f'Extracting the particles from TS {mdObj.tsId}...')
        program = Plugin.getProgram('e2spt_extract.py')
        self.runJob(program, self._genExtractArgs(mdObj), cwd=self._getExtraPath())

    def convertOutputStep(self, mdObj):
        logger.info(f'Converting and unstacking the 3D particles from TS {mdObj.tsId} into MRC...')
        fName = self._getEmanFName(mdObj.tsId)
        # Unstack the 3d particles HDF stack into individual MRC files
        stack3d = join(PARTICLES_3D_DIR, fName)
        self.unstackParticles(stack3d)
        fixVolume(glob.glob(join(self.getStack3dDir(), f"{removeExt(fName)}*.mrc")))
        # # Convert the 2d particles HDF stack into MRC for visualization purposes (compatibility with viewers)
        # stack2d = join(PARTICLES_DIR, fName)
        # inFile = join(PARTICLES_DIR, basename(stack2d))
        # outFile = replaceExt(inFile, 'mrc')
        # self.convertBetweenHdfAndMrc(inFile, outFile)

    def createOutputStep(self, mdObj):
        tsId = mdObj.tsId
        logger.info(f'Creating the outputs corresponding to TS {tsId}...')
        subtomoSet = self.getOutSetOfParticles()
        # Generate the fiducial model (for data visualization purpose)
        fiducialModelGaps = self.getOutSetIfFiducials()

        # Get the projections
        tsProjections = self.getProjections(mdObj)
        landmarkModelGaps = self.createLandmarkModelGaps(mdObj.ts)

        # Because the order is changed when the particles are processed with EMAN, if they have transformation matrix,
        # there will be a data mismatch. To avoid it, The EMAN coordinates will be read from the 3d HDF stack, then the
        # corresponding coordinate will be located among the input set, so the transformation is correctly matched, and
        # finally, the corresponding output unstacked MRC file will be located and added to the output particle.
        emanCoords = self.getEmanCoordsFromHdfStack(tsId)

        _, _, nImgs = mdObj.ts.getDimensions()
        inCoords = mdObj.particles if self.inParticles else mdObj.coords
        subtomoFiles = sorted(glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.mrc' % tsId)))
        stack2d = glob.glob(self._getExtraPath(PARTICLES_DIR, '%s*.hdf' % tsId))[0]
        stack3d = glob.glob(self._getExtraPath(PARTICLES_3D_DIR, '%s*.hdf' % tsId))[0]

        nTotalParticles = mdObj.processingInd
        particleCounter = 0
        for emanCoord, subtomoFile in zip(emanCoords, subtomoFiles):
            self.fillLandmarkModel(landmarkModelGaps, tsProjections, particleCounter, nImgs)
            subtomogram = EmanParticle()
            transform = Transform()
            subtomogram.setFileName(subtomoFile)
            subtomogram.setLocation(particleCounter, subtomoFile)  # The index is stored to be used as
            # the position of the particle in the 3d HDF stack (to handle subsets - LST file gen.)
            subtomogram.setSamplingRate(self.currentSRate)
            subtomogram.setVolName(mdObj.tsId)
            scipionCoord, inCoords = self.getEmanMatchingCoord(emanCoord, inCoords)
            subtomogram.setCoordinate3D(scipionCoord)
            subtomogram.setTsId(scipionCoord.getTomoId())
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
            subtomogram.setInfoJson(mdObj.jsonFile)
            subtomogram.setTsHdf(self._getExtraPath(mdObj.tsHdfName))
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
    def genMdObjDict(self, inTsSet, inCtfSet, coords):
        self.createInitEmanPrjDirs()
        if type(coords) is SetOfSubTomograms:
            coords = coords.getCoordinates3D()
        # Considering the possibility of subsets, let's find the tsIds present in all the sets ob objects introduced,
        # which means the intersection of the tsId lists
        presentCtfTsIds = set(getPresentTsIdsInSet(inCtfSet))
        self.info("TsIds present in the CTF tomo series are: %s" % presentCtfTsIds)
        presentTsSetTsIds = set(getPresentTsIdsInSet(inTsSet))
        self.info("TsIds present in the tilt series are: %s" % presentTsSetTsIds)
        presentTomoSetTsIds = set(coords.getUniqueValues(Coordinate3D.TOMO_ID_ATTR))
        self.info("TsIds present in the coordinates are: %s" % presentTomoSetTsIds)
        presentTsIds = presentCtfTsIds & presentTsSetTsIds & presentTomoSetTsIds
        # The tomograms are obtained as the coordinates precedents. Operating this way, the code considers the case of
        # subsets of coordinates
        presentTomograms = getPresentPrecedents(coords, presentTsIds)

        # Manage the TS, CTF tomo Series and tomograms
        tsIdsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inTsSet if ts.getTsId() in presentTsIds}
        ctfIdsDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inCtfSet if ctf.getTsId() in presentTsIds}
        tomoIdsDict = {tomo.getTsId(): tomo for tomo in presentTomograms}

        # Get the required acquisition data
        self.sphAb = inTsSet.getAcquisition().getSphericalAberration()
        self.voltage = inTsSet.getAcquisition().getVoltage()

        # Split all the data into subunits (EmanMetaData objects) referred to the same tsId
        mdObjDict = {}
        procInd = 0
        for tomoId, ts in sorted(tsIdsDict.items()):
            # Sorting the items before the iteration based on the tsId (key) will ensure the alphabetical order in
            # all the data generated, preventing data mismatching between the alignment files generated by EMAN in
            # posterior refinements and its corresponding version in Scipion format
            tomo = tomoIdsDict[tomoId]
            hdfFileBaseName = f'{tomoId}.hdf'
            iCoords = [coord.clone() for coord in coords.iterCoordinates(volume=tomo)] if coords else None
            iParticles = [particle.clone() for particle in self.inParticles.iterSubtomos(volume=tomo)] if \
                self.inParticles else None
            mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
                                             inTomo=tomo,
                                             tomoHdfName=join(TOMOGRAMS_DIR, hdfFileBaseName),
                                             ts=ts,
                                             tsHdfName=join(TS_DIR, hdfFileBaseName),
                                             ctf=ctfIdsDict[tomoId],
                                             coords=iCoords,
                                             particles=iParticles,
                                             jsonFile=genJsonFileName(self.getInfoDir(), tomoId),
                                             processingInd=procInd)
            procInd += len(iCoords)
        if not mdObjDict:
            raise Exception("There isn't any common tilt series among the coordinates, the CTF tomo series and "
                            "the tilt series chosen.")
        return dict(sorted(mdObjDict.items()))

    def _genExtractArgs(self, mdObj):
        # align2dFile = self.getNewAliFile(is3d=False)
        align3dFile = self.getNewAliFile()
        tomoFiles = join(TOMOGRAMS_DIR, basename(mdObj.inTomo.getFileName()))
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

    def getProjections(self, mdObj):
        """load the corresponding HDF 2d stack file and read from the header the required data to create later a
        landmark model. To do that, the h5py module is used.
        """
        tsId = mdObj.tsId
        stack2dHdf = abspath(self._getExtraPath(PARTICLES_DIR, self._getEmanFName(mdObj.tsId)))
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

