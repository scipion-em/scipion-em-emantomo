# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
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
import copy
import enum
import glob
import logging
import typing
from os.path import abspath, join, basename
from typing import Tuple
import numpy as np
from emantomo import Plugin
from emantomo.constants import PROC_NORMALIZE
from pyworkflow.mapper.sqlite import ID
from pwem.objects import Transform
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, yellowStr, cyanStr
from pyworkflow.utils.path import moveFile, cleanPath, cleanPattern, makePath
from pwem.protocols import EMProtocol
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.constants import BOTTOM_LEFT_CORNER, TR_SCIPION, SCIPION
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfCoordinates3D, SetOfSubTomograms, SubTomogram, TomoAcquisition, Coordinate3D, Tomogram

logger = logging.getLogger(__name__)

SAME_AS_PICKING = 0
OTHER = 1


class OutputExtraction(enum.Enum):
    subtomograms = SetOfSubTomograms


class EmanProtTomoExtraction(EMProtocol, ProtTomoBase):
    """
    Extracts subtomograms directly from tomographic volumes using 3D particle
    coordinates or previously generated subtomograms.

    AI Generated:

    Subtomograms Extraction from Tomogram (EmanProtTomoExtraction) — User Manual
        Overview

        The Subtomograms Extraction from Tomogram protocol is designed to
        generate individual 3D particle volumes directly from reconstructed
        tomograms. Its primary purpose is to isolate regions of interest around
        biological particles so they can later be analyzed through subtomogram
        averaging, classification, structural interpretation, or visualization.

        In cryo-electron tomography workflows, this protocol is typically used
        after particle picking or coordinate generation. Users provide either a
        set of 3D coordinates or an existing set of subtomograms associated with
        those coordinates, and the protocol extracts cubic particle volumes from
        one or more tomograms. The resulting subtomograms preserve the spatial
        relationship with the original tomographic data and are suitable for
        downstream in situ structural biology analysis.

        Biological Context and Typical Applications

        Subtomogram extraction is a fundamental step in cryo-ET studies focused
        on macromolecular complexes inside cells, organelles, membranes, or
        viral assemblies. By isolating individual particles from larger
        tomograms, the protocol enables detailed structural comparison and
        averaging while maintaining the biological context of the specimen.

        Typical applications include ribosome analysis inside cells, membrane
        protein organization, viral capsid studies, cytoskeletal assemblies, and
        characterization of flexible or heterogeneous molecular complexes. The
        protocol is equally suitable for exploratory biological analysis and for
        large-scale high-throughput tomography pipelines.

        Input Coordinates and Tomogram Sources

        The protocol accepts either direct 3D coordinates or previously computed
        subtomograms. When subtomograms are provided, their associated
        coordinates are reused automatically to maintain geometric consistency
        with the original tomographic dataset.

        By default, extraction is performed using the same tomograms associated
        with the picking process. However, the protocol also supports extraction
        from alternative tomograms. This option is especially useful when users
        wish to extract particles from tomograms reconstructed with different
        filtering strategies, denoising methods, alignment refinements, or voxel
        sizes.

        In biological practice, using alternative tomograms can significantly
        improve particle quality. For example, coordinates identified in
        aggressively filtered tomograms may later be re-extracted from cleaner or
        higher-resolution reconstructions for final subtomogram averaging.

        Coordinate and Sampling Consistency

        A critical requirement is that coordinates and tomograms correspond to
        the same biological acquisition and geometric reference frame. The
        protocol automatically evaluates tomogram identifiers to ensure proper
        matching between coordinates and extraction volumes.

        Because coordinates and tomograms may have different sampling rates, the
        protocol accounts for voxel-size scaling to preserve positional accuracy.
        This is particularly important in workflows where coordinates originate
        from binned tomograms while extraction is performed from higher-
        resolution data.

        From a biological perspective, maintaining coordinate precision is
        essential because even small geometric inconsistencies can compromise
        downstream subtomogram alignment and averaging.

        Particle Box Size

        The extraction box size defines the cubic region reconstructed around
        each coordinate. Choosing an appropriate box size is biologically
        important because the box must fully contain the particle together with
        enough surrounding density to avoid edge artifacts during later
        processing.

        Small box sizes may truncate flexible domains or membrane regions,
        whereas excessively large boxes increase computational cost and may
        introduce unrelated surrounding densities. In practice, users often begin
        with moderate box sizes and adjust them according to the expected
        dimensions of the molecular complex.

        Contrast Inversion and Image Conventions

        Different cryo-EM software packages use different intensity conventions.
        The protocol therefore supports optional contrast inversion to ensure
        compatibility with downstream processing tools.

        In many cryo-ET workflows, particles are expected to appear as white
        density against a dark background. Correct contrast orientation is
        important because some reconstruction and alignment algorithms assume a
        particular density convention. Using the wrong convention can negatively
        affect subtomogram alignment and interpretation.

        Subtomogram Normalization

        The protocol optionally normalizes subtomograms before output generation.
        Normalization improves consistency between particles by reducing
        variations in background intensity and overall contrast.

        Two normalization strategies are available. Standard normalization is
        appropriate for many datasets, whereas edge-based normalization is
        particularly useful when particles are surrounded by solvent regions and
        are relatively isolated from neighboring structures.

        Biologically crowded environments, such as intracellular tomograms or
        membrane-embedded systems, may benefit less from edge-based approaches
        because neighboring densities can strongly influence background
        estimation.

        Preservation of Spatial Transformations

        During extraction, the protocol preserves positional and transformation
        information associated with the original coordinates or subtomograms.
        This ensures compatibility with downstream alignment workflows and allows
        reconstructed particles to remain connected to their original biological
        context.

        Coordinate adjustments introduced during extraction are carefully managed
        so that subtomogram orientations and positional relationships remain
        geometrically meaningful. This is especially important in iterative
        subtomogram refinement workflows and in studies involving spatial
        organization within cells.

        Outputs and Downstream Analysis

        The primary output is a set of extracted subtomograms associated with
        their corresponding coordinates and tomographic metadata. These particles
        can be directly used in subtomogram averaging, classification,
        visualization, or structural interpretation pipelines.

        Because the protocol preserves acquisition information and coordinate
        relationships, the resulting subtomograms remain fully traceable to their
        originating tomograms and biological locations.

        Practical Recommendations

        In routine cryo-ET processing, users should first verify that coordinates
        and tomograms correspond to the same acquisition geometry and voxel size.
        Visual inspection of several extracted particles is strongly recommended
        before launching large averaging workflows.

        Moderate normalization is usually beneficial for alignment stability,
        while contrast inversion should match the expectations of downstream
        processing software. When extracting from alternative tomograms, users
        should ensure that coordinate scaling and orientation remain consistent.

        For membrane proteins, crowded intracellular complexes, or highly
        heterogeneous assemblies, careful selection of box size and normalization
        strategy can substantially improve subtomogram quality and downstream
        structural interpretation.

        Final Perspective

        The Subtomograms Extraction from Tomogram protocol provides a robust and
        biologically oriented bridge between tomographic coordinates and
        particle-centered 3D volumes. By preserving geometric consistency,
        acquisition metadata, and particle transformations, the protocol enables
        reliable preparation of subtomograms for modern in situ structural
        biology workflows.

        In cryo-electron tomography, accurate particle extraction is a critical
        prerequisite for meaningful subtomogram averaging and interpretation.
        Careful handling of coordinate precision, tomogram selection, and image
        normalization therefore has a direct impact on the quality of the final
        structural results.
    """
    _label = 'Subtomograms extraction from tomogram'
    _possibleOutputs = OutputExtraction
    stepsExecutionMode = STEPS_PARALLEL
    OUTPUT_PREFIX = _possibleOutputs.subtomograms.name

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomosDict = None
        self.scaleFactor = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputCoordinates', params.PointerParam, label="Coordinates/Subtomograms", important=True,
                      pointerClass=[SetOfCoordinates3D, SetOfSubTomograms],
                      help='Choose coordinates or subtomograms derived from 3d coordinates.')

        form.addParam('tomoSource', params.EnumParam,
                      choices=['same as picking', 'other'],
                      default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted '
                           'from the tomogram used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different tomogram to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'tomogram and coordinates are related ')

        form.addParam('inputTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      condition='tomoSource != %s' % SAME_AS_PICKING,
                      label='Input tomogram',
                      help='Select the tomogram from which to extract.')

        form.addParam('boxSize', params.FloatParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size. '
                           'The wizard will select the box size considering the sampling rate ratio between the '
                           'introduced coordinates and the tomograms that will br used for the extraction.')

        form.addSection(label='Preprocess')
        form.addParam('doInvert', params.BooleanParam, default=True,
                      label='Invert contrast?',
                      help='Invert the contrast if your tomogram is black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize subtomogram?',
                      help='Normalization processor applied to subtomograms before extraction.')

        form.addParam('normproc', params.EnumParam,
                      choices=['normalize', 'normalize.edgemean'],
                      label='Normalize method',
                      condition='doNormalize',
                      default=PROC_NORMALIZE,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Use normalize.edgemean if the particles have a clear solvent background '
                           '(i.e., they are not part of a larger complex or embeded in a membrane)')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tomosDict.keys():
            pId = self._insertFunctionStep(self.writeSetOfCoordinates3D, tsId,
                                           prerequisites=[],
                                           needsGPU=False)
            pId = self._insertFunctionStep(self.extractParticles, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            pId = self._insertFunctionStep(self.convertOutput, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            pId = self._insertFunctionStep(self.createOutputStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            closeSetStepDeps.append(pId)
        self._insertFunctionStep(self.closeOutputSet,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        matchingTsIds, nonMatchingTsIds = self._getMatchingTsIds()
        if len(nonMatchingTsIds) > 0:
            logger.warning(yellowStr(f'Some tsIds do not match: {nonMatchingTsIds}'))
        self.tomosDict = {tsId: tomo.clone() for tomo in self.getInputTomograms()
                          if (tsId := tomo.getTsId()) in matchingTsIds}
        inCoords = self._getSetOfCoordinates()
        # Calculate the scale factor
        firstTomo = self.getInputTomograms().getFirstItem()
        samplingRateInput = inCoords.getSamplingRate()
        samplingRateTomo = firstTomo.getSamplingRate()
        self.scaleFactor = samplingRateInput / samplingRateTomo

    def writeSetOfCoordinates3D(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - Writing the coordinates of tomogram into EMAN format..."))
        inputSet = self._getSetOfCoordinates()
        tomo = self.tomosDict[tsId]

        # Calculate ratio/factor
        samplingRateCoord = inputSet.getSamplingRate()
        samplingRateTomo = self.getInputTomograms().getFirstItem().getSamplingRate()
        scale = samplingRateCoord / samplingRateTomo

        # Iterate in order based on tomogram/Ts id
        with open(self._getCoordsFile(tsId), "w") as emanCoordFile:
            for item in inputSet.iterCoordinates(volume=tomo, orderBy=[Coordinate3D.TOMO_ID_ATTR, ID]):
                coord3D = self._getCoordinateFromItem(item)
                xScaled = coord3D.getX(BOTTOM_LEFT_CORNER) * scale
                yScaled = coord3D.getY(BOTTOM_LEFT_CORNER) * scale
                zScaled = coord3D.getZ(BOTTOM_LEFT_CORNER) * scale
                emanCoordFile.write(f"{xScaled:.3f}\t{yScaled:.3f}\t{zScaled:.3f}\n")

    def extractParticles(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - Extracting the particles from tomogram..."))
        tomo = self.tomosDict[tsId]
        args = '%s ' % abspath(tomo.getFileName())
        args += "--coords %s --boxsize %i" % (f'{tsId}.coords', self.boxSize.get())
        if self.doInvert:
            args += ' --invert'
        if self.doNormalize:
            args += ' --normproc %s' % self.getEnumText('normproc')
        # args += ' --cshrink %i' % (samplingRateTomo / samplingRateCoord)

        program = Plugin.getProgram('e2spt_boxer_old.py')
        self.runJob(program, args, cwd=self._getExtraPath())
        moveFile(self._getExtraPath(join('sptboxer_01', 'basename.hdf')),
                 self._getExtraPath(f'{tsId}.hdf'))
        cleanPath(self._getExtraPath("sptboxer_01"))

    def convertOutput(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - Unstacking the particles extracted from tomogram..."))
        resultsDir = self._getUnstackedResultsDir(tsId)
        makePath(resultsDir)
        coordsFile = self._getCoordsFile(tsId)
        moveFile(coordsFile, join(resultsDir, basename(coordsFile)))
        program = Plugin.getProgram('e2proc3d.py')
        hdfFile = self._getOutHdfCoordsStack(tsId)
        args = ' --unstacking'
        args += ' %s' % hdfFile
        args += ' %s' % join(resultsDir, pwutils.replaceBaseExt(hdfFile, 'mrc'))
        args += ' --apix %.3f' % self.getOutputSamplingRate()
        self.runJob(program, args)
        cleanPattern(hdfFile)

    def createOutputStep(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - Registering the results..."))
        tomo = self.tomosDict[tsId]
        self._registerOutput(tomo)

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, tomo: Tomogram):
        with self._lock:
            outSubtomos = self._createOutputSet()
            inCoords = self._getSetOfCoordinates()
            coordSet = [item.clone() for item in inCoords.iterCoordinates(volume=tomo)]
            self.readSetOfSubTomograms(tomo, outSubtomos, coordSet, self.scaleFactor)

    def closeOutputSet(self):
        self._closeOutputSet()
        outSubtomos = getattr(self, self._possibleOutputs.subtomograms.name, None)
        if outSubtomos.getSize() == 0:
            raise Exception("No particles were extracted!")

    # --------------------------- INFO functions --------------------------------
    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            msg = ("A total of %s subtomograms of size %s were extracted"
                   % (str(self._getSetOfCoordinates().getSize()), self.boxSize.get()))
            if self._tomosOther():
                msg += (" from another set of tomograms: %s"
                        % self.getObjectTag('inputTomogram'))
            msg += " using coordinates from %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")
        if self.doInvert:
            methodsMsgs.append("Inverted contrast on images.")
        if self.doNormalize:
            methodsMsgs.append("Particles were normalised. Using normalization method %s"
                               % self.getEnumText('normproc'))
        return methodsMsgs

    def _summary(self):
        summary = ["Tomogram source: *%s*" % self.getEnumText("tomoSource")]
        if self.getOutputsSize() >= 1:
            summary.append("Particle box size: *%s*" % self.boxSize.get())
            summary.append("Subtomogram extracted: *%s*" % self._getSetOfCoordinates().getSize())
        else:
            summary.append("Output subtomograms not ready yet.")
        if self.doInvert:
            summary.append('*Contrast was inverted.*')
        return summary

    def _validate(self):
        errors = []
        matches, _ = self._getMatchingTsIds()
        if len(matches) == 0:
            errors.append("Cannot relate coordinates tsIds and new tomograms tsIds.")
        return errors

    # --------------------------- UTILS functions ----------------------------------
    def _getSetOfCoordinates(self):
        if self._isInputASetOfSubtomograms():
            return self.inputCoordinates.get().getCoordinates3D()
        else:
            return self.inputCoordinates.get()

    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self._getSetOfCoordinates().getPrecedents()
        else:
            return self.inputTomograms.get()

    def _isInputASetOfSubtomograms(self):
        """ returns true if the input is a set of subtomograms"""
        return isinstance(self.inputCoordinates.get(), SetOfSubTomograms)

    def _tomosOther(self) -> bool:
        """ Return True if other tomograms are used for extract. """
        return self.tomoSource.get() == OTHER

    def _getMatchingTsIds(self) -> Tuple[typing.Set, typing.Set]:
        tomograms = self.getInputTomograms()
        tomoIds = tomograms.getTSIds()
        coords = self._getSetOfCoordinates()
        coordsTomoIds = coords.getTSIds()
        matches = set(tomoIds) & set(coordsTomoIds)
        notMatches = set(tomoIds) ^ set(coordsTomoIds)
        return matches, notMatches

    @staticmethod
    def _getCoordinateFromItem(item) -> Coordinate3D:
        """ Returns the coordinate 3d either because the item is the Coordinate or is a subtomogram"""
        if isinstance(item, Coordinate3D):
            return item
        else:
            return item.getCoordinate3D()

    @staticmethod
    def _getMatrixFromItem(item):
        """ Returns the matrix of the subtomograms otherwise the matrix of the coordinate"""
        if isinstance(item, Coordinate3D):
            return item.getMatrix()
        else:
            return item.getTransform().getMatrix()

    def readSetOfSubTomograms(self,
                              tomo: Tomogram,
                              outputSubTomogramsSet: SetOfSubTomograms,
                              inputSet: list,
                              scaleFactor: int):
        """
        Populates the set of subtomograms

        :param tomo: Tomogram
        :param outputSubTomogramsSet: output set of subtomograms
        :param inputSet: Subtomograms or 3D coordinates set
        :param scaleFactor: factor between the inputSet and the tomogram
        """
        tsId = tomo.getTsId()
        logger.info(cyanStr('Registering the subtomograms from tomogram %s' % tsId))
        subtomoFileList = sorted(glob.glob(join(self._getUnstackedResultsDir(tsId), '*.mrc')))
        for idx, subtomoFile in enumerate(subtomoFileList):
            # logger.info("Registering subtomogram %s - %s" % (counter, subtomoFile))
            subtomogram = SubTomogram()
            transform = Transform()
            subtomogram.setLocation(subtomoFile)
            currentItem = inputSet[idx]
            coord = self._getCoordinateFromItem(currentItem)

            #########################################################################################
            # EMAN work with the coordinates as integers, so we round them and add the decimal part
            # as shifts to the transformation matrix, so all the information is preserved
            origX = coord.getX(SCIPION)
            origY = coord.getY(SCIPION)
            origZ = coord.getZ(SCIPION)
            roundCoordX = round(origX)
            roundCoordY = round(origY)
            roundCoordZ = round(origZ)
            coord.setX(roundCoordX, SCIPION)
            coord.setY(roundCoordY, SCIPION)
            coord.setZ(roundCoordZ, SCIPION)
            subtomogram.setCoordinate3D(coord)
            trMatrix = copy.copy(self._getMatrixFromItem(currentItem))
            shifts = np.array([trMatrix[0, 3], trMatrix[1, 3], trMatrix[2, 3]])
            scaledShifts = scaleFactor * shifts
            trMatrix[0, 3] = scaledShifts[0] + roundCoordX - origX
            trMatrix[1, 3] = scaledShifts[1] + roundCoordY - origY
            trMatrix[2, 3] = scaledShifts[2] + roundCoordY - origY
            #########################################################################################

            # transform.setMatrix(scaleTrMatrixShifts(trMatrix, scaleFactor))
            transform.setMatrix(trMatrix)
            subtomogram.setTransform(transform, convention=TR_SCIPION)
            subtomogram.setVolName(tomo.getFileName())
            outputSubTomogramsSet.append(subtomogram)

    def getOutputSamplingRate(self) -> float:
        return self.getInputTomograms().getSamplingRate()

    def _getCoordsFile(self, tsId: str) -> str:
        return self._getExtraPath(f'{tsId}.coords')

    def _getOutHdfCoordsStack(self, tsId) -> str:
        return self._getExtraPath(f'{tsId}.hdf')

    def _createOutputSet(self) -> SetOfSubTomograms:
        outSubtomos = getattr(self, self._possibleOutputs.subtomograms.name, None)
        if outSubtomos:
            outSubtomos.enableAppend()
        else:
            outSubtomos = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
            outSubtomos.setSamplingRate(self.getOutputSamplingRate())
            outSubtomos.setCoordinates3D(self._getSetOfCoordinates())
            acquisition = TomoAcquisition()
            firstTomo = self.getInputTomograms().getFirstItem()
            acquisition.copyInfo(firstTomo.getAcquisition())
            outSubtomos.setAcquisition(acquisition)
            outSubtomos.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{OutputExtraction.subtomograms.name: outSubtomos})
            self._defineSourceRelation(self._getSetOfCoordinates(), outSubtomos)
        return outSubtomos

    def _getUnstackedResultsDir(self, tsId: str) -> str:
        return self._getExtraPath(tsId)
