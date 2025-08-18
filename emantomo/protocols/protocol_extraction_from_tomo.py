# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
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
from os.path import abspath, join
from typing import Set, Tuple

from emantomo import Plugin
from emantomo.constants import PROC_NORMALIZE
from pyworkflow.mapper.sqlite import ID
from pwem.objects import Transform
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message, yellowStr, cyanStr
from pyworkflow.utils.path import moveFile, cleanPath, cleanPattern
from pwem.protocols import EMProtocol
from tomo.constants import BOTTOM_LEFT_CORNER, TR_SCIPION
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfCoordinates3D, SetOfSubTomograms, SubTomogram, TomoAcquisition, Coordinate3D, Tomogram
# Tomogram type constants for particle extraction
from tomo.utils import scaleTrMatrixShifts

logger = logging.getLogger(__name__)

SAME_AS_PICKING = 0
OTHER = 1


class OutputExtraction(enum.Enum):
    subtomograms = SetOfSubTomograms


class EmanProtTomoExtraction(EMProtocol, ProtTomoBase):
    """ Extraction for Tomo. Uses EMAN2 e2spt_boxer_old.py."""
    _label = 'Subtomograms extraction from tomogram'
    _possibleOutputs = OutputExtraction
    OUTPUT_PREFIX = _possibleOutputs.subtomograms.name

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomosDict = None

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
        form.addParam('doInvert', params.BooleanParam, default=False,
                      label='Invert contrast?',
                      help='Invert the contrast if your tomogram is black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=False,
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

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tomosDict.keys():
            self._insertFunctionStep(self.writeSetOfCoordinates3D, tsId,
                                     needsGPU=False)
            self._insertFunctionStep(self.extractParticles, tsId,
                                     needsGPU=False)
            self._insertFunctionStep(self.convertOutput, tsId,
                                     needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        matchingTsIds, nonMatchingTsIds = self._getMatchingTsIds()
        if len(nonMatchingTsIds) > 0:
            logger.warning(yellowStr(f'Some tsIds do not match: {nonMatchingTsIds}'))
        self.tomosDict = {tomo.getTsId(): tomo.clone() for tomo in self.getInputTomograms()
                          if tomo.getTsId() in matchingTsIds}

    def writeSetOfCoordinates3D(self, tsId: str):
        logger.info(cyanStr("Writing the coordinates of tomogram %s into EMAN format." % tsId))
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
        logger.info(cyanStr("Extracting the particles from tomogram %s." % tsId))
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
        logger.info(cyanStr("Unstacking the particles extracted from tomogram %s." % tsId))
        program = Plugin.getProgram('e2proc3d.py')
        hdfFile = self._getOutHdfCoordsStack(tsId)
        args = ' --unstacking'
        args += ' %s' % hdfFile
        args += ' %s' % self._getExtraPath(pwutils.replaceBaseExt(hdfFile, 'mrc'))
        args += ' --apix %.3f' % self.getOutputSamplingRate()
        self.runJob(program, args)
        cleanPattern(hdfFile)

    def createOutputStep(self):
        logger.info(cyanStr("Registering the results"))

        outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        outputSubTomogramsSet.setSamplingRate(self.getOutputSamplingRate())
        outputSubTomogramsSet.setCoordinates3D(self._getSetOfCoordinates())
        acquisition = TomoAcquisition()

        firstTomo = self.getInputTomograms().getFirstItem()
        acquisition.copyInfo(firstTomo.getAcquisition())
        outputSubTomogramsSet.setAcquisition(acquisition)

        inCoords = self._getSetOfCoordinates()
        samplingRateInput = inCoords.getSamplingRate()
        samplingRateTomo = firstTomo.getSamplingRate()
        factor = samplingRateInput / samplingRateTomo

        for tomo in self.tomosDict.values():
            coordSet = [item.clone() for item in inCoords.iterCoordinates(volume=tomo)]
            self.readSetOfSubTomograms(tomo, outputSubTomogramsSet, coordSet, factor)

        self._defineOutputs(**{OutputExtraction.subtomograms.name: outputSubTomogramsSet})
        self._defineSourceRelation(self._getSetOfCoordinates(), outputSubTomogramsSet)

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

    def _getMatchingTsIds(self) -> Tuple[Set, Set]:
        tomograms = self.getInputTomograms()
        tomoIds = tomograms.getTSIds()
        coords = self._getSetOfCoordinates()
        coordsTomoIds = coords.getTSIds()
        matches = set(tomoIds) & set(coordsTomoIds)
        notMatches = set(tomoIds) ^ set(coordsTomoIds)
        return matches, notMatches

    @staticmethod
    def _getCoordinateFromItem(item):
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
        logger.info(cyanStr('Registegring the subtomograms from tomogram %s' % tsId))
        subtomoFileList = sorted(glob.glob(self._getExtraPath(f'{tsId}*.mrc')))
        for idx, subtomoFile in enumerate(subtomoFileList):
            # logger.info("Registering subtomogram %s - %s" % (counter, subtomoFile))
            subtomogram = SubTomogram()
            transform = Transform()
            subtomogram.setLocation(subtomoFile)
            currentItem = inputSet[idx]
            coord = EmanProtTomoExtraction._getCoordinateFromItem(currentItem)
            subtomogram.setCoordinate3D(coord)
            trMatrix = copy.copy(EmanProtTomoExtraction._getMatrixFromItem(currentItem))
            transform.setMatrix(scaleTrMatrixShifts(trMatrix, scaleFactor))
            subtomogram.setTransform(transform, convention=TR_SCIPION)
            subtomogram.setVolName(tomo.getFileName())
            outputSubTomogramsSet.append(subtomogram)

    def getOutputSamplingRate(self) -> float:
        return self.getInputTomograms().getSamplingRate()

    def _getCoordsFile(self, tsId: str) -> str:
        return self._getExtraPath(f'{tsId}.coords')

    def _getOutHdfCoordsStack(self, tsId) -> str:
        return self._getExtraPath(f'{tsId}.hdf')
