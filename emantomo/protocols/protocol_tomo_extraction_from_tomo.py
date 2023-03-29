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
from os.path import abspath, basename, join
from emantomo import Plugin
from emantomo.constants import PROC_NORMALIZE
from pyworkflow.mapper.sqlite import ID
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import moveFile, cleanPath, cleanPattern
from pwem.protocols import EMProtocol
from tomo.constants import BOTTOM_LEFT_CORNER, TR_EMAN, TR_SCIPION
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfCoordinates3D, SetOfSubTomograms, SubTomogram, TomoAcquisition, Coordinate3D

# Tomogram type constants for particle extraction
from tomo.utils import scaleTrMatrixShifts

SAME_AS_PICKING = 0
OTHER = 1


class OutputExtraction(enum.Enum):
    subtomograms = SetOfSubTomograms


class EmanProtTomoExtraction(EMProtocol, ProtTomoBase):
    """ Extraction for Tomo. Uses EMAN2 e2spt_boxer_old.py."""
    _label = 'Subtomograms extraction from tomogram'
    _devStatus = BETA
    _possibleOutputs = OutputExtraction
    OUTPUT_PREFIX = _possibleOutputs.subtomograms.name
    tomoFiles = []
    lines = []

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Coordinates/Subtomograms", important=True,
                      pointerClass=[SetOfCoordinates3D, SetOfSubTomograms],  help='Choose coordinates or subtomograms derived from 3d coordinates.')

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
        self._insertFunctionStep(self.writeSetOfCoordinates3D)
        self._insertFunctionStep(self.extractParticles)
        self._insertFunctionStep(self.convertOutput)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------

    def _isInputASetOfSubtomograms(self):
        """ returns true if the input is a set of subtomograms"""
        return isinstance(self.inputCoordinates.get(), SetOfSubTomograms)
    def _getSetOfCoordinates(self):
        if self._isInputASetOfSubtomograms():
            return self.inputCoordinates.get().getCoordinates3D()
        else:
            return self.inputCoordinates.get()

    def writeSetOfCoordinates3D(self):

        inputSet = self.inputCoordinates.get()
        coordSet = self._getSetOfCoordinates()

        # Calculate ratio/factor
        samplingRateCoord = coordSet.getSamplingRate()
        samplingRateTomo = self.getInputTomograms().getFirstItem().getSamplingRate()
        scale = samplingRateCoord / samplingRateTomo

        # Store the tomograms to be use
        tomoDict = dict()
        for tomo in self.getInputTomograms():
            tomoDict[tomo.getTsId()] = tomo.clone()

        def onTomogramFinish(coordList):

            self.info("Finishing conversion of tomogram %s." % tomoId)
            if coordList:
                self.lines.append(coordList)
                self.tomoFiles.append(tomo.getFileName())
                emanCoordFile.close()

        # Variables for each tomogram "step"
        tomoId = None
        item_list = []
        emanCoordFile = None

        # Define iterator based on input type
        if self._isInputASetOfSubtomograms():
            iterator = inputSet.iterSubtomos
            orderBy = [SubTomogram.VOL_NAME_FIELD, ID]
        else:
            iterator = inputSet.iterCoordinates
            orderBy = [Coordinate3D.TOMO_ID_ATTR, ID]

        # Iterate in order based on tomogram/Ts id
        for item in iterator(orderBy=orderBy):

            coord3D = self._getCoordinateFromItem(item)
            currentTomoId = coord3D.getTomoId()

            # When changing the tomogram...
            if currentTomoId != tomoId:
                tomoId = currentTomoId

                onTomogramFinish(item_list)

                tomo = tomoDict.get(currentTomoId, None)
                if tomo is None:
                    self.info("Tomogram %s not found in input tomograms set. Coordinates for this tomogram will be skipped." % tomoId)
                    item_list = []
                    continue
                else:
                    # need to open a new file
                    coordFile = self._getExtraPath(pwutils.replaceBaseExt(tomo.getFileName(), 'coords'))
                    emanCoordFile = open(coordFile, "w")
                    item_list = []


            # On each coordinate ... if there is not a tomogram we skip it
            if tomo is None:
                continue

            baseTomoName = basename(tomo.getFileName())
            baseCoordVolume = basename(coord3D.getVolName())

            # Is this check necessary providing the tomoId matches?
            if baseTomoName == baseCoordVolume:
                emanCoordFile.write("%d\t%d\t%d\n" % (coord3D.getX(BOTTOM_LEFT_CORNER) * scale,
                                            coord3D.getY(BOTTOM_LEFT_CORNER) * scale,
                                            coord3D.getZ(BOTTOM_LEFT_CORNER) * scale))
                newItem = item.clone()
                # Do we need the volume?? --> newItem.setVolume(coord3D.getVolume())
                item_list.append(newItem)
            else:
                self.warning("Tomogram name '%s' does not match associated volume in the coordinate '%s'" % (baseTomoName, baseCoordVolume))

        # Last iteration, call onTomogramChange
        onTomogramFinish(item_list)

    def extractParticles(self):
        for tomo in self.tomoFiles:
            args = '%s ' % abspath(tomo)
            args += "--coords %s --boxsize %i" % (pwutils.replaceBaseExt(tomo, 'coords'), self.boxSize.get())
            if self.doInvert:
                args += ' --invert'
            if self.doNormalize:
                args += ' --normproc %s' % self.getEnumText('normproc')
            # args += ' --cshrink %i' % (samplingRateTomo / samplingRateCoord)

            program = Plugin.getProgram('e2spt_boxer_old.py')
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)
            moveFile(self._getExtraPath(join('sptboxer_01', 'basename.hdf')),
                     self._getExtraPath(pwutils.replaceBaseExt(tomo, 'hdf')))
            cleanPath(self._getExtraPath("sptboxer_01"))

    def convertOutput(self):
        program = Plugin.getProgram('e2proc3d.py')
        for hdfFile in glob.glob(self._getExtraPath('*.hdf')):
            args = ' --unstacking'
            args += ' %s' % abspath(hdfFile)
            args += ' %s' % abspath(self._getExtraPath(pwutils.replaceBaseExt(hdfFile, 'mrc')))
            args += ' --apix %.3f' % self.getOutputSamplingRate()
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)
            cleanPattern(hdfFile)

    def createOutputStep(self):

        #Note: using self.lines here prevents the protocol from continuing in case this step code fails!
        outputSet = None
        outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        outputSubTomogramsSet.setSamplingRate(self.getOutputSamplingRate())
        outputSubTomogramsSet.setCoordinates3D(self._getSetOfCoordinates())
        acquisition = TomoAcquisition()

        firstTomo = self.getInputTomograms().getFirstItem()
        acquisition.setAngleMin(firstTomo.getAcquisition().getAngleMin())
        acquisition.setAngleMax(firstTomo.getAcquisition().getAngleMax())
        acquisition.setStep(firstTomo.getAcquisition().getStep())
        outputSubTomogramsSet.setAcquisition(acquisition)

        samplingRateInput = self.inputCoordinates.get().getSamplingRate()
        samplingRateTomo = firstTomo.getSamplingRate()
        factor = samplingRateInput / samplingRateTomo
        counter = 0

        for item in self.getInputTomograms().iterItems():
            for ind, tomoFile in enumerate(self.tomoFiles):
                if basename(tomoFile) == basename(item.getFileName()):
                    coordSet = self.lines[ind]
                    outputSet, counter = self.readSetOfSubTomograms(tomoFile, outputSubTomogramsSet,
                                                                    coordSet, factor, counter)

        self._defineOutputs(**{OutputExtraction.subtomograms.name:outputSet})
        self._defineSourceRelation(self._getSetOfCoordinates(), outputSet)

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
        summary = []
        summary.append("Tomogram source: *%s*"
                       % self.getEnumText("tomoSource"))
        if self.getOutputsSize() >= 1:
            summary.append("Particle box size: *%s*" % self.boxSize.get())
            summary.append("Subtomogram extracted: *%s*" %
                           self._getSetOfCoordinates().getSize())
        else:
            summary.append("Output subtomograms not ready yet.")
        if self.doInvert:
            summary.append('*Contrast was inverted.*')
        return summary

    def _validate(self):
        errors = []
        if self.tomoSource.get() == SAME_AS_PICKING:
            return errors
        tomo_from_coords = self._getSetOfCoordinates().getPrecedents()
        tomoFiles = [pwutils.removeBaseExt(file) for file in self.getInputTomograms().getFiles()]
        coordFiles = [pwutils.removeBaseExt(file) for file in tomo_from_coords.getFiles()]
        numberMatches = len(set(tomoFiles) & set(coordFiles))
        if numberMatches == 0:
            errors.append("Cannot relate Coordinate Tomograms and New Tomograms. In order to stablish a "
                          "relation, the filename of the corresponding Coordinate Tomograms and New Tomogram "
                          "files must be equal. For example, if a coordinate Coordinate Tomogram file is named Tomo_1.mrc, "
                          "then the New Tomogram file to be associated to it should be named Tomo_1.ext "
                          "(being 'ext' any valid extension - '.mrc', '.em'...).\n")
        return errors

    def _warnings(self):
        warnings = []
        if self.tomoSource.get() != SAME_AS_PICKING:
            precedentsSet = self._getSetOfCoordinates().getPrecedents()
            if getattr(precedentsSet.getFirstItem(), '_tsId', None) and \
                    getattr(self.inputTomograms.get().getFirstItem(), '_tsId', None):
                # Match by id (tomoId, tsId)
                tomoIds = [tomo.getTsId() for tomo in self.inputTomograms.get()]
                coordPrecedentsIds = [tomo.getTsId() for tomo in precedentsSet]
                numberMatches = len(set(tomoIds) & set(coordPrecedentsIds))  # Length of the intersection of both lists
                maxNumberFound = max(len(tomoIds), len(coordPrecedentsIds))
                if numberMatches < maxNumberFound:
                    warnings.append("Couldn't find a correspondence between coordinates precedents and the introduced "
                                    "tomograms in which the extraction is desired to be performed. These means that "
                                    "the tsId label is different in both sets of tomograms, at least for some of them.")
                    mismatchesIds = set(coordPrecedentsIds).difference(tomoIds)
                    if mismatchesIds:
                        warnings.append("The following tsIds will not be associated to any New Tomogram (tsIds):")
                        for id in mismatchesIds:
                            warnings.append("\t%s" % id)
                        warnings.append("\n")
            else:
                # Match by filename
                tomoFiles = [pwutils.removeBaseExt(file) for file in self.getInputTomograms().getFiles()]
                coordFiles = [pwutils.removeBaseExt(file) for file in precedentsSet.getFiles()]
                numberMatches = len(set(tomoFiles) & set(coordFiles))
                maxNumberFound = max(len(tomoFiles), len(coordFiles))

                if numberMatches < maxNumberFound:
                    warnings.append("Couldn't find a correspondence between all tomogram files. "
                                    "Association is performed in terms of the file name of the Coordinate Tomograms "
                                    "and the New Tomograms (without the extension). For example, if a Coordinate "
                                    "Tomogram file is named Tomo_1.mrc, then the New Tomogram file file to be "
                                    "associated to it should be named Tomo_1.ext (being 'ext' any valid extension "
                                    "- '.mrc', '.em'...).\n")
                    mismatches_coords = set(coordFiles).difference(tomoFiles)
                    if mismatches_coords:
                        warnings.append("The following Coordinate Tomogram files will not be associated to any New "
                                        "Tomogram (name without extension):")
                        for file in mismatches_coords:
                            warnings.append("\t%s" % file)
                        warnings.append("\n")
                    mismatches_tomos = set(tomoFiles).difference(coordFiles)
                    if mismatches_tomos:
                        warnings.append("The following New Tomogram files will not be associated to any Coordinate Tomogram "
                                        "(name without extension):")
                        for file in mismatches_tomos:
                            warnings.append("\t%s" % file)
                        warnings.append("\n")
        return warnings

    # --------------------------- UTILS functions ----------------------------------

    def _tomosOther(self):
        """ Return True if other tomograms are used for extract. """
        return self.tomoSource == OTHER

    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self._getSetOfCoordinates().getPrecedents()
        else:
            return self.inputTomograms.get()

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

    def readSetOfSubTomograms(self, tomoFile, outputSubTomogramsSet, inputSet, factor, counter):
        """
        Populates the set of subtomograms

        :param tomoFile: tomogram file
        :param outputSubTomogramsSet: output set of subtomograms
        :param inputSet: Subtomograms or 3D coordinates set
        :param factor: factor between the inputSet and the tomogram
        :param counter: counter for eman hdf index
        """

        import time
        time.sleep(10)
        outRegex = self._getExtraPath(pwutils.removeBaseExt(tomoFile) + '-*.mrc')
        subtomoFileList = sorted(glob.glob(outRegex))
        itemList = [item.clone() for item in inputSet] # Get the items (coords or subtomos) in a list)
        for idx, subtomoFile in enumerate(subtomoFileList):
            self.debug("Registering subtomogram %s - %s" % (counter, subtomoFile))
            subtomogram = SubTomogram()
            transform = Transform()
            subtomogram.setLocation(subtomoFile)
            currentItem = itemList[idx]
            coord = EmanProtTomoExtraction._getCoordinateFromItem(currentItem)
            subtomogram.setCoordinate3D(coord)
            trMatrix = copy.copy(EmanProtTomoExtraction._getMatrixFromItem(currentItem))
            transform.setMatrix(scaleTrMatrixShifts(trMatrix, factor))
            subtomogram.setTransform(transform, convention=TR_SCIPION)
            subtomogram.setVolName(tomoFile)
            outputSubTomogramsSet.append(subtomogram)
            counter += 1
        return outputSubTomogramsSet, counter

    def getOutputSamplingRate(self):
        return self.getInputTomograms().getSamplingRate()
