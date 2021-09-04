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


import glob
from os.path import abspath

from pyworkflow import BETA
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import moveFile, cleanPath, cleanPattern

from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition

from emantomo.constants import *

import tomo.constants as const

# Tomogram type constants for particle extraction
from emantomo.convert import setCoords3D2Jsons, tltParams2Json, loadJson, recoverTSFromObj

SAME_AS_PICKING = 0
OTHER = 1


class EmanProtTSExtraction(EMProtocol, ProtTomoBase):
    """ Extraction of subtomograms from tilt serie using EMAN2 e2spt_extract.py."""
    _label = 'extraction from TS'
    _devStatus = BETA
    OUTPUT_PREFIX = 'outputSetOfSubtomogram'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Input Coordinates", important=True,
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')

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
                           'The wizard selects same box size as picking')

        form.addParam('downFactor', params.FloatParam, default=1.0,
                      label='Downsampling factor',
                      help='If 1.0 is used, no downsample is applied. '
                           'Non-integer downsample factors are possible. ')

        form.addParam('tltKeep', params.FloatParam, default=1.0,
                      label='Keep tilt fraction',
                      help='Keep a fraction of tilt images with good score '
                           'determined from tomogram reconstruction')

        form.addParam('rmThr', params.FloatParam, default=-1,
                      label='Contrast threshold',
                      help='Remove 2d particles with high contrast object beyond N '
                           'sigma at 100A. Note that this may result in generating '
                           'fewer particles than selected. Default is -1 (include '
                           'all particles). 0.5 might be a good choice for '
                           'removing gold beads but may need some testing...')

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
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('writeSetOfCoordinates3D')
        self._insertFunctionStep('extractParticles')
        self._insertFunctionStep('convertOutput')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------

    def writeSetOfCoordinates3D(self):
        info_path = self._getExtraPath('info')
        pwutils.makePath(info_path)
        coords = self.inputCoordinates.get()
        tomos = coords.getPrecedents()
        tltSeries = recoverTSFromObj(coords, self)
        self.files = setCoords3D2Jsons(tomos, coords, info_path)
        _ = tltParams2Json(tomos, tltSeries, info_path, mode="a")

    def extractParticles(self):
        for pair in self.files:
            args = os.path.abspath(pair[0])
            args += " --rmbeadthr=%f --shrink=%f --tltkeep=%f --padtwod=1.0  " \
                    "--curves=-1 --curves_overlap=0.5 --compressbits=-1 --boxsz_unbin=%d  " \
                    "--threads=%d" \
                    % (self.rmThr.get(), self.downFactor.get(),
                       self.tltKeep.get(), self.boxSize.get(), self.numberOfThreads.get())
            program = emantomo.Plugin.getProgram('e2spt_extract.py')
            self.runJob(program, args, cwd=self._getExtraPath())

    def convertOutput(self):
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        part_path = self._getExtraPath(os.path.join('particles3d', '*.hdf'))
        for hdfFile in glob.glob(part_path):
            args = ' --unstacking'
            args += ' %s' % abspath(hdfFile)
            args += ' %s' % abspath(self._getExtraPath(pwutils.removeBaseExt(hdfFile) + '.mrc'))
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        self.outputSubTomogramsSet.setSamplingRate(self.getInputTomograms().getSamplingRate() / self.downFactor.get())
        self.outputSubTomogramsSet.setCoordinates3D(self.inputCoordinates)
        acquisition = TomoAcquisition()
        acquisition.setAngleMin(self.getInputTomograms().getFirstItem().getAcquisition().getAngleMin())
        acquisition.setAngleMax(self.getInputTomograms().getFirstItem().getAcquisition().getAngleMax())
        acquisition.setStep(self.getInputTomograms().getFirstItem().getAcquisition().getStep())
        self.outputSubTomogramsSet.setAcquisition(acquisition)
        coordSet = self.inputCoordinates.get()
        for tomo in coordSet.getPrecedents().iterItems():
            tomoFile = tomo.getFileName()
            coordSet = [coord.clone() for coord in coordSet.iterCoordinates(volume=tomo)]
            outputSet = self.readSetOfSubTomograms(tomoFile,
                                                   self.outputSubTomogramsSet,
                                                   coordSet)

        self._defineOutputs(outputSetOfSubtomogram=outputSet)
        self._defineSourceRelation(self.inputCoordinates, outputSet)

    # --------------------------- INFO functions --------------------------------
    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            msg = ("A total of %s subtomograms of size %s were extracted"
                   % (str(self.inputCoordinates.get().getSize()), self.boxSize.get()))
            if self._tomosOther():
                msg += (" from another set of tomograms: %s"
                        % self.getObjectTag('inputTomogram'))
            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")
        if self.doInvert:
            methodsMsgs.append("Inverted contrast on images.")
        if self.downFactor.get() != 1:
            methodsMsgs.append("Subtomograms downsample by factor %d."
                               % self.downFactor.get())
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
                           self.inputCoordinates.get().getSize())
        else:
            summary.append("Output subtomograms not ready yet.")
        if self.doInvert:
            summary.append('*Contrast was inverted.*')
        return summary

    def _validate(self):
        errors = []
        if self.tomoSource.get() == SAME_AS_PICKING:
            return errors
        tomo_from_coords = self.inputCoordinates.get().getPrecedents()
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
        if self.tomoSource.get() == SAME_AS_PICKING:
            return warnings
        tomo_from_coords = self.inputCoordinates.get().getPrecedents()
        tomoFiles = [pwutils.removeBaseExt(file) for file in self.getInputTomograms().getFiles()]
        coordFiles = [pwutils.removeBaseExt(file) for file in tomo_from_coords.getFiles()]
        numberMatches = len(set(tomoFiles) & set(coordFiles))
        if numberMatches < max(len(tomoFiles), len(coordFiles)):
            warnings.append("Couldn't find a correspondence between all tomogram files. "
                            "Association is performed in terms of the file name of the Coordinate Tomograms and the New Tomograms "
                            "(without the extension). For example, if a Coordinate Tomogram file is named Tomo_1.mrc, then the New Tomogram file "
                            "file to be associated to it should be named Tomo_1.ext (being 'ext' any valid extension "
                            "- '.mrc', '.em'...).\n")
            mismatches_coords = set(coordFiles).difference(tomoFiles)
            if mismatches_coords:
                warnings.append("The following Coordinate Tomogram files will not be associated to any New Tomogram "
                                "(name without extension):")
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
            return self.inputCoordinates.get().getPrecedents()
        else:
            return self.inputTomograms.get()

    def readSetOfSubTomograms(self, tomoFile, outputSubTomogramsSet, coordSet):
        if "__" in tomoFile:
            tomoFile = tomoFile.split("__")[0]
        else:
            parentFolder = pwutils.removeBaseExt(os.path.dirname(tomoFile))
            tomoFile = '%s-%s' % (parentFolder, tomoFile)
        outRegex = self._getExtraPath(pwutils.removeBaseExt(tomoFile) + '*.mrc')
        subtomoFileList = sorted(glob.glob(outRegex))
        for counter, subtomoFile in enumerate(subtomoFileList):
            subtomogram = SubTomogram()
            subtomogram.cleanObjId()
            subtomogram.setLocation(subtomoFile)
            subtomogram.setCoordinate3D(coordSet[counter])
            subtomogram.setTransform(coordSet[counter]._eulerMatrix)
            subtomogram.setVolName(tomoFile)
            outputSubTomogramsSet.append(subtomogram)
        return outputSubTomogramsSet
