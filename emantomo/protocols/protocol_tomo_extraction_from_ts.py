# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es) [1]
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
from emantomo.convert import setCoords3D2Jsons, tltParams2Json, loadJson, recoverTSFromObj, jsonFilesFromSet, ctf2Json

SAME_AS_PICKING = 0
OTHER = 1


class EmanProtTSExtraction(EMProtocol, ProtTomoBase):
    """ Extraction of subtomograms from tilt serie using EMAN2 e2spt_extract.py."""
    _label = 'extraction from tilt-series'
    _devStatus = BETA
    OUTPUT_PREFIX = 'outputSetOfSubtomogram'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Input coordinates", important=True,
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('inputCTF', params.PointerParam, label="Input ctf", allowsNull=True,
                      pointerClass='SetOfCTFTomoSeries',
                      help='Optional - Estimated CTF for the tilts series associates to the '
                           'tomograms used to pick the input coordinates. Will be taken into '
                           'account during the extraction if provided.')
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
        self.json_files, self.tomo_files = jsonFilesFromSet(tomos, info_path)
        _ = setCoords3D2Jsons(self.json_files, coords)
        _ = tltParams2Json(self.json_files, tltSeries, mode="a")

        if self.inputCTF.get() is not None:
            _ = ctf2Json(self.json_files, self.inputCTF.get(), mode='a')

    def extractParticles(self):
        for file in self.tomo_files:
            args = os.path.abspath(file)
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
            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")
        if self.downFactor.get() != 1:
            methodsMsgs.append("Subtomograms downsample by factor %d."
                               % self.downFactor.get())
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
        pass

    def _warnings(self):
        pass

    # --------------------------- UTILS functions ----------------------------------
    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        return self.inputCoordinates.get().getPrecedents()

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
            subtomogram.setTransform(coordSet[counter]._eulerMatrix, convention=const.TR_EMAN)
            subtomogram.setVolName(tomoFile)
            outputSubTomogramsSet.append(subtomogram)
        return outputSubTomogramsSet
