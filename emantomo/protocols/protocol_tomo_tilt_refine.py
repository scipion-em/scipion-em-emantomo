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
import re
import numpy as np

from pyworkflow import BETA
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params

from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow.utils import removeBaseExt

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition

from emantomo.constants import *

from tomo.objects import AverageSubTomogram, SetOfAverageSubTomograms

from emantomo.convert import setCoords3D2Jsons, tltParams2Json, recoverTSFromObj, refinement2Json, convertImage, \
    jsonFilesFromSet, ctf2Json

SAME_AS_PICKING = 0
OTHER = 1


class EmanProtRefineTS(EMProtocol, ProtTomoBase):
    """ Sub-tilt refinement from Eman e2spt_tiltrefine.py"""
    _label = 'sub-tilt refinement'
    _devStatus = BETA
    OUTPUT_PREFIX = 'outputSetOfSubtomogram'
    filter_choices = ['local', 'localwiener', 'global']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', params.PointerParam, label="Input subtomograms", important=True,
                      pointerClass='SetOfSubTomograms', help='Select the inpu SetOfSubTomograms.')
        form.addParam('inputAverage', params.PointerParam, label='Input average', important=True,
                      pointerClass='AverageSubTomogram', help='Select the input AverageSubTomogram')
        # form.addParam('halfMaps', params.PointerParam, label='Half maps', important=True,
        #               pointerClass='SetOfAverageSubTomograms',
        #               help='Select the half maps associated to the input AverageSubTomogram')
        form.addParam('inputCTF', params.PointerParam, label="Input ctf", allowsNull=True,
                      pointerClass='SetOfCTFTomoSeries',
                      help='Optional - Estimated CTF for the tilts series associates to the '
                           'tomograms used to extract the input subtomograms. Will be taken into '
                           'account during the refinement process if provided.')
        form.addParam('nIters', params.IntParam, label='Number of iterations', default=4)
        form.addParam('keep', params.FloatParam, label='Proportion of tilts to keep', default=0.8)
        form.addParam('maxAlt', params.FloatParam, label='Maximum volume altitude', default=45.0,
                      help='Maximum altitude to insert to volume')
        form.addParam('topHat', params.EnumParam, label='Filter applied during refinement',
                      choices=self.filter_choices, default=1,
                      help="Instead of imposing a uniform Wiener filter, use a tophat filter "
                           "(*global* similar to Relion). *local* is a local tophat filter, "
                           "*localwiener* is a local Wiener filter")
        form.addParam('resume', params.BooleanParam, label='Continue from previous run?', default=False)
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('writeJsonInfo')
        self._insertFunctionStep('extract2D')
        self._insertFunctionStep('createRefinementProject')
        self._insertFunctionStep('refinementStep')
        # self._insertFunctionStep('convertOutput')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------

    def writeJsonInfo(self):
        info_path = self._getExtraPath('info')
        pwutils.makePath(info_path)
        coords = self.inputSubtomos.get().getCoordinates3D()
        tomos = coords.getPrecedents()
        tltSeries = recoverTSFromObj(coords, self)
        self.json_files, self.tomo_files = jsonFilesFromSet(tomos, info_path)
        _ = setCoords3D2Jsons(self.json_files, coords)
        _ = tltParams2Json(self.json_files, tltSeries, mode="a")

        if self.inputCTF.get() is not None:
            _ = ctf2Json(self.json_files, self.inputCTF.get(), mode='a')

    def extract2D(self):
        boxSize = self.inputSubtomos.get().getCoordinates3D().getBoxSize()
        for file in self.tomo_files:
            args = os.path.abspath(file)
            args += " --rmbeadthr=-1 --shrink=1.0 --tltkeep=1.0 --padtwod=1.0  " \
                    "--curves=-1 --curves_overlap=0.5 --compressbits=-1 --boxsz_unbin=%d  " \
                    "--threads=%d" \
                    % (boxSize, self.numberOfThreads.get())
            if self.resume.get():
                args += " --resume"
            program = emantomo.Plugin.getProgram('e2spt_extract.py')
            self.runJob(program, args, cwd=self._getExtraPath())

    def createRefinementProject(self):
        particles_2d_file = os.path.abspath(glob.glob(self._getExtraPath(os.path.join('particles3d', '*.hdf')))[0])
        project_path = self._getExtraPath('spt_00')
        pwutils.makePath(project_path)
        refinement2Json(self, self.inputSubtomos.get())
        convertImage(self.inputAverage.get().getFileName(), os.path.join(project_path, 'threed_01.hdf'))
        half1, half2 = self.inputAverage.get().getHalfMaps().split(',')
        self.convertHalfToEman(project_path, half1)
        self.convertHalfToEman(project_path, half2)
        args_fsc = "%s %s --calcfsc %s" % \
                   (os.path.join(project_path, 'threed_01_even.hdf'),
                    os.path.join(project_path, 'fsc_masked_01.txt'),
                    os.path.join(project_path, 'threed_01_odd.hdf'))
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        self.runJob(program, args_fsc)
        program = emantomo.Plugin.getProgram('e2proclst.py')
        self.runJob(program, '%s --create %s' % (particles_2d_file,
                                                 os.path.abspath(os.path.join(project_path, 'input_ptcls.lst'))),
                    cwd=self._getExtraPath())

    def refinementStep(self):
        args = "--path=%s --iter=-1 --niters=%d --keep=%f --maxalt=45.0 " \
               "--mask=auto  --threads=%d --parallel=thread:%d --tophat=%s" \
               % ('spt_00/', self.nIters.get(), self.keep.get(), self.numberOfThreads.get(),
                  self.numberOfThreads.get(), self.filter_choices[self.topHat.get()])
        program = emantomo.Plugin.getProgram('e2spt_subtilt_old.py')
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        lastImage_hdf = self.getLastFromOutputPath("threed_\d+.hdf")
        lastImage_mrc = self._getExtraPath("Average_refined.mrc")
        convertImage(lastImage_hdf, lastImage_mrc)
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(lastImage_mrc)
        setOfAverageSubTomograms = self._createSet(SetOfAverageSubTomograms, 'subtomograms%s.sqlite', "")
        setOfAverageSubTomograms.copyInfo(self.inputSubtomos.get())
        setOfAverageSubTomograms.append(averageSubTomogram)
        self._defineOutputs(refinedSubtomogram=setOfAverageSubTomograms)
        self._defineSourceRelation(self.inputSubtomos, setOfAverageSubTomograms)

    # --------------------------- INFO functions --------------------------------
    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("Sub-tilt refinement using e2spt_tiltrefine.")
        return methodsMsgs

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Sub-tilt refinement running...")

        if self.getOutputsSize() >= 1:
            last_iter = self.getLastFromOutputPath("threed_\d+.hdf")
            last_iter = int(re.findall(r'\d+', pwutils.removeBaseExt(last_iter))[0])
            summary.append("Refinement finished succesfully after *%d* iterations."
                           % last_iter)
            summary.append("Refined volume *%s* has been saved succesfully." % self.getObjectTag('refinedSubtomogram'))
        else:
            summary.append("Refined volume not ready yet.")
        return summary

    def _validate(self):
        errorMsg = []
        if not self.inputAverage.get().getHalfMaps():
            errorMsg.append('No halves were detected in the introduced average.')
        return errorMsg

    def _warnings(self):
        pass

    # --------------------------- UTILS functions ----------------------------------
    def getOutputPath(self, *args):
        return os.path.join(self._getExtraPath(os.path.join('spt_00', '_00', *args)))

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob.glob(self.getOutputPath("*"))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, os.path.basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]

    @staticmethod
    def convertHalfToEman(project_path, inHalfName):
        if 'even' in removeBaseExt(inHalfName):
            convertImage(inHalfName, os.path.join(project_path, 'threed_01_even.hdf'))
        else:
            convertImage(inHalfName, os.path.join(project_path, 'threed_01_odd.hdf'))
