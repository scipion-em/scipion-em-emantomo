# coding=utf-8
# **************************************************************************
# *
# * Authors:     David Herreros  (dherreros@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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
import enum
import os
import shutil
from os.path import join

from pyworkflow import BETA

from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam, FloatParam, StringParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.utils import Message, makePath, removeBaseExt, replaceExt
from ..constants import SPT_00_DIR, INPUT_PTCLS_LST, THREED_01, SYMMETRY_HELP_MSG

from ..convert import writeSetOfSubTomograms, refinement2Json
import emantomo

from tomo.protocols import ProtTomoBase
from tomo.objects import AverageSubTomogram


class OutputsAverageSubtomos(enum.Enum):
    averageSubTomos = AverageSubTomogram


class EmanProtSubTomoAverage(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_average.py* EMAN2 program.
    Computes the average a selected subset of a SetOfSubtomograms in the predetermined orientation
    """

    _label = 'average subtomograms'
    _devStatus = BETA
    _possibleOutputs = OutputsAverageSubtomos

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.projectPath = None
        self.volumeFileHdf = None
        self.halfEvenFileHdf = None
        self.halfOddFileHdf = None
        self.hdfSubtomosDir = None
        # self.stepsExecutionMode = STEPS_PARALLEL

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfSubTomogram', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('symmetry', StringParam,
                      label='Symmetry',
                      default='c1',
                      allowsNull=False,
                      help=SYMMETRY_HELP_MSG)
        form.addParam('msWedge', FloatParam,
                      default=3,
                      label="Missing wedge threshold",
                      expertLevel=LEVEL_ADVANCED,
                      help="Threshold for identifying missing data in Fourier"
                           " space in terms of standard deviation of each Fourier"
                           " shell. Default 3.0. If set to 0.0, missing wedge correction"
                           " will be skipped")
        form.addParam('skipPostProc', BooleanParam,
                      default=True,
                      label='Skip post process steps (fsc, mask and filters)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('keepHdfFile', BooleanParam,
                      default=False,
                      label='Keep hdf files?',
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to Yes, the generated files will be saved in both HDF and MRC formats. They are '
                           'generated in HDF and then converted into MRC. The HDF files are deleted by default to '
                           'save storage.')
        form.addParallelSection(threads=4, mpi=0)

    # --------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.computeAverageStep)
        self._insertFunctionStep(self.converOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.projectPath = self._getExtraPath(SPT_00_DIR)
        self.hdfSubtomosDir = self._getExtraPath("subtomograms")
        makePath(self.hdfSubtomosDir)
        makePath(self.projectPath)

    def convertInputStep(self):
        inSubtomos = self.inputSetOfSubTomogram.get()

        volName = inSubtomos.getFirstItem().getVolName()
        newFn = removeBaseExt(volName).split('__ctf')[0] + '.hdf'
        newFn = join(self.hdfSubtomosDir, newFn)
        writeSetOfSubTomograms(inSubtomos, self.hdfSubtomosDir)

        refinement2Json(self, inSubtomos)

        # Generate a virtual stack of particle represented by a .lst file, as expected by EMAN
        args = ' --create %s %s' % (join(self.projectPath, INPUT_PTCLS_LST), newFn)
        program = emantomo.Plugin.getProgram('e2proclst.py')
        self.runJob(program, args)

    def computeAverageStep(self):
        args = " --path=%s --keep 1 --wedgesigma=%f --threads %i " % \
               (self.projectPath, self.msWedge.get(), self.numberOfThreads.get())
        if self.skipPostProc.get():
            args += '--skippostp '
        program = emantomo.Plugin.getProgram('e2spt_average.py')
        self.runJob(program, args)

    def converOutputStep(self):
        # Also fix the sampling rate as it might be set wrong (the value stored in the hdf header may be referred
        # to the original binning, and it will be also in the header of the resulting mrc file
        sRate = self.inputSetOfSubTomogram.get().getSamplingRate()
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        self.volumeFileHdf = self._getExtraPath(SPT_00_DIR, THREED_01 + '.hdf')
        self.halfEvenFileHdf = self._getExtraPath(SPT_00_DIR, THREED_01 + '_even.hdf')
        self.halfOddFileHdf = self._getExtraPath(SPT_00_DIR, THREED_01 + '_odd.hdf')
        filesToConvert = [self.volumeFileHdf, self.halfEvenFileHdf, self.halfOddFileHdf]
        for hdfFile in filesToConvert:
            args = "%s %s --apix %f" % (hdfFile, replaceExt(hdfFile, "mrc"), sRate)
            self.runJob(program, args)
            # Remove the hdf files if requested
            if not self.keepHdfFile.get():
                os.remove(hdfFile)

        # Finally, remove the hdf subtomograms generated in the convert input step, if requested
        if not self.keepHdfFile.get():
            shutil.rmtree(self.hdfSubtomosDir)

    def createOutputStep(self):
        mrcExt = 'mrc'
        inSubtomos = self.inputSetOfSubTomogram.get()
        volume = AverageSubTomogram()
        volume.setFileName(replaceExt(self.volumeFileHdf, mrcExt))
        volume.setHalfMaps([replaceExt(self.halfEvenFileHdf, mrcExt), replaceExt(self.halfOddFileHdf, mrcExt)])
        volume.setSamplingRate(inSubtomos.getSamplingRate())
        volume.fixMRCVolume()

        self._defineOutputs(**{OutputsAverageSubtomos.averageSubTomos.name: volume})
        self._defineSourceRelation(inSubtomos, volume)

    # --------------- INFO functions -------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'averageSubTomos'):
            summary.append("Average not ready yet.")
        else:
            summary.append("Average has been reconstructed.")
        return summary
