# coding=utf-8
# **************************************************************************
# *
# * Authors:     David Herreros  (dherreros@cnb.csic.es) [2]
# *
# * [2] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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

from pyworkflow import utils as pwutils, BETA
import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol
from pyworkflow.protocol import STEPS_PARALLEL

from ..convert import writeSetOfSubTomograms, refinement2Json
import emantomo

from tomo.protocols import ProtTomoBase
from tomo.objects import AverageSubTomogram

class OutputAverage(enum.Enum):
    averageSubTomos = AverageSubTomogram

class EmanProtSubTomoAverage(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_average.py* EMAN2 program.
    It computed the average from a SetOfSubtomograms based on their alignment.
    """

    _label = 'average subtomo'
    _devStatus = BETA
    _possibleOutputs = OutputAverage

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True, label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('msWedge', params.FloatParam, default=3,
                      label="Missing wedge threshold",
                      help="Threshold for identifying missing data in Fourier"
                           " space in terms of standard deviation of each Fourier"
                           " shell. Default 3.0. If set to 0.0, missing wedge correction"
                           " will be skipped")

    #--------------- INSERT steps functions ----------------


    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.computeAverage)
        self._insertFunctionStep(self.createOutputStep)

    #--------------- STEPS functions -----------------------
    def convertInputStep(self):
        storePath = self._getExtraPath("subtomograms")
        pwutils.makePath(storePath)

        volName = self.inputSetOfSubTomogram.get().getFirstItem().getVolName()
        self.newFn = pwutils.removeBaseExt(volName).split('__ctf')[0] + '.hdf'
        self.newFn = pwutils.join(storePath, self.newFn)
        writeSetOfSubTomograms(self.inputSetOfSubTomogram.get(), storePath)

        self.project_path = self._getExtraPath('spt_00')
        pwutils.makePath(self.project_path)

        refinement2Json(self, self.inputSetOfSubTomogram.get())

        program = emantomo.Plugin.getProgram('e2proclst.py')
        self.runJob(program, ' --create %s %s' % (os.path.abspath(os.path.join(self.project_path, 'input_ptcls.lst')),
                                                  os.path.abspath(self.newFn)),
                    cwd=self._getExtraPath())

    def computeAverage(self):
        args = " --path=%s --keep 1 --skippostp --wedgesigma=%f" % \
               (self.project_path, self.msWedge.get())
        program = emantomo.Plugin.getProgram('e2spt_average.py')
        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args)

        # Fix the sampling rate as it might be set wrong
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        volumeFile = self._getExtraPath(os.path.join("spt_00", "threed_01.hdf"))
        args = "--apix %f %s %s" % (self.inputSetOfSubTomogram.get().getSamplingRate(),
                                    volumeFile, pwutils.replaceExt(volumeFile, "mrc"))
        self.runJob(program, args)

    def createOutputStep(self):
        imgSet = self.inputSetOfSubTomogram.get()
        volume = AverageSubTomogram()
        volumeFile = self._getExtraPath(os.path.join("spt_00", "threed_01.mrc"))
        volume.setFileName(volumeFile)
        volume.setSamplingRate(imgSet.getSamplingRate())
        volume.fixMRCVolume()

        self._defineOutputs(**{OutputAverage.averageSubTomos.name:volume})
        self._defineSourceRelation(self.inputSetOfSubTomogram, volume)

    #--------------- INFO functions -------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'averageSubTomos'):
            summary.append("Average not ready yet.")
        else:
            summary.append("Average has been reconstructed.")
        return summary