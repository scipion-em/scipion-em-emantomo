# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *              Arnau Sanchez  (arnau@eyeseetea.com) [1]
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


import os
import glob
import numpy as np

import emantomo

from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms, TomoAcquisition

from emantomo.convert import writeJson, jsonFilesFromSet, tltParams2Json


class EmanProtTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2tomogram.py* EMAN2 program.

    Tomogram reconstruction from aligned tilt series.
    Tomograms are not normally reconstructed at full resolution, generally limited to 1k x 1k or 2k x 2k,
    but the tilt-series are aligned at full resolution. For high resolution subtomogram averaging, the raw
    tilt-series data is used, based on coordinates from particle picking in the downsampled tomograms.
    On a typical workstation reconstruction takes about 4-5 minutes per tomogram.
    """
    _label = 'tomo reconstruction'
    _devStatus = BETA
    choices = ['1k', '2k', '4k']
    resolution = {'1k': 1024, '2k': 2048, '4k': 4096}

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('tiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series", important=True,
                      help='Select the set of tilt series to reconstruct a tomogram')

        form.addSection(label='Optimization')

        form.addParam('outsize', params.EnumParam,
                      display=params.EnumParam.DISPLAY_COMBO,
                      default=0,
                      choices=['1k', '2k', '4k'],
                      label='Size of output tomograms',
                      help='Size of output tomograms:\n\t*1k*: Binning 4\n\t*2k*: Binning 2'
                           '\n\t*4k*: Binning 1')

        form.addParam('clipz', params.IntParam,
                      default=-1,
                      label='Z thickness',
                      help='Z thickness of the final tomogram output. default is -1, (5/16 of tomogram length)')

        form.addParam('bytile', params.BooleanParam,
                      default=True,
                      label='By tiles',
                      help='Make final tomogram by tiles')

        form.addParam('autoclipxy', params.BooleanParam,
                      default=True,
                      condition="bytile==True",
                      label='Non-Square Tomograms?',
                      help='Optimize the x-y shape of the tomogram to fit in the tilt images. '
                           'Only works in bytile reconstruction. '
                           'Useful for non square cameras.')

        form.addParam('filterto', params.FloatParam,
                      default=0.45,
                      label='ABS Filter',
                      help='Filter to abs')

        form.addParam('rmbeadthr', params.FloatParam,
                      default=-1.0,
                      label='Density value threshold for removing beads',
                      help='"Density value threshold (of sigma) for removing beads. High contrast objects beyond this '
                           'value will be removed. Default is -1 for not removing. try 10 for removing fiducials')

        form.addParam('correctrot', params.BooleanParam,
                      default=False,
                      label='Correct rotation',
                      help='Correct for global rotation and position sample flat in tomogram')

        form.addParam('normslice', params.BooleanParam,
                      default=False,
                      label='Normalize slices',
                      help='Normalize each 2D slice')

        form.addParam('extrapad', params.BooleanParam,
                      default=False,
                      label='Extra pad',
                      help='Pad extra for tilted reconstruction. Slower and costs more memory, but reduces boundary '
                           'artifacts when the sample is thick')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createCommandStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def createCommandStep(self):
        command_params = {
            'outsize': self.choices[self.outsize.get()],
            'clipz': self.clipz.get(),
            'filterto': self.filterto.get(),
            'rmbeadthr': self.rmbeadthr.get(),
            'threads': self.numberOfThreads.get()
        }

        args = (' --notmp --outsize=%(outsize)s'
                ' --clipz=%(clipz)d --filterto=%(filterto)f'
                ' --rmbeadthr=%(rmbeadthr)f --load --niter=0,0,0,0'
                )

        if self.bytile.get():
            args += ' --bytile'
        if self.autoclipxy.get() and self.bytile.get():
            args += ' --autoclipxy'
        if self.correctrot.get():
            args += ' --correctrot'
        if self.normslice.get():
            args += ' --normslice'
        if self.extrapad.get():
            args += ' --extrapad'

        args += ' --threads=%(threads)d'
        
        self._processInput()
        for path in self.tlt_files:
            args_file = os.path.abspath(path) + " " + args

            program = emantomo.Plugin.getProgram("e2tomogram.py")
            self._log.info('Launching: ' + program + ' ' + args_file % command_params)
            self.runJob(program, args_file % command_params, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        tilt_series = self.tiltSeries.get()
        size_x = tilt_series.getXDim()

        exponent = np.ceil(np.log2(size_x / self.resolution[self.choices[self.outsize.get()]]).clip(min=0))
        binning = 2 ** exponent
        sampling_rate = binning * tilt_series.getSamplingRate()

        # Output 1: Main tomograms
        tomograms_paths = self._getOutputTomograms()
        tomograms = self._createSetOfTomograms('tomograms%s.sqlite')
        tomograms.copyInfo(tilt_series)
        tomograms.setSamplingRate(sampling_rate)

        for tomogram_path in tomograms_paths:
            self._log.info('Main tomogram: ' + tomogram_path)
            tomogram = Tomogram()
            tomogram.setFileName(tomogram_path)
            for ts in tilt_series:
               if ts.getTsId() in tomogram_path:
                    tomogram.setSamplingRate(ts.getSamplingRate())
                    tomogram.setTsId(ts.getTsId())
                    # Set tomogram acquisition
                    acquisition = TomoAcquisition()
                    acquisition.setAngleMin(ts.getFirstItem().getTiltAngle())
                    acquisition.setAngleMax(ts[ts.getSize()].getTiltAngle())
                    acquisition.setStep(self.getAngleStepFromSeries(ts))
                    tomogram.setAcquisition(acquisition)
            # Set default tomogram origin
            tomogram.setOrigin(newOrigin=False)
            tomograms.append(tomogram)

        self._defineOutputs(tomograms=tomograms)
        self._defineSourceRelation(self.tiltSeries, tomograms)

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def getAngleStepFromSeries(ts):
        """ This method return the average angles step from a series. """
        angleStepAverage = 0
        for i in range(1, ts.getSize()):
            angleStepAverage += abs(ts[i].getTiltAngle()-ts[i+1].getTiltAngle())

        angleStepAverage /= ts.getSize()-1

        return angleStepAverage

    def _processInput(self):
        tilt_series = self.tiltSeries.get()
        info_path = self._getExtraPath('info')
        pwutils.makePath(info_path)
        self.json_files, self.tlt_files = jsonFilesFromSet(tilt_series, info_path)
        _ = tltParams2Json(self.json_files, tilt_series, mode="w")

    def _getRawTiltAngles(self):
        tilt_series = self.tiltSeries.get()
        for tilt in tilt_series.iterItems():
            tilt_angles_dict = tilt.aggregate(["MAX"], "_tiltAngle", ["_tiltAngle"])
            tilt_angles = ["%f" % d['_tiltAngle'] for d in tilt_angles_dict]
            file_tilt_angles = self._getExtraPath(pwutils.removeBaseExt(tilt.getTsId()) + ".txt")
            with open(file_tilt_angles, 'w') as file:
                file.write('\n'.join(tilt_angles))

    def _getOutputTomograms(self):
        pattern = os.path.join(self._getExtraPath("tomograms"), '*.hdf')
        files = glob.glob(pattern)
        assert files, "Output tomogram file not found"
        return [os.path.abspath(path) for path in sorted(files)]

    # --------------------------- INFO functions ----------------------------
    def _methods(self):
        return [
            "From an aligned tilt series: generated a tomogram using e2tomogram.py",
            "Note: Tiltseries must have the correct Apix values in their headers"
        ]

    def _summary(self):
        tilt_series = self.tiltSeries.get()
        return [
            "Input tilt series: {} (size: {})".format(tilt_series.getName(), tilt_series.getSize()),
        ]
