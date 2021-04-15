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

import emantomo

from pyworkflow import BETA
from pyworkflow.protocol import params

from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms


class EmanProtTomoReconstruction(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2tomogram.py* EMAN2 program.

    Alignment of the tilt-series is performed iteratively in conjunction with tomogram reconstruction.
    Tomograms are not normally reconstructed at full resolution, generally limited to 1k x 1k or 2k x 2k,
    but the tilt-series are aligned at full resolution. For high resolution subtomogram averaging, the raw
    tilt-series data is used, based on coordinates from particle picking in the downsampled tomograms.
    On a typical workstation reconstruction takes about 4-5 minutes per tomogram.
    """
    _label = 'tomo reconstruction'
    _devStatus = BETA
    choices = ['1k', '2k', '4k']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('tiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series", important=True,
                      help='Select the set of tilt series to reconstruct a tomogram')

        form.addParam('rawtlt', params.BooleanParam,
                      default=True,
                      label="Use Tilt Angles stored in metadata?")

        form.addParam('tiltStep', params.FloatParam,
                      default=2.0,
                      condition="rawtlt==False",
                      label='Tilt step',
                      help='Step between tilts.')

        form.addParam('zeroid', params.IntParam,
                      default=-1,
                      condition="rawtlt==False",
                      label='Zero ID',
                      help='Index of the center tilt.')

        form.addParam('npk', params.IntParam,
                      default=20,
                      label='Number of landmarks (npk)',
                      help='Number of landmarks to use (such as gold fiducials)')

        form.addParam('bxsz', params.IntParam,
                      default=32,
                      label='Box size for tracking (bxsz)',
                      help='Box size of the particles for tracking. May be helpful to use a larger one for '
                           'fiducial-less cases')

        form.addSection(label='Optimization')

        form.addParam('tltkeep', params.FloatParam,
                      default=0.9,
                      label='Fraction to keep',
                      help='Fraction (0.0 -> 1.0) of tilts to keep in the reconstruction')

        form.addParam('tltax', params.FloatParam,
                      allowsNull=True,
                      default=None,
                      label='Title axis angle',
                      help='Angle of the tilt axis. The program will calculate one if this option is not provided')

        form.addParam('outsize', params.EnumParam,
                      display=params.EnumParam.DISPLAY_COMBO,
                      default=0,
                      choices=['1k', '2k', '4k'],
                      label='Size of output tomograms',
                      help='Size of output tomograms:\n\t*1k*: Binning 4\n\t*2k*: Binning 2'
                           '\n\t*4k*: Binning 1')

        form.addParam('niter', params.StringParam,
                      default='2,1,1,1',
                      label='Number of iterations',
                      help='Number of iterations for bin8, bin4, bin2 images')

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

        form.addParam('load', params.BooleanParam,
                      default=False,
                      label='Load existing parameters',
                      help='Load existing tilt parameters')

        form.addParam('pkMindist', params.FloatParam,
                      default=0.125,
                      label='Min landmarks distance',
                      help='Minimum distance between landmarks, as fraction of micrograph length')

        form.addParam('pkkeep', params.FloatParam,
                      default=0.9,
                      label='Landmarks to keep',
                      help='Fraction of landmarks to keep in the tracking')

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

        form.addParallelSection(threads=12, mpi=0)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createCommandStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def createCommandStep(self):
        command_params = {
            'rawtlt': self.rawtlt.get(),
            'tiltStep': self.tiltStep.get(),
            'zeroid': self.zeroid.get(),
            'npk': self.npk.get(),
            'bxsz': self.bxsz.get(),
            'tltkeep': self.tltkeep.get(),
            'tltax': self.tltax.get(),
            'outsize': self.choices[self.outsize.get()],
            'niter': self.niter.get(),
            'clipz': self.clipz.get(),
            'pk_mindist': self.pkMindist.get(),
            'pkkeep': self.pkkeep.get(),
            'filterto': self.filterto.get(),
            'rmbeadthr': self.rmbeadthr.get(),
            'threads': self.numberOfThreads.get()
        }

        args = " ".join(self._getInputPaths())

        args += (' --notmp --npk=%(npk)d --tltkeep=%(tltkeep)f --outsize=%(outsize)s'
                 ' --bxsz=%(bxsz)d --niter=%(niter)s --clipz=%(clipz)d'
                 ' --pk_mindist=%(pk_mindist)f --pkkeep=%(pkkeep)f'
                 ' --filterto=%(filterto)f --rmbeadthr=%(rmbeadthr)f'
                 )

        if command_params["rawtlt"]:
            rawtlt_file = self._getRawTiltAngles()
            args += ' --rawtlt=%s' % rawtlt_file
        else:
            args += ' --tltstep=%(tiltStep)f --zeroid=%(zeroid)d'
        if command_params["tltax"] is not None:
            args += ' --tltax=%(tltax)f'
        if self.bytile.get():
            args += ' --bytile'
        if self.autoclipxy.get() and self.bytile.get():
            args += ' --autoclipxy'
        if self.load.get():
            args += ' --load'
        if self.correctrot.get():
            args += ' --correctrot'
        if self.normslice.get():
            args += ' --normslice'
        if self.extrapad.get():
            args += ' --extrapad'

        args += ' --threads=%(threads)d'

        program = emantomo.Plugin.getProgram("e2tomogram.py")
        self._log.info('Launching: ' + program + ' ' + args % command_params)
        self.runJob(program, args % command_params, cwd=self._getExtraPath(),
                    numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        tilt_series = self.tiltSeries.get()

        sampling_rate = tilt_series.getSamplingRate()
        if self.outsize.get() == 0:
            sampling_rate *= 4
        elif self.outsize.get() == 1:
            sampling_rate *= 2

        # Output 1: Main tomograms
        tomograms_paths = self._getOutputTomograms()
        tomograms = self._createSet(SetOfTomograms, 'tomograms%s.sqlite', "")
        tomograms.copyInfo(tilt_series)
        tomograms.setSamplingRate(sampling_rate)

        for tomogram_path in tomograms_paths:
            self._log.info('Main tomogram: ' + tomogram_path)
            tomogram = Tomogram()
            tomogram.setFileName(tomogram_path)
            tomogram.copyInfo(tilt_series)
            tomogram.setSamplingRate(sampling_rate)
            tomograms.append(tomogram)

        self._defineOutputs(tomograms=tomograms)
        self._defineSourceRelation(self.tiltSeries, tomograms)

    def _getInputPaths(self):
        tilt_series = self.tiltSeries.get()
        return [path for item in tilt_series for path in item.getFiles()]

    def _getRawTiltAngles(self):
        tilt_series = self.tiltSeries.get()
        tilt_angles = [tilt.aggregate(["MAX"], "_tiltAngle", ["_tiltAngle"]) for tilt in tilt_series.iterItems()]
        tilt_angles = ["%f" % d['_tiltAngle'] for dictionary in tilt_angles for d in dictionary]
        file_tilt_angles = self._getExtraPath("tilt_angles.txt")
        with open(file_tilt_angles, 'w') as file:
            file.write('\n'.join(tilt_angles))
        return os.path.abspath(file_tilt_angles)

    def _getOutputTomograms(self):
        pattern = os.path.join(self._getExtraPath("tomograms"), '*.hdf')
        files = glob.glob(pattern)
        assert files, "Output tomogram file not found"
        return [os.path.abspath(path) for path in sorted(files)]

    def _methods(self):
        return [
            "From an unaligned tilt series: aligned, and generated a tomogram using e2tomogram.py",
            "Note: Tiltseries must have the correct Apix values in their headers"
        ]

    def _summary(self):
        tilt_series = self.tiltSeries.get()
        return [
            "Input tilt series: {} (size: {})".format(tilt_series.getName(), tilt_series.getSize()),
            "Tilt step: {}".format(self.tiltStep.get())
        ]
