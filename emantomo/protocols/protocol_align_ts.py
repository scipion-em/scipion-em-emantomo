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
import math

import emantomo

from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol
from pwem.convert.transformations import euler_matrix
import pwem.objects as data

from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj

from emantomo.convert import loadJson


class EmanProtAlignTs(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2tomogram.py* EMAN2 program.

    Alignment of the tilt-series is performed iteratively.
    """
    _label = 'align tilt series'
    _devStatus = BETA

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

        form.addParam('niter', params.StringParam,
                      default='2,1,1,1',
                      label='Number of iterations',
                      help='Number of iterations for bin8, bin4, bin2 images')

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

        form.addParallelSection(threads=4, mpi=0)

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
            'niter': self.niter.get(),
            'pk_mindist': self.pkMindist.get(),
            'pkkeep': self.pkkeep.get(),
            'filterto': self.filterto.get(),
            'threads': self.numberOfThreads.get()
        }

        args = (' --notmp --npk=%(npk)d --tltkeep=%(tltkeep)f'
                ' --bxsz=%(bxsz)d --niter=%(niter)s'
                ' --pk_mindist=%(pk_mindist)f --pkkeep=%(pkkeep)f'
                ' --filterto=%(filterto)f --dryrun'
                )

        if command_params["rawtlt"]:
            self._getRawTiltAngles()
        if command_params["tltax"] is not None:
            args += ' --tltax=%(tltax)f'

        args += ' --threads=%(threads)d'

        for path in self._getInputPaths():

            tlt_file = os.path.abspath(self._getExtraPath(pwutils.removeBaseExt(path) + ".txt"))

            if os.path.isfile(tlt_file):
                args_file = '--rawtlt=%s ' % tlt_file + args
            else:
                args_file = '--tltstep=%(tiltStep)f --zeroid=%(zeroid)d ' + args

            args_file = path + " " + args_file

            program = emantomo.Plugin.getProgram("e2tomogram.py")
            self._log.info('Launching: ' + program + ' ' + args_file % command_params)
            self.runJob(program, args_file % command_params, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

    def createOutputStep(self):
        tilt_series = self.tiltSeries.get()
        out_tilt_series = self.getOutputSetOfTiltSeries(tilt_series)
        json_paths = self._getJsonPaths()

        for tilt_serie in tilt_series.iterItems():
            out_tilt_serie = tomoObj.TiltSeries(tsId=tilt_serie.getTsId())
            out_tilt_serie.copyInfo(tilt_serie)
            out_tilt_series.append(out_tilt_serie)
            item_basename = os.path.basename(os.path.dirname(tilt_serie.getFirstItem().getFileName())) +\
                            '-' + tilt_serie.getTsId() + '_info'
            json_path = [path for path in json_paths if item_basename == pwutils.removeBaseExt(path)]
            params = loadJson(json_path[0])['tlt_params']
            for idx, tiltImage in enumerate(tilt_serie.iterItems()):
                out_tiltImage = tomoObj.TiltImage()
                out_tiltImage.copyInfo(tiltImage, copyId=True)
                out_tiltImage.setLocation(tiltImage.getLocation())
                matrix = euler_matrix(math.radians(params[idx][2]),
                                      math.radians(params[idx][3]),
                                      math.radians(params[idx][4]),
                                      'szyx')
                matrix[0, 3], matrix[1, 3] = params[idx][0], params[idx][1]
                transform = data.Transform()
                transform.setMatrix(matrix)
                out_tiltImage.setTransform(transform)
                out_tilt_serie.append(out_tiltImage)

        self._defineOutputs(alignedTiltSeries=out_tilt_series)
        self._defineSourceRelation(self.tiltSeries, out_tilt_series)

    def getOutputSetOfTiltSeries(self, tiltSeries):
        outputSetOfTiltSeries = self._createSetOfTiltSeries()
        outputSetOfTiltSeries.copyInfo(tiltSeries)
        outputSetOfTiltSeries.setDim(tiltSeries.getDim())
        return outputSetOfTiltSeries

    def _getJsonPaths(self):
        pattern = os.path.join(self._getExtraPath("info"), '*.json')
        files = glob.glob(pattern)
        assert files, "Json file with aligment information not found"
        return [os.path.abspath(path) for path in sorted(files)]

    def _getInputPaths(self):
        tilt_series = self.tiltSeries.get()
        return [os.path.abspath(path) for item in tilt_series for path in item.getFiles()]

    def _getRawTiltAngles(self):
        tilt_series = self.tiltSeries.get()
        for tilt in tilt_series.iterItems():
            tilt_angles_dict = tilt.aggregate(["MAX"], "_tiltAngle", ["_tiltAngle"])
            tilt_angles = ["%f" % d['_tiltAngle'] for d in tilt_angles_dict]
            file_tilt_angles = self._getExtraPath(pwutils.removeBaseExt(tilt.getTsId()) + ".txt")
            with open(file_tilt_angles, 'w') as file:
                file.write('\n'.join(tilt_angles))

    def _methods(self):
        return [
            "From an unaligned tilt series: aligned titlt series by e2tomogram.py",
            "Note: Tiltseries must have the correct Apix values in their headers"
        ]

    def _summary(self):
        tilt_series = self.tiltSeries.get()
        return [
            "Input tilt series: {} (size: {})".format(tilt_series.getName(), tilt_series.getSize()),
            "Tilt step: {}".format(self.tiltStep.get())
        ]
