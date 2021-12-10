# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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


import os
import numpy as np

import emantomo

from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils
from pyworkflow.object import Set

from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfCTFTomoSeries, CTFTomoSeries, CTFTomo

from emantomo.convert import writeJson, loadJson, tltParams2Json, jsonFilesFromSet


class EmanProtEstimateCTF(EMProtocol, ProtTomoBase):
    """
    Protocol for CTF estimation from tilt series using e2spt_tomoctf.py
    """
    _label = 'ctf estimation'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('tiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series", important=True,
                      help='Select the set of tilt series from which the CTF will be '
                           'estimated')
        lineDefocus = form.addLine('Defocus Range (um)',
                             help="Search range of defocus (start, end, step). Note that "
                                  "they must be introduced in microns.")

        lineDefocus.addParam('minDefocus', params.FloatParam, default=2.0, label='min')
        lineDefocus.addParam('maxDefocus', params.FloatParam, default=7.0, label='max')
        lineDefocus.addParam('stepDefocus', params.FloatParam, default=0.02, label='Step')

        linePhaseShift = form.addLine('Phase shift Range (degrees)',
                                   help="Search range of the phase shift (start, end, step). To avoid"
                                        "the phase shift search use min 0.0, max 1.0, and step 1.0.")
        linePhaseShift.addParam('minPhaseShift', params.FloatParam, default=10.0, label='min')
        linePhaseShift.addParam('maxPhaseShift', params.FloatParam, default=15.0, label='max')
        linePhaseShift.addParam('stepPhaseShift', params.FloatParam, default=5.0, label='Step')

        form.addParam('tilesize', params.IntParam, label='Tile size',
                      default=256,
                      help='Size of tile to calculate FFT. Default is 256.')
        form.addParam('nref', params.IntParam, label='Number of references',
                      default=15,
                      help='Using N tilt images near the center tilt to estimate '
                           'the range of defocus for all images. Default is 15.')
        form.addParam('stepx', params.IntParam, label='Step in X direction',
                      default=20,
                      help='Number of tiles to generate on x-axis (different defocus)')
        form.addParam('stepy', params.IntParam, label='Step in Y direction',
                      default=40,
                      help='Number of tiles to generate on y-axis (same defocus)')
        # form.addParallelSection(threads=4, mpi=0)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('writeJsonInfo')
        self._insertFunctionStep('createCommandStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def writeJsonInfo(self):
        self._processInput()
        self.tlt_files = [os.path.abspath(tlt_file) for tlt_file in self.tlt_files]

    def createCommandStep(self):
        aquisition = self.tiltSeries.get().getAcquisition()
        cs = aquisition.getSphericalAberration()
        voltage = aquisition.getVoltage()

        args = " ".join(self.tlt_files)
        args += " --dfrange=%f,%f,%f --psrange=%f,%f,%f --tilesize=%d --voltage=%d" \
                " --cs=%f --nref=%d --stepx=%d --stepy=%d --threads=%d" \
                % (self.minDefocus.get(), self.maxDefocus.get(), self.stepDefocus.get(),
                  self.minPhaseShift.get(), self.maxPhaseShift.get(), self.stepPhaseShift.get(),
                  self.tilesize.get(), voltage, cs, self.nref.get(),
                  self.stepx.get(), self.stepy.get(), self.numberOfThreads.get())

        program = emantomo.Plugin.getProgram("e2spt_tomoctf.py")
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        set_tilt_series = self.tiltSeries.get()
        info_path = self._getExtraPath('info')
        outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                             template='CTFmodels%s.sqlite')
        outputSetOfCTFTomoSeries.setSetOfTiltSeries(set_tilt_series)
        outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)

        for tilt_serie in set_tilt_series.iterItems(iterate=False):
            newCTFTomoSeries = CTFTomoSeries()
            newCTFTomoSeries.copyInfo(tilt_serie)
            newCTFTomoSeries.setTiltSeries(tilt_serie)
            newCTFTomoSeries.setObjId(tilt_serie.getObjId())
            newCTFTomoSeries.setTsId(tilt_serie.getTsId())
            outputSetOfCTFTomoSeries.append(newCTFTomoSeries)
            json_file = os.path.join(info_path,
                                     os.path.basename(os.path.dirname(tilt_serie.getFirstItem().getFileName())) +
                                     '-' + tilt_serie.getTsId() + '_info.json')
            json_dict = loadJson(json_file)
            defocus = json_dict['defocus']
            phase_shift = json_dict['phase']
            for idx, tiltImage in enumerate(tilt_serie.iterItems()):
                defocusU = defocusV = 10000.0 * defocus[idx]
                newCTFTomo = CTFTomo()
                newCTFTomo.setIndex(idx + 1)
                newCTFTomo.setPhaseShift(phase_shift[idx])
                newCTFTomo.setDefocusU(defocusU)
                newCTFTomo.setDefocusV(defocusV)
                newCTFTomo.setDefocusAngle(0)
                newCTFTomoSeries.append(newCTFTomo)
            outputSetOfCTFTomoSeries.update(newCTFTomoSeries)

        outputSetOfCTFTomoSeries.write()
        self._store()
        outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_CLOSED)
        self._defineOutputs(estimatedCTF=outputSetOfCTFTomoSeries)
        self._defineSourceRelation(self.tiltSeries, outputSetOfCTFTomoSeries)

    def _processInput(self):
        tilt_series = self.tiltSeries.get()
        info_path = self._getExtraPath('info')
        pwutils.makePath(info_path)
        self.json_files, self.tlt_files = jsonFilesFromSet(tilt_series, info_path)
        _ = tltParams2Json(self.json_files, tilt_series, mode="w")

    def _methods(self):
        return [
            "CTF estimation from tilt series through e2spt_tomoctf.py"
        ]

    def _summary(self):
        pass
