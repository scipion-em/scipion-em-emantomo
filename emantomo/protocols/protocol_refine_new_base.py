# coding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import shutil
from enum import Enum
from os.path import exists

from emantomo import Plugin
from emantomo.convert import convertBetweenHdfAndMrc
from emantomo.convert.lstConvert import EmanLstReader
from emantomo.objects import EmanSetOfParticles
from pwem.convert.headers import fixVolume
from pwem.emlib.image import ImageHandler
from pwem.objects import SetOfFSCs
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, StringParam, EnumParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from emantomo.constants import SYMMETRY_HELP_MSG, REFERENCE_NAME, SPT_00_DIR


class EmanProtRefineNewBase(ProtEmantomoBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------- DEFINE param functions ----------------

    # --------------- INSERT steps functions ----------------

    # --------------- STEPS functions -----------------------

    # --------------- UTILS functions -----------------------

    # --------------- INFO functions ------------------------
    def _validate(self):
        errorMsg = []
        refVol = self.getAttrib(REF_VOL)
        if refVol:
            # Check the dimensions
            ih = ImageHandler()
            x, y, z, _ = ih.getDimensions(refVol.getFileName())
            refVolDims = (x, y, z)
            inParticles = self.getAttrib(IN_SUBTOMOS)
            inParticlesDims = inParticles.getBoxSize()
            if refVolDims != inParticlesDims:
                errorMsg.append(f'The dimensions of the reference volume {refVolDims} px and the particles '
                                f'{inParticlesDims} px must be the same')
            # Check the sampling rate
            tol = 1e-03
            inParticlesSRate = inParticles.getSamplingRate()
            refVolSRate = refVol.getSamplingRate()
            if abs(inParticlesSRate - refVolSRate) >= tol:
                errorMsg.append(
                    f'The sampling rate of the input particles [{inParticlesSRate} Å/pix] and the reference volume '
                    f'[{refVolSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')
        return errorMsg
