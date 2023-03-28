# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              David Herreros Calero (dherreros@cnb.csic.es)
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import subprocess

import eman2
import pwem
import pyworkflow.utils as pwutils
from scipion.install.funcs import VOID_TGZ

import emantomo.constants as emanConst

__version__ = "4.0.0"
_logo = "eman2_logo.png"
_references = ['GALAZMONTOYA2015279', 'BELL201625']
_url = "https://blake.bcm.edu/emanwiki/EMAN2/"


class Plugin(eman2.Plugin):

    @classmethod
    def defineBinaries(cls, env):
        pass
