# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)
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

import glob
import os
import threading
from os.path import abspath, join, exists
from emantomo.constants import TOMOGRAMS_DIR
from pyworkflow import utils as pwutils
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.utils.process import runJob
import emantomo
from emantomo.convert import loadJson, jsonFileFromTomoFile, writeViewerCoordsCounterFile, getViewerInfoPath, \
    writeTomoCoordsJson


class EmanDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, **kwargs):
        self.proc = None
        self.tomo = None
        self.path = path
        self.provider = kwargs.get("provider", None)
        self.coords = kwargs.get("coords", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.doubleClickOnTomogram,
                                   allowSelect=False,
                                   lockGui=False,
                                   cancelButton=True,
                                   **kwargs)

    def refresh_gui(self):
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui)
        else:
            outFile = '*%s_info*.json' % pwutils.removeBaseExt(self.tomo.getFileName().split("__")[0])
            jsonPath = os.path.join(getViewerInfoPath(self.path), outFile)
            jsonPath = glob.glob(jsonPath)[0]
            jsonDict = loadJson(jsonPath)
            count = len(jsonDict["boxes_3d"])
            self.tomo.count = count
            writeViewerCoordsCounterFile(self.tomo.getTsId(), self.path, count)
            self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.launchEmanForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def launchEmanForTomogram(self, tomo):
        tomoFile = join(TOMOGRAMS_DIR, tomo.getTsId() + '.hdf')  # PPPT --> use the HDF file generated in the convert
        if not exists(join(self.path, tomoFile)):  # Any other case
            tomoFile = abspath(tomo.getFileName())
        if self.coords:
            jsonFile = jsonFileFromTomoFile(tomoFile, join(self.path, 'info'))
            writeTomoCoordsJson(self.coords, tomo.getTsId(), jsonFile)
        program = emantomo.Plugin.getProgram("e2spt_boxer.py")
        arguments = "%s --box3d" % tomoFile
        runJob(None, program, arguments, env=emantomo.Plugin.getEnviron(), cwd=self.path)
