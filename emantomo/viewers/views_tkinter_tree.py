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

from pyworkflow import utils as pwutils
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.utils.path import moveFile, cleanPath, copyFile
from pyworkflow.utils.process import runJob

import emantomo


class EmanDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, **kwargs):
        self.path = path
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.doubleClickOnTomogram,
                                   **kwargs)

    def refresh_gui(self):
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui)
        else:
            outFile = '*%s_info.json' % pwutils.removeBaseExt(self.tomo.getFileName().split("__")[0])
            pattern = os.path.join(self.path, "info", outFile)
            files = glob.glob(pattern)

            moveFile((files[0]), os.path.join(self.path, os.path.basename(files[0])))
            cleanPath(os.path.join(self.path, "info"))
            self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchEmanForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def lanchEmanForTomogram(self, tomo):
        self._moveCoordsToInfo(tomo)

        program = emantomo.Plugin.getProgram("e2spt_boxer.py")
        arguments = "%s" % os.path.abspath(tomo.getFileName())
        runJob(None, program, arguments, env=emantomo.Plugin.getEnviron(), cwd=self.path)

    def _moveCoordsToInfo(self, tomo):
        fnCoor = '*%s_info.json' % pwutils.removeBaseExt(tomo.getFileName().split("__")[0])
        pattern = os.path.join(self.path, fnCoor)
        files = glob.glob(pattern)

        if files:
            infoDir = pwutils.join(os.path.abspath(self.path), 'info')
            pathCoor = os.path.join(infoDir, os.path.basename(files[0]))
            pwutils.makePath(infoDir)
            copyFile(files[0], pathCoor)
