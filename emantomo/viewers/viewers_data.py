# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
import datetime
import os
from os.path import getmtime

import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils

import pwem.viewers.views as vi
import tomo.objects
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider
from tomo.protocols.protocol_import_coordinates import ProtImportCoordinates3D

from ..convert import setCoords3D2Jsons, jsons2SetCoords3D, jsonFilesFromSet
from .views_tkinter_tree import EmanDialog



class EmanDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        ProtImportCoordinates3D,
        tomo.objects.SetOfCoordinates3D
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.SetOfCoordinates3D):
            outputCoords = obj
            tomos = outputCoords.getPrecedents()

            volIds = outputCoords.aggregate(["MAX", "COUNT"], "_volId", ["_volId"])
            volIds = [(d['_volId'], d["COUNT"]) for d in volIds]

            tomoList = []
            for objId in volIds:
                tomogram = tomos[objId[0]].clone()
                tomogram.count = objId[1]
                tomoList.append(tomogram)

            path = self.protocol._getExtraPath()
            info_path = self.protocol._getExtraPath('info')

            tomoProvider = TomogramsTreeProvider(tomoList, info_path, 'json',)

            if not os.path.exists(info_path):
                pwutils.makePath(info_path)
            json_files, _ = jsonFilesFromSet(tomos, info_path)
            _ = setCoords3D2Jsons(json_files, outputCoords)

            time = datetime.datetime.now()

            EmanDialog(self._tkRoot, path, provider=tomoProvider)

            if hasAnyFileChanged(json_files, time):

                import tkinter as tk
                frame = tk.Frame()
                if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                    jsons2SetCoords3D(self.protocol, outputCoords.getPrecedents(), info_path)

        elif issubclass(cls, ProtImportCoordinates3D):
            if obj.getOutputsSize() >= 1:
                for _, out in obj.iterOutputAttributes(tomo.objects.SetOfCoordinates3D):
                    lastOutput = out
            self._visualize(lastOutput)

        return views

def hasAnyFileChanged(files, time):
    """ Returns true if any of the files in files list has been changed after 'time'"""
    for file in files:
        if hasFileChangedSince(file, time):
            return True

    return False

def hasFileChangedSince(file, time):
    """ Returns true if the file has changed after 'time'"""
    modTime = datetime.datetime.fromtimestamp(getmtime(file))
    return time < modTime

