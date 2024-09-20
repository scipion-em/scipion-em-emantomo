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
import glob
import os
from os.path import getmtime, join, basename

import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.utils import removeBaseExt
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils

import pwem.viewers.views as vi
from tomo.objects import SetOfCoordinates3D, Tomogram, Coordinate3D
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider
from tomo.protocols.protocol_import_coordinates import ProtImportCoordinates3D

from ..convert import jsons2SetCoords3D, jsonFilesFromSet
from .views_tkinter_tree import EmanDialog


class EmanDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _name = 'Open with Eman'
    _targets = [
        ProtImportCoordinates3D,
        SetOfCoordinates3D
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

        if issubclass(cls, SetOfCoordinates3D):
            outputCoords = obj
            tomoIdField = Coordinate3D.TOMO_ID_ATTR
            countField = 'COUNT'
            tomos = outputCoords.getPrecedents()
            tsIdDictList = outputCoords.aggregate(countField, tomoIdField, [tomoIdField])

            tomoList = []
            tomoCountDict = {}
            for tsIdDict in tsIdDictList:
                tsId = tsIdDict[tomoIdField]
                nCoords = tsIdDict[countField]
                tomogram = tomos.getItem(Tomogram.TS_ID_FIELD, tsId).clone()
                tomogram.count = nCoords
                tomoCountDict[tsId] = nCoords
                tomoList.append(tomogram)

            path = self.protocol._getExtraPath()
            info_path = join(path, 'info')

            tomoProvider = TomogramsTreeProvider(tomoList, info_path, 'json')

            if not os.path.exists(info_path):
                pwutils.makePath(info_path)
            json_files, _ = jsonFilesFromSet(tomos, info_path)

            EmanDialog(self._tkRoot, path, provider=tomoProvider)

            if hasAnyFileChanged(tomoCountDict, info_path):

                import tkinter as tk
                frame = tk.Frame()
                if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                    jsons2SetCoords3D(self.protocol, outputCoords.getPrecedents(), info_path)

        elif issubclass(cls, ProtImportCoordinates3D):
            if obj.getOutputsSize() >= 1:
                for _, out in obj.iterOutputAttributes(SetOfCoordinates3D):
                    lastOutput = out
            self._visualize(lastOutput)

        return views


def hasAnyFileChanged(counterDict, infoPath):
    countFileList = glob.glob(os.path.join(infoPath, '*.count'))
    if countFileList:
        lastModFile = max(countFileList, key=os.path.getmtime)
        with open(lastModFile, 'r') as file:
            count = file.read()
            tsId = removeBaseExt(lastModFile)
            if int(count) != counterDict[tsId]:
                return True
            else:
                return False
    else:
        return False


def hasFileChangedSince(file, time):
    """ Returns true if the file has changed after 'time'"""
    modTime = datetime.datetime.fromtimestamp(getmtime(file))
    return time < modTime
