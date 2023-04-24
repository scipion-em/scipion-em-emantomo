# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
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
from os.path import exists
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TS, IN_CTF, IN_TOMOS
from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import PointerParam
from emantomo.convert import jsons2SetCoords3D, ts2Json, ctfTomo2Json, loadJson
from emantomo.viewers.views_tkinter_tree import EmanDialog
from tomo.protocols import ProtTomoPicking
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


class EmanProtTomoBoxing_NEW(ProtEmantomoBase, ProtTomoPicking):
    """ Manual picker for Tomo. Uses EMAN2 e2spt_boxer.py."""
    _label = 'Subtomograms manual picking'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt series with alignment, non-interpolated',
                      # expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      help='Tilt series with alignment (non interpolated) used in the tomograms reconstruction. '
                           'To be deprecated!!')
        form.addParam(IN_CTF, PointerParam,
                      label="CTF tomo series",
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      help='Estimated CTF for the tilt series associates to the tomograms used to pick the input '
                           'coordinates. The corresponding tilt series data will be also accessed through them.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertTsStep, mdObj)
            self._insertFunctionStep(self.convertTomoStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
        self._insertFunctionStep(self.launchBoxingGUIStep, mdObjDict.values(), interactive=True)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inCtfSet = getattr(self, IN_CTF).get()
        inTsSet = self.getTs()
        tomograms = getattr(self, IN_TOMOS).get()
        return self.genMdObjDict(inTsSet, inCtfSet, tomograms=tomograms)

    def writeData2JsonFileStep(self, mdObj):
        mode = 'a' if exists(mdObj.jsonFile) else 'w'
        ts2Json(mdObj, mode=mode)
        ctfTomo2Json(mdObj, self.sphAb, self.voltage, mode='a')

    def launchBoxingGUIStep(self, mdObjList):
        # Manage the coordinate counter per tomogram
        tomoList = []
        for mdObj in mdObjList:
            tomo = mdObj.inTomo
            jsonFile = mdObj.jsonFile
            count = 0
            if exists(jsonFile):
                jsonDict = loadJson(jsonFile)
                coords = jsonDict.get("boxes_3d", [])
                count = len(coords)
            tomo.count = count
            tomoList.append(tomo)

        tomoProvider = TomogramsTreeProvider(tomoList, self.getInfoDir(), "json")
        self.dlg = EmanDialog(None, self._getExtraPath(), provider=tomoProvider)

        # Open dialog to request confirmation to create output
        import tkinter as tk
        frame = tk.Frame()
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
            self._createOutput()

    def _createOutput(self):
        jsons2SetCoords3D(self, self.inputTomograms, self.getInfoDir())

    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.inputTomograms is None:
            return ['Input tomogram not available yet.']

        methodsMsgs.append("Input tomograms imported of dims %s." %(
                              str(self.inputTomograms.get().getDim())))

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                msg = self.getMethods(output)
                methodsMsgs.append("%s: %s" % (self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs
