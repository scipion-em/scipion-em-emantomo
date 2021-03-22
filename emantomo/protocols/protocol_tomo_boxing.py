# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
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

from pyworkflow import BETA
from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import BooleanParam, PointerParam, EnumParam

from emantomo.convert import setCoords3D2Jsons, jsons2SetCoords3D
from emantomo.viewers.views_tkinter_tree import EmanDialog

from tomo.protocols import ProtTomoPicking
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


class EmanProtTomoBoxing(ProtTomoPicking):
    """ Manual picker for Tomo. Uses EMAN2 e2spt_boxer.py.
    """
    _label = 'tomo boxer'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('selection', EnumParam, choices=['Yes', 'No'], default=1,
                      label='Modify previous coordinates?', display=EnumParam.DISPLAY_HLIST,
                      help='This option allows to add and/or remove coordinates to a previous SetOfCoordinates')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates", condition='selection == 0',
                      allowsNull=True, pointerClass='SetOfCoordinates3D',
                      help='Select the previous SetOfCoordinates you want to modify')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # Copy input coordinates to Extra Path
        if self.inputCoordinates.get():
            self._insertFunctionStep('copyInputCoords')

        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    def _createOutput(self):
        jsons2SetCoords3D(self, self.inputTomograms.get(), self._getExtraPath())

    # --------------------------- STEPS functions -----------------------------
    def copyInputCoords(self):
        setCoords3D2Jsons(self.inputTomograms.get(), self.inputCoordinates.get(), self._getExtraPath())

    def launchBoxingGUIStep(self):

        if self.getOutputsSize() >= 1:
            lastOutput = [output for _, output in self.iterOutputAttributes()][-1]
            setCoords3D2Jsons(self.inputTomograms.get(), lastOutput, self._getExtraPath())

        tomoList = [tomo.clone() for tomo in self.inputTomograms.get().iterItems()]

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), "json")

        self.dlg = EmanDialog(None, self._getExtraPath(), provider=tomoProvider,)

        # Open dialog to request confirmation to create output
        import tkinter as tk
        frame = tk.Frame()
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
            self._createOutput()

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
