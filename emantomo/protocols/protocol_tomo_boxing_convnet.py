# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
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
import glob


from pwem.emlib.image import ImageHandler

from pyworkflow import BETA
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import IntParam, BooleanParam, StringParam, USE_GPU, GPU_LIST, LEVEL_ADVANCED

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D
import tomo.constants as const

from emantomo.convert import loadJson, readSetOfCoordinates3D
import emantomo


class EmanProtTomoConvNet(ProtTomoPicking):
    """Eman Deep Learning based picking for Tomography
    """
    _label = 'tomo boxer convnet'
    _devStatus = BETA
    # nn_boxSize = 96

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)
        form.addHidden(USE_GPU, BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                           Select the one you want to use.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")
        form.addParam('boxSize', IntParam, label="Box Size",
                      default=96,
                      help='Final box size for the coordinates')
        form.addParam('groupId', IntParam, label="GroupId", default=1,
                      help="Select a group ID that will be given to the particles. This value is useful to indentify "
                           "different structures in a SetOfCoordinates3D when different sets are joint.")

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    def convertInputStep(self):
        out_path = self._getExtraPath('tomograms')
        info_path = self._getExtraPath('info')
        pwutils.makePath(out_path)
        pwutils.makePath(info_path)
        # program = emantomo.Plugin.getProgram("e2proc3d.py")
        for tomo in self.inputTomograms.get().iterItems():
            tomo_file = tomo.getFileName()
            # tomo_file_hdf = pwutils.removeBaseExt(tomo_file) + ".hdf"
            # dim = tomo.getDimensions()
            # Only rescale Tomomgrams if needed. Otherwise create a symbolic link to save space
            # if self.minBoxSize.get() < self.nn_boxSize:
            #     out_file = os.path.join(out_path, pwutils.removeBaseExt(tomo_file) + ".mrc")
            #     factor = self.nn_boxSize / self.minBoxSize.get()
            #     ImageHandler.scaleSplines(tomo_file + ':mrc', out_file, factor)
            # else:
            #     args = "%s %s --process normalize --clip 927,927,300" % (tomo_file, os.path.join(out_path, tomo_file_hdf))
            #     pwutils.runJob(None, program, args, env=emantomo.Plugin.getEnviron())
            # args = "%s %s --process normalize --clip %d,%d,%d" \
            #        % (tomo_file, os.path.join(out_path, tomo_file_hdf), max(dim), max(dim), dim[2])
            # pwutils.runJob(None, program, args, env=emantomo.Plugin.getEnviron())
            out_file = os.path.join(out_path, pwutils.removeBaseExt(tomo_file))
            pwutils.createLink(tomo_file, out_file)
            self.writeInfoJson(tomo_file, info_path)

    def launchBoxingGUIStep(self):
        program = emantomo.Plugin.getProgram("e2spt_boxer_convnet.py")
        args = "--label particles_00"
        if self.useGpu.get():
            args += " --gpuid %s" % self.getGpuList()[0]
        pwutils.runJob(None, program, args, env=emantomo.Plugin.getEnviron(), cwd=self._getExtraPath())
        self._createOutput()

    def _createOutput(self):
        setTomograms = self.inputTomograms.get()
        outPath = self._getExtraPath("info")
        coord3DSetDict = {}
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coord3DSet = self._createSetOfCoordinates3D(setTomograms, suffix)
        coord3DSet.setName("tomoCoord")
        coord3DSet.setPrecedents(setTomograms)
        coord3DSet.setSamplingRate(setTomograms.getSamplingRate())
        coord3DSet.setBoxSize(self.boxSize.get())
        for tomo in setTomograms.iterItems():
            outFile = '*%s_info.json' % pwutils.removeBaseExt(tomo.getFileName().split("__")[0])
            pattern = os.path.join(outPath, outFile)
            files = glob.glob(pattern)

            if not files or not os.path.isfile(files[0]):
                continue

            jsonFnbase = files[0]
            jsonBoxDict = loadJson(jsonFnbase)

            index = int((list(jsonBoxDict["class_list"].keys()))[0])
            coord3DSetDict[index] = coord3DSet

            # Populate Set of 3D Coordinates with 3D Coordinates
            # factor = self.minBoxSize.get() / self.nn_boxSize if self.minBoxSize.get() is not None else 1
            # FIXME: Correct the scaling factor when there is a mismatch between the sr in the header and in Scipion
            # FIXME: Could be useful in the future?
            # sr = setTomograms.getSamplingRate()
            # if mrcfile.validate(tomo.getFileName()):
            #     with mrcfile.open(tomo.getFileName()) as mrc:
            #         sr_header = mrc.voxel_size.tolist()[0]
            # factor *= sr / sr_header
            readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, tomo.clone(),
                                   origin=const.CENTER_GRAVITY, groupId=self.groupId.get())

        name = self.OUTPUT_PREFIX + suffix
        args = {}
        args[name] = coord3DSet
        self._defineOutputs(**args)
        self._defineSourceRelation(setTomograms, coord3DSet)

        # Update Outputs
        for index, coord3DSet in coord3DSetDict.items():
            self._updateOutputSet(name, coord3DSet, state=coord3DSet.STREAM_CLOSED)

    # --------------------------- UTILS functions -----------------------------
    def writeInfoJson(self, tomo_file, info_path):
        # boxSize = self.minBoxSize.get() if self.minBoxSize.get() else self.nn_boxSize
        contents = '{ "boxes_3d": [], "apix_unbin": %.2f, ' \
                   '"class_list": { "0": { "boxsize": 96, "name": ' \
                   '"particles_00"} } }' % (self.inputTomograms.get().getSamplingRate())
        info_file = os.path.join(info_path, pwutils.removeBaseExt(tomo_file) + "_info.json")
        with open(info_file, 'w') as fid:
            fid.write(contents)

    # --------------------------- INFO functions -----------------------------
    def getInfo(self, output):
        msg = '\tNumber of particles picked: *%d* \n' % output.getSize()
        msg += '\tParticle box size: *%d*' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                msg = self.getInfo(output)
                methodsMsgs.append("%s: \n %s" % (self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def _summary(self):
        summary = []
        if self.getOutputsSize() < 1:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        else:
            for key, output in self.iterOutputAttributes():
                msg = self.getInfo(output)
                summary.append("%s: \n %s" % (self.getObjectTag(output), msg))
        return summary

    def validate(self):
        errors = []
        dim = self.inputTomograms.get().getFirstItem().getDimensions()
        if dim[0] != dim[1] and not emantomo.Plugin.isVersion(emantomo.constants.V_CB):
            errors.append("Error: input tomograms must be square. Please, use a resizing protocol or reconstruct "
                          "your tomograms so X and Y dimensions match.")
        return errors

    def _warnings(self):
        warnings = []
        # if self.minBoxSize.get() < self.nn_boxSize:
        #     warnings.append("Boxsize is smaller than the minimum size allowed by Eman (96). This implies "
        #                     "that temporary rescaled Tomograms will be created so your boxsize corresponds "
        #                     " to a size of 96 to work with Eman. This may occupy a large space in disk.")
        return warnings


