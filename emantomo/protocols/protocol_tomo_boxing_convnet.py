# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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
from os.path import join

from emantomo.constants import TOMOGRAMS_DIR, APIX_UNBIN, TS_FILE, TS_DIR
from emantomo.objects import EmanMetaData
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TOMOS
from emantomo.utils import getPresentTsIdsInSet, genJsonFileName
from pwem.emlib.image import ImageHandler

from pyworkflow import BETA
import pyworkflow.utils as pwutils
from pyworkflow.utils import createLink, getExt
from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import IntParam, BooleanParam, StringParam, USE_GPU, GPU_LIST, LEVEL_ADVANCED

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D
import tomo.constants as const

from emantomo.convert import loadJson, readSetOfCoordinates3D, writeJson
import emantomo
from tomo.utils import getNonInterpolatedTsFromRelations


class EmanProtTomoConvNet(ProtTomoPicking, ProtEmantomoBase):
    """
    Performs deep learning based particle picking on tomographic data using
    EMAN tomography tools. The protocol identifies candidate particles directly
    within reconstructed tomograms and generates a set of 3D coordinates that
    can be used for downstream subtomogram extraction, classification, or
    refinement workflows.

    AI Generated:

    Tomogram ConvNet Picking (EmanProtTomoConvNet) - User Manual
        Overview

        The Tomogram ConvNet Picking protocol provides an interactive deep
        learning approach for detecting particles inside cryo-electron
        tomograms. Its main objective is to accelerate particle identification
        in complex cellular or in situ environments where manual annotation
        would otherwise be extremely time consuming. The protocol integrates
        EMAN neural network based particle detection tools into tomography
        workflows and produces sets of 3D coordinates suitable for subsequent
        subtomogram averaging procedures.

        In biological cryo-electron tomography, particle picking is often one
        of the most challenging steps because tomograms contain low signal,
        missing wedge artifacts, crowded environments, and heterogeneous
        structures. Automated neural network assisted picking can greatly
        improve throughput while still allowing users to supervise and validate
        the detected coordinates interactively.

        Inputs and Biological Context

        The protocol requires a set of tomograms as input. These tomograms are
        expected to represent reconstructed cellular or purified specimen
        volumes where the target particles are visible at least partially.
        Typical use cases include ribosome localization inside cells, membrane
        protein detection, viral particle identification, or extraction of
        macromolecular complexes from crowded cytoplasmic environments.

        The protocol is especially useful when datasets contain many tomograms
        with repetitive structures that can be recognized visually by a neural
        network model. Since tomographic datasets often exhibit anisotropic
        resolution due to the missing wedge, automated detection may still
        require careful biological inspection to avoid false positives or
        missed particles.

        Deep Learning Assisted Picking

        The neural network based strategy attempts to recognize structural
        patterns within tomographic slices and volumes in a manner similar to
        modern image recognition approaches. Rather than relying only on
        thresholding or template matching, the method uses learned structural
        features to identify candidate particles.

        From a biological perspective, this approach is advantageous for
        flexible complexes, crowded intracellular regions, or specimens with
        variable contrast where classical template matching may fail. However,
        successful picking still depends strongly on tomogram quality, contrast
        transfer correction, reconstruction quality, and the visual consistency
        of the target structure.

        Interactive Workflow

        The protocol launches an interactive graphical environment for particle
        annotation and verification. This allows the user to supervise the
        neural network predictions and refine the selection interactively.

        In practical biological workflows, users commonly inspect several
        representative tomograms first to verify that the detected coordinates
        correspond to meaningful macromolecular structures. This step is
        essential because automated detection methods can occasionally identify
        contaminants, fiducials, membrane edges, or reconstruction artifacts as
        particles.

        The interactive strategy provides an effective balance between
        automation and expert supervision. High throughput screening can be
        achieved while preserving biological confidence in the resulting
        coordinate sets.

        GPU Acceleration and Performance

        The protocol supports GPU acceleration for neural network execution.
        This significantly improves processing speed and makes the method
        practical for large tomographic datasets. GPU execution is particularly
        important when working with high resolution tomograms or large numbers
        of volumes.

        Although CPU execution may still be possible in some environments, GPU
        acceleration is generally recommended for routine use. The quality of
        the results is not determined by the GPU itself, but computational
        performance and responsiveness improve substantially during interactive
        analysis.

        Coordinate Generation and Box Size

        The protocol produces a set of 3D particle coordinates associated with
        the original tomograms. These coordinates define the particle centers
        and can later be used for subtomogram extraction.

        The selected box size determines the dimensions of the extracted
        subtomograms in downstream processing. Biologically, the box should be
        large enough to fully contain the macromolecular complex together with
        a reasonable margin of surrounding density. Boxes that are too small
        may truncate important structural regions, whereas excessively large
        boxes increase computational cost and introduce unnecessary background
        noise.

        The protocol also supports assigning group identifiers to coordinates.
        This can be useful when combining multiple particle populations or
        distinguishing biologically distinct structures inside the same
        tomographic dataset.

        Tomogram Requirements and Practical Considerations

        Reliable neural network based picking depends strongly on tomogram
        quality. Well aligned tilt series, appropriate reconstruction methods,
        and adequate contrast restoration generally improve particle detection.
        Poorly reconstructed tomograms or datasets with strong artifacts may
        reduce detection accuracy substantially.

        Square tomograms are recommended for compatibility with the EMAN
        processing environment. Users should also verify that voxel sizes and
        coordinate conventions remain consistent throughout the tomography
        workflow to avoid extraction or alignment problems later.

        In crowded cellular environments, it is often advisable to visually
        inspect representative regions after picking to estimate the false
        positive rate. Some biological targets with strong conformational
        variability or weak contrast may still require partial manual curation.

        Outputs and Their Interpretation

        The primary output is a set of 3D particle coordinates linked to the
        input tomograms. These coordinates preserve the spatial context of the
        detected particles and can be directly used for subtomogram extraction,
        averaging, classification, or visualization.

        The output coordinates should be interpreted as candidate particle
        locations rather than definitive biological assignments. Downstream
        subtomogram classification and averaging remain essential to separate
        true particles from contaminants or incorrect detections.

        Practical Recommendations

        For routine biological workflows, it is advisable to begin with a box
        size slightly larger than the expected particle diameter and visually
        inspect the detected coordinates before large scale extraction. When
        processing heterogeneous cellular tomograms, validating picks on a few
        representative tomograms often prevents propagation of systematic
        detection errors.

        GPU acceleration is recommended whenever available, particularly for
        large datasets or interactive annotation sessions. Users should also
        monitor particle density carefully, since extremely dense coordinate
        distributions may indicate over-picking caused by noise or membrane
        features.

        Final Perspective

        Deep learning based particle picking represents a major advancement for
        cryo-electron tomography workflows because it reduces manual effort and
        enables large scale in situ structural studies. Nevertheless, automated
        detection should always be interpreted within the biological context of
        the sample. Careful validation of coordinates, appropriate box size
        selection, and integration with downstream classification workflows are
        essential for obtaining biologically reliable subtomogram datasets.
    """
    _label = 'tomo boxer convnet'
    _devStatus = BETA
    # nn_boxSize = 96

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)
        form.addHidden(USE_GPU, BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                           Select the one you want to use.")
        form.addParam('boxSize', IntParam, label="Box Size",
                      default=96,
                      help='Final box size for the coordinates')
        form.addParam('groupId', IntParam, label="GroupId", default=1,
                      help="Select a group ID that will be given to the particles. This value is useful to indentify "
                           "different structures in a SetOfCoordinates3D when different sets are joint.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

    # --------------------------- STEPS functions -----------------------------
    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        self._insertFunctionStep(self.convertInputStep, mdObjDict)
        self._insertFunctionStep(self.launchBoxingGUIStep, interactive=True)

    def _initialize(self):
        inTomoSet = self.inputTomograms.get()
        # inTsSet = getNonInterpolatedTsFromRelations(getattr(self, IN_TOMOS), self)
        self.createInitEmanPrjDirs()
        # Manage the tomograms
        presentTsIds = set(getPresentTsIdsInSet(inTomoSet))
        # tsIdsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inTsSet if ts.getTsId() in presentTsIds}
        tomoIdsDict = {tomo.getTsId(): tomo.clone() for tomo in inTomoSet if tomo.getTsId() in presentTsIds}

        # Split all the data into subunits (EmanMetaData objects) referred to the same tsId
        mdObjDict = {}
        for tomoId in presentTsIds:
            mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
                                             # ts=tsIdsDict[tomoId],
                                             tsHdfName=join(TS_DIR, f'{tomoId}.hdf'),
                                             inTomo=tomoIdsDict[tomoId],
                                             tomoHdfName=join(TOMOGRAMS_DIR, f'{tomoId}.hdf'),
                                             jsonFile=genJsonFileName(self.getInfoDir(), tomoId))
        return mdObjDict

    def convertInputStep(self, mdObjDict):
        for mdObj in mdObjDict.values():
            inFile = mdObj.inTomo.getFileName()
            outFile = self._getExtraPath(TOMOGRAMS_DIR, mdObj.tsId + getExt(inFile))
            createLink(inFile, outFile)
            # We have to write each json file with the minimum info required without having to ask for the TS,
            # which is:
            # {
            # "apix_unbin": 3.929999828338623,
            # "tlt_file": "tiltseries/tomoa.hdf"
            # }
            jsonDict = {
                APIX_UNBIN: 3.93,
                TS_FILE: mdObj.tsHdfName
            }
            writeJson(jsonDict, mdObj.jsonFile)

    # def convertInputStep(self):
    #     out_path = self._getExtraPath('tomograms')
    #     info_path = self._getExtraPath('info')
    #     pwutils.makePath(out_path)
    #     pwutils.makePath(info_path)
    #     # program = emantomo.Plugin.getProgram("e2proc3d.py")
    #     for tomo in self.inputTomograms.get().iterItems():
    #         tomo_file = tomo.getFileName()
    #         # tomo_file_hdf = pwutils.removeBaseExt(tomo_file) + ".hdf"
    #         # dim = tomo.getDimensions()
    #         # Only rescale Tomomgrams if needed. Otherwise create a symbolic link to save space
    #         # if self.minBoxSize.get() < self.nn_boxSize:
    #         #     out_file = os.path.join(out_path, pwutils.removeBaseExt(tomo_file) + ".mrc")
    #         #     factor = self.nn_boxSize / self.minBoxSize.get()
    #         #     ImageHandler.scaleSplines(tomo_file + ':mrc', out_file, factor)
    #         # else:
    #         #     args = "%s %s --process normalize --clip 927,927,300" % (tomo_file, os.path.join(out_path, tomo_file_hdf))
    #         #     pwutils.runJob(None, program, args, env=emantomo.Plugin.getEnviron())
    #         # args = "%s %s --process normalize --clip %d,%d,%d" \
    #         #        % (tomo_file, os.path.join(out_path, tomo_file_hdf), max(dim), max(dim), dim[2])
    #         # pwutils.runJob(None, program, args, env=emantomo.Plugin.getEnviron())
    #         out_file = os.path.join(out_path, pwutils.removeBaseExt(tomo_file))
    #         pwutils.createLink(tomo_file, out_file)
    #         self.writeInfoJson(tomo_file, info_path)

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


