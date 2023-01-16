# coding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team  (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
from enum import Enum
from os.path import join, abspath

import emantomo
from emantomo.constants import SYMMETRY_HELP_MSG, INPUT_PTCLS_LST, SUBTOMOGRAMS_DIR, SPT_00, REFS_DIR
from emantomo.convert import writeSetOfSubTomograms, refinement2Json
from pwem import ALIGN_3D
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam, StringParam, LEVEL_ADVANCED, FloatParam, BooleanParam
from pyworkflow.utils import Message, removeBaseExt, makePath
from tomo.objects import SetOfSubTomograms, SetOfClassesSubTomograms
from tomo.protocols import ProtTomoBase


class mraOutputObjects(Enum):
    subtomograms = SetOfSubTomograms
    classes = SetOfClassesSubTomograms


class EmanProtPcaKMeansClassifySubtomos(EMProtocol, ProtTomoBase):
    """This protocol wraps *e2spt_classify.py* EMAN2 program, which performs a
    multi-reference refinement of subtomograms."""

    _label = 'Multi-reference refinement of subtomograms'
    _devStatus = BETA
    _possibleOutputs = mraOutputObjects

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.subtomosDir = None
        self.spt00Dir = None
        self.refsDir = None

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Subtomograms')
        form.addParam('inRefs', PointerParam,
                      label='References',
                      important=True,
                      pointerClass='Subtomogram, SetOfSubtomograms, Volume, SetOfVolumes')
        form.addParam('mask', PointerParam,
                      label='Mask (opt.)',
                      allowsNull=True,
                      pointerClass='VolumeMask')

        form.addSection(label='Optimization')
        form.addParam('nIter', IntParam,
                      default=3,
                      label='Number of iterations')
        form.addParam('sym', StringParam,
                      default='c1',
                      label='Symmetry',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('targetRes', FloatParam,
                      default=20,
                      expertLevel=LEVEL_ADVANCED,
                      label='Target resolution [Å]',
                      help='The specified target resolution to avoid under-filtering artifacts.')
        form.addParam('mass', FloatParam,
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Mass of the particle [kDa]',
                      help='The rough mass of the particle in kilodaltons, used to normalize by mass. '
                           'Due to resolution effects, not always the true mass.')
        form.addParam('doLocalFilter', BooleanParam,
                      default=False,
                      label='Do local filter (top hat)?')

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.mraStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.subtomosDir = self._getExtraPath(SUBTOMOGRAMS_DIR)
        self.spt00Dir = self._getExtraPath(SPT_00)
        self.refsDir = self._getExtraPath(REFS_DIR)
        makePath(*[self.subtomosDir, self.spt00Dir, self.refsDir])

    def convertInputStep(self):
        volName = self.inSubtomos.get().getFirstItem().getVolName()
        stackHdf = join(self.subtomosDir, removeBaseExt(volName).split('__ctf')[0] + '.hdf')

        # Convert the particles to HDF if necessary
        inSubtomos = self.inSubtomos.get()
        writeSetOfSubTomograms(inSubtomos, self.subtomosDir, lignType=ALIGN_3D)

        # Generate the refinement JSON expected to be in the SPT directory
        refinement2Json(self, inSubtomos)

        # Generate the LST file expected to be in the SPT directory
        program = emantomo.Plugin.getProgram('e2proclst.py')
        args = ' --create %s %s' % (join(self.spt00Dir, INPUT_PTCLS_LST), abspath(stackHdf))
        self.runJob(program, args)

        # Convert the references into HDF if necessary
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        sRate = inSubtomos.getSamplingRate()
        for ref in self.inRefs.get():
            args = "--apix %f %s %s" % (sRate, abspath(ref), self._getExtraPath(removeBaseExt(ref) + '.hdf'))
            self.runJob(program, args)

    def mraStep(self):
        pass

    def createOutputStep(self):
        pass

    # --------------- INFO functions -------------------------
    # --------------- UTILS functions ------------------------
    def genMraCmd(self):
        args =
