# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
import glob
from enum import Enum
from os.path import basename, join
from emantomo import Plugin
from emantomo.constants import SYMMETRY_HELP_MSG, REFERENCE_NAME, TOMOGRAMS_DIR
from emantomo.convert import jsons2SetCoords3D
from emantomo.protocols.protocol_base import ProtEmantomoBase, REF_VOL, IN_TOMOS
from pyworkflow.protocol import PointerParam, StringParam, FloatParam, LEVEL_ADVANCED, IntParam, BooleanParam, GE, LE
from pyworkflow.utils import Message
from tomo.objects import SetOfCoordinates3D


class OutputsTemplateMatch(Enum):
    coordinates = SetOfCoordinates3D


class EmanProtTemplateMatching(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_tempmatch.py* EMAN2 program.
    It is a reference-based picking (template matching).
    """

    _label = 'Template matching picking'
    _possibleOutputs = OutputsTemplateMatch
    OUTPUT_PREFIX = _possibleOutputs.coordinates.name

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inTomos = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      important=True,
                      help='Specify tomograms containing reference-like particles to be exctracted. These '
                           'should be "dark contrast".')
        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Reference volume",
                      help="This should be 'light contrast'.")

        form.addSection(label='options')
        form.addParam('nptcl', IntParam,
                      default=500,
                      label='Maximum no. particles picked among the tomograms',
                      help='If a higher number of particles is detected, the program will take the best N, '
                           'being N the value of the current parameter.')
        form.addParam('delta', FloatParam,
                      default=30,
                      label='Anglular sampling to rotate the reference (deg.',
                      help='The lower value, the higher number of orientations that will be checked.')
        form.addParam('dthr', FloatParam,
                      default=-1,
                      label='Distance threshold (Å)',
                      help='Particles closer than the value introduced will be removed. If default =-1, '
                           'it will be considered as half of the box size of the reference.')
        form.addParam('vthr', FloatParam,
                      default=5,
                      label='Template matching threshold (n sigma)',
                      validators=[GE(1), LE(10)],
                      help='Particles with score lower than the introduced value will be removed. Admitted values '
                           'are [1, 10], where 1 means very insensitive picking (pick a lot of particles with a high '
                           'probability of having a a high number of false positives) and 10 means very sensitive '
                           'picking (strict picking, with less particles picked and very low presence of false '
                           'positives, but with high probability of some or even a lot of particles not to be picked.')
        form.addParam('minvol', IntParam,
                      default=-1,
                      label='Minimum peak volume',
                      expertLevel=LEVEL_ADVANCED,
                      help='If default=-1, this filter is not applied.')
        form.addParam('maxvol', IntParam,
                      default=-1,
                      label='Maximum peak volume',
                      expertLevel=LEVEL_ADVANCED,
                      help='If default=-1, this filter is not applied.')
        form.addParam('symmetry', StringParam,
                      default='c1',
                      label='Symmetry of the reference',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('rmedge', BooleanParam,
                      default=True,
                      label='Remove particles on the edge?')
        form.addParam('rmgold', BooleanParam,
                      default=True,
                      label='Remove particles near gold fiducials?')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.prepareEmanPrj)
        self._insertFunctionStep(self.convertRefVolStep)
        self._insertFunctionStep(self.templateMatchingStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inTomos = getattr(self, IN_TOMOS).get()

    def prepareEmanPrj(self):
        super().createInitEmanPrjDirs()
        sRate = self.inTomos.getSamplingRate()
        for tomo in self.inTomos:
            inTomoName = tomo.getFileName()
            # Required to ensure that the sampling rate is correct in the header, as EMAN reads and compares it
            # with the sampling rate of the reference volume to scale the data
            self.convertOrLink(inTomoName, tomo.getTsId(), TOMOGRAMS_DIR, sRate)

    def templateMatchingStep(self):
        program = Plugin.getProgram("e2spt_tempmatch.py")
        self.runJob(program, self._genTempMatchArgs(), cwd=self._getExtraPath())

    def createOutputStep(self):
        jsons2SetCoords3D(self, self.inTomos, self.getInfoDir())

    # --------------------------- UTILS functions -----------------------------
    def _genTempMatchArgs(self):
        args = [f" {self._genTomolist()}",
                f"--ref {REFERENCE_NAME}.hdf ",
                f"--nptcl {self.nptcl.get()}",
                f"--dthr {self.dthr.get():.2f}",
                f"--vthr {self.vthr.get():.2f}",
                f"--minvol {self.minvol.get()}",
                f"--maxvol {self.maxvol.get()}",
                f"--delta {self.delta.get():.2f}",
                f"--sym {self.symmetry.get()}",
                f"--threads {self.numberOfThreads.get()}"]
        if self.rmedge.get():
            args.append("--rmedge")
        if self.rmgold.get():
            args.append('--rmgold')
        return ' '.join(args)

    def _genTomolist(self):
        tomoFileList = [join(TOMOGRAMS_DIR, basename(tomoFile)) for tomoFile
                        in glob.glob(join(self.getTomogramsDir(), '*'))]
        return ' '.join(tomoFileList)

    # -------------------------- INFO functions -------------------------------
