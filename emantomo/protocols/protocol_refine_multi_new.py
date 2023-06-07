# coding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
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
from enum import Enum
from os.path import join

from emantomo import Plugin
from pwem.objects import SetOfFSCs
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, StringParam, EnumParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from emantomo.constants import SYMMETRY_HELP_MSG, SETS_DIR, REFERENCE_NAME

# # 3D maps filtering options
# WIENER = 'wiener'
# GLOBAL = 'global'
# LOCAL_WIENER = 'localwiener'
# LOCAL = 'local'
# filteringKeys = [WIENER, GLOBAL, LOCAL_WIENER, LOCAL]
# mapFilterDict = dict(zip(filteringKeys, range(len(filteringKeys))))


class EmanRefineNewOutputs(Enum):
    subtomograms = SetOfSubTomograms
    subtomogramAverage = AverageSubTomogram
    FSCs = SetOfFSCs


class EmanProtMultiRefinementNew(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_refinemulti_new.py* EMAN2 program.
    Multi-reference classification for the new (2021) SPT refinement protocol.
    """

    _label = 'Multi-reference classification'
    _possibleOutputs = EmanRefineNewOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_SUBTOMOS, PointerParam,
                      pointerClass='EmanSetOfParticles',
                      label='Particles',
                      important=True)
        form.addParam('maskRef', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask the references prior to to classif. (opt.)',
                      allowsNull=True)

        form.addSection(label='Classification')
        form.addParam('nClasses', IntParam,
                      default=1,
                      label="Number of classes",
                      important=True)
        form.addParam('nIter', IntParam,
                      default=5,
                      label='No. iterations')
        form.addParam('symmetry', StringParam,
                      default='c1',
                      label='Sym. to apply to the average',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('breakSym', StringParam,
                      allowsEmpty=True,
                      label='Break specified symmetry',
                      help='If empty, no symmetry will be broken.')

        form.addSection(label='Alignment')
        form.addParam('doAlignment', BooleanParam,
                      default=True,
                      label='Do alignment?')
        form.addParam('maskAlign', PointerParam,
                      pointerClass='VolumeMask',
                      label='Apply mask to the 3D alignment ref. in each iter. (opt.)',
                      allowsNull=True,
                      help="Not applied to the average, which will follow normal EMAN's masking routine.")
        form.addParam('maxRes', FloatParam,
                      default=20,
                      condition='doAlignment',
                      important=True,
                      label='Max. resolution (Å)',
                      help='Maximum resolution (the smaller number) to consider in alignment (in Å).\n'
                           'Since gold-standard validation is not used here, setting this parameter is mandatory.')
        form.addParam('minRes', FloatParam,
                      default=0,
                      condition='doAlignment',
                      label='Min. resolution (Å)',
                      help='Minimum resolution (the larger number) to consider in alignment (in Å).')
        form.addParam('maxAng', IntParam,
                      default=-1,
                      condition='doAlignment',
                      label='Maximum angular diff. (deg.)',
                      help='maximum angle difference for local alignment (in degrees)')
        form.addParam('maxShift', IntParam,
                      default=-1,
                      condition='doAlignment',
                      label='Maximum shift (pix.)',
                      help='If set to -1, it will be estimated as maxShift=boxSize/6.')

        form.addSection(label='Extra params')
        form.addParam('threadsPostProc', IntParam,
                      default=10,
                      label='Threads for post-processing')
        form.addParam('make3dThread', BooleanParam,
                      default=False,
                      label='Do make3d in threading mode with shared memory?',
                      expertLevel=LEVEL_ADVANCED,
                      help='Safer for large boxes.')
        form.addParam('extraParams', StringParam,
                      label="Extra params",
                      help="Here you can add any extra parameters to run Eman's  new subtomogram refinement. "
                           "Parameters should be written in Eman's command line format (--param val)")
        form.addParallelSection(threads=4, mpi=0)

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        self._insertFunctionStep(super().createEmanPrjPostExtractionStep)
        self._insertFunctionStep(super().convertRefVolStep)
        self._insertFunctionStep(super().buildEmanSetsStep)
        self._insertFunctionStep(self.refineStep)

    # --------------- STEPS functions -----------------------
    def refineStep(self):
        program = Plugin.getProgram("e2spt_refine_new.py")
        self.runJob(program, self._genRefineCmd(), cwd=self._getExtraPath())

    # --------------------------- UTILS functions ------------------------------
    def _genRefineCmd(self):
        inParticles = getattr(self, IN_SUBTOMOS).get()
        # Input params
        args = '--ptcls %s ' % self._getLstFile()
        if self.getRefVol():
            args += '--ref %s ' % (REFERENCE_NAME + '.hdf')
        args += '--startres %.2f ' % self.startRes.get()
        # Refinement params
        args += '--iters %s ' % self.iters.get()
        args += '--sym %s ' % self.symmetry.get()
        args += '--keep %.2f ' % self.pkeep.get()
        args += '--tophat %s ' % mapFilterDict[filteringKeys[self.topHat.get()]]
        args += '--maxres %.2f ' % self.maxResAli.get()
        args += '--minres %.2f ' % self.minResAli.get()
        if self._doGoldStandard():
            args += '--goldstandard '
        else:
            args += '--goldcontinue '
            args += '--loadali2d %s ' % inParticles.getAli2d()
            args += '--loadali3d %s ' % inParticles.getAli3d()
        # Local refinement params
        if self.doLocalRefine.get():
            args += '--localrefine '
            args += '--maxang %i ' % self.maxAng.get()
            args += '--maxshift %i ' % self.maxShift.get()
            args += '--smooth %.2f ' % self.sooth.get()
            args += '--smoothN %i ' % self.smoothN.get()
        # Extra params
        args += '--parallel=thread:%i ' % self.numberOfThreads.get()
        args += '--threads %i ' % self.threadsPostProc.get()
        if self.make3dThread.get():
            args += '--m3dthread '
        if self.extraParams.get():
            args += '%s ' % self.extraParams.get()
        args += '--verbose 9 '
        return args

    def _doGoldStandard(self):
        """
        From EMAN doc:
            - goldstandard: "Phase randomize the reference to the starting resolution (--startres) independently for
              the even/odd subsets of particles".
            - goldcontinue: "Continue from previous gold standard refinement. Ues the _even/_odd version of the given
              reference".
        Thus, we'll apply it if it is the first refinement, which means that the set of particles introduced do not
        contain an ali2d nor ali3d files in the corresponding attributes.
        """
        inParticles = getattr(self, IN_SUBTOMOS).get()
        return True if not inParticles.getAli2d() and not inParticles.getAli3d() else False

    # --------------------------- INFO functions --------------------------------
