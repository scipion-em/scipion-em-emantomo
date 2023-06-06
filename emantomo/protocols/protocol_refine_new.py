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
import ast
import json
from enum import Enum
import numpy as np
from emantomo import Plugin
from emantomo.convert import convertBetweenHdfAndMrc
from emantomo.objects import EmanSetOfParticles
from pwem.objects import SetOfFSCs, Transform
from pyworkflow.object import Float
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, StringParam, EnumParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from tomo.constants import TR_EMAN
from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from emantomo.constants import SYMMETRY_HELP_MSG, REFERENCE_NAME, PARTICLE_FILE, SCORE, MATRIX, PARTICLE_IND, \
    DEFOCUS, CLASS, PROJ_MATRIX, PART3D_ID, TILT_ID, ROT_TR_MATRIX, TR_MATRIX, EMAN_SCORE

# 3D maps filtering options
WIENER = 'wiener'
GLOBAL = 'global'
LOCAL_WIENER = 'localwiener'
LOCAL = 'local'
filteringKeys = [WIENER, GLOBAL, LOCAL_WIENER, LOCAL]
mapFilterDict = dict(zip(filteringKeys, range(len(filteringKeys))))


class EmanRefineNewOutputs(Enum):
    subtomograms = SetOfSubTomograms
    subtomogramAverage = AverageSubTomogram
    FSCs = SetOfFSCs


class EmanProtTomoRefinementNew(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_refine_new.py* EMAN2 program.
    This refinement protocol performs subtomogram, subtilt and defocus refinement. The 2D subtilt particles are used
    instead of 3D subvolumes in the subtomogram refinement step. Moreover, this program now can model the localized 2D
    particle motion by considering the motion trajectory of each particle along with its neighbor.
    """

    _label = 'subtomogram refinement new'
    _possibleOutputs = EmanRefineNewOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inParticles = None
        self.avgHdf = None
        self.evenHdf = None
        self.oddFnHdf = None

    # --------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_SUBTOMOS, PointerParam,
                      pointerClass='EmanSetOfParticles',
                      label='Particles',
                      important=True)
        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      allowsNull=True,
                      label="Reference volume (opt.)")
        form.addParam('startRes', FloatParam,
                      default=50,
                      label='Refinement initial resolution (Å)',
                      help='This will be the maximum resolution considered for the first iteration. In later '
                           'iterations, the maximum resolution is calculated from the FSC of the previous iteration '
                           '(unless the parameter max. resolution is specified).')

        form.addSection(label='Refinement')
        form.addParam('iters', StringParam,
                      default='p3,t2,p,t,r,d',
                      label='Iteration information',
                      important=True,
                      help='Input types of refinement separated by comma:\n\n'
                           '\t- *p*: 3d particle translation-rotation.\n'
                           '\t- *t*: subtilt translation.\n'
                           '\t- *r*: subtilt translation-rotation.\n'
                           '\t- *d*: subtilt defocus.\n\n'
                           'Default is p,p,p,t,p,p,t,r,d. Character followed by number is also acceptable. p3 = p,p,p.')
        form.addParam('symmetry', StringParam,
                      default='c1',
                      label='Symmetry',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('pkeep', FloatParam,
                      default=0.95,
                      label='Particle keep',
                      help='Fraction of particles to keep. Note that this is controlled at three separate steps.\n'
                           'When default = 0.95, it removes:\n\n'
                           '\t- The worst 5% 3D particles.\n'
                           '\t- The 5% 2D subtilt with the worst score.\n'
                           '\t- The 5% of subtilt with the largest drift.\n\n'
                           'Also accept comma separated values 0.9,0.5,0.5 to set different keep thresholds for the '
                           'three removal operations described before.')
        form.addParam('topHat', EnumParam,
                      choices=list(mapFilterDict.keys()),
                      display=EnumParam.DISPLAY_COMBO,
                      default=mapFilterDict[WIENER],
                      label='3D map filtering',
                      help='Options to to filter the 3D maps:\n\n'
                           '\t- wiener: wiener filter based on FSC curve. default mode in most programs.\n'
                           '\t- global: tophat filter across the map at the resolution cutoff 0.143 from '
                           'fsc_masked_xx.txt.\n'
                           '\t- localwiener: wiener filter based on the fsc curve of local regions from the even/odd '
                           'maps.\n'
                           '\t- local: tophat filter based on local resolution calculated from the even/odd maps at '
                           '0.143 cutoff.')
        form.addParam('maxResAli', FloatParam,
                      default=0,
                      label='Max. resolution in alignment (Å)',
                      help='The program will determine maximum resolution each round from the FSC of the previous '
                           'round by default.')
        form.addParam('minResAli', FloatParam,
                      default=0,
                      label='Min. resolution in alignment (Å)')

        form.addSection(label='Local refine')
        form.addParam('doLocalRefine', BooleanParam,
                      default=False,
                      label='Do local refine?',
                      help='Perform only local search around the solution from the last iteration.')
        form.addParam('maxAng', IntParam,
                      default=30,
                      condition='doLocalRefine',
                      label='Maximum angular diff. (deg.)',
                      help='maximum angle difference from starting point for local refine (in degrees)')
        form.addParam('maxShift', IntParam,
                      default=-1,
                      condition='doLocalRefine',
                      label='Maximum shift (pix.)',
                      help='If set to -1, it will be estimated as maxShift= boxSize/6.')
        group = form.addGroup('Motion correction', condition='doLocalRefine')
        group.addParam('smooth', IntParam,
                       default=100,
                       label='Smooth motion factor',
                       help='Controls how many of its neighbors are considered to model the local motion. '
                            'Smoother local motion with larger numbers.')
        group.addParam('smoothN', IntParam,
                       default=15,
                       label='No. neighboring particles used for smoothing',
                       help='Used to control how many of its neighbors are considered to model the local motion.')

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
        self._initialize()
        self._insertFunctionStep(super().createEmanPrjPostExtractionStep)
        self._insertFunctionStep(super().convertRefVolStep)
        self._insertFunctionStep(super().buildEmanSetsStep)
        self._insertFunctionStep(self.refineStep)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.inParticles = self.getAttrib(IN_SUBTOMOS)
        self.inSamplingRate = self.inParticles.getSamplingRate()
        noIters = self._getNoIters()
        self.avgHdf = self.getRefinedAverageFn(noIters)
        self.evenHdf = self.getRefineEvenFn(noIters)
        self.oddFnHdf = self.getRefineOddFn(noIters)
        self.ali2dLst = self.getAli2dFile(noIters)
        self.ali3dLst = self.getAli3dFile(self._getNoPIters())

    def refineStep(self):
        program = Plugin.getProgram("e2spt_refine_new.py")
        self.runJob(program, self._genRefineCmd(), cwd=self._getExtraPath())

    def convertOutputStep(self):
        # Average and halves
        inFiles = [self.avgHdf, self.evenHdf, self.oddFnHdf]
        args = '--apix %d' % self.inSamplingRate
        for inFile in inFiles:
            outFile = inFile.replace('.hdf', '.mrc')
            convertBetweenHdfAndMrc(self, inFile, outFile, args)

    def createOutputStep(self):
        HDF = '.hdf'
        MRC = '.mrc'
        # Output 1: average with halves
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(self.avgHdf.replace(HDF, MRC))
        averageSubTomogram.setHalfMaps([self.evenHdf.replace(HDF, MRC), self.oddFnHdf.replace(HDF, MRC)])
        averageSubTomogram.setSamplingRate(self.inParticles.getSamplingRate())
        # # Output 2: aligned particles
        outParticles = EmanSetOfParticles.create(self._getPath(), template='emanParticles%s.sqlite')
        outParticles.copyInfo(self.inParticles)
        align3dData = self.parse3dAlignFile()
        for particle, alignDict in zip(self.inParticles, align3dData):
            outParticle = particle.clone()
            setattr(outParticle, EMAN_SCORE, Float(alignDict[SCORE]))
            outParticle.setTransform(Transform(alignDict[ROT_TR_MATRIX]), convention=TR_EMAN)
            outParticles.append(outParticle)


        # outputSetOfSubTomograms.setCoordinates3D(self.inParticles.getCoordinates3D())
        # updateSetOfSubTomograms(self.inParticles, outputSetOfSubTomograms, particleParams)
        # # Output 3: FSCs
        # fscs = self._createSet(SetOfFSCs, 'fsc%s.sqlite', "")
        # fscMasked = self.getLastFromOutputPath('fsc_masked_\d+.txt')
        # fscUnmasked = self.getLastFromOutputPath('fsc_unmasked_\d+.txt')
        # fscTight = self.getLastFromOutputPath('fsc_maskedtight_\d+.txt')
        #
        # emanFSCsToScipion(fscs, fscMasked, fscUnmasked, fscTight)
        outputs = {self._possibleOutputs.subtomogramAverage.name: averageSubTomogram,
                   self._possibleOutputs.subtomograms.name: outParticles}#,
        #            EmanTomoRefinementOutputs.FSCs.name: fscs}
        #
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inParticles, averageSubTomogram)
        self._defineSourceRelation(self.inParticles, outParticles)

    # --------------------------- UTILS functions ------------------------------
    def _genRefineCmd(self):
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
            args += '--loadali2d %s ' % self.inParticles.getAli2d()
            args += '--loadali3d %s ' % self.inParticles.getAli3d()
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
        return True if not self.inParticles.getAli2d() and not self.inParticles.getAli3d() else False

    def _getNoIters(self):
        """From Eman doc: Default is p,p,p,t,p,p,t,r,d. Character followed by number is also acceptable. p3 = p,p,p."""
        iterStr = self.iters.get()
        itersList = iterStr.split(',')
        noIters = 0
        for i in itersList:
            if len(i) > 1:
                noIters += int(i[1:])
            else:
                noIters += 1
        return noIters

    def _getNoPIters(self):
        """From Eman doc: aliptcls3d files are only produced for 'p' iterations"""
        iterStr = self.iters.get()
        itersList = iterStr.split(',')
        return len(itersList) - itersList[::-1].index("p")

    def parse3dAlignFile(self):
        keys = [PARTICLE_IND, PARTICLE_FILE, SCORE, ROT_TR_MATRIX]
        list_of_lists = []
        lastRow = np.array([0, 0, 0, 1])

        with open(self.ali3dLst, 'r') as f:
            for line in f:
                lineContentsList = line.split('\t')
                if len(lineContentsList) > 1:  # There are some explicative lines at the beginning
                    jsonData = json.loads(lineContentsList[2])
                    matrixAsList = ast.literal_eval(jsonData[ROT_TR_MATRIX][MATRIX])
                    matrix = np.array(matrixAsList).reshape(3, 4)
                    matrix = np.vstack([matrix, lastRow])
                    list_of_lists.append([
                        lineContentsList[0],
                        lineContentsList[1],
                        jsonData[SCORE],
                        matrix])

        # Create a list of dictionaries based on the lists filled before
        return [dict(zip(keys, values)) for values in list_of_lists]

    def parse2dAlignFile(self):
        keys = [PARTICLE_IND, PARTICLE_FILE, SCORE, CLASS, DEFOCUS,
                PART3D_ID, TILT_ID, TR_MATRIX, PROJ_MATRIX]
        list_of_lists = []
        lastRow = np.array([0, 0, 0, 1])

        with open(self.ali2dLst, 'r') as f:
            for line in f:
                lineContentsList = line.split('\t')
                if len(lineContentsList) > 1:
                    jsonData = json.loads(lineContentsList[2])
                    particleInd = lineContentsList[0]
                    particleFile = lineContentsList[1]
                    nClass = jsonData[CLASS]
                    defocus = jsonData[DEFOCUS]
                    trMatrix = ast.literal_eval(jsonData['dxf'][MATRIX])
                    trMatrix = np.array(trMatrix).reshape(3, 4)
                    trMatrix = np.vstack([trMatrix, lastRow])
                    ptcl3dId = jsonData[PART3D_ID]
                    score = jsonData[SCORE]
                    tiltId = jsonData[TILT_ID]
                    projMatrix = ast.literal_eval(jsonData[PROJ_MATRIX][MATRIX])
                    projMatrix = np.array(projMatrix).reshape(3, 4)
                    projMatrix = np.vstack([projMatrix, lastRow])
                    list_of_lists.append(
                        [particleInd, particleFile, score, nClass, defocus, ptcl3dId, tiltId, trMatrix, projMatrix])

        # Create a list of dictionarys based on the lists filled before
        return [dict(zip(keys, values)) for values in list_of_lists]

    # --------------------------- INFO functions --------------------------------
