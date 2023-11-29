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
import shutil
from enum import Enum
from os.path import exists

from emantomo import Plugin
from emantomo.convert import convertBetweenHdfAndMrc
from emantomo.convert.lstConvert import EmanLstReader
from emantomo.objects import EmanSetOfParticles
from pwem.convert.headers import fixVolume
from pwem.objects import SetOfFSCs
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, StringParam, EnumParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from emantomo.constants import SYMMETRY_HELP_MSG, REFERENCE_NAME, SPT_00_DIR

# 3D maps filtering options
WIENER = 'wiener'
GLOBAL = 'global'
LOCAL_WIENER = 'localwiener'
LOCAL = 'local'
filteringKeys = [WIENER, GLOBAL, LOCAL_WIENER, LOCAL]
mapFilterDict = dict(zip(filteringKeys, range(len(filteringKeys))))


class EmanRefineNewOutputs(Enum):
    subtomograms = SetOfSubTomograms
    average = AverageSubTomogram
    FSCs = SetOfFSCs


class EmanProtTomoRefinementNew(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_refine_new.py* EMAN2 program.
    This refinement protocol performs subtomogram, subtilt and defocus refinement. The 2D subtilt particles are used
    instead of 3D subvolumes in the subtomogram refinement step. Moreover, this program now can model the localized 2D
    particle motion by considering the motion trajectory of each particle along with its neighbor.
    """

    _label = 'subtomogram refinement pppt'
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
                      help='Options to filter the 3D maps:\n\n'
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
                      label='Do local refine? (only for p iterations)',
                      help='Perform only local search around the solution from the previous alignment.')
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
        self._insertFunctionStep(self.createEmanPrjPostExtractionStep)
        self._insertFunctionStep(self.convertRefVolStep)
        self._insertFunctionStep(self.refineStep)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.inParticles = self.getAttrib(IN_SUBTOMOS)
        self.inSamplingRate = self.inParticles.getSamplingRate()
        self.noIters = self._getNoIters()
        self.avgHdf = self.getRefinedAverageFn(self.noIters)
        self.evenHdf = self.getRefineEvenFn(self.noIters)
        self.oddFnHdf = self.getRefineOddFn(self.noIters)
        self.ali2dLst = self.getAli2dFile(self.noIters)
        self.ali3dLst = self.getAli3dFile(self._getNoPIters())

    def refineStep(self):
        # In case of continuing from this step, the previous results dir will be removed to avoid EMAN creating one
        # for each execution (one for each continue)
        refineDir = self.getRefineDir()
        if exists(refineDir):
            shutil.rmtree(refineDir)
        self.buildEmanSets(outAliPath=None)
        program = Plugin.getProgram("e2spt_refine_new.py")
        self.runJob(program, self._genRefineCmd(), cwd=self._getExtraPath())

    def convertOutputStep(self):
        # Average and halves
        inFiles = [self.avgHdf, self.evenHdf, self.oddFnHdf]
        args = '--apix %.3f' % self.inSamplingRate
        for inFile in inFiles:
            outFile = inFile.replace('.hdf', '.mrc')
            convertBetweenHdfAndMrc(self, inFile, outFile, args)
            fixVolume(outFile)

    def createOutputStep(self):
        HDF = '.hdf'
        MRC = '.mrc'
        # Output 1: average with halves
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(self.avgHdf.replace(HDF, MRC))
        averageSubTomogram.setHalfMaps([self.evenHdf.replace(HDF, MRC), self.oddFnHdf.replace(HDF, MRC)])
        averageSubTomogram.setSamplingRate(self.inSamplingRate)

        # Output 2: aligned particles
        outParticles = EmanSetOfParticles.create(self._getPath(), template='emanParticles%s.sqlite')
        outParticles.copyInfo(self.inParticles)
        outParticles.setAli2dLstFile(self.ali2dLst)
        outParticles.setAli3dLstFile(self.ali3dLst)
        EmanLstReader.align3dLst2Scipion(self.ali3dLst, self.inParticles, outParticles)

        # Output 3: FSCs
        fscs = self.genFscs(self.noIters)

        # Define outputs and relations
        outputs = {self._possibleOutputs.average.name: averageSubTomogram,
                   self._possibleOutputs.subtomograms.name: outParticles,
                   self._possibleOutputs.FSCs.name: fscs}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inParticles, averageSubTomogram)
        self._defineSourceRelation(self.inParticles, outParticles)
        self._defineSourceRelation(self.inParticles, fscs)

    # --------------------------- UTILS functions ------------------------------
    def _genRefineCmd(self):
        new2dAlignFile = self.getNewAliFile(is3d=False)
        new3dAlignFile = self.getNewAliFile()
        # Input params
        if exists(self._getExtraPath(new3dAlignFile)):
            args = [f'--ptcls={new3dAlignFile}', f'--loadali3d']
        else:
            args = [f'--ptcls={self._getLstFile()}']
        if self.getRefVol():
            args.append(f'--ref={REFERENCE_NAME}.hdf')
        args.append(f'--startres={self.startRes.get():.2f}')
        # Refinement params
        args.append(f'--iters={self.iters.get()}')
        args.append(f'--sym={self.symmetry.get()}')
        args.append(f'--keep={self.pkeep.get():.2f}')
        args.append(f'--tophat={filteringKeys[self.topHat.get()]}')
        args.append(f'--maxres={self.maxResAli.get():.2f}')
        args.append(f'--minres={self.minResAli.get():.2f}')
        args.append(self._getGoldModeCmd(self.inParticles))
        if exists(self._getExtraPath(new2dAlignFile)):
            args.append(f'--loadali2d={new2dAlignFile}')
        # Local refinement params
        if self.doLocalRefine.get():
            args.append('--localrefine')
            args.append(f'--maxang={self.maxAng.get()}')
            args.append(f'--maxshift={self.maxShift.get()}')
            args.append(f'--smooth={self.smooth.get():.2f}')
            args.append(f'--smoothN={self.smoothN.get()}')
        # Extra params
        args.append(f'--parallel=thread:{self.numberOfThreads.get()}')
        args.append(f'--threads={self.threadsPostProc.get()}')
        if self.make3dThread.get():
            args.append('--m3dthread')
        if self.extraParams.get():
            args.append(self.extraParams.get())
        args.append('--verbose=9')

        return ' '.join(args)

    @staticmethod
    def _getGoldModeCmd(inParticles):
        """
        From EMAN doc:
            - goldstandard: "Phase randomize the reference to the starting resolution (--startres) independently for
              the even/odd subsets of particles".
            - goldcontinue: "Continue from previous gold standard refinement. Ues the _even/_odd version of the given
              reference".
        Thus, we'll apply it if it is the first refinement, which means that the set of particles introduced do not
        contain an ali2d nor ali3d files in the corresponding attributes.
        """
        if not inParticles.getAli2dLstFile() and not inParticles.getAli3dLstFile():
            return '--goldstandard'
        else:
            return '--goldcontinue'

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
        compressedList = iterStr.split(',')
        unrolledList = []
        for i in compressedList:
            if len(i) == 1:
                unrolledList.append(i)
            else:
                unrolledList += [i[0]] * int(i[1])
        return len(unrolledList) - unrolledList[::-1].index('p')

    # --------------------------- INFO functions --------------------------------
    # def _validate(self):
    #     errorMsgs = []
    #     inParticles = self.getAttrib(IN_SUBTOMOS)
    #     inRef = self.getAttrib(REF_VOL)
    #     if not self._getGoldModeCmd(inParticles) and inRef:
    #         if not inRef.hasHalfMaps():
    #             errorMsgs.append('If the introduced particles have some kind of alignment or orientation, the '
    #                              'initial volume introduced needs to have the corresponding even and odd halves. '
    #                              'Try to replace the initial volume by the resulting average of a previous refine '
    #                              'protocol.')
    #     return errorMsgs
