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
import itertools
from enum import Enum
from os.path import exists, basename, join, abspath

import numpy as np

from emantomo import Plugin
from emantomo.convert.lstConvert import EmanLstReader
from emantomo.objects import EmanSetOfParticles
from pwem.convert.headers import fixVolume
from pwem.objects import SetOfFSCs
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, StringParam, EnumParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from tomo.objects import AverageSubTomogram, SetOfSubTomograms, SetOfAverageSubTomograms, SetOfClassesSubTomograms
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from emantomo.constants import SYMMETRY_HELP_MSG, THREED, ALI3D_BASENAME, ALI2D_BASENAME, SPTCLS_00_DIR, CLASS


class EmanMultiRefineNewOutputs(Enum):
    subtomograms = SetOfSubTomograms
    classes = SetOfClassesSubTomograms
    FSCs = SetOfFSCs


class EmanProtMultiRefinementNew(ProtEmantomoBase):
    """
    This protocol wraps *e2spt_refinemulti_new.py* EMAN2 program.
    Multi-reference classification for the new (2021) SPT refinement protocol.
    """

    _label = 'Multi-reference classification'
    _possibleOutputs = EmanMultiRefineNewOutputs
    maskRefOutName = 'maskRef'
    maskAlignOutName = 'maskAlign'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.convertedRefDict = None
        self.numClasses = None
        self.refFiles = None

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_SUBTOMOS, PointerParam,
                      pointerClass='EmanSetOfParticles',
                      label='Aligned particles',
                      important=True)
        form.addParam('maskRef', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask the references prior to to classif. (opt.)',
                      allowsNull=True)
        form.addParam('maxRes', FloatParam,
                      default=20,
                      important=True,
                      label='Max. resolution (Å)',
                      help='Maximum resolution (the smaller number) to consider in alignment (in Å).\n'
                           'Since gold-standard validation is not used here, setting this parameter is mandatory.')
        form.addParam('minRes', FloatParam,
                      default=0,
                      label='Min. resolution (Å)',
                      help='Minimum resolution (the larger number) to consider in alignment (in Å).')

        form.addSection(label='Classification')
        group = form.addGroup('Reference management')
        group.addLine('Reference managing cases info --->',
                      important=True,
                      help='1) If no. classes > 0, the initial references will be generated by random classification.\n'
                           '2) If a volume is provided and no. classes > 0, it will be duplicated the first ref N times '
                           'with phase randomization at 2 x max. resolution.\n'
                           '3) If no. classes < 0, the volume or set of volumes provided will be loaded as the references.')
        group.addParam('nClasses', IntParam,
                       default=-1,
                       label="Number of classes")
        group.addParam(REF_VOL, PointerParam,
                       pointerClass='Volume, SetOfVolumes',
                       allowsNull=True,
                       label="Reference volume (opt.)")
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
                      condition='doAlignment',
                      label='Apply mask to the 3D alignment ref. in each iter. (opt.)',
                      allowsNull=True,
                      help="Not applied to the average, which will follow normal EMAN's masking routine.")
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
        self._initialize()
        self._insertFunctionStep(self.createEmanPrjPostExtractionStep)
        self._insertFunctionStep(self.convertReferencesStep)
        self._insertFunctionStep(self.buildEmanSetsStep)
        self._insertFunctionStep(self.refineMultiStep)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.inParticles = self.getAttrib(IN_SUBTOMOS)
        self.inSamplingRate = self.inParticles.getSamplingRate()
        self.numClasses = self.getNoRefs()
        self.refFiles = self.getReferenceFiles()
        # if refFiles:
        #     self.convertedRefDict = {refFile: basename(refFile) for refFile in refFiles}

    def convertReferencesStep(self):
        pass
        # maskRef = self.maskRef.get()
        # maskAlign = self.maskAlign.get()
        # if maskRef:
        #     self.convertOrLink(abspath(maskRef.getFileName()), self.maskRefOutName, self._getExtraPath(),
        #                        self.inSamplingRate)
        # if maskAlign:
        #     self.convertOrLink(abspath(maskAlign.getFileName()), self.maskAlignOutName, self._getExtraPath(),
        #                        self.inSamplingRate)
        # if self.convertedRefDict:
        #     for refFile, baseName in self.convertedRefDict.items():
        #         self.convertOrLink(refFile, baseName, '', self.inSamplingRate)

    def refineMultiStep(self):
        program = Plugin.getProgram("e2spt_refinemulti_new.py")
        self.runJob(program, self._genRefineMultiCmd(), cwd=self._getExtraPath())

    def convertOutputStep(self):
        for classId in range(self.numClasses):
            inFile = join(SPTCLS_00_DIR, self.getOutputThreed(classId, 'hdf', onlyBaseName=True))
            outFile = inFile.replace('.hdf', '.mrc')
            self.convertBetweenHdfAndMrc(inFile, outFile, extraArgs=f'--apix {self.inSamplingRate:.3f}')
            fixVolume(outFile)

    def createOutputStep(self):
        inParticlesPointer = self.getAttrib(IN_SUBTOMOS, getPointer=True)
        outParticles = EmanSetOfParticles.create(self._getPath(), template='emanParticles%s.sqlite')
        outParticles.copyInfo(self.inParticles)
        align2dFiles = []
        align3dFiles = []
        for classId in range(self.numClasses):
            # Get the 2d and 3d align files
            align2dFiles.append(self.getAlign2dFile(classId))
            align3dFiles.append(self.getAlign3dFile(classId))

        # Output 1: aligned particles
        # outParticles.setAli2dLstFile(self.ali2dLst)
        # outParticles.setAli3dLstFile(self.ali3dLst)
        align3dData = EmanLstReader.align3dLst2Scipion(align3dFiles, self.inParticles, outParticles)

        # Output 2: Classes subtomograms
        classes3d = SetOfClassesSubTomograms.create(self._getPath(), template='subtomogramClasses%s.sqlite')
        classes3d.setImages(inParticlesPointer)
        classes3d.classifyItems(updateItemCallback=self._updateParticle,
                                updateClassCallback=self._updateClass,
                                itemDataIterator=iter(align3dData))

        # TODO: check if the FSCs and halves are generated if alignment is selected

        # Define outputs and relations
        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: outParticles,
                               self._possibleOutputs.classes.name: classes3d})
        self._defineSourceRelation(inParticlesPointer, outParticles)
        self._defineSourceRelation(inParticlesPointer, classes3d)

    # --------------------------- UTILS functions ------------------------------
    def _genRefineMultiCmd(self):
        hdf = '.hdf'
        refMask = self.maskRef.get()
        alignMask = self.maskAlign.get()
        breakSym = self.breakSym.get()
        new2dAlignFile = self.getNewAliFile(is3d=False)
        nClasses = self.numClasses
        args = [
            ' '.join(self.refFiles) if self.refFiles else '',
            f'--ptcls {self.getNewAliFile()}',
            f'--niter {self.nIter.get()}',
            f'--sym {self.symmetry.get()}',
            f'--parallel thread:{self.numberOfThreads.get()}',
            f'--threads {self.threadsPostProc.get()}',
            f'--maxres {self.maxRes.get():.2f}',
            f'--minres {self.minRes.get():.2f}',
            '--loadali3d'
        ]
        if nClasses > 0:
            args.append(f'--nref {nClasses}')
        if exists(self._getExtraPath(new2dAlignFile)):
            args.append(f'--loadali2d {new2dAlignFile}')
        if refMask:
            # args.append(f'--maskref {self.maskRefOutName + hdf}')
            args.append(f'--maskref {abspath(refMask.getFileName())}')

        if breakSym:
            args.append(f'--breaksym {breakSym}')
        if self.doAlignment.get():
            args.append(f'--maxang {self.maxAng.get()}')
            args.append(f'--maxshift {self.maxShift.get()}')
            if alignMask:
                # args.append(f'--maskalign {self.maskAlignOutName + hdf}')
                args.append(f'--maskalign {abspath(alignMask.getFileName())}')
        else:
            args.append('--skipali')
        if self.make3dThread.get():
            args.append('--m3dthread')
        return ' '.join(args)

    def getReferenceFiles(self):
        refs = self.getAttrib(REF_VOL)
        refList = None
        if refs:
            refList = [abspath(vol.getFileName()) for vol in refs]
        return refList

    def getNoRefs(self):
        nClasses = self.nClasses.get()
        refs = self.getAttrib(REF_VOL)
        return nClasses if nClasses > 0 else len(refs)

    def getOutputAvgFile(self, classNum):
        return self.getOutputFile(classNum, THREED, 'hdf')

    def getAlign2dFile(self, classNum):
        return self.getOutputFile(classNum, ALI2D_BASENAME, 'lst')

    def getAlign3dFile(self, classNum):
        return self.getOutputFile(classNum, ALI3D_BASENAME, 'lst')

    def getOutputFile(self, classNum, baseName, ext, onlyBaseName=False):
        nIters = self.nIter.get()
        fName = f'{baseName}{nIters:02}_{classNum:02}.{ext}'
        return fName if onlyBaseName else join(self.getMultiRefineDir(), fName)

    def getOutputThreed(self, classNum, ext, onlyBaseName=False):
        return self.getOutputFile(classNum, THREED + '_', ext, onlyBaseName=onlyBaseName)

    @staticmethod
    def _updateParticle(item, row):
        item.setClassId(row[CLASS])

    def _updateClass(self, item):
        item.setAlignment3D()
        representative = item.getRepresentative()
        representative.setLocation(self.getOutputThreed(item.getObjId(), 'mrc'))

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        errorMsg = []
        inTrMatrix = self.getAttrib(IN_SUBTOMOS).getFirstItem().getTransform().getMatrix()
        if np.allclose(inTrMatrix, np.eye(4), atol=1e-4):
            errorMsg.append('No alignment was detected in the introduced particles.')
        return errorMsg
