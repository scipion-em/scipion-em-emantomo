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
import enum
from glob import glob
import re
from os.path import abspath, basename, join
from pwem.objects import SetOfFSCs
from pyworkflow import BETA
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from emantomo.convert import writeSetOfSubTomograms, getLastParticlesParams, updateSetOfSubTomograms, emanFSCsToScipion
import emantomo
import pwem.constants as emcts
from pyworkflow.utils import createLink
from tomo.protocols import ProtTomoBase
from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from ..constants import SUBTOMOGRAMS_DIR, SPT_00_DIR, SYMMETRY_HELP_MSG


class EmanTomoRefinementOutputs(enum.Enum):
    subtomograms = SetOfSubTomograms
    subtomogramAverage = AverageSubTomogram
    FSCs = SetOfFSCs


class EmanProtTomoRefinement(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_refine.py* EMAN2 program.

    Protocol to performs a conventional iterative subtomogram averaging
    using the full set of particles.
    It will take a set of subtomograms (particles) and a subtomogram(reference,
    potentially coming from the initial model protocol)
    and 3D reconstruct a subtomogram.
    It also builds a set of subtomograms that contains the original particles
    plus the score, coverage and align matrix per subtomogram .
    """

    _outputClassName = 'SubTomogramRefinement'
    _label = 'subtomogram refinement'
    _devStatus = BETA
    _possibleOutputs = EmanTomoRefinementOutputs

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.refFileIn = None
        self.alignType = None
        self.newFn = None
        self.refFileOut = None

    # --------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('inputRef', params.PointerParam,
                      pointerClass='Volume',
                      allowsNull=False,
                      default=None,
                      label='Input Ref SubTomogram',
                      help='3D reference for initial model generation.'
                           'No reference is used by default.')

        form.addSection(label='Optimization')
        form.addParam('niter', params.IntParam, default=5,
                      label='Number of iterations',
                      help='The number of iterations to perform.')
        form.addParam('mass', params.FloatParam, default=-1,
                      label='Mass:',
                      help='Mass normalization. Default=-1 Ignores mass')
        form.addParam('pkeep', params.FloatParam, default=0.8,
                      label='Particle keep:',
                      help='Fraction of particles to keep')
        form.addParam('goldstandard', params.IntParam, default=-1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Gold standard:',
                      help='initial resolution for gold standard refinement')
        form.addParam('goldcontinue', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Gold continue',
                      help='continue from an existing gold standard refinement')
        form.addParam('maskFile', params.PointerParam,
                      allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      pointerClass='VolumeMask',
                      label='Mask file',
                      help='Mask file to be applied to initial model')
        form.addParam('setsf', params.PointerParam, allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      pointerClass='VolumeMask', label='Structure factor',
                      help='Select the structure factor')
        form.addParam('sym', params.StringParam, default='c1',
                      label='Symmetry',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('localfilter', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Local filter',
                      help='use tophat local')
        form.addParam('maxtilt', params.FloatParam, default=90.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='maxtilt',
                      help='Explicitly zeroes data beyond specified tilt angle.'
                           'Assumes tilt axis exactly on Y and zero tilt in X-Y'
                           'plane. Default 90 (no limit).')
        form.addParam('useAlign', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Use previous alignments?')
        form.addParam('maxAng', params.FloatParam, default=25,
                      condition='useAlign==True',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum angular change',
                      help='Maximum anglular difference in refine mode (in degrees)')
        form.addParam('extraParams', params.StringParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Extra params",
                      help='Here you can add any extra parameters to run Eman subtomogram refinement. '
                           'Parameters should be written in Eman command line format (--param=val)')

        form.addParallelSection(threads=4, mpi=0)

    # --------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        # TODO: Get the basename.hdf from the inputSetOfSubTomogram
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.refinementSubtomogram)
        # TODO: Set and show the output
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def convertInputStep(self):
        storePath = self._getExtraPath(SUBTOMOGRAMS_DIR)
        pwutils.makePath(storePath)
        if self.useAlign.get():
            self.alignType = emcts.ALIGN_3D
        else:
            self.alignType = emcts.ALIGN_NONE
        writeSetOfSubTomograms(self.inputSetOfSubTomogram.get(), storePath, alignType=self.alignType)
        self.newFn = glob(join(storePath, '*.hdf'))[0]

        # Fix the sampling rate as it might be set wrong
        inSubtomos = self.inputSetOfSubTomogram.get()
        sampling_rate = inSubtomos.getSamplingRate()
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        args = "--apix %f %s %s" % (sampling_rate, self.newFn, self.newFn)
        self.runJob(program, args)
        self.refFileIn = self.inputRef.get().getFileName()
        self.refFileOut = self._getExtraPath('refVol.hdf')
        # The same with the reference volume
        if self.refFileIn.endswith('.mrc'):
            args = "--apix %f %s %s" % (sampling_rate, self.refFileIn, self.refFileOut)
            self.runJob(program, args)
        else:
            createLink(self.refFileIn, self.refFileOut)

    def refinementSubtomogram(self):
        """Run the Subtomogram refinement. """
        args = ' %s' % abspath(self.newFn)
        args += ' --verbose=9'
        args += (' --reference=%s ' % abspath(self.refFileOut))
        args += (' --mass=%f' % self.mass.get())
        args += ' --goldstandard=%d ' % self.goldstandard.get()
        args += ' --pkeep=%f ' % self.pkeep.get()
        args += ' --sym=%s ' % self.sym.get()
        args += ' --maxtilt=%s ' % self.maxtilt
        args += ' --path=%s ' % abspath(self.getOutputPath())
        args += ' --niter=%d' % self.niter.get()
        if self.goldcontinue:
            args += ' --goldcontinue '
        if self.maskFile.get():
            args += ' --mask=%s' % abspath(self.maskFile.get().getFileName())
        if self.localfilter:
            args += ' --localfilter '
        # if self.numberOfMpi > 1:
        #     args += ' --parallel=mpi:%d:%s' % (self.numberOfMpi.get(), SCRATCHDIR)
        # else:
        #     args += ' --parallel=thread:%d' % self.numberOfThreads.get()
        if self.alignType != emcts.ALIGN_NONE:
            args += ' --refine --maxang=%f' % self.maxAng.get()
        if self.extraParams.get():
            args += ' ' + self.extraParams.get()
        args += ' --threads=%d' % self.numberOfThreads.get()

        program = emantomo.Plugin.getProgram('e2spt_refine.py')
        self.runJob(program, args)

    def convertOutputStep(self):
        self.hdfToMrc("threed_\d+.hdf", self.getAverageFn())
        self.hdfToMrc("threed_\d+_even.hdf", self.getEvenFn())
        self.hdfToMrc("threed_\d+_odd.hdf", self.getOddFn())

    def createOutputStep(self):
        inputSetOfSubTomograms = self.inputSetOfSubTomogram.get()
        # Output 1: AverageSubTomogram
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(self.getAverageFn())
        averageSubTomogram.setHalfMaps([self.getEvenFn(), self.getOddFn()])
        averageSubTomogram.setSamplingRate(inputSetOfSubTomograms.getSamplingRate())
        # Output 2: setOfSubTomograms
        particleParams = getLastParticlesParams(self.getOutputPath())
        outputSetOfSubTomograms = self._createSet(SetOfSubTomograms, 'subtomograms%s.sqlite', "particles")
        outputSetOfSubTomograms.copyInfo(inputSetOfSubTomograms)
        outputSetOfSubTomograms.setCoordinates3D(inputSetOfSubTomograms.getCoordinates3D())
        updateSetOfSubTomograms(inputSetOfSubTomograms, outputSetOfSubTomograms, particleParams)
        # Output 3: FSCs
        fscs= self._createSet(SetOfFSCs, 'fsc%s.sqlite', "")
        fscMasked = self.getLastFromOutputPath('fsc_masked_\d+.txt')
        fscUnmasked = self.getLastFromOutputPath('fsc_unmasked_\d+.txt')
        fscTight = self.getLastFromOutputPath('fsc_maskedtight_\d+.txt')

        emanFSCsToScipion(fscs, fscMasked, fscUnmasked, fscTight)
        outputs = {EmanTomoRefinementOutputs.subtomogramAverage.name: averageSubTomogram,
                   EmanTomoRefinementOutputs.subtomograms.name: outputSetOfSubTomograms,
                   EmanTomoRefinementOutputs.FSCs.name: fscs}

        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputSetOfSubTomogram, averageSubTomogram)
        self._defineSourceRelation(self.inputSetOfSubTomogram, outputSetOfSubTomograms)

    # --------------- INFO functions -------------------------
    def _summary(self):
        summary = ["Subtomograms source: %s" % (self.inputSetOfSubTomogram.get().getFileName())]

        if self.inputRef.get() is not None:
            summary.append("Referenced Tomograms source: %s" % (self.inputRef.get().getFileName()))

        if self.getOutputsSize() >= 1:
            summary.append("Subtomogram refinement Completed")
        else:
            summary.append("Subtomogram refinement not ready yet.")

        return summary

    def _methods(self):
        inputSetOfSubTomgrams = self.inputSetOfSubTomogram.get()
        return [
            "Applied refinement using e2spt_refine (stochastic gradient descent)",
            "A total of %d particles of dimensions %s were used"
            % (inputSetOfSubTomgrams.getSize(), inputSetOfSubTomgrams.getDimensions()),
        ]

    # --------------------------- UTILS functions ------------------------------
    def hdfToMrc(self, pattern, mrcName):
        # Fix the sampling rate as it might be set wrong
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        lastImage = self.getLastFromOutputPath(pattern)
        args = "--apix %f %s %s" % (self.inputSetOfSubTomogram.get().getSamplingRate(),
                                    lastImage, mrcName)
        self.runJob(program, args)

    def getAverageFn(self):
        return self._getExtraPath("Average_refined.mrc")

    def getEvenFn(self):
        return self._getExtraPath("even.mrc")

    def getOddFn(self):
        return self._getExtraPath("odd.mrc")

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob(self.getOutputPath("*"))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]

    def getOutputPath(self, *args):
        return join(self._getExtraPath(SPT_00_DIR, *args))

    # @staticmethod
    # def getOutputFile(folderpattern, folder, files, pattern):
    #     pattern = "^" + folderpattern + pattern
    #     outputList = list()
    #     for file in files:
    #         if re.match(pattern, file) is not None:
    #             outputList.append(file.replace(folder, ""))
    #     lastIteration = max(re.findall(r'\d+', ''.join(outputList)))
    #
    #     output = [file for file in outputList if lastIteration in file]
    #     return folder + output.pop()
    #
    # @staticmethod
    # def getLastOutputFolder(files):
    #     folder = "./spt_"
    #     validFolders = [file for file in files if folder in file]
    #     folderSuffix = max(re.findall(r'\d+', ''.join(validFolders)))
    #     folder = folder + folderSuffix
    #     return folder
