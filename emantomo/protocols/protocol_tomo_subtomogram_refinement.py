# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *              Ignacio del Cano  (idelcano@eyeseetea.com) [1]
# *              Arnau Sanchez  (arnau@eyeseetea.com) [1]
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


from glob import glob
import os
import re

from pyworkflow import BETA
from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
# from pyworkflow.protocol import STEPS_PARALLEL

from pwem.protocols import EMProtocol

from emantomo.convert import writeSetOfSubTomograms, getLastParticlesParams, updateSetOfSubTomograms
import emantomo

from tomo.protocols import ProtTomoBase
from tomo.objects import AverageSubTomogram, SetOfSubTomograms, SetOfAverageSubTomograms

from .. import SCRATCHDIR

SAME_AS_PICKING = 0


class EmanProtTomoRefinement(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_refine.py* EMAN2 program.

    Protocol to performs a conventional iterative subtomogram averaging
    using the full set of particles.
    It will take a set of subtomograms (particles) and a subtomogram(reference,
    potentially comming from the initial model protocol)
    and 3D reconstruct a subtomogram.
    It also builds a set of subtomograms that contains the original particles
    plus the score, coverage and align matrix per subtomogram .
    """

    _outputClassName = 'SubTomogramRefinement'
    _label = 'subtomogram refinement'
    _devStatus = BETA
    OUTPUT_PREFIX = 'outputSetOfClassesSubTomograms'
    OUTPUT_DIR = "spt_00"

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True, label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('inputRef', params.PointerParam,
                      pointerClass='Volume', allowsNull=True,
                      default=None, label='Input Ref Tomogram',
                      help='3D reference for initial model generation.'
                           'No reference is used by default.')

        form.addSection(label='Optimization')
        form.addParam('niter', params.IntParam, default=5,
                      label='Number of iterations',
                      help='The number of iterations to perform.')
        form.addParam('mass', params.FloatParam, default=500.0,
                      label='Mass:',
                      help='Default=500.0')
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
        form.addParam('maskFile', params.PointerParam, allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      pointerClass='VolumeMask', label='Mask file',
                      help='Select the mask object')
        form.addParam('setsf', params.PointerParam, allowsNull=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      pointerClass='VolumeMask', label='Structure factor',
                      help='Select the structure factor')
        form.addParam('sym', params.StringParam, default='c1',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Symmetry',
                      help='Symmetry (Default: c1')
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

        form.addParallelSection(threads=4, mpi=1)

    # --------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        # TODO: Get the basename.hdf from the inputSetOfSubTomogram
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('refinementSubtomogram')
        # TODO: Set and show the output
        self._insertFunctionStep('createOutputStep')

    # --------------- STEPS functions -----------------------
    def convertInputStep(self):
        storePath = self._getExtraPath("subtomograms")
        pwutils.makePath(storePath)
        writeSetOfSubTomograms(self.inputSetOfSubTomogram.get(), storePath)
        self.newFn = glob(os.path.join(storePath, '*.hdf'))[0]

    def refinementSubtomogram(self):
        """ Run the Subtomogram refinement. """
        args = ' %s' % self.newFn
        if self.inputRef.get() is not None:
            reference = self.inputRef.get().getFileName()
            reference = reference.split(":")[0]
            args += (' --reference=%s ' % reference)
        args += (' --mass=%f' % self.mass)
        args += ' --goldstandard=%d ' % self.goldstandard
        args += ' --pkeep=%f ' % self.pkeep
        args += ' --sym=%s ' % self.sym
        args += ' --maxtilt=%s ' % self.maxtilt
        args += ' --path=%s ' % self.getOutputPath()
        if self.niter > 1:
            args += ' --niter=%d' % self.niter
        if self.goldcontinue:
            args += ' --goldcontinue '
        if self.localfilter:
            args += ' --localfilter '
        if self.numberOfMpi > 1:
            args += ' --parallel=mpi:%d:%s' % (self.numberOfMpi.get(), SCRATCHDIR)
        else:
            args += ' --parallel=thread:%d' % self.numberOfThreads.get()
        args += ' --threads=%d' % self.numberOfThreads.get()

        program = emantomo.Plugin.getProgram('e2spt_refine.py')
        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args)

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob(self.getOutputPath("*"))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, os.path.basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]

    def createOutputStep(self):
        lastImage = self.getLastFromOutputPath("threed_\d+.hdf")
        inputSetOfSubTomograms = self.inputSetOfSubTomogram.get()

        # Output 1: AverageSubTomogram
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(lastImage)
        setOfAverageSubTomograms = self._createSet(SetOfAverageSubTomograms, 'subtomograms%s.sqlite', "")
        setOfAverageSubTomograms.copyInfo(inputSetOfSubTomograms)
        setOfAverageSubTomograms.append(averageSubTomogram)

        # Output 2: setOfSubTomograms
        particleParams = getLastParticlesParams(self.getOutputPath())
        outputSetOfSubTomograms = self._createSet(SetOfSubTomograms, 'subtomograms%s.sqlite', "particles")
        outputSetOfSubTomograms.copyInfo(inputSetOfSubTomograms)
        outputSetOfSubTomograms.setCoordinates3D(inputSetOfSubTomograms.getCoordinates3D())
        updateSetOfSubTomograms(inputSetOfSubTomograms, outputSetOfSubTomograms, particleParams)

        self._defineOutputs(averageSubTomogram=setOfAverageSubTomograms, outputParticles=outputSetOfSubTomograms)
        self._defineSourceRelation(self.inputSetOfSubTomogram, setOfAverageSubTomograms)
        self._defineSourceRelation(self.inputSetOfSubTomogram, outputSetOfSubTomograms)

    def getOutputPath(self, *args):
        return os.path.join(self._getExtraPath(self.OUTPUT_DIR, *args))

    def getOutputFile(self, folderpattern, folder, files, pattern):
        pattern = "^" + folderpattern + pattern
        outputList = list()
        for file in files:
            if re.match(pattern, file) is not None:
                outputList.append(file.replace(folder, ""))
        lastIteration = max(re.findall(r'\d+', ''.join(outputList)))

        output = [file for file in outputList if lastIteration in file]
        return folder + output.pop()

    def getLastOutputFolder(self, files):
        folder = "./spt_"
        validFolders = [file for file in files if folder in file]
        folderSuffix = max(re.findall(r'\d+', ''.join(validFolders)))
        folder = folder + folderSuffix
        return folder

    # --------------- INFO functions -------------------------
    def _summary(self):
        summary = []
        summary.append("Set Of SubTomograms source: %s" % (self.inputSetOfSubTomogram.get().getFileName()))

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
