# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *              Matias Garcia   (matias@eyeseetea.com) [1]
# *              David Herreros  (dherreros@cnb.csic.es) [2]
# *
# * [1] EyeSeeTea Ltd, London, UK
# * [2] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
import os
from glob import glob

from pyworkflow import utils as pwutils, BETA
import pyworkflow.protocol.params as params

from pwem.protocols import EMProtocol
from pyworkflow.protocol import STEPS_PARALLEL

from ..convert import writeSetOfSubTomograms, refinement2Json, loadJson
import emantomo

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfAverageSubTomograms


class EmanProtTomoClassifySubtomos(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_classify.py* EMAN2 program.
    Protocol to performs a conventional iterative subtomogram averaging
    using the full set of particles.
    It will take a set of subtomograms (particles) and a set of Averages(reference)
    and 3D reconstruct a subtomogram.
    It also builds a set of subtomograms that contains the original particles
    plus the score, coverage and align matrix per subtomogram .
    """

    _outputClassName = 'MultiReferenceRefinement'
    _label = 'classify subtomos'
    OUTPUT_PREFIX = 'outputSetOf3DClassesSubTomograms'
    OUTPUT_DIR = "spt_00"
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True, label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('inputRef', params.PointerParam,
                      pointerClass='AverageSubTomogram',
                      default=None, label='Reference average',
                      help='If not provided, a reference will be created from the input subtomograms.')

        form.addSection(label='Optimization')
        form.addParam('mask', params.PointerParam,
                      label='Mask',
                      allowsNull=True,
                      pointerClass='VolumeMask',
                      help='Select a 3D Mask to be applied to the initial model')

        form.addParam('nClass', params.IntParam, default=2,
                       label='Number of classes')

        form.addParam('sym', params.StringParam, default='c1',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Symmetry',
                      help='Symmetry (Default: c1')

        form.addParam('maxRes', params.FloatParam, default=10,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum resolution',
                      help='Filter particles to this resolution (in Angstroms) before'
                           ' classification')

        form.addParam('nBasis', params.IntParam, default=3,
                      label='Number of basis',
                      help='Number of PCA basis vectors. Default is 3')

        form.addParam('wedgeFill', params.BooleanParam, default=False,
                       label='Fill missing wedge?',
                       help='Determine whether to fill the missing wedge before classification or not')

        form.addParam('clean', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Remove outliers?',
                      help='Remove outliers before PCA')

        form.addParam('shrink', params.FloatParam, allowsNull=True,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Shrink factor',
                       help='Shrink particles before classification')

    #--------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        #TODO: Get the basename.hdf from the inputSetOfSubTomogram
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('pcaClassification')
        self._insertFunctionStep('createOutputStep')

    #--------------- STEPS functions -----------------------
    def convertInputStep(self):
        storePath = self._getExtraPath("subtomograms")
        pwutils.makePath(storePath)

        volName = self.inputSetOfSubTomogram.get().getFirstItem().getVolName()
        self.newFn = pwutils.removeBaseExt(volName).split('__ctf')[0] + '.hdf'
        self.newFn = pwutils.join(storePath, self.newFn)
        writeSetOfSubTomograms(self.inputSetOfSubTomogram.get(), storePath)

        project_path = self._getExtraPath('spt_00')
        pwutils.makePath(project_path)
        refinement2Json(self, self.inputSetOfSubTomogram.get())

        program = emantomo.Plugin.getProgram('e2proclst.py')
        self.runJob(program, ' --create %s %s' % (os.path.abspath(os.path.join(project_path, 'input_ptcls.lst')),
                                                  os.path.abspath(self.newFn)),
                    cwd=self._getExtraPath())


        program = emantomo.Plugin.getProgram('e2proc3d.py')
        args = "%s %s" % (self.inputRef.get().getFileName(), self._getExtraPath(os.path.join('spt_00', 'threed_01.hdf')))
        self.runJob(program, args)

    def pcaClassification(self):
        """ Run the pca classification. """
        args = " --path=%s --iter=1 --nclass=%d --nbasis=%d" % \
               (os.path.abspath(self._getExtraPath('spt_00')), self.nClass.get(), self.nBasis.get())

        if self.mask.get() is None:
            args += " --mask=none"
        else:
            args += " --mask=%s" % self.mask.get().getFileName()

        args += " --sym=%s" % self.sym.get()

        if not self.wedgeFill.get():
            args += " --nowedgefill"

        if self.clean.get():
            args += " --clean"

        if self.shrink.get():
            args += " --shrink=%f" % self.shrink.get()

        program = emantomo.Plugin.getProgram('e2spt_pcasplit.py')
        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        # self._manageGeneratedFiles()
        subtomoSet = self.inputSetOfSubTomogram.get()
        classes3D = self._createSetOfClassesSubTomograms(subtomoSet)
        self._fillClassesFromJsons(classes3D)

        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(subtomoSet, classes3D)

        # Create a SetOfVolumes and define its relations
        volumes = SetOfAverageSubTomograms.create(self._getPath(),
                                                  template='avgSubtomograms%s.sqlite',
                                                  suffix='')
        volumes.setSamplingRate(subtomoSet.getSamplingRate())

        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(subtomoSet, volumes)

    #--------------- UTILS functions ------------------------
    def _fillClassesFromJsons(self, clsSet):
        self.particle_class = {}
        for classID in range(self.nClass.get()):
            keys_class = list(loadJson(self._getExtraPath(os.path.join("sptcls_00", "particle_parms_%02d.json" % (classID + 1)))).keys())
            partID = [int(item.split(", ")[1][:-1]) for item in keys_class]
            self.particle_class[classID] = partID

        # Initialization of variables used during clasiffy, specially in the callbacks: _updateParticle and _updateCLass
        self.particleCounter = 0 # Particle counter that should match the value in particles list input
        self.key_list = list(self.particle_class.keys())
        self.val_list = list(self.particle_class.values())

        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=itertools.count(0))

    def _updateParticle(self, item, row):

        idx = self.particleCounter
        position = [i for i, sublist in enumerate(self.val_list) if idx in sublist][0]
        item.setClassId(self.key_list[position] + 1)
        self.particleCounter += 1

    def _updateClass(self, item):
        classId = item.getObjId()
        item.getRepresentative().setLocation(self._getExtraPath(os.path.join("sptcls_00", "threed_%02d.hdf" % classId)))

    #--------------- INFO functions -------------------------
    def _summary(self):
        pass

    def _methods(self):
        pass
