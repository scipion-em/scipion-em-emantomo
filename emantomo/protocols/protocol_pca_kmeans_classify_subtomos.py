# coding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team  (scipion@cnb.csic.es) [1]
# *              David Herreros  (dherreros@cnb.csic.es) [1]
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
import itertools
from enum import Enum
from os.path import join, abspath
from pwem import ALIGN_3D
from pyworkflow import BETA
from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam, IntParam, StringParam, LEVEL_ADVANCED, FloatParam, \
    BooleanParam
from pyworkflow.utils import makePath, removeBaseExt, removeExt
from ..constants import SYMMETRY_HELP_MSG, SUBTOMOGRAMS_DIR, SPT_00_DIR, INPUT_PTCLS_LST, SPTCLS_00_DIR
from ..convert import writeSetOfSubTomograms, refinement2Json, loadJson
import emantomo
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SetOfClassesSubTomograms, SetOfAverageSubTomograms


class pcaOutputObjects(Enum):
    subtomograms = SetOfSubTomograms
    classes = SetOfClassesSubTomograms
    representatives = SetOfAverageSubTomograms


class EmanProtPcaKMeansClassifySubtomos(EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_pcasplit.py* EMAN2 program.
    Protocol to performs a PCA (Principal Component Analysis) and K-means classification
    using the full set of particles.
    """

    _label = 'PCA-K Means classification of subtomograms'
    _devStatus = BETA
    _possibleOutputs = pcaOutputObjects

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.subtomosDir = None
        self.spt00Dir = None
        self.subtomoClassDict = {}

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Subtomograms')
        form.addParam('mask', PointerParam,
                      label='Mask (opt.)',
                      allowsNull=True,
                      pointerClass='VolumeMask')

        form.addSection(label='Optimization')
        form.addParam('nClass', IntParam,
                      default=2,
                      label='Number of classes')
        form.addParam('sym', StringParam,
                      default='c1',
                      expertLevel=LEVEL_ADVANCED,
                      label='Symmetry',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('maxRes', FloatParam,
                      default=30,
                      expertLevel=LEVEL_ADVANCED,
                      label='Maximum resolution [Å]',
                      help='Filter particles to this resolution (in Angstroms) before'
                           ' classification')
        form.addParam('nBasis', IntParam,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of PCA basis vectors',
                      help='Number of PCA components. If 0 or None, then they will be automatically estimated '
                           'considering the size of the input subtomograms.')
        form.addParam('wedgeFill', BooleanParam,
                      default=False,
                      label='Fill missing wedge?',
                      help='Determine whether to fill the missing wedge before classification or not')
        form.addParam('clean', BooleanParam,
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Remove outliers?',
                      help='Remove outliers before PCA')
        form.addParam('shrink', IntParam,
                      default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Shrink factor',
                      help='Shrink particles before classification.')

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.pcaClassification)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.subtomosDir = self._getExtraPath(SUBTOMOGRAMS_DIR)
        self.spt00Dir = self._getExtraPath(SPT_00_DIR)
        makePath(*[self.subtomosDir, self.spt00Dir])

    def convertInputStep(self):
        volName = self.inSubtomos.get().getFirstItem().getVolName()
        stackHdf = join(self.subtomosDir, removeBaseExt(volName).split('__ctf')[0] + '.hdf')

        # Convert the particles to HDF if necessary
        writeSetOfSubTomograms(self.inSubtomos.get(), self.subtomosDir, lignType=ALIGN_3D)

        # Generate the refinement JSON expected to be in the SPT directory
        refinement2Json(self, self.inSubtomos.get())

        # Generate the LST file expected to be in the SPT directory
        program = emantomo.Plugin.getProgram('e2proclst.py')
        args = ' --create %s %s' % (join(self.spt00Dir, INPUT_PTCLS_LST), abspath(stackHdf))
        self.runJob(program, args)

        # Average the introduced particles to get the refined expected to be in the SPT directory
        args = " --path=%s --keep 1" % self.spt00Dir
        program = emantomo.Plugin.getProgram('e2spt_average.py')
        self.runJob(program, args)

    def pcaClassification(self):
        """ Run the pca classification. """
        args = " --path=%s" % SPT_00_DIR
        args += " --nclass=%i" % self.nClass.get()
        args += " --nbasis=%d" % self._estimatePCAComps()
        args += " --sym=%s" % self.sym.get()
        if self.mask.get() is None:
            args += " --mask=none"
        else:
            args += " --mask=%s" % abspath(self.mask.get().getFileName())
        if not self.wedgeFill.get():
            args += " --nowedgefill"
        if self.clean.get():
            args += " --clean"
        if self.shrink.get():
            args += " --shrink=%i" % self.shrink.get()

        args += ' --iter=1 --verbose=9'
        program = emantomo.Plugin.getProgram('e2spt_pcasplit.py')
        self.runJob(program, args, cwd=self._getExtraPath())

    def convertOutputStep(self):
        """Convert the class representatives and their halves from HDF to MRC"""
        # They follow the syntax, for each class (e. g. class 1), threed_01.hdf, threed_01_even.hdf, threed_01_odd.hdf,
        # located in directory extra/sptcls_00/
        program = emantomo.Plugin.getProgram('e2proc3d.py')
        sRate = self.inSubtomos.get().getSamplingRate()
        for classFile in glob.glob(self._getExtraPath(SPTCLS_00_DIR, 'threed_*.hdf')):
            outMrcFile = self._getExtraPath(SPTCLS_00_DIR, removeBaseExt(classFile) + '.mrc')
            args = "--apix %f %s %s" % (sRate, classFile, outMrcFile)
            self.runJob(program, args)

    def createOutputStep(self):
        # 1) Subtomograms
        inSubtomoSet = self.inSubtomos.get()
        outSubtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
        outSubtomoSet.copyInfo(inSubtomoSet)

        # Generate a dictionary with the classId and the corresponding particles determined by eman
        for classID in range(self.nClass.get()):
            keys_class = list(loadJson(
                self._getExtraPath(join("sptcls_00", "particle_parms_%02d.json" % (classID + 1)))).keys())
            partID = [int(item.split(", ")[1][:-1]) for item in keys_class]
            self.subtomoClassDict[classID] = partID

        # Fill the output set of subtomograms, adding the classId to the input ones
        for i, subtomo in enumerate(inSubtomoSet):
            for classId, idList in self.subtomoClassDict.items():
                if i in idList:
                    subtomo.setClassId(classId + 1)
                    outSubtomoSet.append(subtomo)
                    continue

        # 2) Classes subtomograms
        classes3D = self._createSetOfClassesSubTomograms(inSubtomoSet)
        classes3D.setImages(inSubtomoSet)
        self._fillClassesFromJsons(classes3D)

        # 3) Set of averages with the representative of each class
        # Create a SetOfVolumes and define its relations
        averages = SetOfAverageSubTomograms.create(self._getPath(), template='avgSubtomograms%s.sqlite', suffix='')
        averages.setSamplingRate(inSubtomoSet.getSamplingRate())
        for class3D in classes3D:
            representative = class3D.getRepresentative()
            classId = class3D.getObjId()
            representative.setObjId(classId)  # The class objId is the class number
            representative.setClassId(classId)
            representative.setHalfMaps(self.getHalvesFromClassRepresentative(representative))
            averages.append(representative)

        # Define the outputs and relations
        self._defineOutputs(**{pcaOutputObjects.classes.name: classes3D,
                               pcaOutputObjects.subtomograms.name: outSubtomoSet,
                               pcaOutputObjects.representatives.name: averages})
        self._defineSourceRelation(inSubtomoSet, outSubtomoSet)
        self._defineSourceRelation(inSubtomoSet, classes3D)
        self._defineSourceRelation(inSubtomoSet, averages)

    # --------------- UTILS functions ------------------------
    def _fillClassesFromJsons(self, clsSet):
        # Initialization of variables used during clasiffy, specially in the callbacks: _updateParticle and _updateCLass
        self.particleCounter = 0  # Particle counter that should match the value in particles list input
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=itertools.count(0))

    def _updateParticle(self, item, row):
        idx = self.particleCounter
        for classId, particleIndList in self.subtomoClassDict.items():
            if idx in particleIndList:
                item.setClassId(classId + 1)
                break
        self.particleCounter += 1

    def _updateClass(self, item):
        classId = item.getObjId()
        item.setAlignment3D()
        representative = item.getRepresentative()
        representative.setLocation(self._getExtraPath(join(SPTCLS_00_DIR, "threed_%02d.mrc" % classId)))
        representative.setHalfMaps(self.getHalvesFromClassRepresentative(representative))

    def _estimatePCAComps(self):
        # Antonio Martinez-Sanchez rule... seems to work fine
        pcaComps = self.nBasis.get()
        if not pcaComps:
            x, y, _ = self.inSubtomos.get().getDimensions()
            pcaComps = round(1 / 100 * 0.5 * x * y)

        return pcaComps

    @staticmethod
    def getHalvesFromClassRepresentative(representative):
        ext = '.mrc'
        repPathAndBasename = removeExt(representative.getFileName())
        return [repPathAndBasename + '_even' + ext, repPathAndBasename + '_odd' + ext]

    # --------------- INFO functions -------------------------
    def _summary(self):
        pass

    def _methods(self):
        pass

    def _validate(self):
        errors = []
        nSubtomos = self.inSubtomos.get().getSize()
        nClasses = self.nClass.get()
        nPcaComps = self.nBasis.get()
        if nClasses > nSubtomos:
            errors.append('The number of classes (%i) cannot be greater than the number of subtomograms (%i).' %
                          (nClasses, nSubtomos))
        if nPcaComps:
            if nPcaComps > nSubtomos:
                errors.append('The number of PCA components must be between 0 and min(n_samples, n_features).')
        return errors
