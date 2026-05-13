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
from os.path import exists, join, abspath
from emantomo import Plugin
from emantomo.convert.lstConvert import EmanLstReader
from emantomo.objects import EmanSetOfParticles
from pwem.convert.headers import fixVolume
from pwem.emlib.image import ImageHandler
from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, StringParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from tomo.objects import SetOfSubTomograms, SetOfClassesSubTomograms, SubTomogram, AverageSubTomogram
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL
from emantomo.constants import SYMMETRY_HELP_MSG, THREED, ALI3D_BASENAME, ALI2D_BASENAME, SPTCLS_00_DIR, CLASS, \
    SPT_00_DIR


class EmanMultiRefineNewOutputs(Enum):
    subtomograms = SetOfSubTomograms
    classes = SetOfClassesSubTomograms


class EmanProtMultiRefinementNew(ProtEmantomoBase):
    """
    Performs multi-reference subtomogram classification and refinement using
    the EMAN2 subtomogram averaging framework. The protocol is designed to
    separate heterogeneous particle populations into structurally meaningful
    classes while preserving alignment information and producing refined
    three-dimensional class averages. It is particularly useful in cryo-electron
    tomography workflows where multiple conformational states, compositional
    differences, or flexible assemblies coexist within the same dataset.

    AI Generated:

    Multi-reference Classification (EmanProtMultiRefinementNew) — User Manual
        Overview

        The Multi-reference Classification protocol performs iterative
        classification and refinement of subtomograms using multiple structural
        references. Its main purpose is to identify and separate distinct
        structural states within heterogeneous cryo-electron tomography datasets.
        Instead of producing a single consensus average, the protocol organizes
        particles into biologically interpretable groups that each converge
        toward an independent three-dimensional reconstruction.

        In practical cryo-ET studies, biological samples often contain multiple
        conformations, assembly intermediates, ligand-bound states, or damaged
        particles. Multi-reference refinement allows these populations to be
        separated computationally so that each state can be analyzed
        independently. This is especially important for flexible complexes,
        molecular machines, membrane-associated assemblies, and cellular
        structures displaying continuous or discrete variability.

        Inputs and Classification Strategy

        The protocol requires a set of aligned subtomograms as the primary
        input. These particles are expected to already possess approximate
        orientations and shifts from previous subtomogram refinement steps.
        Optionally, one or more reference volumes may be provided to guide the
        classification process.

        The protocol supports several biologically relevant workflows. When no
        references are supplied, the classification can begin from internally
        generated initial averages, allowing exploratory analysis of unknown
        heterogeneity. When a single reference is supplied together with several
        target classes, the protocol generates multiple competing initial models
        derived from the same structure. This strategy is useful when subtle
        conformational differences are expected but no prior structural states
        are available. Alternatively, users may directly provide multiple
        reference volumes representing known or hypothesized structural states.

        From a biological perspective, the quality and diversity of the initial
        references strongly influence classification stability. References that
        are too similar may lead to poor separation, whereas references that are
        excessively different from the particles may bias the refinement toward
        unrealistic solutions.

        Resolution Control and Biological Interpretation

        The classification procedure operates within a user-defined resolution
        range. Restricting the maximum resolution is particularly important
        because classification is not performed under strict gold-standard
        validation conditions. Conservative resolution limits generally improve
        robustness and reduce the risk of overfitting.

        Lower-resolution classification often captures large conformational
        rearrangements and domain motions more reliably than aggressive
        high-resolution refinement. For many biological studies, beginning with
        moderate resolution limits provides more stable and interpretable class
        separation before proceeding to higher-resolution refinement of selected
        classes.

        The minimum resolution parameter allows the exclusion of very low
        frequency information when necessary. This can help reduce the influence
        of global shape differences or tomographic artifacts that are not
        biologically meaningful.

        Symmetry and Symmetry Breaking

        Symmetry can be applied during averaging to improve signal-to-noise
        ratio when the biological assembly is known to possess rotational or
        point-group symmetry. Proper symmetry application can significantly
        enhance reconstruction quality and improve convergence.

        However, biological systems frequently display partial symmetry,
        asymmetry, or localized deviations from ideal symmetry. The protocol
        therefore supports symmetry breaking strategies that allow asymmetric
        features to emerge during classification. This capability is especially
        valuable for studying ligand occupancy, asymmetric maturation pathways,
        or functional transitions within otherwise symmetric complexes.

        Alignment and Local Refinement

        The protocol optionally performs local alignment refinement during
        classification. When enabled, particles undergo constrained orientation
        and translational optimization around their previous alignment
        parameters. This improves class consistency while preserving the
        continuity of the refinement process.

        Local alignment is particularly effective when the input particles are
        already reasonably aligned from previous subtomogram averaging steps.
        Small angular and translational searches typically provide the best
        balance between stability and refinement precision.

        If alignment refinement is disabled, the protocol focuses exclusively on
        classification using the previously assigned orientations. This mode is
        useful when users wish to preserve an existing alignment solution or
        avoid excessive refinement drift during exploratory heterogeneity
        analysis.

        Masking and Focused Classification

        The protocol supports masks for both reference preparation and alignment
        refinement. These masks are among the most biologically important tools
        for controlling classification behavior because they determine which
        structural regions contribute most strongly to particle assignment.

        Reference masks are applied before classification and can emphasize
        conserved structural regions while suppressing solvent or highly mobile
        domains. Alignment masks focus the local refinement procedure onto
        selected structural regions during iterative optimization.

        Focused classification becomes especially important for flexible
        assemblies, membrane proteins with variable peripheral domains, or
        complexes containing mobile subunits. In such systems, carefully chosen
        masks often determine whether subtle conformational states can be
        separated successfully.

        Practical Considerations

        Reliable classification requires that all particles, references, and
        masks share consistent voxel size and box dimensions. Mismatches in
        sampling or dimensions can introduce instability and reduce biological
        interpretability.

        In practical workflows, users commonly begin with a modest number of
        classes and moderate resolution limits. Overclassification may fragment
        the dataset into noisy subclasses that lack biological significance,
        whereas too few classes may obscure meaningful structural variability.

        Classification results should always be interpreted together with visual
        inspection, particle distribution, and biological context. A class with
        very few particles or poorly resolved structural features may represent
        noise, preferred orientations, or alignment artifacts rather than a true
        molecular state.

        Outputs and Biological Meaning

        The protocol produces refined class averages together with updated
        aligned particle sets and class assignments. Each resulting class
        represents a subset of particles sharing similar structural features and
        orientations.

        Biologically, these outputs can reveal conformational landscapes,
        assembly intermediates, ligand-dependent structural rearrangements, or
        previously hidden molecular states. The aligned particles can also serve
        as inputs for additional rounds of focused refinement or downstream
        structural interpretation.

        Final Perspective

        Multi-reference subtomogram classification is one of the most powerful
        strategies for studying structural heterogeneity in cryo-electron
        tomography. Rather than forcing all particles into a single consensus
        reconstruction, the protocol enables the separation and interpretation
        of distinct biological states. Successful application depends on careful
        reference selection, conservative refinement settings, thoughtful mask
        design, and critical biological interpretation of the resulting classes.
    """

    _label = 'Multi-reference classification pppt'
    _possibleOutputs = EmanMultiRefineNewOutputs
    maskRefOutName = 'maskRef'
    maskAlignOutName = 'maskAlign'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

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
        self._addBinThreads(form)

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
                      label='Do alignment?',
                      help='If set to No (default), it will skip the alignment entirely when aligned particles are '
                           'provided. Otherwise a local orientation search will still be performed.')
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
                      label='Maximum shift (px)',
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

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.createEmanPrjPostExtractionStep)
        self._insertFunctionStep(self.refineMultiStep)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    def _initialize(self):
        self.inParticles = self.getAttrib(IN_SUBTOMOS)
        self.inSamplingRate = self.inParticles.getSamplingRate()
        self.numClasses = self.getNoRefs()
        self.refFiles = self.getReferenceFiles()

    def refineMultiStep(self):
        # In case of continuing from this step, the previous results dir will be removed to avoid EMAN creating one
        # for each execution (one for each continue)
        multiRefineDir = self.getMultiRefineDir()
        if exists(multiRefineDir):
            shutil.rmtree(multiRefineDir)
        self.buildEmanSets()
        program = Plugin.getProgram("e2spt_refinemulti_new.py")
        self.runJob(program, self._genRefineMultiCmd(), cwd=self._getExtraPath())

    def convertOutputStep(self):
        for classId in range(self.numClasses):
            inFile = join(SPTCLS_00_DIR, self.getOutputThreed(classId, 'hdf', onlyBaseName=True))
            outFile = inFile.replace('.hdf', '.mrc')
            self.convertBetweenHdfAndMrc(inFile, outFile, extraArgs=f'--apix {self.inSamplingRate:.3f}')
            fixVolume(self._getExtraPath(outFile))

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

        # Define outputs and relations
        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: outParticles,
                               self._possibleOutputs.classes.name: classes3d})
        self._defineSourceRelation(inParticlesPointer, outParticles)
        self._defineSourceRelation(inParticlesPointer, classes3d)

    # --------------------------- UTILS functions ------------------------------
    def _genRefineMultiCmd(self):
        refMask = self.maskRef.get()
        alignMask = self.maskAlign.get()
        breakSym = self.breakSym.get()
        new2dAlignFile = self.getNewAliFile(is3d=False, outPath=SPT_00_DIR)
        nClasses = self.nClasses.get()
        args = [
            ' '.join(self.refFiles) if self.refFiles else '',
            f'--ptcls {self.getNewAliFile(outPath=SPT_00_DIR)}',
            f'--niter {self.nIter.get()}',
            f'--sym {self.symmetry.get()}',
            f'--parallel thread:{self.binThreads.get()}',
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
            if isinstance(refs, SubTomogram):
                refList = [abspath(refs.getFileName())]
            else:
                refList = [abspath(vol.getFileName()) for vol in refs]
        return refList

    def getNoRefs(self):
        nClasses = self.nClasses.get()
        refs = self.getAttrib(REF_VOL)
        if nClasses > 0:
            return nClasses
        else:
            return 1 if isinstance(refs, SubTomogram) else len(refs)

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
        tol = 1e-03
        ih = ImageHandler()
        nClasses = self.nClasses.get()
        refs = self.getAttrib(REF_VOL)
        inParticles = self.getAttrib(IN_SUBTOMOS)
        inParticlesDims = inParticles.getBoxSize()
        inParticlesSRate = inParticles.getSamplingRate()
        if nClasses <= 0 and not refs:
            errorMsg.append('At least a number of classes or a set of references is required.')
        if refs:
            # Dimensions
            ref = refs if type(refs) is AverageSubTomogram else refs.getFirstItem()
            x, y, z = ref.getDim()
            refsDims = (x, y, z)
            if refsDims != inParticlesDims:
                errorMsg.append(f'The dimensions of the referece/s {refsDims} px and the particles '
                                f'{inParticlesDims} px must be the same')
            # Sampling rate
            refsSRate = refs.getSamplingRate()
            if abs(inParticlesSRate - refsSRate) >= tol:
                errorMsg.append(
                    f'The sampling rate of the input particles [{inParticlesSRate} Å/pix] and the reference volume/s '
                    f'[{refsSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')
        maskRef = self.maskRef.get()
        if maskRef:
            # Dimensions
            x, y, z, _ = ih.getDimensions(maskRef.getFileName())
            maskDims = (x, y, z)
            maskSRate = maskRef.getSamplingRate()
            if inParticlesDims != maskDims:
                errorMsg.append(f'The dimensions of the introduced mask {maskDims} px and the particles '
                                f'{inParticlesDims} px must be the same')
            if abs(inParticlesSRate - maskSRate) >= tol:
                errorMsg.append(
                    f'The sampling rate of the input particles [{inParticlesSRate} Å/pix] and the mask '
                    f'[{maskSRate} Å/pix] are not equal within the specified tolerance [{tol} Å/pix].')
        return errorMsg
