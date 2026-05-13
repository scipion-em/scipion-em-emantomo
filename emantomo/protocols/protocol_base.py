# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]             
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import logging
import re
from os.path import join, abspath, basename

import numpy as np

from emantomo import Plugin
from emantomo.constants import INFO_DIR, TOMOGRAMS_DIR, TS_DIR, SETS_DIR, PARTICLES_DIR, PARTICLES_3D_DIR, \
    REFERENCE_NAME, TOMOBOX, SPT_00_DIR, THREED, ALI3D_BASENAME, ALI2D_BASENAME, FSC_MASKED_BNAME, FSC_UNMASKED_BNAME, \
    FSC_MASKED_TIGHT_BNAME, LST_LINE, PART3D_ID, INIT_MODEL_DIR, SPT_CLS_00_FIR
from emantomo.convert import emanFSCsToScipion
from emantomo.convert.lstConvert import EmanLstReader, EmanLstWriter
from emantomo.objects import EmanParticle, EmanSetOfParticles, EmanMetaData
from pwem.objects import SetOfFSCs
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer, String
from pyworkflow.protocol import IntParam
from pyworkflow.utils import makePath, createLink, redStr, cyanStr
from tomo.protocols import ProtTomoBase


logger = logging.getLogger(__name__)

IN_TS = 'inputTS'
IN_COORDS = 'inputCoordinates'
IN_CTF = 'inputCTF'
IN_SUBTOMOS = 'inputSubtomos'
IN_BOXSIZE = 'boxSize'
IN_TOMOS = 'inputTomograms'
REF_VOL = 'refVol'


class ProtEmantomoBase(EMProtocol, ProtTomoBase):
    """
    Base protocol for EMAN tomography workflows. It provides the common
    infrastructure required to prepare, organize, convert, and manage
    tomographic datasets and subtomogram data within EMAN-based cryo-ET
    processing pipelines.

    AI Generated:

    EMAN Tomography Base Protocol (ProtEmantomoBase) — User Manual
        Overview

        The EMAN Tomography Base protocol acts as the foundational layer for
        tomography workflows that rely on EMAN software within Scipion. Its
        main purpose is to standardize how tomographic datasets, subtomograms,
        references, metadata, and refinement resources are prepared and managed
        before higher-level processing tasks are performed.

        In cryo-electron tomography workflows, many downstream operations such
        as subtomogram alignment, averaging, classification, refinement, or CTF
        estimation depend on a consistent project structure and carefully
        organized metadata. This protocol provides that common environment so
        that derived tomography protocols can operate reliably and consistently
        across different biological datasets.

        Biological Context

        Cryo-electron tomography experiments frequently involve large and
        heterogeneous datasets containing tilt series, reconstructed tomograms,
        extracted subtomograms, and multiple refinement iterations. Managing
        these datasets correctly is essential because inconsistencies in file
        organization, sampling rates, coordinate systems, or metadata handling
        can propagate errors throughout the entire structural analysis workflow.

        The protocol is designed to simplify this process by automatically
        organizing the information required by EMAN tomography procedures. This
        allows biological users to focus more on the interpretation of
        structures and less on technical preparation steps.

        Data Preparation and Conversion

        One of the central responsibilities of the protocol is preparing input
        data in formats compatible with EMAN tomography applications. Cryo-ET
        datasets are commonly generated in different formats depending on the
        acquisition system and reconstruction software. The protocol harmonizes
        these datasets so they can be processed uniformly.

        During this preparation stage, tilt series, tomograms, particle stacks,
        and reference maps may be converted or linked into the working project
        structure. Sampling rate consistency is especially important because
        downstream alignment and refinement procedures assume that voxel sizes
        are physically meaningful and internally coherent.

        For biological interpretation, preserving the correct sampling rate is
        critical. Incorrect voxel scaling may lead to inaccurate structural
        dimensions, unreliable fitting, or incorrect interpretation of molecular
        conformations.

        Project Organization

        The protocol establishes a standardized EMAN tomography project
        hierarchy containing directories for tilt series, tomograms, metadata,
        particle stacks, and refinement outputs. This organization ensures that
        all subsequent tomography operations can locate the required resources
        consistently.

        In practical cryo-ET workflows, this structured organization becomes
        particularly important when working with large collections of tomograms
        or iterative subtomogram refinement projects. Proper organization also
        facilitates reproducibility, collaborative work, and long-term project
        maintenance.

        Subtomogram and Particle Management

        The protocol provides utilities for managing subtomogram particles and
        their associated alignment information. This includes maintaining
        relationships between particle stacks, geometric orientations, and
        refinement metadata throughout iterative processing stages.

        From a biological perspective, preserving particle orientation and
        alignment information is essential because these transformations encode
        the spatial relationship between molecular complexes and the original
        tomographic volume. Reliable propagation of this information enables
        meaningful averaging and structural interpretation.

        The protocol also supports workflows in which particles originate from
        oriented picking procedures or previous refinement iterations. This is
        especially important in advanced subtomogram averaging studies where
        orientation priors contribute significantly to refinement stability.

        Reference Volume Handling

        Many tomography workflows rely on one or more reference volumes for
        alignment or refinement. The protocol supports preparation and
        organization of these references, including half-map handling for
        independent validation workflows.

        In biological practice, the choice of reference can strongly influence
        refinement behavior. High-quality references generally improve alignment
        stability, whereas biased or structurally incorrect references may lead
        to reference-driven artifacts. Proper handling of independent half maps
        is particularly important for reliable resolution estimation and
        validation.

        Refinement and Iterative Processing

        The protocol supports iterative refinement workflows by maintaining
        refinement directories, alignment files, and averaged reconstructions
        generated during successive processing rounds. This organization enables
        users to monitor refinement progression and evaluate reconstruction
        quality over time.

        In practical cryo-ET studies, iterative refinement is often necessary
        to improve structural resolution and recover weak biological signals.
        The ability to track intermediate iterations helps users identify
        convergence problems, overfitting, or structural heterogeneity.

        FSC and Validation Support

        The protocol includes support for organizing and interpreting Fourier
        Shell Correlation data generated during refinement workflows. FSC
        analysis is central to cryo-EM and cryo-ET validation because it
        provides an estimate of reconstruction reproducibility and effective
        resolution.

        Biological users should interpret FSC curves carefully. While higher
        resolution values are often desirable, they are only meaningful if the
        reconstruction remains biologically plausible and free from overfitting.
        Independent half-map validation remains one of the most important
        safeguards in modern cryo-EM analysis.

        Parallel Execution and Resource Management

        The protocol allows users to define the number of EMAN processing
        threads used during execution. This enables efficient processing of
        large tomography datasets while adapting resource usage to the available
        computational infrastructure.

        In practical environments such as institutional cryo-EM facilities or
        shared computing clusters, balancing processing speed and memory usage
        becomes essential. Large tomograms and particle stacks may require
        substantial memory resources, especially during parallel processing.

        Practical Recommendations

        In routine cryo-ET workflows, users should ensure that sampling rates,
        box sizes, and coordinate systems remain consistent before starting
        downstream processing. Careful organization of metadata and reference
        files significantly reduces the risk of refinement instability or
        reconstruction artifacts.

        When working with large datasets, it is advisable to monitor storage
        usage and computational resources closely because tomography workflows
        can rapidly generate large intermediate files and multiple refinement
        iterations.

        For subtomogram averaging studies, maintaining accurate particle
        orientation information and validation procedures is essential for
        obtaining biologically meaningful structures.

        Final Perspective

        The EMAN Tomography Base protocol serves as the operational backbone of
        EMAN-based cryo-electron tomography workflows within Scipion. Although
        it primarily manages infrastructure and data organization, its role is
        biologically important because reliable structural interpretation
        depends fundamentally on consistent metadata handling, accurate sampling
        relationships, and reproducible refinement organization.

        For most cryo-ET users, careful preparation and management of the
        tomography environment are as important as the reconstruction algorithms
        themselves. Well-organized workflows improve reproducibility, reduce
        processing errors, and support more reliable biological conclusions.
    """

    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inParticles = None
        self.inSamplingRate = -1.0
        self.scaleFactor = 1.0
        self.voltage = 300.0
        self.sphAb = 2.7
        self.failedItems = []

    @staticmethod
    def _addBinThreads(form):
        form.addParam('binThreads', IntParam,
                      label='Emantomo threads',
                      default=6,
                      important=True,
                      help='Number of threads used by EMAN each time it is called in the protocol execution. For '
                           'example, if 2 Scipion threads and 3 Emantomo threads are set, the tomograms will be '
                           'processed in groups of 2 at the same time with a call of EMAN with 3 threads each, so '
                           '6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of tomograms specified by Scipion threads.')

    # --------------------------- STEPS functions -----------------------------
    def convertTsStep(self, mdObj: EmanMetaData):
        # The converted TS must be unbinned, because EMAN will read the sampling rate from its header. This is why
        # the TS associated to the CTF is the one considered first. Later, when generating the json, the TS alignment
        # parameters are read from the introduced TS and the shifts are scaled to at the unbinned scale
        tsId = mdObj.tsId
        try:
            logger.info(cyanStr(f'Converting TS {tsId} into HDF...'))
            inTsFName = mdObj.ts.getFirstItem().getFileName()
            sRate = mdObj.ts.getSamplingRate()
            self.convertOrLink(inTsFName, mdObj.tsId, TS_DIR, sRate)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))

    def createEmanPrjPostExtractionStep(self):
        inSubtomos = getattr(self, IN_SUBTOMOS).get()
        # Create project dir structure
        self.createInitEmanPrjDirs()
        infoDir = self.getInfoDir()
        tsDir = self.getTsDir()
        # tomoDir = self.getTomogramsDir()
        stack2dDir = self.getStack2dDir()
        stack3dDir = self.getStack3dDir()
        makePath(self.getSetsDir(), stack2dDir, stack3dDir, self.getRefineDir())
        # Get the unique values of the files to be linked
        dataDict = inSubtomos.getUniqueValues([EmanParticle.INFO_JSON,
                                               EmanParticle.TS_HDF,
                                               # EmanParticle.TOMO_HDF,
                                               EmanParticle.STACK_2D_HDF,
                                               EmanParticle.STACK_3D_HDF])
        # Link the files
        for infoJson, tsFile, stack2d, stack3d in zip(dataDict[EmanParticle.INFO_JSON],
                                                      dataDict[EmanParticle.TS_HDF],
                                                      # dataDict[EmanParticle.TOMO_HDF],
                                                      dataDict[EmanParticle.STACK_2D_HDF],
                                                      dataDict[EmanParticle.STACK_3D_HDF]):
            createLink(infoJson, join(infoDir, basename(infoJson)))
            createLink(tsFile, join(tsDir, basename(tsFile)))
            # createLink(tomoFile, join(tomoDir, basename(tomoFile)))
            createLink(stack2d, join(stack2dDir, basename(stack2d)))
            createLink(stack3d, join(stack3dDir, basename(stack3d)))

    def convertRefVolStep(self):
        inRef = self.getRefVol()
        if inRef:
            sRate = inRef.getSamplingRate()
            self.convertOrLink(inRef.getFileName(), REFERENCE_NAME, '', sRate)
            if inRef.hasHalfMaps():
                halves = inRef.getHalfMaps(asList=True)
                self.convertOrLink(halves[0], f'{REFERENCE_NAME}_even', '', sRate)
                self.convertOrLink(halves[1], f'{REFERENCE_NAME}_odd', '', sRate)

    # --------------------------- UTILS functions ----------------------------------
    def buildEmanSets(self, outAliPath=SPT_00_DIR):
        # LST with the particles
        EmanLstWriter.writeSimpleLst(self.inParticles, self.getLstEmanRelPath())
        # LST with the 2d/3d particles and the corresponding alignments
        align3dFile = getattr(self.inParticles, EmanSetOfParticles.ALI_3D, String()).get()
        align2dFile = getattr(self.inParticles, EmanSetOfParticles.ALI_2D, String()).get()
        transformMatrix = self.inParticles.getFirstItem().getTransform(convention=None).getMatrix()
        particlesFromOrientedPicking = np.any(abs(transformMatrix - np.eye(4)) > 1e-3)
        if align3dFile or particlesFromOrientedPicking:
            new3dAlignFile = self._getExtraPath(self.getNewAliFile(outPath=outAliPath))
            EmanLstWriter.writeAlign3dLst(self.inParticles, new3dAlignFile)
        if align2dFile:
            new2dAlignFile = self._getExtraPath(self.getNewAliFile(outPath=outAliPath, is3d=False))
            dataDictList = EmanLstReader.read2dParticles(align2dFile)
            self.write2dLst(new2dAlignFile, dataDictList)

    def getObjByName(self, name):
        """Return an object, from a protocol, named 'name' instead of a pointer."""
        obj = getattr(self, name, None)
        if obj and type(obj) == Pointer:
            return obj.get()
        else:
            return obj

    def getBoxSize(self):
        """If the user provides a box size value. If not (expected behaviour), read it from the set of subtomograms
        introduced."""
        boxSizeFromForm = self.getAttrib(IN_BOXSIZE)
        return boxSizeFromForm if boxSizeFromForm else self.getAttrib(IN_SUBTOMOS).getCoordinates3D().getBoxSize()

    def getRefVol(self):
        return self.getAttrib(REF_VOL)

    def createInitEmanPrjDirs(self):
        """Create in the current protocol's extra path the initial directory structure of an EMAN tomo project: info,
        tiltseries and tomograms dirs"""
        makePath(self.getInfoDir(), self.getTomogramsDir(), self.getTsDir())

    def getInfoDir(self):
        return self._getExtraPath(INFO_DIR)

    def getTomogramsDir(self):
        return self._getExtraPath(TOMOGRAMS_DIR)

    def getTsDir(self):
        return self._getExtraPath(TS_DIR)

    def getStack2dDir(self):
        return self._getExtraPath(PARTICLES_DIR)

    def getStack3dDir(self):
        return self._getExtraPath(PARTICLES_3D_DIR)

    def getSetsDir(self):
        return self._getExtraPath(SETS_DIR)

    def getInitModelDir(self):
        return self._getExtraPath(INIT_MODEL_DIR)

    def getRefineDir(self):
        return self._getExtraPath(SPT_00_DIR)

    def getMultiRefineDir(self):
        return self._getExtraPath(SPT_CLS_00_FIR)

    @staticmethod
    def getNewAliFile(is3d=True, outPath=None):
        aliFile = 'particle_info_3d.lst' if is3d else 'particle_info_2d.lst'  # Names are hardcoded in some parts of EMAN's native code
        return join(outPath, aliFile) if outPath else aliFile

    def getAttrib(self, attribName, getPointer=False):
        attribPointer = getattr(self, attribName)
        return attribPointer if getPointer else attribPointer.get()

    def convertOrLink(self, inFile, tsId, outDir, sRate):
        """Fill the simulated EMAN project directories with the expected data at this point of the pipeline.
        Also convert the precedent tomograms into HDF files if they are not. The converted filename will be the tsId,
        avoiding added suffixes this way"""
        hdf = '.hdf'
        outFile = join(outDir, tsId + hdf)
        if inFile.endswith(hdf):
            createLink(inFile, self._getExtraPath(outFile))
        else:
            self.convertBetweenHdfAndMrc(abspath(inFile), outFile, extraArgs=f'--apix {sRate:.3f}')
        return outFile

    def convertBetweenHdfAndMrc(self, inFile, outFile, extraArgs=''):
        program = Plugin.getProgram("e2proc3d.py")
        args = '%s %s ' % (inFile, outFile)
        self.runJob(program, args + extraArgs, cwd=self._getExtraPath())

    def _getLstFile(self):
        lstFile = glob.glob(join(self.getSetsDir(), TOMOBOX + '*.lst'))[0]
        return join(SETS_DIR, basename(lstFile))

    def getLstEmanRelPath(self):
        return join(self.getSetsDir(), f'{TOMOBOX}.lst')

    def getRefinedAverageFn(self, iterNumber, ext='hdf', half=None):
        # Example, for 9 iterations, the resulting file would be called threed_09.hdf
        pattern = self._getExtraPath(SPT_00_DIR, THREED + f'_{iterNumber:02d}')
        return pattern + f'_{half}.{ext}' if half else pattern + f'.{ext}'

    def getRefineEvenFn(self, iterNumber):
        # Example, for 9 iterations, the resulting file would be called threed_09_even.hdf
        return self.getRefinedAverageFn(iterNumber, half='even')

    def getRefineOddFn(self, iterNumber):
        # Example, for 9 iterations, the resulting file would be called threed_09_odd.hdf
        return self.getRefinedAverageFn(iterNumber, half='odd')

    def getAli3dFile(self, iterNum):
        return self._getExtraPath(SPT_00_DIR, f'{ALI3D_BASENAME}{iterNum:02d}.lst')

    def getAli2dFile(self, iterNum):
        return self._getExtraPath(SPT_00_DIR, f'{ALI2D_BASENAME}{iterNum:02d}.lst')

    def genFscs(self, iterNum):
        sptPath = self._getExtraPath(SPT_00_DIR)
        fscs = SetOfFSCs.create(self._getPath(), template='fscs%s.sqlite')
        fscMasked = join(sptPath, f'{FSC_MASKED_BNAME}{iterNum:02d}.txt')
        fscUnmasked = join(sptPath, f'{FSC_UNMASKED_BNAME}{iterNum:02d}.txt')
        fscTight = join(sptPath, f'{FSC_MASKED_TIGHT_BNAME}{iterNum:02d}.txt')
        emanFSCsToScipion(fscs, fscMasked, fscUnmasked, fscTight)
        return fscs

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob.glob(join(self.getRefineDir(), '*'))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]

    def write2dLst(self, new2dAlignFile, dataDictList):
        presentAbsIndices = [particle.getAbsIndex() for particle in self.inParticles]
        lines = []
        for dataDict in dataDictList:
            line = dataDict[LST_LINE]
            particle3dInd = dataDict[PART3D_ID]
            if particle3dInd in presentAbsIndices:
                lines.append(line.replace('\n', ''))
        EmanLstWriter.lines2LstFile(lines, new2dAlignFile)
