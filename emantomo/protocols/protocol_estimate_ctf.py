# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *              Scipion Team (scipion@cnb.csic.es) [1]
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

from enum import Enum
from os.path import exists, join
from emantomo import Plugin
from emantomo.constants import TS_DIR
from emantomo.objects import EmanMetaData
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TS
from emantomo.utils import getPresentTsIdsInSet, genJsonFileName
from pyworkflow.protocol import PointerParam, FloatParam, IntParam, BooleanParam
from pyworkflow.utils import Message
from tomo.objects import SetOfCTFTomoSeries, CTFTomoSeries, CTFTomo
from emantomo.convert import loadJson, ts2Json


class EstimateCtfOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class EmanProtEstimateCTF(ProtEmantomoBase):
    """
    Estimates the contrast transfer function parameters for tilt-series images
    in cryo-electron tomography experiments.

    AI Generated:

    CTF Estimation (EmanProtEstimateCTF) — User Manual
        Overview

        The CTF Estimation protocol determines the contrast transfer function
        parameters associated with tilt-series images in cryo-electron
        tomography workflows. Its main purpose is to estimate defocus and,
        optionally, phase shift values for each tilt image so that downstream
        tomographic reconstruction and subtomogram analysis can properly account
        for microscope-induced image distortions.

        In cryo-ET experiments, accurate CTF estimation is essential because the
        electron microscope modifies image contrast in a way that strongly
        affects structural interpretation and achievable resolution. Correctly
        estimated CTF parameters are therefore a prerequisite for high-quality
        tomogram reconstruction, subtomogram averaging, and quantitative
        structural analysis.

        Biological Importance of CTF Estimation

        Biological macromolecules embedded in vitrified samples are typically
        imaged under defocus conditions to enhance contrast. Although this
        improves particle visibility, it also introduces oscillatory distortions
        in the recorded signal. The CTF estimation process attempts to measure
        these distortions so they can later be corrected computationally.

        In practical biological workflows, accurate CTF estimation improves the
        interpretability of cellular structures, membrane complexes, viral
        assemblies, ribosomes, and other macromolecular systems observed in situ.
        Poor CTF estimation may lead to reduced resolution, reconstruction
        artifacts, or incorrect interpretation of structural features.

        Inputs and Experimental Requirements

        The protocol requires a set of aligned tilt series as input. These tilt
        series should already contain the acquisition geometry and metadata
        associated with the microscope settings used during data collection.

        Because the estimation relies on microscope acquisition parameters such
        as acceleration voltage and spherical aberration, the tilt-series
        metadata must accurately reflect the imaging conditions. Inconsistent or
        incomplete acquisition information can compromise the reliability of the
        resulting CTF models.

        Defocus Search Range

        One of the central parameters of the protocol is the defocus search
        range. Users define the minimum and maximum expected defocus values
        together with the search increment. The protocol then evaluates possible
        solutions within this interval.

        Choosing an appropriate defocus range is biologically important because
        unrealistic values can lead to unstable estimations or convergence to
        incorrect solutions. Narrow search ranges generally improve speed and
        robustness when approximate acquisition conditions are already known,
        whereas wider ranges may be required for exploratory datasets or poorly
        characterized acquisitions.

        In most cryo-ET experiments, moderate underfocus values are commonly
        used to balance image contrast and preservation of high-resolution
        information.

        Phase Shift Estimation

        The protocol optionally supports phase shift estimation, which is
        particularly relevant for datasets acquired using phase plates. In these
        experiments, the microscope intentionally introduces an additional phase
        modification to improve image contrast at low spatial frequencies.

        When phase shift estimation is enabled, the protocol searches for
        physically meaningful phase values together with the defocus estimation.
        This can substantially improve reconstruction quality for phase-plate
        datasets but may increase computational complexity.

        For standard cryo-ET acquisitions without phase plates, phase shift
        estimation is typically unnecessary and may be safely disabled.

        Tile-Based Frequency Analysis

        CTF estimation is performed using tiled regions extracted from the tilt
        images. These local regions are analyzed in Fourier space to identify
        characteristic oscillatory patterns associated with the microscope CTF.

        The tile size parameter determines the balance between spatial locality
        and frequency precision. Larger tiles generally provide more reliable
        frequency information but may include structural heterogeneity or sample
        contamination. Smaller tiles improve local specificity but may reduce
        estimation stability in noisy datasets.

        In biological practice, moderate tile sizes often provide the best
        compromise for cellular tomography data.

        Reference Tilts and Sampling Strategy

        The protocol uses a subset of central tilt images as references to
        estimate the general defocus behavior across the tilt series. Central
        tilts are usually preferred because they suffer less from projection
        distortions, increased sample thickness, and signal degradation compared
        to highly tilted images.

        Additional sampling parameters control the spatial distribution of tiles
        across the image plane. Broader sampling improves robustness against
        local contamination, ice gradients, or uneven illumination, although it
        also increases computational cost.

        Biological datasets containing thick cellular material, fiducial beads,
        or strongly heterogeneous regions may benefit from more extensive
        sampling strategies.

        Outputs and Interpretation

        The protocol produces a set of CTF models associated with the input tilt
        series. Each tilt image receives estimated defocus values and, when
        applicable, phase shift information.

        These outputs can subsequently be used in tomographic reconstruction,
        CTF correction, subtomogram extraction, and averaging workflows. The
        estimated parameters are especially important for achieving improved
        resolution and reducing systematic imaging artifacts.

        Biologically, reliable CTF estimation contributes directly to the
        interpretability of macromolecular structures and to the reproducibility
        of cryo-ET analyses.

        Practical Recommendations

        In routine cryo-electron tomography workflows, users should begin with
        defocus ranges consistent with the microscope acquisition settings and
        verify that the estimated values vary smoothly across the tilt series.
        Extremely irregular estimates may indicate poor signal quality,
        contamination, incorrect metadata, or unsuitable search ranges.

        For conventional cryo-ET datasets, disabling phase shift estimation is
        usually appropriate unless phase plates were explicitly used during data
        acquisition.

        Central tilt images generally provide the most reliable signal for CTF
        estimation, while highly tilted projections may contain stronger noise
        and sample-thickness artifacts. Visual inspection of the resulting CTF
        behavior across tilts is therefore strongly recommended before
        proceeding to reconstruction or subtomogram averaging.

        Final Perspective

        The CTF Estimation protocol provides a critical foundation for accurate
        cryo-electron tomography processing by characterizing the optical effects
        introduced during image acquisition. Through systematic estimation of
        defocus and optional phase shifts, the protocol enables downstream
        tomographic reconstruction and subtomogram analysis to achieve improved
        structural fidelity and biological interpretability.

        In modern cryo-ET workflows, careful CTF estimation is not merely a
        technical preprocessing step but a key determinant of the final
        structural quality and reliability of biological conclusions.
    """
    _label = 'ctf estimation'
    _possibleOutputs = EstimateCtfOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series",
                      important=True)
        lineDefocus = form.addLine('Defocus range (μm)',
                                   help="Search range of defocus (start, end, step). Note that "
                                        "they must be introduced in microns.")

        lineDefocus.addParam('minDefocus', FloatParam, default=2.0, label='min')
        lineDefocus.addParam('maxDefocus', FloatParam, default=7.0, label='max')
        lineDefocus.addParam('stepDefocus', FloatParam, default=0.02, label='Step')
        form.addParam('doPhaseShiftSearch', BooleanParam,
                      label='Do phase shift search?',
                      default=False)
        linePhaseShift = form.addLine('Phase shift range (deg.)',
                                      condition='doPhaseShiftSearch',
                                      help="Search range of the phase shift (start, end, step). To avoid"
                                           "the phase shift search use min 0.0, max 1.0, and step 1.0.")
        linePhaseShift.addParam('minPhaseShift', FloatParam, default=0, label='min', condition='doPhaseShiftSearch')
        linePhaseShift.addParam('maxPhaseShift', FloatParam, default=1, label='max', condition='doPhaseShiftSearch')
        linePhaseShift.addParam('stepPhaseShift', FloatParam, default=1, label='Step', condition='doPhaseShiftSearch')

        form.addParam('tilesize', IntParam,
                      label='Size of tile to calculate FFT',
                      default=256)
        form.addParam('nref', IntParam,
                      label='Number of references',
                      default=15,
                      help='Using N tilt images near the center tilt to estimate the range of defocus for all images.')
        form.addParam('stepx', IntParam,
                      label='Step in X direction',
                      default=20,
                      help='Number of tiles to generate on x-axis (different defocus)')
        form.addParam('stepy', IntParam,
                      label='Step in Y direction',
                      default=40,
                      help='Number of tiles to generate on y-axis (same defocus)')

        self._addBinThreads(form)
        form.addParallelSection(threads=1, mpi=0)

        # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertTsStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
            self._insertFunctionStep(self.estimateCtfStep, mdObj)
        self._insertFunctionStep(self.createOutputStep, mdObjDict)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTsSet = self.getAttrib(IN_TS)
        self.createInitEmanPrjDirs()
        # Manage the TS
        presentTsIds = set(getPresentTsIdsInSet(inTsSet))
        tsIdsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inTsSet if ts.getTsId() in presentTsIds}
        # Get the required acquisition data
        self.sphAb = inTsSet.getAcquisition().getSphericalAberration()
        self.voltage = inTsSet.getAcquisition().getVoltage()
        # Split all the data into subunits (EmanMetaData objects) referred to the same tsId
        mdObjDict = {}
        for tomoId, ts in tsIdsDict.items():
            mdObjDict[tomoId] = EmanMetaData(tsId=tomoId,
                                             ts=ts,
                                             tsHdfName=join(TS_DIR, f'{tomoId}.hdf'),
                                             jsonFile=genJsonFileName(self.getInfoDir(), tomoId))
        return mdObjDict

    @staticmethod
    def writeData2JsonFileStep(mdObj):
        mode = 'a' if exists(mdObj.jsonFile) else 'w'
        ts2Json(mdObj, mode=mode)

    def estimateCtfStep(self, mdObj):
        program = Plugin.getProgram('e2spt_tomoctf.py')
        self.runJob(program, self._genCtfEstimationArgs(mdObj), cwd=self._getExtraPath())

    def createOutputStep(self, mdObjDict):
        inTsSet = self.getAttrib(IN_TS)
        outCtfSet = SetOfCTFTomoSeries.create(self._getPath(), template='CTFmodels%s.sqlite')
        outCtfSet.setSetOfTiltSeries(inTsSet)

        for tsId, mdObj in mdObjDict.items():
            ts = mdObj.ts
            newCTFTomoSeries = CTFTomoSeries()
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            newCTFTomoSeries.setObjId(ts.getObjId())
            newCTFTomoSeries.setTsId(tsId)
            outCtfSet.append(newCTFTomoSeries)
            jsonDict = loadJson(mdObj.jsonFile)
            defocus = jsonDict['defocus']
            phase_shift = jsonDict['phase']
            for idx, tiltImage in enumerate(ts.iterItems()):
                defocusU = defocusV = 10000.0 * defocus[idx]
                newCTFTomo = CTFTomo()
                newCTFTomo.setIndex(idx + 1)
                newCTFTomo.setAcquisitionOrder(tiltImage.getAcquisitionOrder())
                if phase_shift[idx] != 0:
                    newCTFTomo.setPhaseShift(phase_shift[idx])
                newCTFTomo.setDefocusU(defocusU)
                newCTFTomo.setDefocusV(defocusV)
                newCTFTomo.setDefocusAngle(0)
                newCTFTomoSeries.append(newCTFTomo)

        self._defineOutputs(**{self._possibleOutputs.CTFs.name: outCtfSet})
        self._defineSourceRelation(getattr(self, IN_TS), outCtfSet)

    # --------------------------- UTILS functions -----------------------------
    def _genCtfEstimationArgs(self, mdObj):
        args = '%s ' % mdObj.tsHdfName
        args += '--dfrange %s ' % ','.join([str(self.minDefocus.get()),
                                            str(self.maxDefocus.get()),
                                            str(self.stepDefocus.get())])
        args += '--psrange %s ' % ','.join([str(self.minPhaseShift.get()),
                                            str(self.maxPhaseShift.get()),
                                            str(self.stepPhaseShift.get())])
        args += '--tilesize %i ' % self.tilesize.get()
        args += '--voltage %i ' % self.voltage
        args += '--cs %.2f ' % self.sphAb
        args += '--nref %i ' % self.nref.get()
        args += '--stepx %i ' % self.stepx.get()
        args += '--stepy %i ' % self.stepy.get()
        args += '--threads %i ' % self.binThreads.get()
        args += '--verbose 9 '
        return args


