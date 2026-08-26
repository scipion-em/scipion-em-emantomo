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
from pyworkflow.protocol import PointerParam
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_SUBTOMOS, REF_VOL


class EmanProtRefineNewBase(ProtEmantomoBase):
    """
    Provides a common foundation for subtomogram refinement workflows within EMAN-based cryo-electron tomography processing pipelines. The protocol establishes the shared input structure required for iterative refinement procedures, allowing users to combine subtomogram particle datasets with optional reference volumes in order to improve structural consistency, orientation accuracy, and reconstruction quality during downstream refinement stages.

    AI Generated:

    Refinement Base Protocol (EmanProtRefineNewBase) — User Manual
        Overview

        The Refinement Base protocol serves as a shared framework for tomography refinement workflows based on
        EMAN subtomogram processing strategies. Its purpose is to standardize the biological and computational
        inputs required for iterative refinement procedures while providing a consistent environment for
        advanced subtomogram alignment and reconstruction tasks.

        In cryo-electron tomography, refinement is a critical stage where initially extracted subtomograms are
        progressively aligned and averaged in order to improve signal quality and recover higher-resolution
        structural information. This protocol defines the common organization used by multiple refinement
        workflows and ensures that subtomogram datasets and reference structures are handled consistently.

        Biological Context and Typical Applications

        Subtomogram refinement is commonly applied when studying macromolecular complexes directly within
        cellular or native environments. Typical biological applications include membrane proteins, ribosomes,
        viral assemblies, cytoskeletal structures, and flexible molecular machines observed in situ.

        In these workflows, particles are first identified and extracted from tomograms, often producing noisy
        and heterogeneous datasets. Refinement protocols then iteratively improve particle orientations and
        positions relative to a common reference, gradually increasing structural consistency across the
        dataset.

        The refinement process becomes especially important when aiming to distinguish conformational states,
        improve interpretability of density maps, or prepare structures for downstream biological analysis and
        modeling.

        Inputs and General Workflow

        The protocol is centered around two primary biological inputs: a set of subtomogram particles and an
        optional reference volume. The particle dataset contains the extracted three-dimensional observations
        of the molecular complexes under study, while the reference volume provides an initial structural
        target that guides alignment and refinement.

        The reference volume is optional because some refinement strategies can begin from low-resolution or
        internally generated initial models. However, when a biologically meaningful reference is available,
        it often improves convergence stability and accelerates refinement.

        In practical workflows, the reference should ideally resemble the expected biological structure and
        share compatible voxel size and box dimensions with the subtomograms. Excessively inaccurate references
        may bias refinement or reduce alignment reliability.

        Refinement Philosophy

        Iterative refinement aims to improve structural coherence across all particles while preserving
        biologically meaningful variability. During refinement, particles are repeatedly aligned against the
        evolving reference, producing progressively improved reconstructions.

        From a biological perspective, refinement is not merely a computational optimization procedure. The
        process directly influences the interpretability of structural features such as secondary structure
        elements, membrane organization, ligand binding states, or conformational variability.

        Careful refinement therefore becomes essential for extracting reliable biological conclusions from
        noisy tomographic data.

        Computational Considerations

        Tomography refinement workflows are often computationally demanding because they involve repeated
        three-dimensional alignment operations across large particle datasets. The protocol supports execution
        environments where computational resources such as processing threads and binning strategies can be
        adjusted to balance performance and accuracy.

        Lower-resolution exploratory refinements are often useful during the early stages of analysis because
        they allow rapid validation of alignment behavior. Higher-resolution refinements can later be applied
        once particle quality and structural homogeneity have been confirmed.

        Outputs and Biological Interpretation

        Refinement workflows built upon this protocol typically generate improved particle alignments,
        transformation parameters, and progressively enhanced average reconstructions. These outputs serve as
        the structural basis for downstream interpretation, classification, and visualization.

        Biologically, improved refinement quality often leads to clearer structural features, reduced noise,
        and better discrimination between molecular states. However, users should remain aware that excessive
        refinement against inaccurate references or heterogeneous datasets may artificially suppress meaningful
        variability.

        Practical Recommendations

        In routine biological applications, it is advisable to begin refinement using carefully curated
        subtomogram datasets with obvious artifacts and damaged particles removed beforehand. Initial
        references should be biologically plausible but not overly detailed in order to minimize reference
        bias.

        Iterative inspection of refinement results is strongly recommended. Users should evaluate whether
        structural features become progressively clearer while ensuring that biologically meaningful
        heterogeneity is not lost during averaging.

        When working with flexible or compositionally diverse complexes, classification before or during
        refinement may substantially improve structural quality.

        Final Perspective

        Refinement is one of the central stages of subtomogram averaging because it transforms noisy particle
        observations into biologically interpretable structural reconstructions. The Refinement Base protocol
        provides the organizational framework necessary to support these workflows consistently across EMAN
        tomography refinement strategies.

        Careful dataset preparation, appropriate reference selection, and iterative biological validation are
        the key elements for obtaining reliable and meaningful refinement results in cryo-electron tomography.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------- DEFINE param functions ----------------
    def _addCommonInputParams(self, form):
        form.addParam(IN_SUBTOMOS, PointerParam,
                      pointerClass='EmanSetOfParticles',
                      label='Particles',
                      important=True,
                      help='Select the set of subtomograms to build an initial model')
        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      allowsNull=True,
                      label="Reference volume (opt.)")
        self._addBinThreads(form)

    # --------------- INSERT steps functions ----------------

    # --------------- STEPS functions -----------------------

    # --------------- UTILS functions -----------------------

    # --------------- INFO functions ------------------------
