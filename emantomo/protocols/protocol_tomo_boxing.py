# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
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

from pyworkflow import BETA
from pyworkflow.utils.properties import Message
from pyworkflow import utils as pwutils
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import BooleanParam, PointerParam, EnumParam

from emantomo.convert import setCoords3D2Jsons, jsons2SetCoords3D, jsonFilesFromSet
from emantomo.viewers.views_tkinter_tree import EmanDialog

from tomo.protocols import ProtTomoPicking
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


class EmanProtTomoBoxing(ProtTomoPicking):
    """
    Provides an interactive environment for manual particle picking in cryo-electron tomography using the EMAN2 tomography boxing tools. The protocol allows users to visually inspect tomograms and define three-dimensional particle coordinates directly within the tomographic volume, supporting the creation of curated coordinate datasets for subtomogram averaging, particle extraction, or downstream structural analysis.

    AI Generated:

    Manual Tomogram Picking (EmanProtTomoBoxing) — User Manual
        Overview

        The Manual Tomogram Picking protocol is designed for interactive identification of particles within
        cryo-electron tomography datasets. Its primary purpose is to allow users to manually define
        biologically relevant coordinates directly inside reconstructed tomograms. This process is especially
        important in tomography projects where automatic particle detection may fail because of low contrast,
        structural heterogeneity, crowded cellular environments, or uncommon particle shapes.

        In practical biological workflows, manual picking is commonly used during exploratory studies,
        validation of automated picking strategies, preparation of high-quality training datasets, or
        generation of carefully curated particle sets for subtomogram averaging. The protocol provides a
        visual workflow where users inspect tomographic slices and place coordinates interactively, ensuring
        that selected particles correspond to meaningful biological structures.

        Biological Context and Typical Applications

        Cryo-electron tomography often involves highly heterogeneous samples such as organelles, membrane
        complexes, ribosomes inside cells, viral assemblies, or flexible macromolecular structures embedded
        in native environments. In these situations, fully automated approaches may incorrectly identify
        contaminants, noise, or unrelated densities as particles.

        Manual inspection remains one of the most reliable strategies when biological interpretation is
        critical. Users can visually discriminate intact complexes from damaged particles, recognize
        preferred conformations, and avoid regions affected by reconstruction artifacts or contamination.
        This makes manual picking especially valuable during the early stages of a project or when studying
        rare structural states.

        Inputs and Coordinate Management

        The protocol requires a set of tomograms that will serve as the visualization source for particle
        selection. Each tomogram is inspected independently, allowing users to navigate through slices and
        identify particles in three dimensions.

        Optionally, previously generated coordinates can be loaded and modified. This feature is especially
        useful in iterative workflows where users want to refine earlier annotations, remove false positives,
        add newly identified particles, or correct particle positions after additional biological inspection.
        Existing coordinate sets therefore become editable working datasets rather than immutable results.

        This capability is particularly important in collaborative environments where several rounds of
        annotation and validation may occur before producing a final curated particle collection.

        Interactive Picking Workflow

        The protocol launches an interactive graphical environment specialized for tomographic particle
        annotation. Users can browse individual tomograms, inspect orthogonal slice views, and place
        particle coordinates directly on relevant densities.

        In biological practice, users typically focus on identifying recognizable structural signatures such
        as membrane-associated densities, repeating assemblies, viral capsids, ribosomes, or filamentous
        structures. Because tomographic contrast may vary substantially between datasets, manual inspection
        provides flexibility that is difficult to achieve with rigid automated criteria.

        The workflow also supports revisiting tomograms multiple times. This is useful when users initially
        perform broad exploratory picking and later refine selections after visualizing averages or obtaining
        additional structural information.

        Coordinate Interpretation and Organization

        The resulting output consists of three-dimensional coordinate sets associated with the original
        tomograms. These coordinates define the spatial locations of particles and can subsequently be used
        for subtomogram extraction, alignment, averaging, or classification workflows.

        Particle box size plays an important biological role because it determines the amount of surrounding
        structural context included during extraction. Small box sizes may truncate flexible domains or
        membrane regions, whereas excessively large boxes increase computational cost and may introduce
        unnecessary background noise.

        Careful selection of coordinate positions is equally important. Coordinates should ideally be centered
        on structurally meaningful regions so that downstream averaging procedures converge reliably.

        Updating Existing Coordinate Sets

        One particularly useful aspect of the protocol is the ability to continue working from previously
        generated coordinate annotations. In many tomography projects, particle selection evolves over time
        as users gain additional biological insight or improve reconstruction quality.

        For example, a user may initially select all visible ribosome-like densities and later decide to
        exclude membrane-associated subpopulations or damaged particles after preliminary averaging. The
        protocol supports this iterative curation process naturally.

        This approach is especially valuable in cellular tomography studies where structural complexity makes
        one-pass annotation unrealistic.

        Practical Recommendations

        In routine biological applications, it is often advisable to begin with a relatively conservative
        picking strategy, selecting only the clearest and most interpretable particles. High-confidence
        particles typically produce better initial averages and can later guide broader particle selection.

        Users should inspect tomograms at multiple slice depths to confirm that candidate particles represent
        coherent three-dimensional structures rather than isolated noisy features. Careful attention should
        also be paid to crowded cellular regions, missing wedge artifacts, and reconstruction boundaries,
        which can complicate interpretation.

        When working with highly heterogeneous samples, maintaining separate coordinate sets for different
        structural classes may simplify downstream analysis and classification.

        Outputs and Biological Interpretation

        The protocol produces curated three-dimensional coordinate datasets linked to the original tomograms.
        These outputs preserve the spatial context of each selected particle and form the foundation for
        subsequent subtomogram analysis workflows.

        From a biological perspective, the quality of the coordinate set strongly influences all downstream
        processing stages. Accurate particle localization improves alignment stability, enhances average
        resolution, and reduces the inclusion of unrelated densities.

        Consequently, manual particle picking should not be viewed merely as a preparatory task, but as a
        biologically meaningful interpretation step that directly affects the structural conclusions obtained
        from tomography experiments.

        Final Perspective

        Manual picking remains one of the most important validation and curation tools in cryo-electron
        tomography. Although automated approaches continue to improve, direct user inspection provides a
        level of biological interpretation and contextual awareness that is often essential for challenging
        datasets.

        Careful coordinate placement, iterative refinement, and thoughtful inspection of tomographic features
        are key elements for generating reliable particle datasets suitable for high-quality subtomogram
        averaging and structural analysis.
    """
    _label = 'Manual picking'

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)
        self.info_path = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('selection', EnumParam, choices=['Yes', 'No'], default=1,
                      label='Modify previous coordinates?', display=EnumParam.DISPLAY_HLIST,
                      help='This option allows to add and/or remove coordinates to a previous SetOfCoordinates')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates", condition='selection == 0',
                      allowsNull=True, pointerClass='SetOfCoordinates3D',
                      help='Select the previous SetOfCoordinates you want to modify')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        # Copy input coordinates to Extra Path
        if self.inputCoordinates.get():
            self._insertFunctionStep(self.copyInputCoords)
        # Launch Boxing GUI
        self._insertFunctionStep(self.launchBoxingGUIStep, interactive=True)

    def _createOutput(self):
        jsons2SetCoords3D(self, self.inputTomograms.get(), self.info_path)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.info_path = self._getExtraPath('info')

    def copyInputCoords(self):
        pwutils.makePath(self.info_path)
        self.json_files, self.tomo_files = jsonFilesFromSet(self.inputTomograms.get(), self.info_path)
        setCoords3D2Jsons(self.json_files, self.inputCoordinates.get())

    def launchBoxingGUIStep(self):
        inTomos = self.inputTomograms.get()
        lastOutput = None
        if self.getOutputsSize() >= 1:
            pwutils.makePath(self.info_path)
            self.json_files, self.tomo_files = jsonFilesFromSet(inTomos, self.info_path)
            lastOutput = [output for _, output in self.iterOutputAttributes()][-1]
            setCoords3D2Jsons(self.json_files, lastOutput)

        if lastOutput is not None:
            volIds = lastOutput.aggregate(["MAX", "COUNT"], "_volId", ["_volId"])
            volIds = dict([(d['_volId'], d["COUNT"]) for d in volIds])
        else:
            volIds = dict()

        tomoList = []
        for tomo in inTomos.iterItems():
            tomogram = tomo.clone()
            if tomo.getObjId() in volIds:
                tomogram.count = volIds[tomo.getObjId()]
            else:
                tomogram.count = 0
            tomoList.append(tomogram)

        tomoProvider = TomogramsTreeProvider(tomoList, self.info_path, "json")

        self.dlg = EmanDialog(None, self._getExtraPath(), provider=tomoProvider,)

        # Open dialog to request confirmation to create output
        import tkinter as tk
        frame = tk.Frame()
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
            self._createOutput()

    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.inputTomograms is None:
            return ['Input tomogram not available yet.']

        methodsMsgs.append("Input tomograms imported of dims %s." %(
                              str(self.inputTomograms.get().getDim())))

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                msg = self.getMethods(output)
                methodsMsgs.append("%s: %s" % (self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs
