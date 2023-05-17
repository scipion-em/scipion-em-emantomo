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
# Extra attributes for refinement
EMAN_SCORE = 'eman_score'
EMAN_COVERAGE = 'eman_coverage'

# Default label for picking with eman
TOMOBOX = 'tomobox'

# SQLITE labels to manage the present objects in the input sets
TS_ID = "_tsId"
GROUP_ID = "_groupId"

# Directories and files of an EMAN project
INFO_DIR = 'info'
TS_DIR = 'tiltseries'
TLT_DIR = 'tlt'
INTERP_TS = 'tiltseries_ali'
TOMOGRAMS_DIR = 'tomograms'
PARTICLES_DIR = 'particles'
PARTICLES_3D_DIR = 'particles3d'
SUBTOMOGRAMS_DIR = 'subtomograms'
SETS_DIR = 'sets'
INIT_MODEL_DIR = 'sptsgd_00'
INIT_MODEL_NAME = 'output_cls0.hdf'
INIT_MODEL_MRC = 'initialModel.mrc'
INIR_MODEL_NAME_OLD = 'output.hdf'
SPT_00_DIR = 'spt_00'
REFS_DIR = 'references'
REFERENCE_NAME = 'reference'
PARTICLES_LST_FILE = '%s.lst' % TOMOBOX
INPUT_PTCLS_LST = 'input_ptcls.lst'
THREED_01 = 'threed_01'
SPTCLS_00_DIR = 'sptcls_00'

SYMMETRY_HELP_MSG = 'Supported symmetries are :\n\n' \
                    'Cn - single rotational n-fold symmetry axis\n' \
                    'Dn - dihedral rotational n-fold, eg - GroEL is ~D7\n' \
                    'icos - icosahedral, (5,3,2)\n' \
                    'oct - octahedral, symmetry of a cube (4,3,2)\n' \
                    'tet - tetrahedral, (3,2)\n' \
                    'h - helical\n\n' \
                    'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for a detailed description of symmetry in Eman.'


