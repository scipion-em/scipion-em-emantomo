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
# Supported versions
VERSIONS = ['2.99.52']
EMAN_DEFAULT_VER_NUM = VERSIONS[-1]

# Extra attributes for refinement
EMAN_SCORE = '_eman_score'
EMAN_COVERAGE = '_eman_coverage'
EMAN_ALI_LOSS = '_eman_ali_loss'
EMAN_OFF_TILT_AXIS = '_eman_off_tilt_axis'

# Norm options
PROC_NORMALIZE = 0
PROC_EDGEMEAN = 1

# Info json fields
ALI_LOSS = 'ali_loss'
APIX_UNBIN = 'apix_unbin'
TS_FILE = 'tlt_file'
TLT_PARAMS = 'tlt_params'

# Default label for picking with eman
TOMOBOX = 'tomobox'

# SQLITE labels to manage the present objects in the input sets
TS_ID = "_tsId"
GROUP_ID = "_groupId"

# EMAN alignment LST files data
PARTICLE_IND = 'particleInd'
PARTICLE_FILE = 'particleFile'
SCORE = 'score'
MATRIX = 'matrix'
CLASS = 'class'
DEFOCUS = 'defocus'
PART3D_ID = 'ptcl3d_id'
TILT_ID = 'tilt_id'
PROJ_MATRIX = 'xform.projection'
ROT_TR_MATRIX = 'xform.align3d'
TR_MATRIX = 'translationMatrix'
LST_LINE = 'lstLine'

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
INIT_MODEL_NAME = 'output_cls%i.hdf'
INIT_MODEL_MRC = 'initialModel%i.mrc'
INIR_MODEL_NAME_OLD = 'output.hdf'
SPT_00_DIR = 'spt_00'
SPT_CLS_00_FIR = 'sptcls_00'
REFS_DIR = 'references'
REFERENCE_NAME = 'reference'
INPUT_PTCLS_LST = 'input_ptcls.lst'
THREED = 'threed'
THREED_01 = 'threed_01'
SPTCLS_00_DIR = 'sptcls_00'
ALI2D_BASENAME = 'aliptcls2d_'
ALI3D_BASENAME = 'aliptcls3d_'
FSC_MASKED_BNAME = 'fsc_masked_'
FSC_UNMASKED_BNAME = 'fsc_unmasked_'
FSC_MASKED_TIGHT_BNAME = 'fsc_maskedtight_'

SYMMETRY_HELP_MSG = 'Supported symmetries are :\n\n' \
                    'Cn - single rotational n-fold symmetry axis\n' \
                    'Dn - dihedral rotational n-fold, eg - GroEL is ~D7\n' \
                    'icos - icosahedral, (5,3,2)\n' \
                    'oct - octahedral, symmetry of a cube (4,3,2)\n' \
                    'tet - tetrahedral, (3,2)\n' \
                    'h - helical\n\n' \
                    'See http://blake.bcm.edu/emanwiki/EMAN2/Symmetry for a detailed description of symmetry in Eman.'


