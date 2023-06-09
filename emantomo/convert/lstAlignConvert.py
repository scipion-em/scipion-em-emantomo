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
import ast
import json

import numpy as np

from emantomo.constants import PARTICLE_IND, PARTICLE_FILE, ROT_TR_MATRIX, SCORE, MATRIX, EMAN_SCORE, CLASS, DEFOCUS, \
    PART3D_ID, TR_MATRIX, PROJ_MATRIX, TILT_ID, LST_LINE
from pwem.objects import Transform
from pyworkflow.object import Float
from tomo.constants import TR_EMAN


class EmanLstReader:

    @staticmethod
    def read(lstFileName):
        """Reads an existing 3d alignment LST file. Example of the file contents:

        #LSX
        # This file is in fast LST format. All lines after the next line have exactly the number of characters shown on the next line. This MUST be preserved if editing.
        # 236
        0	particles3d/ts_25__tomobox.hdf	{"score":-0.27148754218956683,"xform.align3d":{"__class__":"Transform","matrix":"[0.117509,-0.107427,0.987244,-0.532006,0.992805,-0.0103203,-0.119294,-7.39777,0.023004,0.994159,0.105441,-1.86965]"}}
        1	particles3d/ts_25__tomobox.hdf	{"score":-0.22619236232343745,"xform.align3d":{"__class__":"Transform","matrix":"[0.53198,0.111643,-0.839364,-0.574314,0.838809,-0.204981,0.504364,-7.58973,-0.115745,-0.972378,-0.202693,2.56157]"}}
        [...]

        :param lstFileName: path of the LST file.
        :returns: a list of dictionaries based on the lists filled before, one for each line (3d particle) following
        each element this structure:
                {
                PARTICLE_IND: particle index in the corresponding 3d HDF file stack,
                PARTICLE_FILE: path of the 3d HDF stack containing the current particle,
                SCORE: EMAN's calculated score,
                ROT_TR_MATRIX: 3d transformation matrix calculated by EMAN.
                LST_LINE: current line as in the original file. Useful for sub-setting.
                }
        """
        keys = [PARTICLE_IND, PARTICLE_FILE, SCORE, ROT_TR_MATRIX, LST_LINE]
        list_of_lists = []
        lastRow = np.array([0, 0, 0, 1])

        with open(lstFileName, 'r') as f:
            for line in f:
                lineContentsList = line.split('\t')
                if len(lineContentsList) > 1:  # There are some explicative lines at the beginning
                    jsonData = json.loads(lineContentsList[2])
                    matrixAsList = ast.literal_eval(jsonData[ROT_TR_MATRIX][MATRIX])
                    matrix = np.array(matrixAsList).reshape(3, 4)
                    matrix = np.vstack([matrix, lastRow])
                    list_of_lists.append([
                        lineContentsList[0],
                        lineContentsList[1],
                        jsonData[SCORE],
                        matrix,
                        line
                    ])
                else:  # Header lines (required to replicate file when sub-setting)
                    list_of_lists.append([
                        None,
                        None,
                        None,
                        None,
                        line
                    ])

        return [dict(zip(keys, values)) for values in list_of_lists]

    @staticmethod
    def align3dLst2Scipion(lstFileName, inParticles, outParticles):
        """Converts the data from an existing EMAN's 3d align LST file into a Scipion EmanSetOfParticles object.
        :param lstFileName: path of the LST file.
        :param inParticles: input EmanSetOfParticles.
        :param outParticles: output EmanSetOfParticles, expected to contain the set info. It will be filled here with
        the updated particles.
        """
        align3dData = EmanLstReader.read(lstFileName)
        for particle, alignDict in zip(inParticles, align3dData):
            outParticle = particle.clone()
            setattr(outParticle, EMAN_SCORE, Float(alignDict[SCORE]))
            outParticle.setTransform(Transform(alignDict[ROT_TR_MATRIX]), convention=TR_EMAN)
            outParticles.append(outParticle)

    @staticmethod
    def align2dLst2Scipion(lstFileName):
        """Reads an existing 2d alignment LST file. Example of the file contents:

        #LSX
        # This file is in fast LST format. All lines after the next line have exactly the number of characters shown on the next line. This MUST be preserved if editing.
        # 382
        0	particles/ts_25__tomobox.hdf	{"class":0,"defocus":0,"dxf":{"__class__":"Transform","matrix":"[1,0,0,-0.345638,-0,1,0,0.0120683,0,-0,1,0]"},"ptcl3d_id":0,"score":-0.727735565517016,"tilt_id":0,"xform.projection":{"__class__":"Transform","matrix":"[0.994593,-0.0657139,-0.0804149,1.01098,-0.0332251,-0.934989,0.353117,-3.19245,-0.0983923,-0.348536,-0.932117,5.92621]"}}
        1	particles/ts_25__tomobox.hdf	{"class":0,"defocus":0,"dxf":{"__class__":"Transform","matrix":"[1,0,0,-1.24626,-0,1,0,0.130769,0,-0,1,0]"},"ptcl3d_id":0,"score":-0.5650576706552629,"tilt_id":1,"xform.projection":{"__class__":"Transform","matrix":"[0.995015,-0.0643083,-0.0762159,0.701939,-0.0284292,-0.9155,0.401314,-3.81944,-0.095584,-0.397146,-0.912765,5.72374]"}}
        [...]

        :returns: a list of dictionaries based on the lists filled before, one for each line (3d particle) following
        each element this structure:
                {
                PARTICLE_IND: particle index in the corresponding 2d HDF file stack,
                PARTICLE_FILE: path of the 3d HDF stack containing the current particle,
                SCORE: EMAN's calculated score,
                CLASS: class number,
                DEFOCUS: defocus value,
                PART3D_ID: index of the corresponding 3d particle in the HDF stack file
                TILT_ID: index of the tilt image in the corresponding stack,
                TR_MATRIX: translation matrix,
                PROJ_MATRIX: projection matrix.
                LST_LINE: current line as in the original file. Useful for sub-setting.
                }
        """
        keys = [PARTICLE_IND, PARTICLE_FILE, SCORE, CLASS, DEFOCUS, PART3D_ID, TILT_ID, TR_MATRIX,
                PROJ_MATRIX, LST_LINE]
        list_of_lists = []
        lastRow = np.array([0, 0, 0, 1])

        with open(lstFileName, 'r') as f:
            for line in f:
                lineContentsList = line.split('\t')
                if len(lineContentsList) > 1:
                    jsonData = json.loads(lineContentsList[2])
                    particleInd = lineContentsList[0]
                    particleFile = lineContentsList[1]
                    nClass = jsonData[CLASS]
                    defocus = jsonData[DEFOCUS]
                    trMatrix = ast.literal_eval(jsonData['dxf'][MATRIX])
                    trMatrix = np.array(trMatrix).reshape(3, 4)
                    trMatrix = np.vstack([trMatrix, lastRow])
                    ptcl3dId = jsonData[PART3D_ID]
                    score = jsonData[SCORE]
                    tiltId = jsonData[TILT_ID]
                    projMatrix = ast.literal_eval(jsonData[PROJ_MATRIX][MATRIX])
                    projMatrix = np.array(projMatrix).reshape(3, 4)
                    projMatrix = np.vstack([projMatrix, lastRow])
                    list_of_lists.append(
                        [particleInd, particleFile, score, nClass, defocus, ptcl3dId, tiltId,
                         trMatrix, projMatrix, line])

        return [dict(zip(keys, values)) for values in list_of_lists]


class EmanLstWriter:

    def scipion2Align3dLst(self):
        pass
