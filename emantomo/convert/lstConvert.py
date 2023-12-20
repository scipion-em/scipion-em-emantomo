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
from os.path import basename, join
import numpy as np
from emantomo.constants import PARTICLE_IND, PARTICLE_FILE, ROT_TR_MATRIX, SCORE, MATRIX, CLASS, DEFOCUS, \
    PART3D_ID, TR_MATRIX, PROJ_MATRIX, TILT_ID, LST_LINE, PARTICLES_3D_DIR
from emantomo.utils import getPresentPrecedents
from pwem.objects import Transform
from pyworkflow.utils import removeBaseExt
from tomo.constants import TR_EMAN
from tomo.objects import Coordinate3D


class EmanLstReader:

    @staticmethod
    def _getClassIdFromFileName(lstFileName):
        """The name structure is name_iter_class.lst  --> aliptcls2d_05_00.lst"""
        lstBaseName = removeBaseExt(lstFileName)
        return int(lstBaseName.split('_')[-1])

    @staticmethod
    def read3dParticles(lstFileName):
        """Reads an existing 3d alignment LST file [name_iter_class.lst  --> aliptcls2d_05_00.lst]. Example of the file contents:

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
                CLASS: class id read from the filename,
                ROT_TR_MATRIX: 3d transformation matrix calculated by EMAN,
                LST_LINE: current line as in the original file. Useful for sub-setting.
                }
        """
        keys = [PARTICLE_IND, PARTICLE_FILE, SCORE, CLASS, ROT_TR_MATRIX, LST_LINE]
        list_of_lists = []
        lastRow = np.array([0, 0, 0, 1])

        with open(lstFileName, 'r') as f:
            classId = EmanLstReader._getClassIdFromFileName(lstFileName)
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
                        classId,
                        matrix,
                        line
                    ])

        return [dict(zip(keys, values)) for values in list_of_lists]

    @staticmethod
    def align3dLst2Scipion(lstFileNames, inParticles, outParticles):
        """Converts the data from an existing EMAN's 3d align LST file into a Scipion EmanSetOfParticles object.
        :param lstFileNames: path of the LST file. It can also be a list of list of LST files.
        :param inParticles: input EmanSetOfParticles.
        :param outParticles: output EmanSetOfParticles, expected to contain the set info. It will be filled here with
        the updated particles.
        """
        if isinstance(lstFileNames, str):
            align3dData = EmanLstReader.read3dParticles(lstFileNames)
            for particle, alignDict in zip(inParticles, align3dData):
                if alignDict[PARTICLE_IND]:  # Will be None for the header lines
                    outParticle = particle.clone()
                    outParticle.setTransform(Transform(alignDict[ROT_TR_MATRIX]), convention=TR_EMAN)
                    outParticle.setEmanScore(alignDict[SCORE])
                    outParticles.append(outParticle)

        else:  # Multiple classification case. To avoid a DDBB access for each particle, the whole set is read once and
            # structured into a list of dictionaries to carry out the particle matching
            align3dData = []
            IN_PARTICLE_OBJ = 'inParticleObj'
            for lstFile in lstFileNames:
                align3dData.extend(EmanLstReader.read3dParticles(lstFile))
            inParticlesDicts = [{IN_PARTICLE_OBJ: particle.clone(),
                                 PARTICLE_IND: particle.getIndex(),
                                 PARTICLE_FILE: join(PARTICLES_3D_DIR, basename(particle.getStack3dHdf()))}
                                for particle in inParticles]

            for classId, alignDict in enumerate(align3dData):
                # The filter() function is more efficient than using a list comprehension because it doesn’t create a
                # new list in memory. Instead, it returns an iterator that generates the filtered values on-the-fly
                matchingInPartDict = list(filter(lambda d:
                                                 d.get(PARTICLE_FILE) == alignDict.get(PARTICLE_FILE) and
                                                 d.get(PARTICLE_IND) == int(alignDict.get(PARTICLE_IND)),
                                                 inParticlesDicts))
                outParticle = matchingInPartDict[0].get(IN_PARTICLE_OBJ)
                outParticle.setTransform(Transform(alignDict[ROT_TR_MATRIX]), convention=TR_EMAN)
                outParticle.setEmanScore(alignDict[SCORE])
                outParticle.setClassId(alignDict[CLASS])
                outParticles.append(outParticle)

        return align3dData

    @staticmethod
    def read2dParticles(lstFileName):
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

    @staticmethod
    def writeSimpleLst(emanParticles, outLstFileName):
        """Writes a simple list of particle files in LST format. Example:

        #LSX
        # This file is in fast LST format. All lines after the next line have exactly the number of characters shown on the next line. This MUST be preserved if editing.
        # 37
        0	particles3d/tomoa__tomobox.hdf
        1	particles3d/tomoa__tomobox.hdf
        0	particles3d/tomoc__tomobox.hdf
        1	particles3d/tomoc__tomobox.hdf
        2	particles3d/tomoc__tomobox.hdf
        3	particles3d/tomoc__tomobox.hdf

        Thus, only the stack file and the particle index in the corresponding stack are required.

        """
        # Get the 3d stacks present in the introduced set
        coords = emanParticles.getCoordinates3D()
        presentTsIds = coords.getUniqueValues(Coordinate3D.TOMO_ID_ATTR)
        presentPrecedents = getPresentPrecedents(coords, presentTsIds)
        # Prepare contents
        lines = []
        for tomo in sorted(presentPrecedents, key=lambda t: t.getTsId()):
            for emanParticle in emanParticles.iterSubtomos(tomo):
                lines.append(f'{emanParticle.getIndex()}\t'
                             f'{join(PARTICLES_3D_DIR, basename(emanParticle.getStack3dHdf()))}')

        # Write the LST file
        extraPad = 1 + 1 + 4  # index, tab, 4 blank spaces added by EMAN
        EmanLstWriter.lines2LstFile(lines, outLstFileName, extraPad=extraPad)

    @staticmethod
    def writeAlign3dLst(emanParticles, outLstFileName, isExtractingOrientedPicking=False):
        # Get the 3d stacks present in the introduced set
        coords = emanParticles.getCoordinates3D()
        presentTsIds = coords.getUniqueValues(Coordinate3D.TOMO_ID_ATTR)
        presentPrecedents = getPresentPrecedents(coords, presentTsIds)
        presentPrecedents = sorted(presentPrecedents, key=lambda t: t.getTsId())
        # Prepare contents
        lines = []
        for tomo in presentPrecedents:
            for emanParticle in emanParticles.iterSubtomos(tomo):
                score = emanParticle.getEmanScore()
                if isExtractingOrientedPicking:  # To preserve the oriented picking, the particles should be loaded as
                    # 3d oriented particles
                    matrix = emanParticle.getCoordinate3D().getMatrix(convention=TR_EMAN)
                else:
                    matrix = emanParticle.getTransform(convention=TR_EMAN).getMatrix()
                # Remove the last row and transform into a flatten list to write it in a json string, as expected
                matrix = matrix[:3, :].flatten().tolist()
                particleDict = {
                    SCORE: score,
                    ROT_TR_MATRIX: {
                        "__class__": "Transform",
                        MATRIX: str(matrix)
                    }
                }
                lines.append(f'{emanParticle.getIndex()}\t'
                             f'{join(PARTICLES_3D_DIR, basename(emanParticle.getStack3dHdf()))}\t'
                             f'{json.dumps(particleDict)}'.replace(' ', ''))

        # Write the LST file
        extraPad = 1 + 1 + 4  # index, tab, 4 blank spaces added by EMAN
        EmanLstWriter.lines2LstFile(lines, outLstFileName, extraPad=extraPad)

    @staticmethod
    def lines2LstFile(lines, outLstFileName, extraPad=0):
        linesLen = EmanLstWriter.getLongestLineLen(lines, extraPad=extraPad)
        with open(outLstFileName, 'w') as outLst:
            EmanLstWriter.writeHeaderLines(outLst, linesLen)
            outLst.writelines([f'{line}'.ljust(linesLen) + '\n' for line in lines])

    @staticmethod
    def getLongestLineLen(lines, extraPad=0):
        longestBaseName = max(lines, key=len)
        return len(join(PARTICLES_3D_DIR, longestBaseName)) + extraPad

    @staticmethod
    def writeHeaderLines(fileId, linesLen):
        fileId.write('#LSX\n')
        fileId.write('# This file is in fast LST format. All lines after the next line have exactly the number '
                     'of characters shown on the next line. This MUST be preserved if editing.\n')
        fileId.write(f'# {linesLen + 1}\n')  # + 1 is for the last \n added later when writing the file


