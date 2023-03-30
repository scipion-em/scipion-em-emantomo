# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Laura del Cano (ldelcano@cnb.csic.es) [1]
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import json
from os.path import join

import numpy as np

import emantomo
import numpy
import os
from ast import literal_eval

import pwem.constants as emcts
from pwem.objects import FSC
import pyworkflow.utils as pwutils
from pwem.objects.data import Transform
from pyworkflow.object import Float, RELATION_SOURCE, OBJECT_PARENT_ID, Pointer

import tomo.constants as const
from tomo.objects import SetOfTiltSeries, SetOfTomograms
from tomo.constants import TR_EMAN

from .. import Plugin
from emantomo.constants import EMAN_SCORE, EMAN_COVERAGE


def loadJson(jsonFn):
    """ This function loads the Json dictionary into memory """
    with open(jsonFn) as jsonFile:
        jsonDict = json.load(jsonFile)

    return jsonDict


def writeJson(jsonDict, jsonFn, indent=None):
    """ This function write a Json dictionary """
    with open(jsonFn, 'w') as outfile:
        json.dump(jsonDict, outfile, indent=indent)


def appendJson(jsonDict, jsonFn, indent=None):
    """ Append a new dictionary to a already existing Json file"""
    with open(jsonFn, 'r+') as outfile:
        data = json.load(outfile)
        data.update(jsonDict)
        outfile.seek(0)
        json.dump(data, outfile, indent=indent)


def readCTFModel(ctfModel, filename):
    jsonDict = loadJson(filename)
    keyPos = None
    ctfPhaseShift = 0.0

    if 'ctf_frame' in jsonDict:
        keyPos = jsonDict['ctf_frame'][1]
    elif 'ctf' in jsonDict:
        keyPos = jsonDict['ctf'][0]
    else:
        setWrongDefocus(ctfModel)

    if keyPos:
        defocus = float(keyPos['defocus'])
        defocusAngle = float(keyPos['dfang'])
        dfdiff = float(keyPos['dfdiff'])
        ampcont = float(keyPos['ampcont'])
        defocusU = 10000.0 * defocus + 5000.0 * dfdiff
        defocusV = 20000.0 * defocus - defocusU
        ctfPhaseShift = calculatePhaseShift(ampcont)

        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
        if 'ctf_im2d' in jsonDict:
            # psdFile = jsonDict['ctf_im2d']['__image__'][0]
            fnBase = pwutils.removeExt(filename) + '_jsonimg'
            psdFile = "1@%s.hdf" % fnBase
            if pwutils.exists(psdFile):
                ctfModel.setPsdFile(psdFile)
    ctfModel.setPhaseShift(float(ctfPhaseShift))


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)


def writeCTFModel(ctfObj, filename):
    """ Write a CTFModel object as Xmipp .ctfparam"""
    pass


def jsonToCtfModel(ctfJsonFn, ctfModel):
    """ Create a CTFModel from a json file """
    mdFn = str(ctfJsonFn).replace('particles', 'info')
    mdFn = mdFn.split('__ctf_flip')[0] + '_info.json'
    if pwutils.exists(mdFn):
        readCTFModel(ctfModel, mdFn)


# def readSetOfCoordinates(workDir, micSet, coordSet, invertY=False, newBoxer=False):
#     """ Read from Eman .json files.
#     Params:
#         workDir: where the Eman boxer output files are located.
#         micSet: the SetOfMicrographs to associate the .json, which
#             name should be the same of the micrographs.
#         coordSet: the SetOfCoordinates that will be populated.
#     """
#     if newBoxer:
#         # read boxSize from info/project.json
#         jsonFnbase = pwutils.join(workDir, 'info', 'project.json')
#         jsonBoxDict = loadJson(jsonFnbase)
#         size = int(jsonBoxDict["global.boxsize"])
#     else:
#         # read boxSize from e2boxercache/base.json
#         jsonFnbase = pwutils.join(workDir, 'e2boxercache', 'base.json')
#         jsonBoxDict = loadJson(jsonFnbase)
#         size = int(jsonBoxDict["box_size"])
#
#     jsonFninfo = pwutils.join(workDir, 'info/')
#
#     for mic in micSet:
#         micBase = pwutils.removeBaseExt(mic.getFileName())
#         micPosFn = ''.join(glob.glob(jsonFninfo + '*' + micBase + '_info.json'))
#         readCoordinates(mic, micPosFn, coordSet, invertY)
#     coordSet.setBoxSize(size)


def readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, inputTomo, updateItem=None,
                           origin=const.BOTTOM_LEFT_CORNER, scale=1, groupId=None):
    if "boxes_3d" in jsonBoxDict.keys():
        boxes = jsonBoxDict["boxes_3d"]

        for box in boxes:
            classKey = box[5]
            coord3DSet = coord3DSetDict[classKey]
            coord3DSet.enableAppend()

            newCoord = readCoordinate3D(box, inputTomo, origin=origin, scale=scale)
            if groupId is None:
                newCoord.setGroupId(classKey)
            else:
                newCoord.setGroupId(groupId)

            # Execute Callback
            if updateItem:
                updateItem(newCoord)

            coord3DSet.append(newCoord)


# def readCoordinates(mic, fileName, coordsSet, invertY=False):
#     if pwutils.exists(fileName):
#         jsonPosDict = loadJson(fileName)
#
#         if "boxes" in jsonPosDict:
#             boxes = jsonPosDict["boxes"]
#
#             for box in boxes:
#                 x, y = box[:2]
#
#                 if invertY:
#                     y = mic.getYDim() - y
#
#                 coord = Coordinate()
#                 coord.setPosition(x, y)
#                 coord.setMicrograph(mic)
#                 coordsSet.append(coord)


def readCoordinate3D(box, inputTomo, origin=const.BOTTOM_LEFT_CORNER, scale=1):
    from tomo.objects import Coordinate3D
    x, y, z = scale * numpy.asarray(box[:3])
    coord = Coordinate3D()
    coord.setVolume(inputTomo)
    coord.setPosition(x, y, z, origin)
    return coord


def writeSetOfSubTomograms(subtomogramSet, path, **kwargs):
    """ Convert the imgSet particles to .hdf files as expected by Eman.
        This function should be called from a current dir where
        the images in the set are available.
        """
    firstItem = subtomogramSet.getFirstItem()
    ext = pwutils.getExt(firstItem.getFileName())[1:]
    if ext == 'hdf':
        # create links if input has hdf format
        for fn in subtomogramSet.getFiles():
            newFn = pwutils.removeBaseExt(fn).split('__ctf')[0] + '.hdf'
            newFn = pwutils.join(path, newFn)
            pwutils.createLink(fn, newFn)
            print("   %s -> %s" % (fn, newFn))
    else:
        fileName = ""
        a = 0
        proc = Plugin.createEmanProcess(args='write')

        for subtomo in subtomogramSet.iterItems(orderBy=['_volId', 'id'], direction='ASC'):
            volName, _ = os.path.splitext(os.path.basename(subtomo.getVolName()))

            objDict = subtomo.getObjDict()

            suffix = kwargs.get('suffix', '')

            objDict['hdfFn'] = pwutils.join(path, "%s%s.hdf" % (volName, suffix))

            alignType = kwargs.get('alignType')

            if alignType != emcts.ALIGN_NONE:
                shift, angles = alignmentToRow(subtomo.getTransform(convention=TR_EMAN), alignType)
                # json cannot encode arrays so I convert them to lists
                # json fail if has -0 as value
                objDict['_shifts'] = shift.tolist()
                objDict['_angles'] = angles.tolist()
            objDict['_itemId'] = subtomo.getObjId()

            objFn = objDict['_filename']
            # the index in EMAN begins with 0
            if fileName != objFn:
                fileName = objFn.split(":")[0]
                objDict['_filename'] = fileName
                if objDict['_index'] == 0:
                    a = 0
                else:
                    a = 1
            objDict['_index'] = int(objDict['_index'] - a)

            # Write the e2converter.py process from where to read the image
            print(json.dumps(objDict), file=proc.stdin, flush=True)
            proc.stdout.readline()
        proc.kill()


# def writeSetOfMicrographs(micSet, filename):
#     """ Simplified function borrowed from xmipp. """
#     mdata = md.MetaData()
#
#     for img in micSet:
#         objId = mdata.addObject()
#         imgRow = md.Row()
#         imgRow.setValue(md.MDL_ITEM_ID, objId)
#
#         index, fname = img.getLocation()
#         fn = ImageHandler.locationToXmipp((index, fname))
#         imgRow.setValue(md.MDL_MICROGRAPH, fn)
#
#         if img.isEnabled():
#             enabled = 1
#         else:
#             enabled = -1
#         imgRow.setValue(md.MDL_ENABLED, enabled)
#         imgRow.writeToMd(mdata, objId)
#
#     mdata.write('Micrographs@%s' % filename)


# def readSetOfParticles(lstFile, partSet, copyOrLink, direc):
#     for index, fn in iterLstFile(lstFile):
#         item = Particle()
#         # set full path to particles stack file
#         abspath = os.path.abspath(lstFile)
#         fn = abspath.replace('sets/%s' % os.path.basename(lstFile), '') + fn
#         newFn = pwutils.join(direc, os.path.basename(fn))
#         if not pwutils.exists(newFn):
#             copyOrLink(fn, newFn)
#
#         item.setLocation(index, newFn)
#         partSet.append(item)
#
#
# def writeSetOfParticles(partSet, path, **kwargs):
#     """ Convert the imgSet particles to .hdf files as expected by Eman.
#     This function should be called from a current dir where
#     the images in the set are available.
#     """
#     ext = pwutils.getExt(partSet.getFirstItem().getFileName())[1:]
#     if ext == 'hdf':
#         # create links if input has hdf format
#         for fn in partSet.getFiles():
#             newFn = pwutils.removeBaseExt(fn).split('__ctf')[0] + '.hdf'
#             newFn = pwutils.join(path, newFn)
#             pwutils.createLink(fn, newFn)
#             print("   %s -> %s" % (fn, newFn))
#     else:
#         firstCoord = partSet.getFirstItem().getCoordinate() or None
#         hasMicName = False
#         if firstCoord:
#             hasMicName = firstCoord.getMicName() or False
#
#         fileName = ""
#         a = 0
#         proc = Plugin.createEmanProcess(args='write')
#
#         for i, part in iterParticlesByMic(partSet):
#             micName = micId = part.getMicId()
#             if hasMicName:
#                 micName = pwutils.removeBaseExt(part.getCoordinate().getMicName())
#             objDict = part.getObjDict()
#
#             if not micId:
#                 micId = 0
#
#             suffix = kwargs.get('suffix', '')
#             if hasMicName and (micName != str(micId)):
#                 objDict['hdfFn'] = pwutils.join(path,
#                                                 "%s%s.hdf" % (micName, suffix))
#             else:
#                 objDict['hdfFn'] = pwutils.join(path,
#                                                 "mic_%06d%s.hdf" % (micId, suffix))
#
#             alignType = kwargs.get('alignType')
#
#             if alignType != emcts.ALIGN_NONE:
#                 shift, angles = alignmentToRow(part.getTransform(), alignType)
#                 # json cannot encode arrays so I convert them to lists
#                 # json fail if has -0 as value
#                 objDict['_shifts'] = shift.tolist()
#                 objDict['_angles'] = angles.tolist()
#             objDict['_itemId'] = part.getObjId()
#
#             # the index in EMAN begins with 0
#             if fileName != objDict['_filename']:
#                 fileName = objDict['_filename']
#                 if objDict['_index'] == 0:  # TODO: Index appears to be the problem (when not given it works ok)
#                     a = 0
#                 else:
#                     a = 1
#             objDict['_index'] = int(objDict['_index'] - a)
#             # Write the e2converter.py process from where to read the image
#             print(json.dumps(objDict), file=proc.stdin, flush=True)
#             proc.stdout.readline()
#         proc.kill()


def getImageDimensions(imageFile):
    """ This function will allow us to use EMAN2 to read some formats
     not currently supported by the native image library (Xmipp).
     Underneath, it will call a script to do the job.
    """
    proc = Plugin.createEmanProcess('e2ih.py', args=imageFile)
    return tuple(map(int, proc.stdout.readline().split()))


def convertImage(inputLoc, outputLoc):
    """ This function will allow us to use EMAN2 to write some formats
     not currently supported by the native image library (Xmipp).
     Underneath, it will call an script to do the job.
    """

    def _getFn(loc):
        """ Use similar naming convention as in Xmipp.
        This does not works for EMAN out of here.
        """
        if isinstance(loc, tuple):
            if loc[0] != emcts.NO_INDEX:
                return "%06d@%s" % loc
            return loc[1]
        else:
            return loc

    proc = Plugin.createEmanProcess('e2ih.py', args='%s %s' % (_getFn(inputLoc),
                                                               _getFn(outputLoc)))
    proc.wait()


def iterLstFile(filename):
    with open(filename) as f:
        for line in f:
            if '#' not in line:
                # Decompose Eman filename
                index, filename = int(line.split()[0]) + 1, line.split()[1]
                yield index, filename


def geometryFromMatrix(matrix, inverseTransform, axes='szyz'):
    """ Convert the transformation matrix to shifts and angles."""
    from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -numpy.rad2deg(euler_from_matrix(matrix, axes=axes))
    return shifts, angles


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from given
    2D shifts in X and Y and the 3 euler angles."""
    from pwem.convert.transformations import euler_matrix
    from numpy import deg2rad
    radAngles = -deg2rad(angles)

    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        from numpy.linalg import inv
        M[:3, 3] = -shifts[:3]
        M = inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def alignmentToRow(alignment, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    #     is2D = alignType == em.ALIGN_2D
    #     inverseTransform = alignType == em.ALIGN_PROJ

    # transformation matrix is processed here because
    # it uses routines available through scipion python
    matrix = alignment.getMatrix()
    return geometryFromMatrix(matrix, True)


def rowToAlignment(alignmentList, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
        """
    # use all angles in 2D since we might have mirrors
    # is2D = alignType == em.ALIGN_2D
    inverseTransform = alignType == emcts.ALIGN_PROJ

    alignment = Transform()
    angles = numpy.zeros(3)
    shifts = numpy.zeros(3)
    shifts[0] = alignmentList[3]
    shifts[1] = alignmentList[4]
    shifts[2] = 0
    angles[0] = alignmentList[0]
    angles[1] = alignmentList[1]
    angles[2] = alignmentList[2]

    matrix = matrixFromGeometry(shifts, angles, inverseTransform)
    alignment.setMatrix(matrix)

    return alignment


# def iterParticlesByMic(partSet):
#     """ Iterate the particles ordered by micrograph """
#     for i, part in enumerate(partSet.iterItems(orderBy=['_micId', 'id'],
#                                                direction='ASC')):
#         yield i, part


def iterSubtomogramsByVol(subtomogramSet):
    """ Iterate subtomograms ordered by tomogram """
    items = [subtomo.clone() for subtomo in subtomogramSet.iterItems(orderBy=['_volId', 'id'], direction='ASC')]
    for i, part in enumerate(items):
        yield i, part


def convertReferences(refSet, outputFn):
    """ Simplified version of writeSetOfParticles function.
    Writes out an hdf stack.
    """
    fileName = ""
    a = 0
    proc = Plugin.createEmanProcess(args='write')

    for part in refSet:
        objDict = part.getObjDict()
        objDict['hdfFn'] = outputFn
        objDict['_itemId'] = part.getObjId()

        # the index in EMAN begins with 0
        if fileName != objDict['_filename']:
            fileName = objDict['_filename']
            if objDict['_index'] == 0:
                a = 0
            else:
                a = 1
        objDict['_index'] = int(objDict['_index'] - a)

        # Write the e2converter.py process from where to read the image
        print(json.dumps(objDict), file=proc.stdin)
        proc.stdin.flush()
        proc.stdout.readline()
    proc.kill()


def calculatePhaseShift(ampcont):
    # calculate phase shift as in EMAN2 ctf.cpp
    if -100.0 < ampcont <= 100.0:
        PhaseShift = numpy.arcsin(ampcont / 100.0)
    elif ampcont > 100.0:
        PhaseShift = numpy.pi - numpy.arcsin(2.0 - ampcont / 100.0)
    else:
        PhaseShift = -numpy.pi - numpy.arcsin(-2.0 - ampcont / 100.0)
    ctfPhaseShift = numpy.rad2deg(PhaseShift)

    return ctfPhaseShift


def getLastParticlesParams(directory):
    """
    Return a dictionary containing the params values of the last iteration.

    Key: Particle index (int)
    Value: Dict[{coverage: float, score: float, alignMatrix: list[float]}]
    """
    # JSON files with particles params: path/to/particle_parms_NN.json
    particleParamsPaths = glob.glob(os.path.join(directory, 'particle_parms_*.json'))
    if not particleParamsPaths:
        raise Exception("Particle params files not found")

    lastParticleParamsPath = sorted(particleParamsPaths)[-1]
    particlesParams = json.load(open(lastParticleParamsPath))
    output = {}

    for key, values in particlesParams.items():
        # key: '(path/to/particles/basename.hdf', nParticle)'
        # values: '{"coverage": 1.0, "score": 2.0, "xform.align3d": {"matrix": [...]}}'
        import re
        match = re.search(r'(\d+)\)$', key)
        if not match:
            continue
        particleIndex = int(match.group(1))
        coverage = values.get("coverage")
        score = values.get("score")
        alignMatrix = values.get("xform.align3d", {}).get("matrix")

        if emantomo.Plugin.isVersion(emantomo.constants.V_CB):
            alignMatrix = literal_eval(alignMatrix)

        if coverage and score and alignMatrix:
            customParticleParams = dict(
                coverage=coverage,
                score=score,
                alignMatrix=alignMatrix
            )
            output[particleIndex] = customParticleParams

    return output


def updateSetOfSubTomograms(inputSetOfSubTomograms, outputSetOfSubTomograms, particlesParams):
    """Update a set of subtomograms from a template and copy attributes coverage/score/transform"""

    def updateSubTomogram(subTomogram, index):
        particleParams = particlesParams.get(index)
        if not particleParams:
            print("Could not get params for particle %d" % index)
            setattr(subTomogram, "_appendItem", False)
        else:
            setattr(subTomogram, EMAN_COVERAGE, Float(particleParams["coverage"]))
            setattr(subTomogram, EMAN_SCORE, Float(particleParams["score"]))
            # Create 4x4 matrix from 4x3 e2spt_sgd align matrix and append row [0,0,0,1]
            am = numpy.array(particleParams["alignMatrix"])
            # angles = numpy.array([am[0:3], am[4:7], am[8:11], [0, 0, 0]])
            # shift = numpy.array([am[3], am[7], am[11], 1])
            # matrix = numpy.column_stack((angles, shift.T))
            homogeneous = numpy.array([0, 0, 0, 1])
            matrix = numpy.row_stack((am.reshape(3, 4), homogeneous))
            subTomogram.setTransform(Transform(matrix), convention=TR_EMAN)

    outputSetOfSubTomograms.copyItems(inputSetOfSubTomograms,
                                      updateItemCallback=updateSubTomogram,
                                      itemDataIterator=itertools.count(0))


def jsonFilesFromSet(setScipion, path):
    json_files = []
    if isinstance(setScipion, SetOfTomograms):
        tomo_files = []
        for file in setScipion.getFiles():
            fileBasename = pwutils.removeBaseExt(file)
            if "__" in fileBasename:
                fnInputCoor = '%s_info.json' % fileBasename.split("__")[0]
            else:
                parentFolder = pwutils.removeBaseExt(os.path.dirname(file))
                fnInputCoor = '%s-%s_info.json' % (parentFolder, fileBasename)
            pathInputCoor = pwutils.join(path, fnInputCoor)
            json_files.append(pathInputCoor)
            tomo_files.append(file)
        return json_files, tomo_files
    elif isinstance(setScipion, SetOfTiltSeries):
        tlt_files = []
        for tilt_serie in setScipion.iterItems(iterate=False):
            json_file = os.path.join(path,
                                     os.path.basename(os.path.dirname(tilt_serie.getFirstItem().getFileName())) +
                                     '-' + tilt_serie.getTsId() + '_info.json')
            json_files.append(json_file)
            tlt_files.append(tilt_serie.getFirstItem().getFileName())
        return json_files, tlt_files


def setCoords3D2Jsons(json_files, setCoords, mode="w"):
    paths = []
    groupIds = setCoords.aggregate(["MAX", "COUNT"], "_groupId", ["_groupId"])
    groupIds = set([d['_groupId'] for d in groupIds])
    emanIds = list(range(len(groupIds)))
    dict_eman = dict(zip(groupIds, emanIds))
    for json_file in json_files:
        coords = []
        for coor in setCoords.iterCoordinates():
            tomoName = pwutils.removeBaseExt(coor.getVolume().getFileName())
            if "__" in tomoName:
                tomoName = '%s_info' % tomoName.split("__")[0]
            else:
                tomoName += "_info"
            if tomoName in json_file:
                coords.append([coor.getX(const.BOTTOM_LEFT_CORNER),
                               coor.getY(const.BOTTOM_LEFT_CORNER),
                               coor.getZ(const.BOTTOM_LEFT_CORNER),
                               "manual", 0.0, dict_eman[coor.getGroupId()]])

        if coords:
            coordDict = {"boxes_3d": coords,
                         "class_list": {}
                         }
            for groupId in groupIds:
                coordDict["class_list"]["%s" % dict_eman[groupId]] = {"boxsize": setCoords.getBoxSize(),
                                                                      "name": "particles_%02d" % dict_eman[groupId]}
            if mode == "w":
                writeJson(coordDict, json_file)
                paths.append(json_file)
            elif mode == "a":
                appendJson(coordDict, json_file)
                paths.append(json_file)
    return paths


def jsons2SetCoords3D(protocol, setTomograms, outPath):
    from tomo.objects import SetOfCoordinates3D
    if isinstance(setTomograms, Pointer):
        setTomograms = setTomograms.get()
    coord3DSetDict = {}

    # Subsets do not have this
    if hasattr(protocol, "_getOutputSuffix"):
        suffix = protocol._getOutputSuffix(SetOfCoordinates3D)
        outputname = protocol.OUTPUT_PREFIX + suffix
    else:
        for count in range(100):
            outputname = "coordinates%s" % count
            if not hasattr(protocol, outputname):
                suffix = "user%s" % count
                break

    coord3DSet = SetOfCoordinates3D.create(protocol._getPath(),
                                           prefix=outputname)
    coord3DSet.setPrecedents(setTomograms)
    coord3DSet.setSamplingRate(setTomograms.getSamplingRate())
    first = True
    for tomo in setTomograms.iterItems():
        outFile = '*%s_info.json' % pwutils.removeBaseExt(tomo.getFileName().split("__")[0])
        pattern = os.path.join(outPath, outFile)
        files = glob.glob(pattern)

        if not files or not os.path.isfile(files[0]):
            continue

        jsonFnbase = files[0]
        jsonBoxDict = loadJson(jsonFnbase)

        if first:
            coord3DSet.setBoxSize(int(jsonBoxDict["class_list"]["0"]["boxsize"]))
            first = False

        for key in list(jsonBoxDict["class_list"].keys()):
            coord3DSetDict[int(key)] = coord3DSet

        # Populate Set of 3D Coordinates with 3D Coordinates
        readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, tomo.clone())

    protocol._defineOutputs(**{outputname: coord3DSet})
    protocol._defineSourceRelation(setTomograms, coord3DSet)

    # # Update Outputs
    # for index, coord3DSet in coord3DSetDict.items():
    #     protocol._updateOutputSet(name, coord3DSet, state=coord3DSet.STREAM_CLOSED)


def tltParams2Json(json_files, tltSeries, mode="w"):
    paths = []
    sr = tltSeries.getSamplingRate()
    for idj, json_file in enumerate(json_files):
        tilt_serie = tltSeries[idj + 1]
        tlt_params = []
        for idx, tiltImage in enumerate(tilt_serie.iterItems()):
            paths.append(os.path.abspath(tiltImage.getFileName()))
            tr_matrix = tiltImage.getTransform().getMatrix() if tiltImage.getTransform() is not None else numpy.eye(3)
            a1 = numpy.rad2deg(numpy.arccos(tr_matrix[0, 0]))
            a2 = tiltImage.getTiltAngle()
            a3 = tiltImage.tiltAngleAxis.get() if hasattr(tiltImage, 'tiltAngleAxis') else 0.0
            s1, s2 = tr_matrix[0, 2], tr_matrix[1, 2]
            tlt_params.append([s1, s2, a1, a2, a3])
        tlt_files = os.path.abspath(tilt_serie[1].getFileName())
        if tlt_params:
            tlt_dict = {"apix_unbin": sr,
                        "tlt_file": tlt_files,
                        "tlt_params": tlt_params
                        }
            if mode == "w":
                writeJson(tlt_dict, json_file)
                paths.append(json_file)
            elif mode == "a":
                appendJson(tlt_dict, json_file)
                paths.append(json_file)
    return paths


def ctf2Json(json_files, ctfSeries, mode='w'):
    paths = []
    aquisition = ctfSeries.getSetOfTiltSeries().getAcquisition()
    cs = aquisition.getSphericalAberration()
    voltage = aquisition.getVoltage()
    for idj, json_file in enumerate(json_files):
        ctf_serie = ctfSeries[idj + 1]
        defocus = []
        phase = []
        for idx, ctfTomo in enumerate(ctf_serie.iterItems()):
            defocus_eman = (ctfTomo.getDefocusU() + ctfTomo.getDefocusV()) / 20000.0
            phase.append(ctfTomo.getPhaseShift())
            defocus.append(defocus_eman)
        if defocus and phase:
            ctf_dict = {"cs": cs,
                        "voltage": voltage,
                        "defocus": defocus,
                        "phase": phase
                        }
            if mode == "w":
                writeJson(ctf_dict, json_file)
                paths.append(json_file)
            elif mode == "a":
                appendJson(ctf_dict, json_file)
                paths.append(json_file)
    return paths


def refinement2Json(protocol, subTomos, mode='w'):
    lst_file = protocol._getExtraPath(join('spt_00', 'input_ptcls.lst'))
    json_name = protocol._getExtraPath(join('spt_00', 'particle_parms_01.json'))
    parms_dict = {}

    count = 0
    for subTomo in subTomos.iterSubtomos():
        key = "('%s', %d)" % (os.path.abspath(lst_file), count)
        count += 1
        coverage = getattr(subTomo, EMAN_COVERAGE, Float(0.0)).get()
        score = getattr(subTomo, EMAN_SCORE, Float(-0.0)).get()
        matrix_st = subTomo.getTransform(convention=TR_EMAN).getMatrix()

        if subTomo.hasCoordinate3D():
            matrix_c = subTomo.getCoordinate3D().getMatrix(convention=TR_EMAN)
        else:
            matrix_c = np.eye(4)

        am_st, am_c = [0] * 12, [0] * 12
        am_st[0:3], am_st[4:7], am_st[8:11] = matrix_st[0, :3], matrix_st[1, :3], matrix_st[2, :3]
        am_c[0:3], am_c[4:7], am_c[8:11] = matrix_c[0, :3], matrix_c[1, :3], matrix_c[2, :3]
        am_st[3], am_st[7], am_st[11] = matrix_st[0, 3], matrix_st[1, 3], matrix_st[2, 3]
        am_c[3], am_c[7], am_c[11] = matrix_c[0, 3], matrix_c[1, 3], matrix_c[2, 3]

        if emantomo.Plugin.isVersion(emantomo.constants.V_CB):
            am_c = "[" + ",".join(str(a) for a in am_c) + "]"
            am_st = "[" + ",".join(str(a) for a in am_st) + "]"

        parms_dict[key] = {"coverage": coverage, "score": score,
                           "xform.align3d": {"__class__": "Transform",
                                             "matrix": am_st},
                           "xform.start": {"__class__": "Transform",
                                           "matrix": am_c},
                           }
    if mode == "w":
        writeJson(parms_dict, json_name, indent=1)
    elif mode == "a":
        appendJson(parms_dict, json_name, indent=1)


def recoverTSFromObj(child_obj, protocol):
    p = protocol.getProject()
    graph = p.getSourceGraph(False)
    relations = p.mapper.getRelationsByName(RELATION_SOURCE)
    n = graph.getNode(child_obj.strId())
    connection = []
    while n is not None and not n.getParent().isRoot():
        n = n.getParent()
        connection.append(n.pointer.getUniqueId())
        connection.append(n.pointer.get().strId())
    for rel in relations:
        pObj = p.getObject(rel[OBJECT_PARENT_ID])
        pExt = rel['object_parent_extended']
        pp = Pointer(pObj, extended=pExt)
        if pp.getUniqueId() in connection:
            if isinstance(pp.get(), SetOfTiltSeries) and pp.get().getFirstItem().getFirstItem().hasTransform():
                return pp.get()
    raise ValueError('Could not find any SetOfTiltSeries associated to %s.' % type(child_obj))


def emanFSCsToScipion(fscSet, *fscFiles):
    def _getFscValues(fscFn):
        resolution_inv, frc = [], []
        with open(fscFn) as f1:
            for l in f1:
                resolution_inv.append(float(l.split()[0]))
                frc.append(float(l.split()[1]))

        return resolution_inv, frc

    for fscFile in fscFiles:
        res_inv, frc = _getFscValues(fscFile)
        fsc = FSC(objLabel=os.path.basename(fscFile))
        fsc.setData(res_inv, frc)
        fscSet.append(fsc)
