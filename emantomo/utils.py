# *
# * Authors:     Scipion Team
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************

from os.path import join, abspath
from emantomo.constants import TOMO_ID, TS_ID
from pyworkflow.object import Pointer
from pyworkflow.utils import removeBaseExt, getParentFolder
from tomo.objects import SetOfCoordinates3D, SetOfSubTomograms


def getPresentPrecedents(coordSet, tomoIdsList):
    """There can be precedents without coordinates or not present when matching because of, for
    example, making a subset."""
    return [tomo.clone() for tomo in coordSet.getPrecedents() if tomo.getTsId() in tomoIdsList]


def genTomoJsonFileName(tomoFileName, jsonOutDirName):
    fileBasename = removeBaseExt(tomoFileName)
    if "__" in fileBasename:
        fileBasename = fileBasename.split("__")[0]

    parentDir = removeBaseExt(getParentFolder(tomoFileName))
    return join(jsonOutDirName, '%s-%s_info.json' % (parentDir, fileBasename))


def getBoxSize(prot):
    boxSizeFromForm = getattr(prot, 'boxSize', None)
    return boxSizeFromForm.get() if boxSizeFromForm else prot.inputSubtomos.get().getCoordinates3D().getBoxSize()


def getFromPresentObjects(inSet, labelList):
    idLabel = TOMO_ID if type(inSet) in [SetOfCoordinates3D, SetOfSubTomograms] else TS_ID
    matches = inSet.aggregate(['COUNT'], idLabel, labelList)
    outDict = {}
    for label in labelList:
        uniqueValues = set([d[label] for d in matches])
        outDict[label] = list(uniqueValues)

    return outDict


def genEmanGrouping(groupIds):
    """Generate a matching dictionary between the incoming group ids and the ones expected by EMAN."""
    emanIds = list(range(len(groupIds)))
    emanDict = dict(zip(groupIds, emanIds))
    return emanDict


def genExtract2dCmd(tomoFile, boxSize, nThreads):
    return "%s --rmbeadthr=-1 --shrink=1.0 --tltkeep=1.0 --padtwod=1.0  " \
           "--curves=-1 --curves_overlap=0.5 --compressbits=-1 --boxsz_unbin=%i  " \
           "--threads=%i" \
           % (abspath(tomoFile), boxSize, nThreads)


def genJsonFileName(infoDir, tomoId):
    return join(infoDir, '%s_info.json' % tomoId)


def getObjByName(prot, name):
    """Return an object, from a protocol, named 'name' instead of a pointer."""
    obj = getattr(prot, name, None)
    if obj and type(obj) == Pointer:
        return obj.get()
    else:
        return obj


def genEmanTomoPrj(destDir):
    pass

