# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os

from pwem.wizards.wizard import EmWizard
# from pwem.viewers import CoordinatesObjectView
# from pyworkflow.utils import makePath, cleanPath, readProperties
#
# from . import Plugin
# from .convert import writeSetOfMicrographs
# from .protocols import SparxGaussianProtPicking, EmanProtTomoExtraction, EmanProtTomoTempMatch
from .protocols import EmanProtTomoExtraction, EmanProtTomoTempMatch


# class SparxGaussianPickerWizard(EmWizard):
#     _targets = [(SparxGaussianProtPicking, ['boxSize',
#                                             'lowerThreshold', 'higherThreshold',
#                                             'gaussWidth'])]
#
#     def show(self, form):
#         autopickProt = form.protocol
#         micSet = autopickProt.getInputMicrographs()
#         if not micSet:
#             print('You must specify input micrographs')
#             return
#
#         # ensuring a valid boxSize
#         if autopickProt.boxSize.get() is None:
#             autopickProt.boxSize.set(100)
#
#         project = autopickProt.getProject()
#         coordsDir = project.getTmpPath(micSet.getName())
#
#         cleanPath(coordsDir)
#         makePath(coordsDir)
#
#         micMdFn = os.path.join(coordsDir, "micrographs.xmd")
#         writeSetOfMicrographs(micSet, micMdFn)
#
#         pickerProps = os.path.join(coordsDir, 'picker.conf')
#         params = ['boxSize', 'lowerThreshold', 'higherThreshold', 'gaussWidth']
#         program = Plugin.getBoxerCommand(boxerVersion='old')
#
#         extraParams = "invert_contrast=%s:use_variance=%s:%s" % (
#             autopickProt.doInvert,
#             autopickProt.useVarImg,
#             autopickProt.extraParams)
#
#         args = {
#             "params": ','.join(params),
#             "preprocess": "emprogram %s" % Plugin.getProgram('sxprocess.py'),
#             "picker": "emprogram %s" % program,
#             "convert": "emconvert",
#             'coordsDir': coordsDir,
#             'micsSqlite': micSet.getFileName(),
#             "boxSize": autopickProt.boxSize,
#             "lowerThreshold": autopickProt.lowerThreshold,
#             "higherThreshold": autopickProt.higherThreshold,
#             "gaussWidth": autopickProt.gaussWidth,
#             "extraParams": extraParams
#         }
#
#         with open(pickerProps, "w") as f:
#             f.write("""
#             parameters = %(params)s
#             boxSize.value = %(boxSize)s
#             boxSize.label = Box Size
#             boxSize.help = Box size in pixels
#             lowerThreshold.value =  %(lowerThreshold)s
#             lowerThreshold.label = Lower Threshold
#             lowerThreshold.help = Lower Threshold
#             higherThreshold.help = Higher Threshold
#             higherThreshold.value =  %(higherThreshold)s
#             higherThreshold.label = Higher Threshold
#             gaussWidth.help = Width of the Gaussian kernel used
#             gaussWidth.value =  %(gaussWidth)s
#             gaussWidth.label = Gauss Width
#             runDir = %(coordsDir)s
#             preprocessCommand = %(preprocess)s demoparms --makedb=thr_low=%%(lowerThreshold):thr_hi=%%(higherThreshold):boxsize=%%(boxSize):gauss_width=%%(gaussWidth):%(extraParams)s
#             autopickCommand = %(picker)s --gauss_autoboxer=demoparms --write_dbbox --boxsize=%%(boxSize) --norm=normalize.ramp.normvar %%(micrograph)
#             convertCommand = %(convert)s --coordinates --from emantomo --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s
#             """ % args)
#         process = CoordinatesObjectView(project, micMdFn, coordsDir, autopickProt,
#                                         mode=CoordinatesObjectView.MODE_AUTOMATIC,
#                                         pickerProps=pickerProps).show()
#         process.wait()
#         myprops = readProperties(pickerProps)
#
#         if myprops['applyChanges'] == 'true':
#             for param in params:
#                 form.setVar(param, myprops[param + '.value'])


class EmanTomoExtractionWizard(EmWizard):
    _targets = [(EmanProtTomoExtraction, ['boxSize'])]

    def show(self, form):
        tomoExtractProt = form.protocol
        inputCoordinates = tomoExtractProt.inputCoordinates.get()
        if not inputCoordinates:
            print('You must specify input coordinates')
            return

        boxSize = inputCoordinates.getBoxSize()
        if not boxSize:
            print('These coordinates do not have box size. Please, enter box size manually.')
            return

        if tomoExtractProt.downFactor.get() != 1:
            boxSize = float(boxSize/tomoExtractProt.downFactor.get())

        form.setVar('boxSize', boxSize)


class EmanTomoTempMatchWizard(EmWizard):
    _targets = [(EmanProtTomoTempMatch, ['boxSize'])]

    def show(self, form):
        tomoExtractProt = form.protocol
        inputReference = tomoExtractProt.ref.get()
        if not inputReference:
            print('You must specify input reference volume')
            return

        boxSize = inputReference.getDim()[0]

        form.setVar('boxSize', boxSize)
