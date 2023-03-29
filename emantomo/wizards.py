# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
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

from pwem.wizards.wizard import EmWizard
from pyworkflow.gui import showInfo
from .protocols import EmanProtTomoExtraction


class EmanTomoExtractionWizard(EmWizard):
    _targets = [(EmanProtTomoExtraction, ['boxSize'])]  #, (EmanProtTSExtraction, ['boxSize'])]

    def _getProtocol(self, form)-> EmanProtTomoExtraction:
        return form.protocol
    def show(self, form):
        tomoExtractProt = self._getProtocol(form)
        inputCoordinates = tomoExtractProt._getSetOfCoordinates()
        if not inputCoordinates:
            showInfo('Box wizard', 'You must specify input coordinates', form.root)
            return

        boxSize = inputCoordinates.getBoxSize()
        if not boxSize:
            showInfo('Box wizard', 'These coordinates do not have box size. Please, enter box size manually.', form.root)
            return

        # Picking tomograms from coordinates are different than the introduced tomogrmas for the extraction
        inputTomos = tomoExtractProt.getInputTomograms()
        samplingRateCoord = inputCoordinates.getSamplingRate()
        samplingRateTomo = inputTomos.getSamplingRate()
        scale = samplingRateCoord / samplingRateTomo
        boxSize = float(boxSize * scale)

        form.setVar('boxSize', boxSize)


# class EmanTomoTempMatchWizard(EmWizard):
#     _targets = [(EmanProtTomoTempMatch, ['boxSize'])]
#
#     def show(self, form):
#         tomoExtractProt = form.protocol
#         inputReference = tomoExtractProt.ref.get()
#         if not inputReference:
#             print('You must specify input reference volume')
#             return
#
#         boxSize = inputReference.getDim()[0]
#
#         form.setVar('boxSize', boxSize)


# class EmanTomoResizeWizard(EmWizard):
#     _targets = [(EmanProtTomoResize, ['xDim', 'yDim', 'zDim'])]
#
#     def show(self, form):
#         tomoResizeProt = form.protocol
#         inputTomos = tomoResizeProt.inputTomograms.get()
#         if not inputTomos:
#             print('You must specify input tomograms first')
#             return
#
#         dim = inputTomos.getFirstItem().getDim()
#
#         form.setVar('xDim', dim[0])
#         form.setVar('yDim', dim[1])
#         form.setVar('zDim', dim[2])
