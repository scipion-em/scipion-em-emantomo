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

from pyworkflow.utils import removeBaseExt
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


class EmantomoTomoProvider(TomogramsTreeProvider):
    
    def __init__(self, tomoList, path, mode):
        super().__init__(tomoList, path, mode)

    def getObjectInfo(self, tomo):
        tomogramName = removeBaseExt(tomo)
        if tomo.count == 0:
            return {'key': tomogramName, 'parent': None,
                    'text': tomogramName, 'values': (tomo.count, "TODO"),
                    'tags': "pending"}
        else:
            return {'key': tomogramName, 'parent': None,
                    'text': tomogramName, 'values': (tomo.count, "DONE"),
                    'tags': "done"}
