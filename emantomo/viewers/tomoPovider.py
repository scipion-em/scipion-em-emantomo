from os.path import basename

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
