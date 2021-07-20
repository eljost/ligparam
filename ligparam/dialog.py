import os

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, uic
from parmed.topologyobjects import BondType, AngleType, DihedralType


THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


class TermTable(pg.TableWidget):
    dtypes = {
        BondType: [("k", float), ("req", float)],
        AngleType: [("k", float), ("theteq", float)],
        # scee and scnb are omitted for now
        DihedralType: [("phi_k", float), ("per", int), ("phase", float)],
    }

    def set_terms(self, terms):
        self.terms = terms
        type_ = type(self.terms[0])
        dtype = self.dtypes[type_]
        fields = [field for field, _ in dtype]

        data = list()
        for term in terms:
            d = [getattr(term, field) for field in fields]
            data.append(d)
        data = np.core.records.fromarrays(np.array(data).T, dtype=dtype)
        self.setData(data)


class TermDialog(QtGui.QDialog):
    def __init__(self, nodes, types, terms, *args, **kwargs):
        super().__init__(*args, **kwargs)
        uic.loadUi(THIS_DIR / "dialog.ui", self)

        self.nodes = nodes
        self.types = types
        self.term_table.set_terms(terms)
