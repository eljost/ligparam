import os

from pathlib import Path

import numpy as np
from parmed.topologyobjects import BondType, AngleType, DihedralType
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, uic
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT
from pysisyphus.run import run_scan, get_calc_closure
from pysisyphus.constants import AU2KJPERMOL


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


def run_scan_wrapper(geom, calc_getter, type_, indices, symmetric, steps, step_size):
    scan_kwargs = {
        "type": type_,
        "indices": indices,
        "symmetric": symmetric,
        "steps": steps,
        "step_size": step_size,
        "opt": {
            "type": "rfo",
            "thresh": "gau_loose",
        },
    }
    return run_scan(geom, calc_getter, scan_kwargs)


class TermDialog(QtGui.QDialog):
    def __init__(self, nodes, types, terms, indices, **kwargs):
        super().__init__(**kwargs)
        uic.loadUi(THIS_DIR / "dialog.ui", self)

        self.nodes = nodes
        self.types = types
        self.terms = terms
        self.indices = indices

        dbl_validator = QtGui.QDoubleValidator(self)
        self.step_size.setValidator(dbl_validator)
        node_str = "-".join(nodes)
        self.setWindowTitle(f"Edit parameters for {node_str}")
        self.term_table.set_terms(terms)
        self.run_qm.clicked.connect(self.run_qm_scan)
        self.run_ff.clicked.connect(self.run_ff_scan)

        prim_types = {
            2: PT.BOND,
            3: PT.BEND,
            4: PT.PROPER_DIHEDRAL,
        }
        self.prim_type = prim_types[len(nodes)]
        self.type_.setText(str(self.prim_type))
        self.prim_indices_le.setText(str(indices))

        self.plot.addLegend()

    def update_plot(self, vals, ens, name, **plot_kwargs):
        ens = ens.copy()
        ens -= ens.min()
        ens *= AU2KJPERMOL
        self.plot.plot(vals, ens, name=name, **plot_kwargs)

    def get_scan_kwargs(self):
        steps = self.steps.value()
        step_size = float(self.step_size.text())
        symmetric = self.symmetric.isChecked()
        return steps, step_size, symmetric

    def run_qm_scan(self):
        # geom = geom_loader("azb_mp2.crd", coord_type="redund")
        geom = geom_loader("azb_xtb_opt.xyz", coord_type="redund")

        calc_kwargs = {
            "keywords": "hf sto-3g",
            "pal": 6,
            "mem": 1000,
        }
        # calc_getter = get_calc_closure("ligparam_scan", "orca", kwargs)
        calc_getter = get_calc_closure("ligparam_scan", "xtb", {})

        steps, step_size, symmetric = self.get_scan_kwargs()
        geoms, vals, ens = run_scan_wrapper(
            geom,
            calc_getter,
            self.prim_type,
            self.indices,
            symmetric,
            steps,
            step_size,
        )
        pen = pg.mkPen((255, 0, 0))
        self.update_plot(vals, ens, "QM", pen=pen, symbol="x")

    def run_ff_scan(self):
        geom = geom_loader("azb_gfnff_opt.xyz", coord_type="redund")
        calc_getter = get_calc_closure("ligparam_scan", "xtb", {"gfn": "ff"})
        steps, step_size, symmetric = self.get_scan_kwargs()
        geoms, vals, ens = run_scan_wrapper(
            geom,
            calc_getter,
            self.prim_type,
            self.indices,
            symmetric,
            steps,
            step_size,
        )
        pen = pg.mkPen((0, 255, 0))
        self.update_plot(vals, ens, "FF", pen=pen, symbol="o")
