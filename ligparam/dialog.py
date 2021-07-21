import os

from pathlib import Path

import numpy as np
from parmed.topologyobjects import BondType, AngleType, DihedralType
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, uic
from pysisyphus.calculators.OpenMM import OpenMM
from pysisyphus.constants import AU2KJPERMOL, BOHR2ANG
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT
from pysisyphus.run import run_scan, get_calc_closure
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from simtk.openmm import app


THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


class TermTable(pg.TableWidget):
    dtypes = {
        BondType: [("k", float), ("req", float)],
        AngleType: [("k", float), ("theteq", float)],
        # scee and scnb are omitted for now
        DihedralType: [("phi_k", float), ("per", int), ("phase", float)],
    }

    def __init__(self, *args, **kwargs):
        kwargs["editable"] = True
        super().__init__(*args, **kwargs)

    def set_terms(self, terms):
        self.terms = terms
        # BondType, AngleType or DihedralType
        self.type_ = type(self.terms[0])
        self.dtype = self.dtypes[self.type_]
        self.fields = [field for field, _ in self.dtype]

        data = list()
        for term in terms:
            d = [getattr(term, field) for field in self.fields]
            data.append(d)
        data = np.core.records.fromarrays(np.array(data).T, dtype=self.dtype)
        self.data_backup = data.copy()
        self.setData(data)

    def get_terms(self):
        for i, term in enumerate(self.terms):
            for j, (attr, conv) in enumerate(self.dtype):
                data = self.item(i, j).text()
                setattr(self.terms[i], attr, conv(data))


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
    def __init__(self, nodes, types, terms, indices, top, params, qm_geom, ff_geom, **kwargs):
        super().__init__(**kwargs)
        uic.loadUi(THIS_DIR / "dialog.ui", self)

        self.nodes = nodes
        self.types = types
        self.terms = terms
        self.indices = indices
        self.top = top
        self.params = params
        self.qm_geom = qm_geom
        self.ff_geom = ff_geom

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

        self.ff_scans = 0
        if self.prim_type == PT.BOND:
            step_size = 0.05
            unit = "Å"
        else:
            step_size = 5
            unit = "deg"
        self.step_size.setText(str(step_size))
        self.step_unit.setText(unit)
        self.plot.addLegend()

    def update_plot(self, vals, ens, name, **plot_kwargs):
        vals = vals.copy()
        if self.prim_type == PT.BOND:
            vals *= BOHR2ANG
        else:
            vals = np.rad2deg(vals)
        ens = ens.copy()
        ens -= ens.min()
        ens *= AU2KJPERMOL
        self.plot.plot(vals, ens, name=name, **plot_kwargs)

    def get_scan_kwargs(self):
        steps = self.steps.value()
        step_size = float(self.step_size.text())
        if self.prim_type == PT.BOND:
            step_size /= BOHR2ANG
        else:
            step_size = np.deg2rad(step_size)
        symmetric = self.symmetric.isChecked()
        return steps, step_size, symmetric

    def run_qm_scan(self):
        geom = self.qm_geom.copy()
        calc_kwargs = {
            "keywords": "ri-mp2 6-31G* cc-pvdz/C",
            "pal": 6,
            "mem": 1000,
        }
        calc_getter = get_calc_closure("ligparam_scan", "orca", calc_kwargs)
        # calc_getter = get_calc_closure("ligparam_scan", "xtb", {})

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
        geom = self.ff_geom.copy()
        self.top.coordinates = geom.coords3d * BOHR2ANG

        # Update params with current values from TableWidget
        self.term_table.get_terms()
        print("Using current terms:")
        for i, term in enumerate(self.terms):
            print(f"\t{i:02d}: {term}")

        def calc_getter(**kwargs):
            calc = OpenMM(self.top, self.params, **kwargs)
            return calc

        # Obtain minimum at the current coordinates
        geom.set_calculator(calc_getter())
        opt = RFOptimizer(geom, thresh="gau", max_cycles=150)
        opt.run()

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
        self.setStatusTip("Relaxed scan finished")
        pen = pg.mkPen((0, 255, 0))
        ff_label = self.ff_label.text()
        if ff_label == "default":
            ff_label = f"FF_{self.ff_scans}"

        # self.update_plot(vals, ens, f"FF_{self.ff_scans}", pen=pen, symbol="o")
        self.update_plot(vals, ens, ff_label, pen=pen, symbol="o")
        self.ff_scans += 1
