import os
from pathlib import Path

import numpy as np
from parmed.topologyobjects import BondType, AngleType, DihedralType
import psutil
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, uic
from pysisyphus.calculators.OpenMM import OpenMM
from pysisyphus.constants import AU2KJPERMOL, BOHR2ANG
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT
from pysisyphus.run import run_scan, get_calc_closure
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from ligparam.config import Config
from ligparam.db import select_or_insert_prim_coord, insert_scan, get_scan_data
from ligparam.helpers import log


THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))
CALCULATORS = Config["calculators"].copy()


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

        self.itemEntered.connect(self.item_entered)
        self.itemPressed.connect(self.item_pressed)
        self.itemSelectionChanged.connect(self.item_selection_changed)

    def item_selection_changed(self):
        cur_item = self.currentItem()
        self.item_text_backup = cur_item.text()
        print(f"saved {cur_item} text: {self.item_text_backup}")

    def item_entered(self):
        print("item entered")

    def item_pressed(self):
        print("item pressed")

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
    def __init__(
        self, nodes, types, terms, indices, top, params, qm_geom, ff_geom, **kwargs
    ):
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
        self.run_ff.setDefault(True)
        self.clear_ff.clicked.connect(self.clear_ff_plot)
        self.clear_all.clicked.connect(self.plot.clear)

        # self.calcs = Config["calculators"].copy()
        self.calcs = CALCULATORS.copy()
        calc_levels = list()
        for key, (type_, _) in self.calcs.items():
            calc_levels.append(type_)
            self.calc_level.addItem(str(key))
        self.calc_level.setCurrentIndex(calc_levels.index("psi4"))

        self.ff_label_default = "default"
        self.ff_label.setText(self.ff_label_default)

        prim_types = {
            2: PT.BOND,
            3: PT.BEND,
            4: PT.PROPER_DIHEDRAL,
        }
        self.prim_type = prim_types[len(nodes)]
        self.type_.setText(str(self.prim_type))
        self.prim_indices_le.setText(str(indices))
        self.typed_prim = (self.prim_type, *indices)

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
        self.ff_lines = list()

        # Make PrimCoord available
        select_or_insert_prim_coord(self.typed_prim)
        _, calculator, _ = self.get_calc_data()
        vals_, ens_ = get_scan_data(self.typed_prim, calculator)
        if vals_.size > 0:
            self.plot_qm(vals_, ens_, calculator)

        self.change_counter = 0
        self.term_table.itemChanged.connect(self.item_changed)

    def item_changed(self, item):
        row = item.row()
        col = item.column()
        field = self.term_table.fields[col]
        old_val = self.term_table.item_text_backup
        new_val = item.text()
        text = (
            f"{self.change_counter:03d}: "
            f"Updated {field} in row {row} from {old_val} to {new_val}"
        )
        self.term_history.appendPlainText(text)
        self.change_counter += 1

    def update_plot(self, vals, ens, name, **plot_kwargs):
        vals = vals.copy()
        if self.prim_type == PT.BOND:
            vals *= BOHR2ANG
        else:
            vals = np.rad2deg(vals)
        ens = ens.copy()
        ens -= ens.min()
        ens *= AU2KJPERMOL
        is_ff = plot_kwargs.pop("is_ff", False)
        line = self.plot.plot(vals, ens, name=name, **plot_kwargs)
        if is_ff:
            self.ff_lines.append(line)

    def clear_ff_plot(self):
        for line in self.ff_lines:
            log(line)
            line.clear()
        self.ff_lines = list()

    def get_scan_kwargs(self):
        steps = self.steps.value()
        step_size = float(self.step_size.text())
        if self.prim_type == PT.BOND:
            step_size /= BOHR2ANG
        else:
            step_size = np.deg2rad(step_size)
        symmetric = self.symmetric.isChecked()
        return steps, step_size, symmetric

    def get_calc_data(self):
        calc_key = self.calc_level.currentText()
        calc_type, calc_kwargs = self.calcs[calc_key]
        return calc_key, calc_type, calc_kwargs

    def run_qm_scan(self):
        geom = self.qm_geom.copy()
        calc_key, calc_type, calc_kwargs = self.get_calc_data()

        calc_kwargs = calc_kwargs.copy()
        calc_kwargs["pal"] = psutil.cpu_count(logical=False)
        calc_kwargs["mem"] = 1500
        calc_kwargs["out_dir"] = "qm_calcs"
        calc_getter = get_calc_closure("ligparam_scan", calc_type, calc_kwargs)

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
        # Insert into DB
        insert_scan(self.typed_prim, calc_type, vals, ens)
        self.plot_qm(vals, ens, calc_key)

    def plot_qm(self, vals, ens, calc_key):
        pen = pg.mkPen((255, 0, 0))
        qm_label = f"QM_{calc_key}"
        self.update_plot(vals, ens, qm_label, pen=pen, symbol="x")

    def run_ff_scan(self):
        geom = self.ff_geom.copy()
        self.top.coordinates = geom.coords3d * BOHR2ANG

        # Update params with current values from TableWidget
        self.term_table.get_terms()
        log("Using current terms:")
        for i, term in enumerate(self.terms):
            log(f"\t{i:02d}: {term}")

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
        if ff_label == self.ff_label_default:
            ff_label = f"FF_{self.ff_scans}"

        self.update_plot(vals, ens, ff_label, pen=pen, symbol="o", is_ff=True)
        self.ff_scans += 1
        self.ff_label.setText(self.ff_label_default)
