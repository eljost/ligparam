import argparse
import pickle
import sys

import networkx as nx
import numpy as np
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from pysisyphus.helpers import geom_loader
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, mkQApp

from ligparam.dialog import TermDialog
from ligparam.helpers import log


SHIFT = QtCore.Qt.Key_Shift
SHIFT_MOD = QtCore.Qt.ShiftModifier
PG_DOWN = QtCore.Qt.Key_PageDown
PG_UP = QtCore.Qt.Key_PageUp
BLUE_GREY = pg.mkBrush((101, 120, 149))


class Graph(pg.GraphItem):
    sigKeyRelease = QtCore.pyqtSignal(object)

    def __init__(self):
        self.dragPoint = None
        self.dragOffset = None
        self.textItems = []
        pg.GraphItem.__init__(self)
        self.scatter.sigClicked.connect(self.clicked)
        self.sigKeyRelease.connect(self.keyReleaseEvent)

        self.gather = False
        self.stack = list()
        self._data_list = None
        self.text_counter = 0

    def setData(self, **kwds):
        self.text = kwds.pop("text", [])
        self.data = kwds
        if "pos" in self.data:
            self.npts = self.data["pos"].shape[0]
            self.data["data"] = np.empty(self.npts, dtype=[("index", int)])
            self.data["data"]["index"] = np.arange(self.npts)
        self.setTexts(self.text)
        self.updateGraph()

    def get_text(self):
        return self.texts[self.text_counter]

    def updateData(self, **kwargs):
        self.data.update(kwargs)

    def setTexts(self, text):
        for i in self.textItems:
            i.scene().removeItem(i)
        self.textItems = []
        for t in text:
            item = pg.TextItem(t)
            self.textItems.append(item)
            item.setParentItem(self)

    def updateGraph(self):
        pg.GraphItem.setData(self, **self.data)
        for i, item in enumerate(self.textItems):
            item.setPos(*self.data["pos"][i])

    def get_symbol_brushes(self, selected):
        symbol_brushes = [BLUE_GREY] * self.npts
        for i in selected:
            symbol_brushes[i] = pg.mkBrush(color=[255, 0, 0])
        return symbol_brushes

    def clicked(self, scatter, pts):
        data_list = self.scatter.data.tolist()
        if self.gather:
            point = [tup for tup in data_list if pts[0] in tup][0]
            index = data_list.index(point)
            if (index not in self.stack) and len(self.stack) < 4:
                self.stack.append(index)
                symbol_brushes = self.get_symbol_brushes(self.stack)
                self.updateData(symbolBrush=symbol_brushes)
                self.updateGraph()

    def mousePressEvent(self, event):
        mod = event.modifiers()
        self.gather = mod == SHIFT_MOD
        super().mousePressEvent(event)

    def keyReleaseEvent(self, event):
        key = event.key()
        updated = key in (PG_DOWN, PG_UP)
        if event.key() == PG_DOWN:
            self.text_counter -= 1
        elif event.key() == PG_UP:
            self.text_counter += 1

        if updated:
            self.text_counter %= len(self.texts)
            self.setTexts(self.get_text())
            self.updateGraph()


def get_graph(fn):
    with open(fn, "rb") as handle:
        bonds = pickle.load(handle)
    G = nx.from_dict_of_lists(bonds)
    return G, bonds


class Main(pg.GraphicsLayoutWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Parameter Editor GUI")
        self.show()

        self.vb = self.addViewBox()
        self.vb.setAspectLocked()
        self.graph = Graph()
        self.vb.addItem(self.graph)

    def set_charmm_resi(self, resi, top, params):
        self.resi = resi
        self.top = top
        self.params = params

        G = self.resi.to_networkx()
        nodes = [n for n in G.nodes]
        edges = [e[:-1] for e in nx.to_edgelist(G)]
        pos = nx.kamada_kawai_layout(G)
        pos_arr = np.array([pos[n] for n in nodes])
        self.node_map = {node: i for i, node in enumerate(nodes)}
        self.inv_node_map = {i: node for node, i in self.node_map.items()}
        atoms = self.resi.atoms
        self.type_map = {atom.name: atom.type for atom in atoms}

        # Adjaceny array
        adj = list()
        for from_, to_ in edges:
            adj.append((self.node_map[from_], self.node_map[to_]))
        adj = np.array(adj, dtype=int)

        symbols = ["o" for n in nodes]
        charges = [atom.charge for atom in atoms]
        types = [atom.type for atom in atoms]

        self.graph.texts = (
            nodes,
            [f"{charge:.4f}" for charge in charges],
            types,
        )

        symbol_brushes = [BLUE_GREY] * len(nodes)
        self.graph.setData(
            pos=pos_arr,
            adj=adj,
            size=0.1,
            symbol=symbols,
            pxMode=False,
            text=self.graph.get_text(),
            symbolBrush=symbol_brushes,
        )

    def set_geoms(self, qm_geom, ff_geom):
        self.qm_geom = qm_geom
        self.ff_geom = ff_geom

    def show_term_dialog(self, nodes):
        funcs = {
            2: self.params.bond_types,
            3: self.params.angle_types,
            4: self.params.dihedral_types,
        }
        types = tuple([self.type_map[node] for node in nodes])
        func = funcs[len(types)]
        try:
            terms = func[types]
        except KeyError:
            terms = func[types[::-1]]

        log(types)
        if not isinstance(terms, list):
            terms = [
                terms,
            ]
        for i, term in enumerate(terms):
            log(f"\t{i:02d}: {term}")

        indices = [self.node_map[node] for node in nodes]
        self.td = TermDialog(
            nodes,
            types,
            terms,
            indices,
            self.top,
            self.params,
            self.qm_geom,
            self.ff_geom,
        )
        self.td.exec_()
        log()

    def keyReleaseEvent(self, event):
        super().keyPressEvent(event)
        self.graph.sigKeyRelease.emit(event)

        stack = self.graph.stack
        if event.key() == SHIFT:
            self.graph.stack = list()
            if len(stack) > 1:
                nodes = [self.inv_node_map[i] for i in stack]
                self.show_term_dialog(nodes)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("str")
    parser.add_argument("psf")
    parser.add_argument("resi")
    parser.add_argument("qm_geom_fn")
    parser.add_argument("--ff_geom_fn", default=None)

    return parser.parse_args(args)


def dump_params(top, params, prm_fn):
    top.load_parameters(params)
    resi_params = params.from_structure(top)
    resi_params.write(par=prm_fn)


def run():
    args = parse_args(sys.argv[1:])
    str_ = args.str
    psf = args.psf
    resi_ = args.resi
    qm_geom_fn = args.qm_geom_fn

    qm_geom = geom_loader(qm_geom_fn, coord_type="redund")
    log(f"Loaded QM geometry from '{qm_geom_fn}'")
    if args.ff_geom_fn is None:
        ff_geom_fn = qm_geom_fn
    ff_geom = geom_loader(ff_geom_fn, coord_type="redund")
    log(f"Loaded FF geometry from '{ff_geom_fn}'")

    param_fns = (
        "par_all36_cgenff.prm",
        "top_all36_cgenff.rtf",
        str_,
    )
    log("Loading files using ParmEd:")
    for i, fn in enumerate(param_fns):
        log(f"\t{i:02d}: {fn}")
    params = CharmmParameterSet(*param_fns)
    resis = params.residues
    log(f"Found {len(resis)} residue(s)")
    resi = params.residues[resi_]
    log(f"Chose residue '{resi}'")
    top = CharmmPsfFile(psf)
    log(f"Loaded {psf}")

    mkQApp("Parameter Editor")
    main = Main()
    main.set_charmm_resi(resi, top, params)
    main.set_geoms(qm_geom, ff_geom)
    try:
        pg.exec()
    except Exception as err:
        log(err)
    prm_fn = f"{resi_.lower()}_optimized.prm"
    dump_params(top, params, prm_fn)
    log(f"Dumped optimized parameters to '{prm_fn}'")


if __name__ == "__main__":
    run()
