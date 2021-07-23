import argparse
import sys

import networkx as nx
import numpy as np
from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from pysisyphus.helpers import geom_loader
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, mkQApp

from ligparam.dialog import TermDialog
from ligparam.helpers import log
from ligparam.Graph import Graph


SHIFT = QtCore.Qt.Key_Shift


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

        self.graph.setData(
            pos=pos_arr,
            adj=adj,
            size=0.1,
            symbol=symbols,
            pxMode=False,
            text=self.graph.get_text(),
            atoms=self.qm_geom.atoms,
        )
        self.graph.updateGraph()

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
        log("Term(s) after closing dialog:")
        for i, term in enumerate(terms):
            log(f"\t@@@ {i}: {' '.join([type_ for type_ in types])} {term}")

    def keyReleaseEvent(self, event):
        super().keyPressEvent(event)
        self.graph.sigKeyRelease.emit(event)

        stack = self.graph.stack
        if event.key() == SHIFT:
            self.graph.stack = list()
            print("updated graph")
            self.graph.updateGraph()
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
    main.set_geoms(qm_geom, ff_geom)
    main.set_charmm_resi(resi, top, params)
    try:
        pg.exec()
    except Exception as err:
        log(err)
    prm_fn = f"{resi_.lower()}_optimized.prm"
    dump_params(top, params, prm_fn)
    log(f"Dumped optimized parameters to '{prm_fn}'")


if __name__ == "__main__":
    run()
