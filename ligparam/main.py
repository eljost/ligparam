import argparse
import pickle
import sys

import networkx as nx
import numpy as np
from parmed.charmm import CharmmParameterSet
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, mkQApp

from ligparam.dialog import TermDialog


SHIFT = QtCore.Qt.Key_Shift
SHIFT_MOD = QtCore.Qt.ShiftModifier
BLUE_GREY = pg.mkBrush((101, 120, 149))


class Graph(pg.GraphItem):
    def __init__(self):
        self.dragPoint = None
        self.dragOffset = None
        self.textItems = []
        pg.GraphItem.__init__(self)
        self.scatter.sigClicked.connect(self.clicked)

        self.gather = False
        self.stack = list()
        self._data_list = None

    def setData(self, **kwds):
        self.text = kwds.pop("text", [])
        self.data = kwds
        if "pos" in self.data:
            self.npts = self.data["pos"].shape[0]
            self.data["data"] = np.empty(self.npts, dtype=[("index", int)])
            self.data["data"]["index"] = np.arange(self.npts)
        self.setTexts(self.text)
        self.updateGraph()

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

    def clicked(self, scatter, pts):
        data_list = self.scatter.data.tolist()
        if self.gather:
            point = [tup for tup in data_list if pts[0] in tup][0]
            index = data_list.index(point)
            if (index not in self.stack) and len(self.stack) < 4:
                self.stack.append(index)

                symbol_brushes = [BLUE_GREY] * self.npts
                for i in self.stack:
                    symbol_brushes[i] = pg.mkBrush(color=[255, 0, 0])
                self.updateData(symbolBrush=symbol_brushes)
                self.updateGraph()

    def mousePressEvent(self, event):
        mod = event.modifiers()
        self.gather = mod == SHIFT_MOD
        super().mousePressEvent(event)


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

    def set_charmm_resi(self, resi, params):
        G = resi.to_networkx()
        nodes = [n for n in G.nodes]
        edges = [e[:-1] for e in nx.to_edgelist(G)]
        pos = nx.kamada_kawai_layout(G)
        pos_arr = np.array([pos[n] for n in nodes])
        self.node_map = {node: i for i, node in enumerate(nodes)}
        self.inv_node_map = {i: node for node, i in self.node_map.items()}
        self.type_map = {atom.name: atom.type for atom in resi.atoms}
        self.params = params

        # Adjaceny array
        adj = list()
        for from_, to_ in edges:
            adj.append((self.node_map[from_], self.node_map[to_]))
        adj = np.array(adj, dtype=int)

        symbols = ["o" for n in nodes]

        symbol_brushes = [BLUE_GREY] * len(nodes)
        self.graph.setData(
            pos=pos_arr,
            adj=adj,
            size=0.1,
            symbol=symbols,
            pxMode=False,
            text=nodes,
            symbolBrush=symbol_brushes,
        )

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

        print(types)
        if not isinstance(terms, list):
            terms = [
                terms,
            ]
        for i, term in enumerate(terms):
            print(f"\t{i:02d}: {term}")

        self.td = TermDialog(nodes, types, terms)
        self.td.exec_()
        print()

    def keyReleaseEvent(self, event):
        stack = self.graph.stack

        if event.key() == SHIFT:
            self.graph.stack = list()
            print(f"Shift released: stack={stack}")
            if len(stack) > 1:
                nodes = [self.inv_node_map[i] for i in stack]
                self.show_term_dialog(nodes)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("inp")
    parser.add_argument("resi")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    inp = args.inp
    resi_ = args.resi

    params = CharmmParameterSet(inp)
    print(f"Loaded '{inp}'")
    resis = params.residues
    print(f"Found {len(resis)} residue(s)")
    for i, resi in enumerate(params.residues):
        print(f"\t{i:03d}: {resi}")
    resi = params.residues[resi_]
    print(f"Chose residue '{resi}'")

    mkQApp("Parameter Editor")
    main = Main()
    main.set_charmm_resi(resi, params)
    pg.exec()


if __name__ == "__main__":
    run()
