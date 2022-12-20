import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore

from ligparam.helpers import color


SHIFT_MOD = QtCore.Qt.ShiftModifier
PG_DOWN = QtCore.Qt.Key_PageDown
PG_UP = QtCore.Qt.Key_PageUp


ATOM_COLORS = {
    "h": color((255, 255, 255)),
    "c": color((133, 133, 133)),
    "o": color((255, 0, 0)),
    "n": color((0, 0, 255)),
    "s": color((255, 255, 0)),
    "cl": color((0, 255, 0)),
    "x": color((87, 72, 53)),
}
MAGENTA = color((255, 0, 255))


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
        try:
            self.data["symbol_brushes"] = self.get_symbol_brushes()
        except KeyError:
            pass
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

    def get_symbol_brushes(self):
        symbol_brushes = list()
        for atom in self.data["atoms"]:
            try:
                brush = ATOM_COLORS[atom.lower()]
            except KeyError:
                brush = ATOM_COLORS["x"]
            symbol_brushes.append(brush)
        # symbol_brushes = [ATOM_COLORS[atom.lower()] for atom in self.data["atoms"]]
        for i in self.stack:
            symbol_brushes[i] = MAGENTA
        return symbol_brushes

    def clicked(self, scatter, pts):
        data_list = self.scatter.data.tolist()
        if self.gather:
            point = [tup for tup in data_list if pts[0] in tup][0]
            index = data_list.index(point)
            if (index not in self.stack) and len(self.stack) < 4:
                self.stack.append(index)
                symbol_brushes = self.get_symbol_brushes()
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
