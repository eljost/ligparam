import sys

from parmed.charmm import CharmmParameterSet
from pyqtgraph.Qt import mkQApp
import pyqtgraph as pg

from ligparam import __version__
from ligparam.dialog import TermDialog


def test_version():
    assert __version__ == '0.1.0'


def test_dialog():
    mkQApp("dialog test")
    inp = "azb.str"
    resi = "AZB"
    params = CharmmParameterSet(inp)
    print(f"Loaded '{inp}'")
    resis = params.residues
    print(f"Found {len(resis)} residue(s)")
    for i, resi in enumerate(params.residues):
        print(f"\t{i:03d}: {resi}")
    resi = params.residues[resi]
    print(f"Chose residue '{resi}'")

    btype = ("NG2D1", "NG2D1")
    bterms = [params.bond_types[btype]]
    bnodes = ("N1", "N2")
    binds = (6, 7)

    atype = ('CG2R61', 'CG2R61', 'CG2R61')
    aterms = [params.angle_types[atype]]
    anodes = ("C5", "C4", "C3")

    dtype = ('CG2R61', 'CG2R61', 'NG2D1', 'NG2D1')
    dterms = params.dihedral_types[dtype]
    dnodes = ("C4", "C5", "N1", "N2")

    td = TermDialog(bnodes, btype, bterms, binds)
    # td = TermDialog(anodes, atype, aterms)
    # td = TermDialog(dnodes, dtype, dterms)
    td.show()
    pg.exec()
