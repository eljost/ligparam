import sys

from parmed.charmm import CharmmParameterSet
from PyQt5.QtWidgets import QApplication, QMainWindow

from ligparam import __version__
from ligparam.dialog import TermDialog


def test_version():
    assert __version__ == '0.1.0'


def test_dialog():
    app = QApplication(sys.argv)
    window = QMainWindow()
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
    window.show()

    btype = ('CG2R61', 'CG2R61')
    bterms = [params.bond_types[btype]]
    bnodes = ("C1", "C2")

    atype = ('CG2R61', 'CG2R61', 'CG2R61')
    aterms = [params.angle_types[atype]]
    anodes = ("C5", "C4", "C3")

    dtype = ('CG2R61', 'CG2R61', 'NG2D1', 'NG2D1')
    dterms = params.dihedral_types[dtype]
    dnodes = ("C4", "C5", "N1", "N2")

    # td = TermDialog(btype, bterms)
    # td = TermDialog(atype, aterms)
    # td = TermDialog(dnodes, dtype, dterms)
    td = TermDialog(dnodes, dtype, dterms)
    td.exec_()
    app.exec()
