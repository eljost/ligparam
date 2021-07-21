import sys

from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from pyqtgraph.Qt import mkQApp
import pyqtgraph as pg
from pysisyphus.helpers import geom_loader

from ligparam import __version__
from ligparam.dialog import TermDialog


def test_dialog():
    mkQApp("dialog test")
    inp = "azb.str"
    qm_geom = geom_loader("azb_mp2.crd", coord_type="redund")
    ff_geom = qm_geom.copy()
    param_fns = (
        "par_all36_cgenff.prm",
        "top_all36_cgenff.rtf",
        inp,
    )
    resi = "AZB"
    params = CharmmParameterSet(*param_fns)
    print(f"Loaded '{inp}'")
    resis = params.residues
    print(f"Found {len(resis)} residue(s)")
    resi = params.residues[resi]
    print(f"Chose residue '{resi}'")
    top = CharmmPsfFile("azb.psf")

    btype = ("NG2D1", "NG2D1")
    bterms = [params.bond_types[btype]]
    bnodes = ("N1", "N2")
    binds = (6, 7)

    atype = ("CG2R61", "CG2R61", "CG2R61")
    aterms = [params.angle_types[atype]]
    anodes = ("C5", "C4", "C3")

    dtype = ("CG2R61", "CG2R61", "NG2D1", "NG2D1")
    dterms = params.dihedral_types[dtype]
    dnodes = ("C4", "C5", "N1", "N2")

    td = TermDialog(bnodes, btype, bterms, binds, top, params, qm_geom, ff_geom)
    # td = TermDialog(anodes, atype, aterms)
    # td = TermDialog(dnodes, dtype, dterms)
    td.show()
    pg.exec()
