import numpy as np
import peewee
from peewee import (
    CharField,
    FloatField,
    DecimalField,
    fn,
    ForeignKeyField,
    IntegerField,
    Model,
)

from ligparam.config import Config, get_db


# import logging
# logger = logging.getLogger("peewee")
# logger.addHandler(logging.StreamHandler())
# logger.setLevel(logging.DEBUG)


db = get_db()
VAL_PREC = Config["db"]["val_prec"]
VAL_THRESH = 10.0**(-VAL_PREC)


class BaseModel(Model):
    class Meta:
        database = db


class Calculator(BaseModel):

    name = CharField(unique=True)


class PrimCoord(BaseModel):

    type_ = IntegerField()


class AtomIndex(BaseModel):

    prim = ForeignKeyField(PrimCoord, backref="indices")
    ord_ = IntegerField()  # Order
    ind = IntegerField()


class ScanPoint(BaseModel):

    prim = ForeignKeyField(PrimCoord, backref="scanpoint")
    calc = ForeignKeyField(Calculator, backref="scanpoint")
    value = DecimalField(max_digits=10, decimal_places=VAL_PREC)
    energy = FloatField()

    class Meta:
        database = db
        indexes = (
            (("prim_id", "calc_id", "value"), True),
        )


def get_prim_coord(typed_prim):
    """
    Pure SQL implementation.

    This does not take into account the order of the atom indices, which
    should be a problem for cyclic systems, where multiple bends from the
    same set of atoms can be defined.

      B
     / \
    A---C

    A-B-C
    B-C-A
    C-A-B
    """
    pt, *indices = typed_prim
    pt_int = pt.value

    prim_coords = (
        PrimCoord.select(PrimCoord, AtomIndex)
        .join(AtomIndex)
        # Filter by atom indices and correct type of primitive internal
        .where(AtomIndex.ind << indices, PrimCoord.type_ == pt_int)
        .group_by(PrimCoord)
        .having(fn.COUNT(PrimCoord.id) == len(indices))
    )

    if len(prim_coords) > 1:
        raise Exception("Implement handling of index order!")
    else:
        prim_coord = prim_coords.first()

    return prim_coord


def get_prim_coord_(typed_prim):
    """
    Partial SQL implementation.
    """
    pt, *indices = typed_prim
    pcs = PrimCoord.select()
    ais = AtomIndex.select().where(AtomIndex.ind in indices).order_by(AtomIndex.ord_)
    models = peewee.prefetch(pcs, ais)
    for i, model in enumerate(models, 1):
        minds = [ai.ind for ai in model.indices]
        if minds == indices or minds == indices[::-1]:
            break
    else:
        return None
    prim_coord = PrimCoord.get(PrimCoord.id == i)
    return prim_coord


def select_or_insert_prim_coord(typed_prim):
    prim_coord = get_prim_coord(typed_prim)

    if prim_coord is None:
        pt, *indices = typed_prim
        prim_coord = PrimCoord.create(type_=pt.value)
        for j, ind in enumerate(indices):
            AtomIndex.create(prim=prim_coord, ord_=j, ind=ind)
    return prim_coord


def get_scan_data(typed_prim, calculator):
    prim_coord = get_prim_coord(typed_prim)
    calc = Calculator.get(Calculator.name == calculator)
    # scan_points = (
    # ScanPoint
    # .select()
    # .where(ScanPoint.prim == prim_coord, ScanPoint.calc == calc)
    # )
    scan_points = (
        ScanPoint.select(ScanPoint, Calculator, PrimCoord)
        .join(Calculator)
        .switch(ScanPoint)
        .join(PrimCoord)
        .where(ScanPoint.prim == prim_coord, ScanPoint.calc == calc)
    )
    values = list()
    energies = list()

    if scan_points:
        values, energies = zip(*[(sp.value, sp.energy) for sp in scan_points])
    values = np.array(values, dtype=float)
    energies = np.array(energies)
    return values, energies


def insert_scan(typed_prim, calculator, values, energies):
    calc = Calculator.get(Calculator.name == calculator)
    prim_coord = get_prim_coord(typed_prim)

    # Delete existing values
    present_scan_points = (
        ScanPoint.select()
        .where(ScanPoint.prim==prim_coord, ScanPoint.calc == calc)
    )
    to_delete = [
        psp.id for psp in present_scan_points
        if (np.abs(values - float(psp.value)) < VAL_THRESH).any()
    ]
    if to_delete:
        ScanPoint.delete().where(ScanPoint.id << to_delete).execute()

    for val, en in zip(values, energies):
        ScanPoint.create(prim=prim_coord, calc=calc, value=val, energy=en)


def init_db(db):
    # Initialize Tables
    db.create_tables([PrimCoord, AtomIndex, ScanPoint, Calculator])

    # Initialize Calculator table
    calcs = list()
    for key, (type_, _) in Config["calculators"].items():
        calcs.append(type_)
    preset_calcs = set([calc.name for calc in Calculator.select()])
    missing_calcs = set(calcs) - preset_calcs
    for calc in missing_calcs:
        Calculator.create(name=calc)


init_db(db)
