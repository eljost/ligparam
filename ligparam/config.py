import os
from pathlib import Path

from peewee import SqliteDatabase
from xdg import xdg_config_home
import yaml


def get_config_dir():
    conf_home = xdg_config_home()
    conf_dir = conf_home / "ligparam"

    if not conf_dir.exists():
        os.mkdir(conf_dir)

    return conf_dir


def get_config(path=None):
    if path is None:
        conf_dir = get_config_dir()
        path = conf_dir / "config.yaml"

    with open(path) as handle:
        conf = yaml.load(handle, Loader=yaml.SafeLoader)

    return conf


def get_db(fn=None):
    if fn is None:
        fn = Config["db"]["fn"]
    db = SqliteDatabase(fn, pragmas={"foreign_keys": 1})
    return db


TOPPAR_FNS = (
    "top_all36_prot.rtf",
    "par_all36m_prot.prm",
    "top_all36_carb.rtf",
    "par_all36_carb.prm",
    "top_all36_na.rtf",
    "par_all36_na.prm",
    "toppar_water_ions.str",
    "top_all36_cgenff.rtf",
    "par_all36_cgenff.prm",
)


def get_toppar():
    # First, try to load from TOPPAR variable
    try:
        toppar_path = Path(os.environ["TOPPAR"])
        param_fns = [
            str(tp_fn) for fn in TOPPAR_FNS if (tp_fn := toppar_path / fn).exists()
        ]
    except KeyError:
        param_fns = Config["toppar"]
    return param_fns


Config = get_config()
