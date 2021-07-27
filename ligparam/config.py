import os

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


def get_toppar():
    return Config["toppar"]


Config = get_config()
