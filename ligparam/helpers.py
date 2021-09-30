import logging
from pathlib import Path
import re

import pyqtgraph as pg


__logger = logging.getLogger("ligparam")


def log(msg="", level=logging.INFO):
    __logger.log(level, msg)


def dbg(msg):
    log(msg, level=logging.DEBUG)


def color(clr):
    return pg.mkBrush(clr)


def inc_fn(fn, pattern):
    fn = str(fn)
    num_pattern = r"(\d+)" + pattern
    regex = re.compile(num_pattern)
    mobj = regex.search(fn)
    try:
        counter = int(mobj.group(1))
        sub = True
    except AttributeError:
        counter = 0
        sub = False

    fmt = "03d"
    while True:
        if sub:
            fn_inc = regex.sub(f"{counter:{fmt}}{pattern}", fn)
        else:
            fn_inc = fn + f".{counter:{fmt}}"
        counter += 1
        if not Path(fn_inc).exists():
            break
    return fn_inc
