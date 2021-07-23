import logging

import pyqtgraph as pg


__logger = logging.getLogger("ligparam")


def log(msg="", level=logging.INFO):
    __logger.log(level, msg)


def dbg(msg):
    log(msg, level=logging.DEBUG)


def color(clr):
    return pg.mkBrush(clr)
