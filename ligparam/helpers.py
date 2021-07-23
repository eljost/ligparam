import logging


__logger = logging.getLogger("ligparam")


def log(msg="", level=logging.INFO):
    __logger.log(level, msg)


def dbg(msg):
    log(msg, level=logging.DEBUG)
