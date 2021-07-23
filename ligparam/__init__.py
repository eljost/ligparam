import logging
import sys


logger = logging.getLogger("ligparam")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler("ligparam.log", mode="w", delay=True)
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(message)s")
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
logger.addHandler(stdout_handler)
