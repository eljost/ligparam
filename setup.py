#!/usr/bin/env python

import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name="ligparam",
        version="0.0.1",
        # networkx
        # pysisyphus
        # numpy
        # pyqtgraph
        # pyqt5
        # peewee
        # ParmEd
        entry_points={
            "console_scripts": [
                "ligparam = ligparam.main:run",
            ]
        },
    )
