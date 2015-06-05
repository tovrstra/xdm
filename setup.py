#!/usr/bin/env python

import numpy as np
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name='xdm',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    package_dir = {'xdm': 'xdm'},
    packages=['xdm'],
    scripts=glob("scripts/*.py"),
    cmdclass = {
        'build_ext': build_ext,
    },
    package_data={
        'xdm': ['*.pxd'],
    },
    headers=glob('xdm/*.h'),
    ext_modules=[
        Extension("xdm.cext",
            sources=['xdm/cext.pyx', 'xdm/brhole.cpp'],
            depends=['xdm/brhole.h', 'xdm/brhole.pxd'],
            include_dirs=[np.get_include(), '.'],
            language="c++",
        )
    ],
)
