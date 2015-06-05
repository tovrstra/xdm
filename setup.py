#!/usr/bin/env python

import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {
        'build_ext': build_ext,
    },
    ext_modules=[
        Extension("cext",
            sources=['cext.pyx', 'brhole.cpp'],
            depends=['brhole.h', 'brhole.pxd'],
            include_dirs=[np.get_include(), '.'],
            language="c++",
        )
    ],
)
