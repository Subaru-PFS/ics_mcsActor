#cython: language_level=3
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np
import os

# os.environ["CC"]="gcc"

setup(
    name='centroid',
    ext_modules=[
        Extension('centroid',
                  sources=['centroid_wrapper.pyx', "centroid_win.c", "newCentroid.c"],
                  include_dirs=[np.get_include()])
    ],
    cmdclass={'build_ext': build_ext}
)
