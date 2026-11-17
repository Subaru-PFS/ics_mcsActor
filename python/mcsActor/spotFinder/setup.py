#cython: language_level=3
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

setup(
    name='spotFinder',
    ext_modules=[
        Extension('spotFinder',
                  sources=['spotFinder_wrapper.pyx', 'spotFinder.c'],
                  include_dirs=[np.get_include()],
                  libraries=['m', 'pthread'])
    ],
    cmdclass={'build_ext': build_ext}
)
