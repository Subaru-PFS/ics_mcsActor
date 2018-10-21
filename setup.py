import os
import sdss3tools

from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

centroid = Extension('mcsActor/mpfitCentroid/centroid',
                     sources=[os.path.join('python/mcsActor/mpfitCentroid', f)
                              for f in ("centroid_wrapper.pyx",
                                        "mpfit.c",
                                        "centroid_tidy.c",
                                        "centroid_routines.c",
                                        "fibreid_all.c",
                                        "calling_wrappers.c")],
                     include_dirs=[np.get_include()])

sdss3tools.setup(
    description = "Toy SDSS-3 actor.",
    ext_modules = cythonize(centroid),
    data_dirs = ['source']
)

