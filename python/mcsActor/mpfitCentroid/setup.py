from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np
import os

#os.environ["CC"]="gcc"

setup(
  name = 'centroid',
  ext_modules=[
    Extension('centroid',
              sources=['centroid_wrapper.pyx',"mpfit.c","centroid_tidy.c","centroid_routines.c","fibreid_all.c","calling_wrappers.c"],
              include_dirs=[np.get_include(),
                            "libcentroid"])
    ],
  cmdclass = {'build_ext': build_ext}
)
