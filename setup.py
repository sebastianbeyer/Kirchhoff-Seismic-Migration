from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

ext = Extension("mod_cyMigrate", ["mod_cyMigrate.pyx"],
    include_dirs = [np.get_include()])
                
setup(ext_modules=[ext],
      cmdclass = {'build_ext': build_ext})
