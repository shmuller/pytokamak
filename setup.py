from distutils.core import setup, Extension
import numpy as np

mag_fit = Extension('pytokamak.LP.mag_fit',
                     sources = ['mag_fit/python/mag_fitmodule.c', 
                                'mag_fit/mag_fit.c',
                                'mag_fit/mag_doppel2.c'])

fitfun = Extension('pytokamak.LP.fitfun',
                    sources = ['fitfun/python/fitfunmodule.c'],
                    libraries = ['minpack'])

setup(name = 'pytokamak',
      version = '1.0',
      description = 'Tokamak Data Analysis Package',
      author = 'Dr. Stefan H. Muller',
      include_dirs = [np.get_include()],
      packages = ['pytokamak.utils', 'pytokamak.tokamak', 'pytokamak.LP'],
      package_data = {'pytokamak.tokamak': ['data_aug/*']},
      py_modules = ['pytokamak.__init__'],
      ext_modules = [mag_fit, fitfun])

