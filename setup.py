from distutils.core import setup, Extension
import numpy as np

mag_fit = Extension('LP.mag_fit',
                     sources = ['mag_fit/python/mag_fitmodule.c', 
                                'mag_fit/mag_fit.c',
                                'mag_fit/mag_doppel2.c'])

fitfun = Extension('LP.fitfun',
                    sources = ['fitfun/python/fitfunmodule.c'],
                    libraries = ['minpack'])

setup(name = 'LP',
      version = '1.0',
      description = 'Langmuir Probe Analysis Package',
      include_dirs = [np.get_include()],
      packages = ['utils', 'tokamak', 'LP'],
      package_data = {'tokamak': ['data_aug/*']},
      ext_modules = [mag_fit, fitfun])

