from distutils.core import setup, Extension
import numpy as np
import os

def is_pyfile(s):
    return s.endswith('.py') and not s.startswith('_')

def pyfiles(subdir):
    return [subdir + '.' + s[:-3] for s in os.listdir(subdir) if is_pyfile(s)]

py_modules = pyfiles('LP') + pyfiles('tokamak') + pyfiles('utils')

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
      py_modules = py_modules,
      ext_modules = [mag_fit, fitfun])

