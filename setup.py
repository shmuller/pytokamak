from distutils.core import setup, Extension
import os

is_pyfile = lambda s: s.endswith('.py') and not s.startswith('_')

py_modules = ['LP.' + s[:-3] for s in os.listdir('LP') if is_pyfile(s)]

mag_fit = Extension('LP.mag_fit',
                     sources = ['mag_fit/python/mag_fitmodule.c', 'mag_fit/mag_fit.c'])

setup(name = 'LP',
      version = '1.0',
      description = 'Langmuir Probe Analysis Package',
      py_modules = py_modules,
      ext_modules = [mag_fit])

