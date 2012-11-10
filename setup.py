from distutils.core import setup
import os

is_pyfile = lambda s: s.endswith('.py') and not s.startswith('_')

py_modules = ['LP.' + s[:-3] for s in os.listdir('LP') if is_pyfile(s)]

setup(name = 'LP',
      version = '1.0',
      description = 'Langmuir Probe Analysis Package',
      py_modules = py_modules)

