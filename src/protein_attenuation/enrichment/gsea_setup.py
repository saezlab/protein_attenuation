#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from distutils.core import setup
from Cython.Build import cythonize

setup(
  name='SS GSEA',
  ext_modules=cythonize('gsea.pyx'),
)
