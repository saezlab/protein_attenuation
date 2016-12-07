#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

from setuptools import setup

setup(
    name='protein_attenuation',
    version='1.0',
    author='Emanuel Goncalves',
    author_email='emanuel@ebi.ac.uk',
    description='Proteogenomics analysis of copy-number alterations in protein complexes in tumours ',
    license='GPLv3',
    keywords='protein copy-number tumour',
    url='https://github.com/saezlab/protein_attenuation',
    packages=['protein_attenuation'],
    install_requires=[
        'numpy>=1.11.1',
        'scipy>=0.17.1',
        'pandas>=0.18.1',
        'pydot>=1.0.2',
        'matplotlib>=1.4.3',
        'seaborn>=0.7.1',
        'sklearn>=0.18.1',
        'statsmodels>=0.6.1'
    ],
    long_description=open('README.md').read()
)
