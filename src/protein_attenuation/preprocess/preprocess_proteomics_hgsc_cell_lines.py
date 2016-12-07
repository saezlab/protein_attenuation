#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

from pandas import read_csv, DataFrame
from protein_attenuation.utils import gkn

# -- Import and preprocess proteomics
proteomics = read_csv('./data/hgsc_cell_lines_proteomics.csv')

proteomics = proteomics[[len(i.split(';')) == 1 for i in proteomics['Majority protein IDs']]]

proteomics = proteomics.groupby('Gene names').mean()

# Proteins measured in at least 50% of the samples
proteomics = proteomics[proteomics.count(1) > proteomics.shape[1] * .5]

# Normalise
proteomics = DataFrame({p: gkn(proteomics.ix[p].dropna()) for p in proteomics.index}).T

# -- Export
proteomics.to_csv('./data/hgsc_cell_lines_proteomics_preprocessed.csv')
print '[INFO] HGSC cell lines proteomics: ', './data/hgsc_cell_lines_proteomics_preprocessed.csv'
