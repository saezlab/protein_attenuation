#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
from protein_attenuation.utils import gkn

# -- Protein map
p_map = read_csv('./data/nci60_cell_lines_proteomics_pmap.csv', index_col=0)['Symbol']
p_map = p_map[p_map != ' '].to_dict()

# -- Import and preprocess proteomics
proteomics = read_csv('./data/nci60_cell_lines_proteomics.csv')

proteomics = proteomics[[len(i.split(';')) == 1 for i in proteomics['Protein IDs']]]
proteomics = proteomics[[i in p_map for i in proteomics['Accession']]]
proteomics['gene'] = [p_map[i].split('(')[0].strip() for i in proteomics['Accession']]

proteomics = proteomics.groupby('gene').mean()
proteomics = proteomics.drop(['kda_prot', '10-Sep', '06-Sep', '05-Mar', '01-Sep'])

proteomics = np.log2(proteomics.replace(0, np.nan))

proteomics = proteomics[[i for i in proteomics if i.startswith('iBAQ')]]

# Proteins measured in at least 50% of the samples
proteomics = proteomics[proteomics.count(1) > proteomics.shape[1] * .5]

# Normalise
proteomics = DataFrame({p: gkn(proteomics.ix[p].dropna()) for p in proteomics.index}).T

# Samples names
proteomics.columns = [i.split('_')[2] for i in proteomics]

# -- Export
proteomics.to_csv('./data/nci60_cell_lines_proteomics_preprocessed.csv')
print '[INFO] Exported'
