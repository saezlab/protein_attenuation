#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr
from pandas.stats.misc import zscore
from pandas import DataFrame, read_csv, concat, Series, pivot_table
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# # -- Achilles_v3.3.8
# ess = read_csv('./data/Achilles_v3.3.8.Gs.txt', sep='\t').drop('Description', axis=1).dropna()
# ess['gene'] = [i.split('_')[0] for i in ess['Name']]
# ess = ess.set_index(['gene', 'Name'])
# print ess.shape

# -- Achilles_v2.4.3
ess = read_csv('./data/Achilles_QC_v2.4.3.rnai.Gs.txt', sep='\t').drop('Description', axis=1)
ess['gene'] = [i.split('_')[0] for i in ess['Name']]
ess = ess.set_index(['gene', 'Name'])
print ess.shape


# -- Protein complexes
uniprot = read_uniprot_genename()

corum = {(uniprot[px][0], uniprot[py][0]) for px, py in get_complexes_pairs() if px in uniprot and py in uniprot}
corum = {(px, py) for px, py in corum if px in ess.index.levels[0] and py in ess.index.levels[0]}

c_proteins = {p for px, py in corum for p in [px, py]}
print 'corum', len(corum)


# --
ess_corr = ess.ix[c_proteins].T.corr(method='pearson')
ess_corr.values[np.tril_indices(ess_corr.shape[0], 0)] = np.nan
ess_corr = ess_corr.unstack().dropna()
print ess_corr.shape
