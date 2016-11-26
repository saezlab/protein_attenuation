#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas.stats.misc import zscore
from pandas import DataFrame, read_csv, concat, Series, pivot_table

# -- Achilles_v3.3.8
ess = read_csv('./data/Achilles_v3.3.8.Gs.txt', sep='\t').drop('Description', axis=1)
ess['gene'] = [i.split('_')[0] for i in ess['Name']]
ess = ess.set_index(['gene', 'Name'])
print ess.shape

# # -- Achilles_v2.4.3
# ess = read_csv('./data/Achilles_QC_v2.4.3.rnai.Gs.txt', sep='\t').drop('Description', axis=1)
# ess['gene'] = [i.split('_')[0] for i in ess['Name']]
# ess = ess.set_index(['gene', 'Name'])
# print ess.shape

# -- Attenuated proteins
cors = read_csv('./tables/proteins_correlations.csv', index_col=0)
att_proteins = set(cors[(cors['cluster'] == 1) & (cors['cluster'] == 1)].index).intersection(ess.index.levels[0])


# --
ess.ix[att_proteins].mean().sort_values()