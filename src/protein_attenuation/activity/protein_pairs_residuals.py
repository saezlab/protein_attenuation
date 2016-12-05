#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from protein_attenuation import wd
from sklearn.linear_model import LinearRegression
from pandas import read_csv, DataFrame, concat, Series, pivot_table


# -- Correlated protein pairs
p_pairs = read_csv('%s/tables/top_correlated_protein_pairs_proteomics.csv' % wd)
p_pairs = {(s, t) for p1, p2 in p_pairs[['p1', 'p2']].values for s, t in [(p1, p2), (p2, p1)]}
print len(p_pairs)

# -- Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
proteins = set(proteomics.index)
print proteomics


# -- Estimate residuals
# p1, p2 = 'MSH2', 'MSH6'
def protein_pair_residuals(p1, p2):
    if p1 in proteins and p2 in proteins:
        samples = set(proteomics.ix[[p1, p2]].T.dropna().index)

        x, y = proteomics.ix[[p1], samples].T, proteomics.ix[p2, samples]

        lm = LinearRegression().fit(x, y)

        y_pred = lm.predict(x)

        resid = y - y_pred

        return resid

p_residuals = DataFrame({'%s_%s' % (p1, p2): protein_pair_residuals(p1, p2) for p1, p2 in p_pairs}).T
p_residuals.round(6).to_csv('%s/tables/protein_pairs_residuals.csv' % wd)
print p_residuals
