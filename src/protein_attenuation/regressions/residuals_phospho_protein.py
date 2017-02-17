#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from sklearn.linear_model import LinearRegression
from pandas import DataFrame, read_csv

# -- Import
phoshoproteomics = read_csv('./data/cptac_phosphoproteomics_corrected_normalised.csv', index_col=0)
proteomics = read_csv('./data/cptac_proteomics_tmt_corrected_normalised.csv', index_col=0)


# -- Overlap
samples = set(proteomics).intersection(phoshoproteomics)

proteins = set(proteomics.index)
psites = {p for p in phoshoproteomics.index if p.split('_')[0] in proteins}


# -- Residuals
# p = 'EPS8L1_T187'
def phospho_residual(p):
    print p

    df = DataFrame({'phospho': phoshoproteomics.ix[p], 'prot': proteomics.ix[p.split('_')[0]]}).dropna()

    if df.shape[0] > (len(samples) * .5):
        y, x = df['phospho'], df[['prot']]

        lm = LinearRegression().fit(x, y)
        y_ = y - lm.coef_[0] * df['prot'] - lm.intercept_

        return y_
residuals = {p: phospho_residual(p) for p in psites}
residuals = DataFrame({p: residuals[p] for p in residuals if residuals[p] is not None}).T
residuals.to_csv('./data/residuals_phospho_protein.csv')
print '[INFO] Phosphoproteomics ~ Protein residuals: ', './data/residuals_phospho_protein.csv'
