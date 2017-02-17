#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from sklearn.linear_model import LinearRegression
from pandas import DataFrame, read_csv

# -- Import
transcriptomics = read_csv('./data/tcga_rnaseq_tmt_corrected_normalised.csv', index_col=0)
proteomics = read_csv('./data/cptac_proteomics_tmt_corrected_normalised.csv', index_col=0)


# -- Overlap
genes = set(proteomics.index).intersection(transcriptomics.index)
samples = set(proteomics).intersection(transcriptomics)


# -- Residuals
# g = 'ZW10'
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in genes}).T
residuals.to_csv('./data/residuals_protein_transcript_tmt.csv')
print '[INFO] Proteomics ~ Transcriptomics residuals: ', './data/residuals_protein_transcript_tmt.csv'
