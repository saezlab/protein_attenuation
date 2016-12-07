#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, Series, DataFrame
from scipy.stats.stats import ttest_ind
from sklearn.metrics.regression import mean_squared_error, mean_absolute_error
from sklearn.linear_model.base import LinearRegression


# -- Residuals
residuals = read_csv('./data/residuals_protein_transcript.csv', index_col=0)
print 'residuals', residuals.shape


# -- Copy-number
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape


# -- Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape


# -- Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- Regulatory interactions
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[ppairs_cnv['beta'] > 0]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

associations = {(px, py) for px, py in ppairs_cnv[['px', 'py']].values}
associations_px = {i[0] for i in associations}
associations_py = {i[1] for i in associations}
print len(associations)


# -- Protein correlations
samples = set(cnv).intersection(transcriptomics).intersection(proteomics)
genes = set(cnv.index).intersection(transcriptomics.index).intersection(proteomics.index)

cors = {}
for g in genes:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': proteomics.ix[g, samples]}).dropna().corr()

    cors[g] = {
        'cnv_tran': df.ix['CNV', 'Transcriptomics'],
        'cnv_prot': df.ix['CNV', 'Proteomics'],
        'trans_prot': df.ix['Transcriptomics', 'Proteomics']
    }

cors = DataFrame(cors).T
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
cors['type'] = ['px' if i in associations_px else 'py' for i in cors.index]
print cors

sns.boxplot(x='type', y='diff', data=cors, notch=True)

df = DataFrame([{'px': px, 'py': py, 'rme': np.mean(proteomics.ix[px, samples] - proteomics.ix[py, samples])} for px, py in associations if px in residuals.index])
print df.sort('rme')


# --
genes = list(transcriptomics.index)
samples = set(residuals).intersection(cnv)

# px, py = 'COG3', 'COG2'
coef = []
for px, py in associations:
    py = residuals.ix[py, samples].dropna()
    py = py.apply(lambda x: 1 if x < py.quantile(.25) else 0)

    px = cnv.ix[px, py.index].replace([2, 1], 0).replace([-2, -1], 1)
    px = np.logical_and(px, py).astype(int)

    for g in genes:
        pz = transcriptomics.T.ix[px.index, g]
        t, pval = ttest_ind(pz.ix[np.logical_not(px).values.astype(bool)], pz.ix[px.values.astype(bool)])

        coef.append({'px': px, 'py': py, 'gene': g, 't': t, 'pval': pval})

coef = DataFrame(coef)
print coef.sort('pval')
