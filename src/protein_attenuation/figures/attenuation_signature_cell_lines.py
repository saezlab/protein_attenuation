#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, concat, pivot_table
from scipy.stats.stats import ttest_ind, spearmanr
from sklearn.metrics.ranking import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests


# -- Cell lines: protoemics
brca_proteomics = read_csv('./data/brca_cell_lines_proteomics_preprocessed.csv', index_col=0)
hgsc_proteomics = read_csv('./data/hgsc_cell_lines_proteomics_preprocessed.csv', index_col=0)


# -- Cell lines: transcriptomics
transcriptomics = read_csv('./data/sanger_gene_experssion_rma.tsv', sep='\t')
transcriptomics = pivot_table(transcriptomics, index='GENE_NAME', columns='SAMPLE_NAME', values='Z_SCORE', fill_value=np.nan, aggfunc=np.mean)
print transcriptomics.shape


# -- Cell lines: Copy-number
cnv = read_csv('./data/sanger_copy_number.tsv', sep='\t')
cnv['value'] = [1 if i == 'gain' else (-1 if i == 'low' else 0) for i in cnv['MUT_TYPE']]
cnv['gene'] = [i.split('_')[0] for i in cnv['gene_name']]

cnv = cnv.groupby(['SAMPLE_NAME', 'gene'])['value'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='gene', columns='SAMPLE_NAME', values='value', fill_value=0)
print cnv


# --
proteomics = hgsc_proteomics.copy()

samples = set(proteomics).intersection(transcriptomics).intersection(cnv)
print len(samples)

cors = []
for s in samples:
    df = DataFrame({
        'cnv': cnv[s], 'trans': transcriptomics[s], 'prot': proteomics[s]
    }).dropna().corr()

    cors.append({'sample': s, 'cnv_tran': df.ix['cnv', 'trans'], 'cnv_prot': df.ix['cnv', 'prot']})
cors = DataFrame(cors).dropna().set_index('sample')
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
print cors.sort('diff')


p_attenuation_cor = []
for p in transcriptomics.index:
    df = concat([cors['diff'], transcriptomics.ix[p]], axis=1).dropna()
    p_attenuation_cor.append({'gene': p, 'cor': df.corr().ix[0, 1]})

p_attenuation_cor = DataFrame(p_attenuation_cor).set_index('gene')
print p_attenuation_cor.sort('cor')
