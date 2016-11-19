#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv, concat, Series, pivot_table


# -- Achilles_v3.3.8
crispr = read_csv('./data/Achilles_v3.3.8.Gs.txt', sep='\t', index_col=0).drop('Description', axis=1)
print crispr.shape


# -- Achilles_v2.4.3
shrna = read_csv('./data/Achilles_QC_v2.4.3.rnai.Gs.txt', sep='\t', index_col=0).drop('Description', axis=1)
print shrna.shape


# -- CNV
cnv = read_csv('./data/CCLE_copynumber_byGene_2013-12-03.txt', sep='\t', index_col=1).drop(['EGID', 'CHR', 'CHRLOC', 'CHRLOCEND'], axis=1)
print cnv.shape


# -- Protein attenuation
p_cor = read_csv('./tables/proteins_correlations.csv', index_col=0)
cor_diff = (p_cor['CNV_Transcriptomics'] - p_cor['CNV_Proteomics']).to_dict()
print p_cor


# -- Regulatory interactions
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

px, py = set(ppairs_cnv['px']), set(ppairs_cnv['py'])
print len(px), len(py)


#
c = Series(dict(zip(*(np.unique([i.split('_')[0] for i in crispr.index], return_counts=True))))).sort_values()

crispr.index = [i.split('_')[0] for i in crispr.index]

{g: crispr.ix[g].T.corr().ix[0].values[1:] for g in c[c > 1].index}
