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


# --
plot_df = crispr.unstack().reset_index()
plot_df.columns = ['cell_line', 'gene', 'score']
plot_df['diff'] = [cor_diff[i.split('_')[0]] if i.split('_')[0] in cor_diff else np.nan for i in plot_df['gene']]
plot_df['diff_b'] = ['not attenuated' if i < .1 else ('attenuated' if i > .5 else 'all') for i in plot_df['diff']]
plot_df['type'] = ['px' if i.split('_')[0] in px else ('py' if i.split('_')[0] in py else 'all') for i in plot_df['gene']]
plot_df['tissue'] = ['_'.join(i.split('_')[1:]).lower() for i in plot_df['cell_line']]

#
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='score', y='type', data=plot_df, col='tissue', col_wrap=5, kind='box', fliersize=2, linewidth=0.3, size=2, orient='h', sharex=False)
g.map(plt.axvline, x=0, ls='--', lw=0.3, c='black', alpha=.5)
g.despine()
plt.savefig('./reports/crispr_regulators_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/crispr_regulators_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

#
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='score', y='diff_b', data=plot_df.dropna(), col='tissue', col_wrap=5, kind='box', fliersize=2, linewidth=0.3, size=2, orient='h', sharex=False)
g.map(plt.axvline, x=0, ls='--', lw=0.3, c='black', alpha=.5)
g.despine()
plt.savefig('./reports/crispr_attenuated_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/crispr_attenuated_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


# --
px_cor = {}
samples = set(crispr).intersection(cnv)

for g in px:
    x = crispr.ix[[i for i in crispr.index if g in i.split('_')], samples]
    y = cnv.ix[g, samples]

    px_cor[g] = x.T.corrwith(y).mean()

px_cor = Series(px_cor).dropna()


py_cor = {}
samples = set(crispr).intersection(cnv)

for g in py:
    x = crispr.ix[[i for i in crispr.index if g in i.split('_')], samples]
    y = cnv.ix[g, samples]

    py_cor[g] = x.T.corrwith(y).mean()

py_cor = Series(py_cor).dropna()

#
df = crispr.unstack().reset_index()
df['cnv'] = [cnv.ix[g.split('_')[0], c] if g.split('_')[0] in cnv.index else np.nan for c, g in df[['level_0', 'Name']].values]
df.columns = ['cell_line', 'gene', 'score', 'cnv']

plt.scatter(df['score'], df['cnv'])
