#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from cptac import palette, palette_dbs
from scipy.stats.stats import ttest_ind
from pandas import read_csv, DataFrame, Series
from pymist.enrichment.gsea import gsea
from scipy.stats.stats import spearmanr
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
genes = set(cnv.index).intersection(transcriptomics.index).intersection(proteomics.index)
print 'samples', 'genes', len(samples), len(genes)


# -- CORUM
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {uniprot[g][0] for p in corum for g in p if g in uniprot}.intersection(genes)

corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[g][0] for g in corum_dict[k] if g in uniprot and uniprot[g][0] in genes} for k in corum_dict}

corum_name = get_complexes_name()
print 'corum', len(corum)


# -- Correlations
res = {}
for g in genes:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': proteomics.ix[g, samples]}).dropna().corr()

    res[g] = {
        'CNV_Transcriptomics': df.ix['CNV', 'Transcriptomics'],
        'CNV_Proteomics': df.ix['CNV', 'Proteomics'],
        'Transcriptomics_Proteomics': df.ix['Transcriptomics', 'Proteomics']
    }

    print g

res = DataFrame(res).T
res['Interaction'] = ['Complex' if i in corum else 'All' for i in res.index]
print res

res_dict = {}
for c in corum_dict:
    p_corr = proteomics.ix[corum_dict[c].intersection(genes), samples].T.corr()
    p_corr.values[np.tril_indices(p_corr.shape[0], 0)] = np.nan
    p_corr = p_corr.unstack().dropna()

    c_corr = [cnv.ix[g, samples].corr(proteomics.ix[g, samples]) for g in corum_dict[c] if g in genes]

    res_dict[c] = {'Proteomics': list(p_corr), 'CNV': c_corr}

    print corum_name[c]

# -- Plot
# Scatter of correlations
ax_min, ax_max = np.min([res['CNV_Transcriptomics'].min() * 1.10, res['CNV_Proteomics'].min() * 1.10]), np.min([res['CNV_Transcriptomics'].max() * 1.10, res['CNV_Proteomics'].max() * 1.10])

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'CNV_Transcriptomics', 'CNV_Proteomics', res, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
    space=0, s=20, edgecolor='w', linewidth=.5, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.5
)
# g.plot_marginals(sns.kdeplot, shade=True, color='#595959')

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)

g.x = res[res['Interaction'] == 'All']['CNV_Transcriptomics']
g.y = res[res['Interaction'] == 'All']['CNV_Proteomics']
# g.plot_joint(sns.kdeplot, cmap=sns.light_palette(palette_dbs['All'], as_cmap=True), legend=False, shade=True, shade_lowest=False, n_levels=15, alpha=.5)
g.plot_marginals(sns.kdeplot, color=palette_dbs['All'], shade=True, legend=False)

g.x = res[res['Interaction'] == 'Complex']['CNV_Transcriptomics']
g.y = res[res['Interaction'] == 'Complex']['CNV_Proteomics']
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(palette_dbs['CORUM'], as_cmap=True), legend=False, shade=True, shade_lowest=False, n_levels=6, alpha=.8)
g.plot_marginals(sns.kdeplot, color=palette_dbs['CORUM'], shade=True, legend=False)

# cax = g.fig.add_axes([.98, .4, .01, .2])
# cax.axis('off')
handles = [mlines.Line2D([], [], color=palette_dbs['CORUM' if s == 'Complex' else s], linestyle='-', markersize=15, label=s) for s in ['All', 'Complex']]
g.ax_joint.legend(loc=2, handles=handles)

plt.gcf().set_size_inches(5, 5)

g.set_axis_labels('CNV ~ Transcriptomics', 'CNV ~ Proteomics')
plt.savefig('./reports/correlation_difference_lmplot_corr.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'


# Scatter marginals
plot_df = DataFrame([(np.mean(res_dict[c]['Proteomics']), np.mean(res_dict[c]['CNV'])) for c in res_dict]).dropna()
plot_df.columns = ['Proteomics', 'CNV']

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'CNV', 'Proteomics', plot_df, 'reg', color=palette_dbs['CORUM'], joint_kws={'scatter_kws': {'s': 20, 'edgecolor': 'w', 'linewidth': .5}},
    marginal_kws={'hist': False, 'rug': False}, annot_kws={'template': 'Pearson: {val:.2g}, p-value: {p:.1e}'}, space=0,
    xlim=[-.05, .5], ylim=[-.4, 1.]
)

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)

g.x = plot_df['CNV']
g.y = plot_df['Proteomics']
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(palette_dbs['CORUM'], as_cmap=True), legend=False, shade=True, shade_lowest=False, n_levels=10, alpha=.7)

g.plot_marginals(sns.kdeplot, color=palette_dbs['CORUM'], shade=True, legend=False)

plt.gcf().set_size_inches(5, 5)

g.set_axis_labels('Mean correlation of proteomics with CNV', 'Mean correlation of protein complex proteomics')
plt.savefig('./reports/correlation_difference_lmplot_marginals.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'
