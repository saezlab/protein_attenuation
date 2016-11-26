#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import pickle
import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from pymist.enrichment.gsea import gsea
from sklearn.mixture.gaussian_mixture import GaussianMixture
from statsmodels.stats.multitest import multipletests
from cptac.slapenrich import slapenrich
from cptac.utils import read_gmt
from cptac import palette, palette_cnv_number, default_color
from matplotlib.gridspec import GridSpec
from scipy.stats.stats import pearsonr, ttest_ind
from matplotlib_venn import venn2, venn2_circles
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- CORUM
uniprot = read_uniprot_genename()
with open('./tables/corum_dict_non_redundant.pickle', 'rb') as handle:
    corum_dict = pickle.load(handle)

corum_n = get_complexes_name()

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}
print 'corum', len(corum_proteins)


# -- Import data-sets
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# Samplesheet
samplesheet = Series.from_csv('./data/samplesheet.csv')


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
print 'samples', len(samples)


# -- Sample attenutaion score
cors = []
for s in samples:
    df = DataFrame({
        'cnv': cnv[s], 'trans': transcriptomics[s], 'prot': proteomics[s]
    }).dropna().corr()

    cors.append({'sample': s, 'cnv_tran': df.ix['cnv', 'trans'], 'cnv_prot': df.ix['cnv', 'prot']})
cors = DataFrame(cors).dropna().set_index('sample')
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
print cors.sort('diff')


# -- GMM attenuation
gmm = GaussianMixture(n_components=2).fit(cors[['diff']])

s_type = Series(dict(zip(*(cors[['diff']].index, gmm.predict(cors[['diff']])))))
clusters = Series(dict(zip(*(range(2), gmm.means_[:, 0]))))

cors['cluster'] = [s_type[i] for i in cors.index]
cors.sort(['cluster', 'diff'], ascending=False).to_csv('./tables/samples_correlations.csv')
cors['type'] = [samplesheet.ix[i] for i in cors.index]
cors.to_csv('./tables/samples_correlations.csv')
# cors = read_csv('./tables/samples_correlations.csv', index_col=0)
print cors.sort(['cluster', 'diff'], ascending=False)


# -- Correlations scatter plot
ax_min, ax_max = np.min([cors['cnv_tran'].min() * 1.10, cors['cnv_prot'].min() * 1.10]), np.max([cors['cnv_tran'].max() * 1.10, cors['cnv_prot'].max() * 1.10])

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
g = sns.jointplot(
    'cnv_tran', 'cnv_prot', cors, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
    space=0, s=15, edgecolor='w', linewidth=.1, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.3
)
g.plot_marginals(sns.kdeplot, shade=True, color='#99A3A4', lw=.3)

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)

g.x = cors['cnv_tran']
g.y = cors['cnv_prot']
g.plot_joint(sns.kdeplot, cmap=sns.light_palette('#99A3A4', as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)

plt.gcf().set_size_inches(3, 3)

g.set_axis_labels('Copy-number ~ Transcriptomics\n(Pearson)', 'Copy-number ~ Proteomics\n(Pearson)')
plt.savefig('./reports/samples_correlation_difference_lmplot_corr.png', bbox_inches='tight', dpi=600)
plt.savefig('./reports/samples_correlation_difference_lmplot_corr.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Sample correlation representative cases
plot_df = ['TCGA-AA-A00N', 'TCGA-36-1580', 'TCGA-24-1430']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
gs, pos = GridSpec(3, 2, hspace=.5, wspace=.5), 0
for s in plot_df:
    df = concat([cnv[s], transcriptomics[s], proteomics[s]], axis=1).dropna()
    df.columns = ['cnv', 'trans', 'prot']

    #
    ax = plt.subplot(gs[pos])

    sns.boxplot(x='cnv', y='trans', data=df, ax=ax, palette=palette_cnv_number, sym='', linewidth=.3)
    sns.stripplot(x='cnv', y='trans', data=df, ax=ax, palette=palette_cnv_number, jitter=True, size=3, linewidth=.3, edgecolor='white')
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('Copy-number')
    ax.set_ylabel('Transcriptomics')
    ax.set_title('%s (attenuation = %.2f)\nPearson = %.2f' % (s, cors.ix[s, 'diff'], cors.ix[s, 'cnv_tran']))
    ax.set_ylim(df['trans'].min() * 1.05, df['trans'].max() * 1.05)

    #
    ax = plt.subplot(gs[pos + 1])

    sns.boxplot(x='cnv', y='prot', data=df, ax=ax, palette=palette_cnv_number, sym='', linewidth=.3)
    sns.stripplot(x='cnv', y='prot', data=df, ax=ax, palette=palette_cnv_number, jitter=True, size=3, linewidth=.3, edgecolor='white')
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('Copy-number')
    ax.set_ylabel('Proteomics')
    ax.set_title('%s (attenuation = %.2f)\nPearson = %.2f' % (s, cors.ix[s, 'diff'], cors.ix[s, 'cnv_prot']))
    ax.set_ylim(df['trans'].min() * 1.05, df['trans'].max() * 1.05)

    pos += 2

plt.gcf().set_size_inches(6, 9)
plt.savefig('./reports/samples_correlation_difference_scatter_corr.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/samples_correlation_difference_scatter_corr.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Samples attenuation GMM Scatter plot
pal = {clusters.argmin(): '#2980B9', clusters.argmax(): '#E74C3C'}
ax_min, ax_max = np.min([cors['cnv_tran'].min() * 1.10, cors['cnv_prot'].min() * 1.10]), np.max([cors['cnv_tran'].max() * 1.10, cors['cnv_prot'].max() * 1.10])

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
g = sns.jointplot(
    'cnv_tran', 'cnv_prot', cors, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
    space=0, s=15, edgecolor='w', linewidth=.1, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.1
)

g.x = cors.loc[cors['cluster'] == clusters.argmax(), 'cnv_tran']
g.y = cors.loc[cors['cluster'] == clusters.argmax(), 'cnv_prot']
g.plot_joint(sns.regplot, color=pal[clusters.argmax()], fit_reg=False)
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(pal[clusters.argmax()], as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
g.plot_marginals(sns.kdeplot, color=pal[clusters.argmax()], shade=True, legend=False)

g.x = cors.loc[cors['cluster'] == clusters.argmin(), 'cnv_tran']
g.y = cors.loc[cors['cluster'] == clusters.argmin(), 'cnv_prot']
g.plot_joint(sns.regplot, color=pal[clusters.argmin()], fit_reg=False)
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(pal[clusters.argmin()], as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
g.plot_marginals(sns.kdeplot, color=pal[clusters.argmin()], shade=True, legend=False)

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)

handles = [mpatches.Circle([.5, .5], .5, facecolor=pal[s], label='Not attenuated' if s else 'Attenuated') for s in pal]
plt.legend(loc='top left', handles=handles, title='Samples')

plt.gcf().set_size_inches(3, 3)

g.set_axis_labels('Copy-number ~ Transcriptomics\n(Pearson)', 'Copy-number ~ Proteomics\n(Pearson)')
plt.savefig('./reports/samples_correlation_difference_lmplot_corr_attenuated.png', bbox_inches='tight', dpi=600)
plt.savefig('./reports/samples_correlation_difference_lmplot_corr_attenuated.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Slapenrich complexes amplifications
m_matrix = (cnv.loc[:, cors[cors['cluster'] == clusters.argmax()].index] > 1).astype(int)
m_matrix = m_matrix.loc[:, m_matrix.sum() > 0]
print m_matrix.shape

res, pp = slapenrich(m_matrix, corum_dict, set(cnv.index))
res['name'] = [corum_n[int(i.split(':')[0])] for i in res.index]
res.sort('fdr').to_csv('./tables/slapenrich_tumours.csv')
print res[res['fdr'] < .05].sort('fdr')

print Series(dict(zip(*(np.unique([p for i in res[res['fdr'] < .05].index for p in corum_dict[i]], return_counts=True))))).sort_values()

# Plot
plot_df = res[res['fdr'] < .05].copy().sort('logoddratios', ascending=False)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='logoddratios', y='name', data=plot_df, kind='bar', ci=None, legend_out=False, lw=.3, orient='h', color=palette['Clinical'])
g.despine(trim=True)
plt.xlabel('Log odd ratios')
plt.ylabel('')
plt.gcf().set_size_inches(8, 2)
plt.title('Complex amplification enrichment (FDR < 5%)')
plt.savefig('./reports/samples_correlation_difference_complex.pdf', bbox_inches='tight')
plt.savefig('./reports/samples_correlation_difference_complex.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

# Plot
plot_df = DataFrame([{'protein': p, 'cnv': m_matrix.ix[p].sum()} for c in res[res['fdr'] < .05].index for p in corum_dict[c]])
plot_df = plot_df[plot_df['cnv'] != 0]
plot_df = plot_df.groupby('protein').first().reset_index().sort('cnv', ascending=False)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='cnv', y='protein', data=plot_df, kind='bar', ci=None, legend_out=False, lw=.3, orient='h', color=palette['Clinical'])
g.despine(trim=True)
plt.xlabel('Number of amplifications')
plt.ylabel('')
plt.gcf().set_size_inches(2.5, 2.5)
plt.title('Amplified genes in enriched complexes')
plt.savefig('./reports/samples_correlation_difference_complex_amplifications.pdf', bbox_inches='tight')
plt.savefig('./reports/samples_correlation_difference_complex_amplifications.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


# -- Samples attenuation gene signature
p_attenuation_cor = []
for p in transcriptomics.index:
    df = concat([cors['diff'], transcriptomics.ix[p]], axis=1).dropna()
    p_attenuation_cor.append({'gene': p, 'cor': df.corr().ix[0, 1]})

p_attenuation_cor = DataFrame(p_attenuation_cor).set_index('gene')
p_attenuation_cor.to_csv('./tables/samples_attenuated_gene_signature.csv')
print p_attenuation_cor.ix[corum_proteins].dropna().sort('cor')
