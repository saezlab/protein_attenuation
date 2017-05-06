#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from protein_attenuation import palette
from protein_attenuation.enrichment.gsea import gsea
from sklearn.mixture import GaussianMixture
from statsmodels.stats.multitest import multipletests
from pandas import read_csv, DataFrame, Series, concat
from protein_attenuation.utils import read_gmt, gkn, get_complexes_pairs, read_uniprot_genename


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
genes = set(cnv.index).intersection(transcriptomics.index).intersection(proteomics.index)


# -- CORUM
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {uniprot[g][0] for p in corum for g in p if g in uniprot}.intersection(genes)


# -- GO terms
msigdb_go_bp = read_gmt('./files/c5.bp.v5.1.symbols.gmt')
msigdb_go_cc = read_gmt('./files/c5.cc.v5.1.symbols.gmt')
msigdb_go_mf = read_gmt('./files/c5.mf.v5.1.symbols.gmt')

# -- Uniprot PTMs lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('./files/uniprot_ptms_proteins/%s' % f, sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('./files/uniprot_ptms_proteins/') if f.startswith('uniprot_')
}

# -- Correlations
cors = {}
for g in genes:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': proteomics.ix[g, samples]}).dropna().corr()

    cors[g] = {
        'cnv_tran': df.ix['CNV', 'Transcriptomics'],
        'cnv_prot': df.ix['CNV', 'Proteomics']
    }

cors = DataFrame(cors).T
cors['attenuation'] = cors['cnv_tran'] - cors['cnv_prot']

# -- Attenuation GMM
gmm = GaussianMixture(n_components=2).fit(cors[['attenuation']])

s_type = Series(dict(zip(*(cors[['attenuation']].index, gmm.predict(cors[['attenuation']])))))
clusters = Series(dict(zip(*(range(2), gmm.means_[:, 0]))))

cors['cluster'] = [s_type[i] for i in cors.index]
cors['attenuation_potential'] = ['High' if i == clusters.argmax() else 'Low' for i in cors['cluster']]
cors.sort(['cluster', 'attenuation'], ascending=False).to_csv('./tables/protein_attenuation_table.csv')
# cors = read_csv('./tables/protein_attenuation_table.csv', index_col=0)
# print cors.sort(['cluster', 'attenuation'], ascending=False)
print '[INFO] Protein attenuation potential table: ', './tables/protein_attenuation_table.csv'


# -- Plot scatter of correlations
pal = {'Low': palette['Clinical'], 'High': palette['Transcriptomics']}
ax_min, ax_max = np.min([cors['cnv_tran'].min() * 1.10, cors['cnv_prot'].min() * 1.10]), np.max([cors['cnv_tran'].max() * 1.10, cors['cnv_prot'].max() * 1.10])

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
g = sns.jointplot(
    'cnv_tran', 'cnv_prot', cors, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
    space=0, s=5, edgecolor='w', linewidth=.1, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.1
)

g.x = cors.loc[cors['attenuation_potential'] == 'High', 'cnv_tran']
g.y = cors.loc[cors['attenuation_potential'] == 'High', 'cnv_prot']
g.plot_joint(sns.regplot, color=pal['High'], fit_reg=False, scatter_kws={'s': 5, 'alpha': .5})
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(pal['High'], as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
g.plot_marginals(sns.kdeplot, color=pal['High'], shade=True, legend=False)

g.x = cors.loc[cors['attenuation_potential'] == 'Low', 'cnv_tran']
g.y = cors.loc[cors['attenuation_potential'] == 'Low', 'cnv_prot']
g.plot_joint(sns.regplot, color=pal['Low'], fit_reg=False, scatter_kws={'s': 5, 'alpha': .5})
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(pal['Low'], as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
g.plot_marginals(sns.kdeplot, color=pal['Low'], shade=True, legend=False)

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)

handles = [mpatches.Circle([0, 0], .25, facecolor=pal[s], label=s) for s in pal]
g.ax_joint.legend(loc='upper left', handles=handles, title='Protein\nattenuation')
plt.gcf().set_size_inches(3, 3)

g.set_axis_labels('Copy-number ~ Transcriptomics\n(Pearson)', 'Copy-number ~ Proteomics\n(Pearson)')
plt.savefig('./reports/correlation_difference_lmplot_corr.png', bbox_inches='tight', dpi=600)
plt.close('all')
print '[INFO] Protein attenuation scatter: ', './reports/correlation_difference_lmplot_corr.png'


# -- Enrichment
dataset = gkn(cors['cnv_tran'] - cors['cnv_prot']).to_dict()
signatures = {'PTM': ptms, 'BP': msigdb_go_bp, 'CC': msigdb_go_cc}

df_enrichment = [(t, sig, len(db[sig].intersection(dataset)), gsea(dataset, db[sig], 1000)) for t, db in signatures.items() for sig in db]
df_enrichment = DataFrame([{'type': t, 'signature': s, 'length': l, 'escore': es, 'pvalue': pval} for t, s, l, (es, pval) in df_enrichment]).dropna()
df_enrichment.sort(['pvalue', 'escore']).to_csv('./tables/protein_attenuation_enrichment.csv', index=False)
# df_enrichment = read_csv('./tables/protein_attenuation_enrichment.csv')
print '[INFO] Protein attenuation potential table: ', './tables/protein_attenuation_enrichment.csv'

# Plot - Striplot
plot_df = df_enrichment[(df_enrichment['length'] > 5) & (df_enrichment['type'] != 'MF')].copy()
plot_df['fdr'] = multipletests(plot_df['pvalue'], method='fdr_bh')[1]
plot_df = plot_df[plot_df['fdr'].abs() < .05].sort('escore')
plot_df['signature'] = [i.replace('_', ' ').lower() for i in plot_df['signature']]
plot_df = concat([plot_df.head(20), plot_df.tail(20)])

pal = dict(zip(*(set(plot_df['type']), sns.color_palette('Set1', n_colors=4).as_hex())))

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.stripplot(y='signature', x='escore', hue='type', data=plot_df, palette=pal, s=8)
plt.axvline(0, ls='-', lw=0.3, c='black', alpha=.5)
sns.despine(trim=True)
plt.xlabel('Enrichment score (GSEA)')
plt.ylabel('')
plt.title('Copy-number correlation attenuation enrichment\n(Copy-number~Transcriptomics - Copy-number~Proteomics)')
plt.gcf().set_size_inches(2, 10)
plt.savefig('./reports/protein_correlation_difference_enrichment.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_correlation_difference_enrichment.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Protein attenuation complexes: ', './reports/protein_correlation_difference_enrichment.pdf'

# Boxplot
plot_df = df_enrichment[(df_enrichment['length'] > 5) & (df_enrichment['type'] != 'MF')].copy()
plot_df['escore'] = gkn(plot_df['escore'])
plot_df['type'] = ['Complex/subunit' if 'COMPLEX' in i or 'SUBUNIT' in i else 'Other' for i in plot_df['signature']]

pal = {'Other': palette['Clinical'], 'Complex/subunit': palette['Transcriptomics']}
order = ['Complex/subunit', 'Other']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=1, aspect=2)
g = g.map_dataframe(sns.stripplot, 'escore', 'type', orient='h', size=4, jitter=.2, alpha=.2, linewidth=.1, edgecolor='white', palette=pal, order=order)
g = g.map_dataframe(sns.boxplot, 'escore', 'type', orient='h', linewidth=.3, sym='', palette=pal, notch=True, order=order)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Normalised GSEA enrichment score')
g.despine(trim=True)
plt.savefig('./reports/protein_correlation_difference_enrichment_boxplot.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_correlation_difference_enrichment_boxplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Protein attenuation complexes boxplot: ', './reports/protein_correlation_difference_enrichment_boxplot.pdf'
