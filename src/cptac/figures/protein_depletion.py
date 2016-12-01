#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from statsmodels.stats.multitest import multipletests
from cptac import palette, palette_dbs
from cptac.utils import read_gmt, gkn
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import LinearRegression
from scipy.stats.stats import ttest_ind
from pandas import read_csv, DataFrame, Series, concat
from pymist.enrichment.gsea import gsea
from sklearn.metrics.ranking import roc_curve, auc, roc_auc_score
from scipy.stats.stats import spearmanr
from pymist.utils.stringdb import get_stringdb
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

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}

corum_name = get_complexes_name()
print 'corum', len(corum)


# -- GO terms
msigdb_go_bp = read_gmt('./files/c5.bp.v5.1.symbols.gmt')
msigdb_go_cc = read_gmt('./files/c5.cc.v5.1.symbols.gmt')
msigdb_go_mf = read_gmt('./files/c5.mf.v5.1.symbols.gmt')
print 'msigdb_go_mf', 'msigdb_go_cc', 'msigdb_go_bp', len(msigdb_go_mf), len(msigdb_go_cc), len(msigdb_go_bp)

# -- Ubiquitination
ubq = read_csv('./tables/ubiquitin_site_dataset.csv')
ubq = ubq[ubq['Regulation'] == 'Up_regulated']
ubq = {uniprot[i][0] for i in ubq['uniprot_accession'] if i in uniprot}

# -- Uniprot PTMs lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('./files/uniprot_ptms_proteins/%s' % f, sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('./files/uniprot_ptms_proteins/') if f.startswith('uniprot_')
}
print 'ptms', len(ptms)


# -- Correlations
cors = {}
for g in genes:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': proteomics.ix[g, samples]}).dropna().corr()

    cors[g] = {
        'cnv_tran': df.ix['CNV', 'Transcriptomics'],
        'cnv_prot': df.ix['CNV', 'Proteomics']
    }

cors = DataFrame(cors).T
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
print cors


# -- Attenuation GMM
gmm = GaussianMixture(n_components=2).fit(cors[['diff']])

s_type = Series(dict(zip(*(cors[['diff']].index, gmm.predict(cors[['diff']])))))
clusters = Series(dict(zip(*(range(2), gmm.means_[:, 0]))))

cors['cluster'] = [s_type[i] for i in cors.index]
cors.sort(['cluster', 'diff'], ascending=False).to_csv('./tables/proteins_correlations.csv')
# cors = read_csv('./tables/proteins_correlations.csv', index_col=0)
print cors.sort(['cluster', 'diff'], ascending=False)

t, pval = ttest_ind(cors['cnv_prot'], cors['cnv_tran'], equal_var=False)
print t, pval

# Plot scatter of correlations
pal = {clusters.argmin(): palette['Clinical'], clusters.argmax(): palette['Transcriptomics']}
ax_min, ax_max = np.min([cors['cnv_tran'].min() * 1.10, cors['cnv_prot'].min() * 1.10]), np.max([cors['cnv_tran'].max() * 1.10, cors['cnv_prot'].max() * 1.10])

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
g = sns.jointplot(
    'cnv_tran', 'cnv_prot', cors, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
    space=0, s=5, edgecolor='w', linewidth=.1, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.1
)

g.x = cors.loc[cors['cluster'] == clusters.argmax(), 'cnv_tran']
g.y = cors.loc[cors['cluster'] == clusters.argmax(), 'cnv_prot']
g.plot_joint(sns.regplot, color=pal[clusters.argmax()], fit_reg=False, scatter_kws={'s': 5, 'alpha': .5})
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(pal[clusters.argmax()], as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
g.plot_marginals(sns.kdeplot, color=pal[clusters.argmax()], shade=True, legend=False)

g.x = cors.loc[cors['cluster'] == clusters.argmin(), 'cnv_tran']
g.y = cors.loc[cors['cluster'] == clusters.argmin(), 'cnv_prot']
g.plot_joint(sns.regplot, color=pal[clusters.argmin()], fit_reg=False, scatter_kws={'s': 5, 'alpha': .5})
g.plot_joint(sns.kdeplot, cmap=sns.light_palette(pal[clusters.argmin()], as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)
g.plot_marginals(sns.kdeplot, color=pal[clusters.argmin()], shade=True, legend=False)

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)

handles = [mpatches.Circle([.5, .5], .5, facecolor=pal[s], label='Attenuated' if s else 'Not attenuated') for s in [1, 0]]
plt.legend(loc='top left', handles=handles, title='Protein')
plt.gcf().set_size_inches(3, 3)

g.set_axis_labels('Copy-number ~ Transcriptomics\n(Pearson)', 'Copy-number ~ Proteomics\n(Pearson)')
plt.savefig('./reports/correlation_difference_lmplot_corr.png', bbox_inches='tight', dpi=600)
plt.close('all')
print '[INFO] Plot done'


# -- Enrichment
dataset = gkn(cors['cnv_tran'] - cors['cnv_prot']).to_dict()
signatures = {'PTM': ptms, 'BP': msigdb_go_bp, 'CC': msigdb_go_cc}

df_enrichment = [(t, sig, len(db[sig].intersection(dataset)), gsea(dataset, db[sig], 1000)) for t, db in signatures.items() for sig in db]
df_enrichment = DataFrame([{'type': t, 'signature': s, 'length': l, 'escore': es, 'pvalue': pval} for t, s, l, (es, pval) in df_enrichment]).dropna()
df_enrichment.sort(['pvalue', 'escore']).to_csv('./tables/protein_depletion_enrichment.csv', index=False)
# df_enrichment = read_csv('./tables/protein_depletion_enrichment.csv')
print df_enrichment[df_enrichment['length'] > 5].sort('escore')

# Plot - Striplot
plot_df = df_enrichment[(df_enrichment['length'] > 5) & (df_enrichment['type'] != 'MF')].copy()
plot_df['fdr'] = multipletests(plot_df['pvalue'], method='fdr_bh')[1]
plot_df = plot_df[plot_df['fdr'].abs() < .05].sort('escore')
plot_df['signature'] = [i.replace('_', ' ').lower() for i in plot_df['signature']]
plot_df = concat([plot_df.head(20), plot_df.tail(20)])
print plot_df.sort('fdr')

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
print '[INFO] Plot done'


sigs = {
    'INTEGRATOR_COMPLEX': msigdb_go_cc['INTEGRATOR_COMPLEX'],
    'BIOGENIC_AMINE_METABOLIC_PROCESS': msigdb_go_bp['BIOGENIC_AMINE_METABOLIC_PROCESS'],
    'sulfation': ptms['sulfation']
}

[gsea(dataset, sigs[k], 1000, './reports/protein_correlation_difference_enrichment_gsea_%s.png' % k, plot_title=k.replace('_', ' ').lower(), y2_label='Correlation difference\n(centered)') for k in sigs]

# Boxplot
plot_df = df_enrichment[(df_enrichment['length'] > 5) & (df_enrichment['type'] != 'MF')].copy()
plot_df['escore'] = gkn(plot_df['escore'])
plot_df['type'] = ['Complex/subunit' if 'COMPLEX' in i or 'SUBUNIT' in i else 'Other' for i in plot_df['signature']]

pal = {'Complex/subunit': palette['Transcriptomics'], 'Other': palette['Clinical']}

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=1, aspect=2)
g = g.map_dataframe(sns.stripplot, 'escore', 'type', orient='h', size=4, jitter=.2, alpha=.2, linewidth=.1, edgecolor='white', palette=pal)
g = g.map_dataframe(sns.boxplot, 'escore', 'type', orient='h', linewidth=.3, sym='', notch=True, palette=pal)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Normalised GSEA enrichment score')
g.despine(trim=True)
plt.savefig('./reports/protein_correlation_difference_enrichment_boxplot.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_correlation_difference_enrichment_boxplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
