#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.colors
import statsmodels.api as sm
from pandas.stats.misc import zscore
from scipy.stats.stats import pearsonr
from matplotlib_venn import venn2, venn2_circles
from cptac.utils import hypergeom_test, read_gmt
from cptac import palette
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import protein-pairs
uniprot = read_uniprot_genename()
uniprot_proteins = {v[0] for v in uniprot.values()}

# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum)

# Overlap
p_pairs = corum
p_pairs_proteins = {p for p1, p2 in p_pairs for p in [p1, p2]}
print 'p_pairs', len(p_pairs)


# -- Import data-sets
# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# Overlap
proteins = set(transcriptomics.index).intersection(proteomics.index)
samples = set(transcriptomics).intersection(proteomics)
print len(proteins), len(samples)


# -- Import gene-sets
# Uniprot PTMs lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('./files/uniprot_ptms_proteins/%s' % f, sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('./files/uniprot_ptms_proteins/') if f.startswith('uniprot_')
}
print 'ptms', len(ptms)

# GO terms
msigdb_go_bp = read_gmt('./files/c5.bp.v5.1.symbols.gmt')
msigdb_go_cc = read_gmt('./files/c5.cc.v5.1.symbols.gmt')
msigdb_go_mf = read_gmt('./files/c5.mf.v5.1.symbols.gmt')
print 'msigdb_go_mf', 'msigdb_go_cc', 'msigdb_go_bp', len(msigdb_go_mf), len(msigdb_go_cc), len(msigdb_go_bp)

# Pathways
msigdb_cp = read_gmt('./files/c2.cp.v5.1.symbols.gmt')
msigdb_kegg = read_gmt('./files/c2.cp.kegg.v5.1.symbols.gmt')
msigdb_cgp = read_gmt('./files/c2.cgp.v5.1.symbols.gmt')
print 'msigdb_cp', 'msigdb_kegg', 'msigdb_cgp', len(msigdb_cp), len(msigdb_kegg), len(msigdb_cgp)


# -- Protein pairs correlation
# px, py = 'SDHA', 'SDHB'
def protein_correlation(px, py):
    if px in proteins and py in proteins:
        samples = list(concat([proteomics.ix[[px, py]].T, transcriptomics.ix[[px, py]].T], axis=1).dropna().index)

        p_cor, p_pval = pearsonr(proteomics.ix[px, samples], proteomics.ix[py, samples])
        t_cor, t_pval = pearsonr(transcriptomics.ix[px, samples], transcriptomics.ix[py, samples])

        return {
            'px': px, 'py': py,
            'p_cor': p_cor, 'p_pval': p_pval,
            't_cor': t_cor, 't_pval': t_pval,
            'len': len(samples)
        }

cor_df = {'%s_%s' % (px, py): protein_correlation(px, py) for px, py in p_pairs}
cor_df = DataFrame({i: cor_df[i] for i in cor_df if cor_df[i]}).T

cor_df['diff'] = [p_cor - t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]
cor_df['sum'] = [p_cor + t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]

cor_df['t_fdr'] = multipletests(cor_df['t_pval'], method='fdr_bh')[1]
cor_df['p_fdr'] = multipletests(cor_df['p_pval'], method='fdr_bh')[1]
print cor_df.sort('p_fdr')


# -- Enrichment analysis
bkg = set(p_pairs_proteins).intersection(proteins)
p_pairs = {p for pp in cor_df[(cor_df['p_fdr'] < .05) & (cor_df['p_cor'] > .5) & (cor_df['t_cor'] < .5) & (cor_df['t_fdr'] > .05)].index for p in pp.split('_')}.intersection(bkg)
t_pairs = {p for pp in cor_df[(cor_df['t_fdr'] < .05) & (cor_df['t_cor'] > .5) & (cor_df['p_cor'] < .5) & (cor_df['p_fdr'] > .05)].index for p in pp.split('_')}.intersection(bkg)
print 'bkg', 'p_pairs', 't_pairs', len(bkg), len(p_pairs), len(t_pairs)


def enrichment_hypergeom(signature_type, signature):
    # Hypergeometric
    p_pval, p_len = hypergeom_test(signature.intersection(bkg), bkg, p_pairs)
    t_pval, t_len = hypergeom_test(signature.intersection(bkg), bkg, t_pairs)

    s_len = len(signature.intersection(bkg))

    return {
        'type': signature_type,
        'p_pval': p_pval, 'p_len': p_len, 'p_perc': float(p_len) / s_len,
        't_pval': t_pval, 't_len': t_len, 't_perc': float(t_len) / s_len,
        's_len': s_len
    }

hyper = DataFrame({
    signature_name: enrichment_hypergeom(signature_type, signatures[signature_name]) for signature_type, signatures in [('cc', msigdb_go_cc), ('bp', msigdb_go_bp)]
    for signature_name in signatures if len(signatures[signature_name].intersection(bkg)) > 0
}).T.dropna()
hyper['p_fdr'] = multipletests(hyper['p_pval'], method='fdr_bh')[1]
hyper['t_fdr'] = multipletests(hyper['t_pval'], method='fdr_bh')[1]
print hyper.sort(['p_len', 'p_perc'], ascending=False)[['p_len', 'p_perc', 'p_fdr']]

# Plot
plot_df = hyper[(hyper['p_fdr'] < .05) | (hyper['t_fdr'] < .05)]
plot_df = plot_df[(plot_df['p_len'] > 1) | (plot_df['t_len'] > 1)]

plot_df = plot_df[['p_perc', 't_perc', 'p_fdr', 't_fdr']]
plot_df[['p_fdr', 't_fdr']] = ~(plot_df[['p_fdr', 't_fdr']] < .05)
plot_df.index = [i.lower().replace('_', ' ') for i in plot_df.index]
plot_df = plot_df.sort(['p_perc', 't_perc'], ascending=False)

mask = plot_df[['p_fdr', 't_fdr']].copy()
mask.columns = ['Proteomics', 'Transcriptomics']

plot_df = plot_df[['p_perc', 't_perc']].copy().astype(np.float)
plot_df.columns = ['Proteomics', 'Transcriptomics']

sns.set(font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.clustermap(plot_df * 100, mask=mask, linewidths=.2, annot=True, fmt='.0f', cmap=sns.light_palette('#680018', as_cmap=True), row_cluster=False, col_cluster=False, col_colors=Series(palette).rename('Data'))
plt.title('Overlap')
plt.gcf().set_size_inches(1, 10)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/protein_pairs_goterms_enrichment.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
