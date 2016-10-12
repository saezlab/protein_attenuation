#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from scipy.stats.stats import pearsonr, ttest_ind
from matplotlib.gridspec import GridSpec
from statsmodels.stats.weightstats import ztest
from cptac import palette, default_color, palette_cnv_number
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta


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
print 'samples', len(samples)


# -- Protein complexes interactions
uniprot = read_uniprot_genename()
uniprot_fasta = read_fasta()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum)


# -- E3 and DUB ligases
e3_proteins = {(p1, p2) for p1, p2 in read_csv('./tables/E3_target_pairs.csv')[['Source', 'Target']].values}
corum.update(e3_proteins)
print 'e3_proteins', len(e3_proteins)

dup_proteins = {(p1, p2) for p1, p2 in read_csv('./tables/DUB_target_pairs.csv')[['Source', 'Target']].values}
corum.update(dup_proteins)
print 'dup_proteins', len(dup_proteins)

print 'corum', len(corum)


# --
# g = 'ZW10'
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in set(transcriptomics.index).intersection(proteomics.index)}).T
print 'residuals', residuals.shape


# -- Regressions: Py Residuals ~ Px CNV
# px, py = 'ARID1A', 'DPF2'
def regressions(px, py):
    if py in residuals.index and px in cnv.index:
        # Protein measurements
        y = residuals.ix[py, samples].dropna()
        x = cnv.ix[[px], y.index].T

        # Fit models
        lm = LinearRegression().fit(x, y)

        # Predict
        y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))

        # Log likelihood
        l_lm = log_likelihood(y_true, y_pred)

        # F-statistic
        f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])

        # R-squared
        r = r_squared(y_true, y_pred)

        res = {
            'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
        }

        print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res

ppairs = [regressions(px, py) for px, py in corum]
ppairs = DataFrame([i for i in ppairs if i])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs.to_csv('./tables/ppairs_cnv_regulation_all.csv', index=False)
# ppairs = read_csv('./tables/ppairs_cnv_regulation_all.csv')
print ppairs.sort('fdr')


# -- Create network
network_i = igraph.Graph(directed=True)

# Initialise network lists
edges = [(px, py) for px, py in ppairs[ppairs['fdr'] < .05][['px', 'py']].values]
vertices = list({p for px, py in edges for p in (px, py)})

# Add nodes
network_i.add_vertices(vertices)
print network_i.summary()

# Add edges
network_i.add_edges(edges)
print network_i.summary()

# Draw network
graph = pydot.Dot(graph_type='digraph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3', fillcolor='#CCCCCC', )
graph.set_edge_defaults(color='#CCCCCC', arrowhead='vee')

for edge in network_i.es:
    source_id, target_id = network_i.vs[[edge.source, edge.target]]['name']

    source = pydot.Node(source_id, style='filled', shape='ellipse', penwidth='0')
    source.set_fillcolor(palette['CNV'])

    target = pydot.Node(target_id, style='filled', shape='ellipse', penwidth='0')
    target.set_fillcolor(palette['Proteomics'])

    graph.add_node(source)
    graph.add_node(target)

    edge = pydot.Edge(source, target)
    graph.add_edge(edge)

graph.write_pdf('./reports/ppairs_cnv_regulation_network.pdf')
print '[INFO] Network PDF exported'


# -- Scatter
def ppair_correlation(px, py):
    x, y = zip(*proteomics.ix[[px, py]].T.dropna().values)
    return pearsonr(x, y)

ppairs_signif = ppairs[ppairs['fdr'] < .05].sort('fdr')
ppairs_signif['cor'] = [ppair_correlation(px, py)[0] for px, py in ppairs_signif[['px', 'py']].values]
ppairs_signif.to_csv('./tables/ppairs_cnv_regulation.csv', index=False)

plot_df_short = ppairs_signif.sort('fdr')[:3]
plot_df_short = ppairs_signif.ix[[32758]]

# px, py = 'COG3', 'COG2'
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs, pos = GridSpec(len(plot_df_short), 2, hspace=.5, wspace=.3), 0
for px, py in plot_df_short[['px', 'py']].values:
    #
    ax = plt.subplot(gs[pos])

    y = residuals.ix[py, samples].dropna()
    x = cnv.ix[px, y.index]
    plot_df = concat([x, y], axis=1).dropna().astype(float)

    # sns.regplot(plot_df[px], plot_df[py], ax=ax, color=default_color, fit_reg=True, scatter=False, truncate=False)
    # for c in [0, -1, 1, -2, 2]:
    #     sns.regplot(plot_df[plot_df[px] == c][px], plot_df[plot_df[px] == c][py], ax=ax, color=palette_cnv_number[c], fit_reg=False, truncate=True)
    sns.boxplot(x=px, y=py, data=plot_df, ax=ax, palette=palette_cnv_number, sym='', linewidth=.3)
    sns.stripplot(x=px, y=py, data=plot_df, ax=ax, palette=palette_cnv_number, jitter=True, size=3, linewidth=.3, edgecolor='white')
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    # ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('%s (copy number)' % px)
    ax.set_ylabel('%s (residuals)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(x, y))
    ax.set_ylim(plot_df[py].min() * 1.05, plot_df[py].max() * 1.05)

    #
    ax = plt.subplot(gs[pos + 1])

    y = proteomics.ix[py, samples].dropna()
    x = proteomics.ix[px, y.index]
    plot_df = concat([x, y, cnv.ix[px, y.index].rename('cnv')], axis=1).dropna()

    sns.regplot(plot_df[px], plot_df[py], ax=ax, color=default_color, fit_reg=True, scatter=True, truncate=True)
    for c in [0, -1, 1, -2, 2]:
        sns.regplot(plot_df[plot_df['cnv'] == c][px], plot_df[plot_df['cnv'] == c][py], ax=ax, color=palette_cnv_number[c], fit_reg=False, truncate=True)
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('%s (proteomics)' % px)
    ax.set_ylabel('%s (proteomics)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(plot_df[px], plot_df[py]))
    ax.set_ylim(plot_df[py].min() * 1.05, plot_df[py].max() * 1.05)

    pos += 2

plt.gcf().set_size_inches(5, 2 * len(plot_df_short))
# plt.savefig('%s/reports/ppairs_cnv_regulation_scatter.png' % wd, bbox_inches='tight', dpi=150)
plt.savefig('./reports/ppairs_cnv_regulation_scatter.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Protein sequence length
ppairs_signif = ppairs_signif[ppairs_signif['cor'] > 0]

p_length = DataFrame([{'uniprot': p, 'name': uniprot[p][0], 'length': len(uniprot_fasta[p])} for p in uniprot_fasta if p in uniprot]).groupby('name')['length'].max().to_dict()

plot_df = DataFrame([{'type': t, 'protein': p, 'length': p_length[p]} for px, py in ppairs_signif[['px', 'py']].values for t, p in [('px', px), ('py', py)] if p in p_length])
print plot_df

ttest, pval = ttest_ind(plot_df[plot_df['type'] == 'px']['length'], plot_df[plot_df['type'] == 'py']['length'])
print 'ttest', 'pval', ttest, pval

# Boxplot
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.violinplot(x='type', y='length', data=plot_df, linewidth=.3, cut=0, inner='quartile', split=False, color=palette_cnv_number[0])
sns.stripplot(x='type', y='length', data=plot_df, linewidth=.3, jitter=True, edgecolor='white', split=False, color=default_color)
plt.ylim(0)
sns.despine(trim=True)
plt.ylabel('Protein sequence length (number of AA)')
plt.title('Px (CNV) ~ Py (Residuals)\nT-test: %.2f, p-value: %.2e' % (ttest, pval))
plt.gcf().set_size_inches(2, 4)
plt.savefig('./reports/protein_pairs_protein_info_boxplots.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'

# Histogram
plot_df = DataFrame([{'px': px, 'py': py, 'len_px': p_length[px], 'len_py': p_length[py], 'signif': int(fdr < .05)} for px, py, fdr in ppairs[['px', 'py', 'fdr']].values if px in p_length and py in p_length])
plot_df['diff'] = plot_df['len_px'] - plot_df['len_py']
print plot_df.sort('diff')

z, zpval = ttest_ind(plot_df.loc[plot_df['signif'] == 1, 'diff'].values, plot_df.loc[plot_df['signif'] == 0, 'diff'].values, equal_var=False)

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
for i in [0, 1]:
    values = plot_df[plot_df['signif'] == i]['diff']
    sns.distplot(values, hist=False, kde_kws={'shade': True}, color=palette_cnv_number[i], label='%s (mean = %.2f)' % ('Significant' if i else 'All', np.mean(values)))
plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Protein length difference\nT-test: %.2f, %.2e' % (z, zpval))
plt.xlabel('len(Px) - len(Py)')
sns.despine(trim=True)
plt.legend()
plt.gcf().set_size_inches(4, 2)
plt.savefig('./reports/protein_pairs_protein_info_histogram.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
