#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import pickle
import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pymist.enrichment.gsea import gsea
from statsmodels.stats.multitest import multipletests
from cptac.slapenrich import slapenrich
from cptac.utils import read_gmt
from cptac import palette, palette_cnv_number, default_color
from matplotlib.gridspec import GridSpec
from scipy.stats.stats import pearsonr
from matplotlib_venn import venn2, venn2_circles
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


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

# Residuals
residuals = read_csv('./data/residuals_protein_transcript.csv', index_col=0)
print 'residuals', residuals.shape


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
print 'samples', len(samples)


# -- Improt regression results
ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
print ppairs_cnv.sort('fdr')

ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
print ppairs_trans.sort('fdr')

px_highlight = ['EIF3A', 'RPA2', 'COG3', 'COG6', 'SMARCA2']


# -- E3 and DUB ligases
e3_proteins = {(p1, p2) for p1, p2 in read_csv('./tables/E3_target_pairs.csv')[['Source', 'Target']].values}
print 'e3_proteins', len(e3_proteins)

dup_proteins = {(p1, p2) for p1, p2 in read_csv('./tables/DUB_target_pairs.csv')[['Source', 'Target']].values}
print 'dup_proteins', len(dup_proteins)


# -- Venn: overlap between transcriptomics and CNV
sns.set(style='white', font_scale=1., rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

associations = {
    'Copy-number variation': {(px, py) for px, py, fdr in ppairs_cnv[['px', 'py', 'fdr']].values if fdr < .05},
    'Transcriptomics': {(px, py) for px, py, fdr in ppairs_trans[['px', 'py', 'fdr']].values if fdr < .05}
}

venn2(associations.values(), set_labels=associations.keys(), set_colors=[palette[k] for k in associations])
venn2_circles(associations.values(), linestyle='solid', color='white')

plt.savefig('./reports/regressions_associations_venn.pdf', bbox_inches='tight')
plt.savefig('./reports/regressions_associations_venn.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


# -- Volcano: Transcriptomics
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

plot_df = ppairs_trans.copy()

plt.scatter(
    x=plot_df['beta'], y=-np.log10(plot_df['fdr']),
    s=25, c=[palette['Transcriptomics'] if f < .05 else sns.light_palette(palette['Transcriptomics']).as_hex()[1] for f in plot_df['fdr']],
    linewidths=[.5 if (px, py) in associations['Transcriptomics'] else 0 for px, py in plot_df[['px', 'py']].values],
    alpha=.7, edgecolor=[palette['Copy-number variation'] if (px, py) in associations['Copy-number variation'] else '#FFFFFF' for px, py in plot_df[['px', 'py']].values]
)

for fdr, beta, px, py in plot_df[['fdr', 'beta', 'px', 'py']].values:
    if fdr < .05 and px in px_highlight:
        plt.text(beta, -np.log10(fdr), '%s ~ %s' % (px, py), fontsize=6)

plt.axhline(-np.log10(0.01), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axhline(-np.log10(0.05), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axvline(0, c='#99A3A4', ls='-', lw=.3, alpha=.7)

plt.ylim(0)

sns.despine()
plt.ylabel('Adj. p-value (-log10)')
plt.xlabel('Beta')
plt.gcf().set_size_inches(4.5, 7)
plt.savefig('./reports/regressions_associations_volcano_transcriptomics.png', bbox_inches='tight', dpi=600)
# plt.savefig('./reports/regressions_associations_volcano.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Volcano: Transcriptomics
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

plot_df = ppairs_cnv.copy()
# plot_df = ppairs_cnv[ppairs_cnv['fdr'] < .05]

plt.scatter(
    x=plot_df['beta'], y=-np.log10(plot_df['fdr']),
    s=25, c=[palette['Copy-number variation'] if f < .05 else sns.light_palette(palette['Copy-number variation']).as_hex()[1] for f in plot_df['fdr']],
    linewidths=[.5 if (px, py) in associations['Transcriptomics'] else 0 for px, py in plot_df[['px', 'py']].values],
    alpha=.7, edgecolor=[palette['Transcriptomics'] if (px, py) in associations['Transcriptomics'] else '#FFFFFF' for px, py in plot_df[['px', 'py']].values]
)

for fdr, beta, px, py in plot_df[['fdr', 'beta', 'px', 'py']].values:
    if fdr < .05 and px in px_highlight:
        plt.text(beta, -np.log10(fdr), '%s ~ %s' % (px, py), fontsize=6)

plt.axhline(-np.log10(0.01), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axhline(-np.log10(0.05), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axvline(0, c='#99A3A4', ls='-', lw=.3, alpha=.7)

plt.ylim(0)

sns.despine()
plt.ylabel('Adj. p-value (-log10)')
plt.xlabel('Beta')
plt.gcf().set_size_inches(4.5, 7)
plt.savefig('./reports/regressions_associations_volcano.png', bbox_inches='tight', dpi=600)
# plt.savefig('./reports/regressions_associations_volcano.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
plot_df = DataFrame([
    ('Copy-number variation', 'Positive', ppairs_cnv[(ppairs_cnv['fdr'] < .05) & (ppairs_cnv['beta'] > 0)].shape[0]),
    ('Copy-number variation', 'Negative', ppairs_cnv[(ppairs_cnv['fdr'] < .05) & (ppairs_cnv['beta'] < 0)].shape[0]),
    ('Transcriptomics', 'Positive', ppairs_trans[(ppairs_trans['fdr'] < .05) & (ppairs_trans['beta'] > 0)].shape[0]),
    ('Transcriptomics', 'Negative', ppairs_trans[(ppairs_trans['fdr'] < .05) & (ppairs_trans['beta'] < 0)].shape[0])
])
plot_df.columns = ['Data', 'Type', 'Counts']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='Type', y='Counts', data=plot_df, row='Data', sharey=False, kind='bar', ci=None, legend_out=False, lw=.3, color='#99A3A4', aspect=1, size=2)
g.despine()
g.set_ylabels('Count')
g.set_titles('{row_name}')
plt.gcf().set_size_inches(1.5, 4)
plt.savefig('./reports/regressions_associations_count.pdf', bbox_inches='tight')
plt.savefig('./reports/regressions_associations_count.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


# -- Export regulations table
reg_table = ppairs_cnv[ppairs_cnv['fdr'] < .05].copy()
reg_table = reg_table[[(px, py) in associations['Transcriptomics'] for px, py in reg_table[['px', 'py']].values]]

reg_table_counts = Series(dict(zip(*(np.unique(reg_table['px'], return_counts=True)))))
reg_table['counts'] = [1. / reg_table_counts.ix[px] for px in reg_table['px']]

reg_table.sort(['counts', 'f_pval']).drop('counts', axis=1)[['px', 'py', 'beta', 'fdr']].to_csv('./tables/complex_regulators.txt', index=False)
print '[INFO] Export done'


# -- Create network
network_i = igraph.Graph(directed=True)

# Initialise network lists
edges = [(px, py) for px, py in ppairs_cnv[ppairs_cnv['fdr'] < .05][['px', 'py']].values if (px, py) in associations['Transcriptomics']]
vertices = list({p for px, py in edges for p in (px, py)})

# Add nodes
network_i.add_vertices(vertices)
print network_i.summary()

# Add edges
network_i.add_edges(edges)
print network_i.summary()

# Draw network
graph = pydot.Dot(graph_type='digraph')

graph.set_graph_defaults(packMode='clust', pack='true')
graph.set_node_defaults(fontcolor='white', penwidth='5', fillcolor='#CCCCCC', width='1', height='1', fontsize='20', fontname='sans-serif')
graph.set_edge_defaults(color='#CCCCCC', arrowhead='vee', penwidth='2.')

edges = DataFrame([{'index': i.index, 'degree': network_i.degree(i.index, mode='OUT'), 'fdr': -np.log10(ppairs_cnv.loc[ppairs_cnv['px'] == i['name'], 'fdr'].min())} for i in network_i.vs])
edges = edges[edges['degree'] != 0]
edges = edges.sort(['degree', 'fdr'], ascending=False)

for i in edges.index:
    for edge in network_i.es[network_i.incident(i)]:
        source_id, target_id = network_i.vs[[edge.source, edge.target]]['name']

        source = pydot.Node(source_id, style='filled', shape='ellipse', penwidth='0')
        source.set_fillcolor(palette['Copy-number variation'])

        target = pydot.Node(target_id, style='filled', shape='ellipse', penwidth='0')
        target.set_fillcolor(palette['Proteomics'])

        graph.add_node(source)
        graph.add_node(target)

        edge = pydot.Edge(source, target)
        graph.add_edge(edge)

graph.write_png('./reports/regressions_associations_network.png')
graph.write_pdf('./reports/regressions_associations_network.pdf')
print '[INFO] Network PDF exported'


# -- Scatter
def ppair_correlation(px, py):
    x, y = zip(*proteomics.ix[[px, py]].T.dropna().values)
    return pearsonr(x, y)

plot_df = ppairs_cnv[(ppairs_cnv['px'].isin(['COG3', 'SMARCA2'])) & (ppairs_cnv['fdr'] < .05)]
plot_df['cor'] = [ppair_correlation(px, py)[0] for px, py in plot_df[['px', 'py']].values]

# px, py = 'COG3', 'COG2'
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
gs, pos = GridSpec(len(plot_df), 2, hspace=.5, wspace=.5), 0
for px, py in plot_df.sort('cor', ascending=False)[['px', 'py']].values:
    #
    ax = plt.subplot(gs[pos])

    y = residuals.ix[py, samples].dropna()
    x = cnv.ix[px, y.index]
    df = concat([x, y], axis=1).dropna().astype(float)

    sns.boxplot(x=px, y=py, data=df, ax=ax, palette=palette_cnv_number, sym='', linewidth=.3)
    sns.stripplot(x=px, y=py, data=df, ax=ax, palette=palette_cnv_number, jitter=True, size=3, linewidth=.3, edgecolor='white')
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('%s (copy number)' % px)
    ax.set_ylabel('%s (residuals)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(x, y))
    ax.set_ylim(df[py].min() * 1.05, df[py].max() * 1.05)

    #
    ax = plt.subplot(gs[pos + 1])

    y = proteomics.ix[py, samples].dropna()
    x = proteomics.ix[px, y.index]
    df = concat([x, y, cnv.ix[px, y.index].rename('cnv')], axis=1).dropna()

    sns.regplot(df[px], df[py], ax=ax, color=default_color, fit_reg=True, scatter=True, truncate=True, line_kws={'linewidth': .3})
    for c in [0, -1, 1, -2, 2]:
        sns.regplot(df[df['cnv'] == c][px], df[df['cnv'] == c][py], ax=ax, color=palette_cnv_number[c], fit_reg=False, truncate=True)
    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('%s (proteomics)' % px)
    ax.set_ylabel('%s (proteomics)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(df[px], df[py]))
    ax.set_ylim(df[py].min() * 1.05, df[py].max() * 1.05)

    pos += 2

plt.gcf().set_size_inches(6, 3 * len(plot_df))
plt.savefig('./reports/regressions_associations_scatter.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/regressions_associations_scatter.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Px and Py attenuation scores histograms
types = DataFrame([{'Px': px, 'Py': py} for px, py in associations['Copy-number variation'].intersection(associations['Transcriptomics'])])

plot_df = read_csv('./tables/proteins_correlations.csv', index_col=0)
plot_df = DataFrame(plot_df['CNV_Transcriptomics'] - plot_df['CNV_Proteomics'])
plot_df['type'] = ['Px' if i in set(types['Px']) else ('Py' if i in set(types['Py']) else 'Other') for i in plot_df.index]
plot_df.columns = ['diff', 'type']

pal = {'Px': '#2980B9', 'Py': '#E74C3C', 'Other': '#99A3A4'}

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
for t in ['Px', 'Py', 'Other']:
    sns.distplot(plot_df.ix[plot_df['type'] == t, 'diff'], hist=False, kde_kws={'shade': True, 'linewidth': .3}, label=t, color=pal[t])

plt.axvline(0, ls='-', lw=0.3, c='black', alpha=.5)
sns.despine(trim=True)
plt.xlabel('Copy-number correlation attenuation')
plt.ylabel('Density')
plt.gcf().set_size_inches(3, 2)
plt.savefig('./reports/protein_attenuation_histogram.pdf', bbox_inches='tight')
plt.savefig('./reports/protein_attenuation_histogram.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
