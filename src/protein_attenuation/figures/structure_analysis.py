#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, Series, DataFrame, pivot_table
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name


# -- Protein gene names
uniprot = read_uniprot_genename()
uniprot_fasta = read_fasta()


# -- CORUM
corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[g][0] for g in corum_dict[k] if g in uniprot} for k in corum_dict}

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}
print 'corum', len(corum_proteins)


# -- Improt regression results
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans_interactions = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans_interactions for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

px, py = set(ppairs_cnv['px']), set(ppairs_cnv['py'])


# -- Structural information
structure = read_csv('./files/interfaces.tab', sep='\t')

structure = structure[[s in uniprot for s in structure['PROTEIN']]]
structure = structure[[s in uniprot for s in structure['PARTNER']]]

structure['PROTEIN_name'] = [uniprot[i][0] for i in structure['PROTEIN']]
structure['PARTNER_name'] = [uniprot[i][0] for i in structure['PARTNER']]
print structure.head()


# -- Count percentage of interfaces
protein_len = Series({uniprot[k][0]: len(uniprot_fasta[k]) for k in uniprot_fasta if k in uniprot})

interfaces = structure.groupby(['PROTEIN_name', 'POS']).first().reset_index()
interfaces_counts = Series(dict(zip(*(np.unique(interfaces['PROTEIN_name'], return_counts=True)))))

interactions = structure.groupby(['PROTEIN_name', 'PARTNER_name']).first().reset_index()
interactions_counts = Series(dict(zip(*(np.unique(interactions['PROTEIN_name'], return_counts=True)))))

structure_df = DataFrame({p: {'length': protein_len.ix[p], 'interface': interfaces_counts.ix[p], 'interactions': interactions_counts.ix[p]} for p in set(structure['PROTEIN_name'])}).T
structure_df['type'] = ['Px' if i in px else ('Py' if i in py else ('Complex' if i in corum_proteins else 'Other')) for i in structure_df.index]

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(structure_df, size=1, aspect=2)
g = g.map_dataframe(sns.stripplot, x='interface', y='type', orient='h', size=2, jitter=.2, alpha=.35, linewidth=.1, edgecolor='white', color='#99A3A4')
g = g.map_dataframe(sns.boxplot, x='interface', y='type', orient='h', linewidth=.3, sym='', color='#99A3A4', notch=True)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Number of residues on interfaces')
g.set_titles('{row_name}')
g.despine(trim=True)
g.set_ylabels('')
plt.savefig('./reports/structure_interfaces_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/structure_interfaces_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(structure_df, size=1, aspect=2)
g = g.map_dataframe(sns.stripplot, x='interactions', y='type', orient='h', size=2, jitter=.2, alpha=.35, linewidth=.1, edgecolor='white', color='#99A3A4')
g = g.map_dataframe(sns.boxplot, x='interactions', y='type', orient='h', linewidth=.3, sym='', color='#99A3A4', notch=True)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Number interacting proteins')
g.set_titles('{row_name}')
g.despine(trim=True)
g.set_ylabels('')
plt.savefig('./reports/structure_interactions_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/structure_interactions_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


# --
plot_df = structure.groupby(['PROTEIN_name', 'POS'])['AA'].agg(lambda x: set(x)).reset_index()
plot_df = plot_df[[len(i) == 1 for i in plot_df['AA']]]
plot_df['AA'] = [list(i)[0] for i in plot_df['AA']]

plot_df = plot_df.groupby(['PROTEIN_name', 'AA']).count().reset_index()
plot_df = pivot_table(plot_df, index='PROTEIN_name', columns='AA', values='POS', fill_value=0)
plot_df = plot_df.astype(float).divide(plot_df.sum(1), axis=0)

plot_df = plot_df.unstack().reset_index()
plot_df.columns = ['aa', 'gene', 'freq']
plot_df['type'] = ['Px' if i in px else ('Py' if i in py else ('Complex' if i in corum_proteins else 'Other')) for i in plot_df['gene']]

order = Series({aa: plot_df[(plot_df['aa'] == aa) & (plot_df['type'] == 'Px')]['freq'].mean() - plot_df[(plot_df['aa'] == aa) & (plot_df['type'] == 'Py')]['freq'].mean() for aa in set(plot_df['aa'])}).sort_values()

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df[[i in ['Px', 'Py'] for i in plot_df['type']]], size=4, aspect=1, legend_out=True)
# g = g.map_dataframe(sns.stripplot, x='interactions', y='type', orient='h', size=2, jitter=.2, alpha=.35, linewidth=.1, edgecolor='white', color='#99A3A4')
g = g.map_dataframe(sns.boxplot, x='freq', y='aa', hue='type', orient='h', linewidth=.3, fliersize=2, notch=True, order=order.index)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Amino acid frequency in the protein interfaces')
g.set_titles('{row_name}')
g.despine(trim=True)
g.add_legend()
g.set_ylabels('')
plt.savefig('./reports/structure_aa_freq_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/structure_aa_freq_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
