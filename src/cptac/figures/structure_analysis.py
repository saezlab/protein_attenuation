#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, Series, DataFrame
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
structure_df['percentage'] = structure_df['interface'].astype(float) / structure_df['length'] * 100
structure_df['type'] = ['Px' if i in px else ('Py' if i in py else ('Complex' if i in corum_proteins else 'Other')) for i in structure_df.index]

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='percentage', y='type', data=structure_df, kind='box', linewidth=.3, color='#99A3A4', notch=True, aspect=2, size=1, order=['Other', 'Complex', 'Py', 'Px'], orient='h', fliersize=2)
g.despine()
g.set_xlabels('Percentage of residues in interfaces')
g.set_ylabels('')
plt.gcf().set_size_inches(3, 1.5)
plt.savefig('./reports/structure_interfaces_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/structure_interfaces_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='interactions', y='type', data=structure_df, kind='box', linewidth=.3, color='#99A3A4', notch=True, aspect=2, size=1, order=['Other', 'Complex', 'Py', 'Px'], orient='h', fliersize=2)
g.despine()
g.set_xlabels('Number of interacting proteins with interfaces')
g.set_ylabels('')
plt.gcf().set_size_inches(3, 1.5)
plt.savefig('./reports/structure_interactions_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/structure_interactions_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
