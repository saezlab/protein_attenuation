#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pickle
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from cptac.slapenrich import slapenrich
from scipy.stats import chi2
from scipy.stats.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from sklearn.metrics.ranking import roc_auc_score
from cptac.utils import log_likelihood, f_statistic, r_squared, randomise_matrix
from pymist.utils.corumdb import get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pandas import read_csv, pivot_table, concat, Series, DataFrame


# -- Protein attenuation
p_cor = read_csv('./tables/proteins_correlations.csv', index_col=0)
p_att = p_cor['CNV_Transcriptomics'] - p_cor['CNV_Proteomics']


# -- CORUM
uniprot = read_uniprot_genename()
with open('./tables/corum_dict_non_redundant.pickle', 'rb') as handle:
    corum_dict = pickle.load(handle)

corum_n = get_complexes_name()

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}
print 'corum', len(corum_proteins)


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


# -- Import data
# Samplesheets
samplesheet = read_csv('./data/sanger_samplesheet.csv', index_col=0).dropna(subset=['TCGA'])
c_samplesheet = read_csv('./data/sanger_drug_samplesheet.csv')
d_targets = c_samplesheet.groupby('Name')['Putative_target'].agg(lambda x: set(x))

# Drug response
drug = read_csv('./data/sanger_drug_response_auc.csv', index_col=1).drop('cosmic', axis=1).T
# drug = read_csv('./data/sanger_drug_response_ic50.csv', index_col=1).drop('cosmic', axis=1).T

# Copy-number
cnv = read_csv('./data/sanger_copy_number.tsv', sep='\t')
cnv['value'] = [1 if i == 'gain' else (-1 if i == 'low' else 0) for i in cnv['MUT_TYPE']]
cnv['gene'] = [i.split('_')[0] for i in cnv['gene_name']]

cnv = cnv.groupby(['SAMPLE_NAME', 'gene'])['value'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='gene', columns='SAMPLE_NAME', values='value', fill_value=0)
print cnv


# --
mutation_m = cnv.copy()
background = set(mutation_m.index)

pathways = {
    'attenuated': set(p_att[p_att > p_att.quantile(.85)].index).intersection(corum_proteins).intersection(background)
}

res, pp = slapenrich(mutation_m, pathways, background)
print res

burden = pp['attenuated'].replace(0, np.nan).dropna()
print burden.sort_values()


# --
# drugs = ['Bortezomib', 'MG-132', 'AUY922', 'SNX-2112', '17-AAG', 'Elesclomol', 'CCT018159']
drugs = ['Bortezomib']

plot_df = DataFrame([{'drug': d, 'cell': c, 'auc': drug.ix[d, c], 'cnv': burden[c]} for d in drugs for c in burden.index if c in drug.columns]).dropna()
plot_df['burden'] = ['High' if i < .05 else 'Low' for i in plot_df['cnv']]
plot_df['counts'] = [cnv.ix[corum_proteins, i].sum() for i in plot_df['cell']]
print plot_df.sort('auc')


sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=3, aspect=1, legend_out=False)
g = g.map_dataframe(sns.stripplot, x='burden', y='auc', orient='v', split=True, size=2, jitter=.2, alpha=.6, linewidth=.1, edgecolor='white')
g = g.map_dataframe(sns.boxplot, x='burden', y='auc', orient='v', linewidth=.3, sym='', notch=True)
g.despine(trim=True)
g.set_axis_labels('Copy-number burden', 'Drug response (AUC)')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/drug_response_proteasome_boxplot.pdf', bbox_inches='tight')
plt.savefig('./reports/drug_response_proteasome_boxplot.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'


sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.jointplot(x='auc', y='cnv', data=plot_df, space=0)
g.set_axis_labels('Drug response (AUC)', 'Copy-number burden')
sns.despine(trim=True)
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/drug_response_proteasome.pdf', bbox_inches='tight')
plt.savefig('./reports/drug_response_proteasome.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Done'
