#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pickle
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from cptac.slapenrich import slapenrich
from scipy.stats import chi2
from scipy.stats.stats import pearsonr, spearmanr, ttest_ind
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from sklearn.metrics.ranking import roc_auc_score
from cptac.utils import log_likelihood, f_statistic, r_squared, randomise_matrix
from pymist.utils.corumdb import get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pandas import read_csv, pivot_table, concat, Series, DataFrame


# -- CORUM
uniprot = read_uniprot_genename()
with open('./tables/corum_dict_non_redundant.pickle', 'rb') as handle:
    corum_dict = pickle.load(handle)

corum_n = get_complexes_name()

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}
print 'corum', len(corum_proteins)

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

# Gene expression
trans = read_csv('./data/sanger_gene_experssion_rma.tsv', sep='\t')
trans['value'] = [1 if i == 'over' else (-1 if i == 'under' else 0) for i in trans['REGULATION']]

# # trans = trans.groupby(['SAMPLE_NAME', 'GENE_NAME'])['value'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
# # trans = trans.reset_index()
# # trans = trans.dropna()
# # trans = pivot_table(trans, index='GENE_NAME', columns='SAMPLE_NAME', values='value', fill_value=0)
# # print trans

trans = pivot_table(trans, index='GENE_NAME', columns='SAMPLE_NAME', values='Z_SCORE', fill_value=np.nan, aggfunc=np.mean)
print trans

# --
sig = read_csv('./tables/samples_attenuated_gene_signature.csv', index_col=0)['cor']

burden = trans.ix[sig.index].corrwith(sig).sort_values()
burden.name = 'burden'
print burden.sort_values()


# --
drugs = ['Bortezomib', 'MG-132', 'AUY922', 'SNX-2112', '17-AAG', 'Elesclomol', 'CCT018159']
# drugs = ['Bortezomib', 'MG-132', ]
# drugs = list(drug.index)

plot_df = DataFrame([{'drug': d, 'cell': c, 'auc': drug.ix[d, c], 'cnv': burden[c]} for d in drugs for c in burden.index if c in drug.columns]).dropna()
plot_df['burden'] = ['High' if i > .1 else ('Low' if i < -.1 else 'ND') for i in plot_df['cnv']]
print plot_df.sort('auc')

t, pval = ttest_ind(plot_df[plot_df['burden'] == 'Low']['auc'], plot_df[plot_df['burden'] == 'High']['auc'])
print 'T-test: %.2f, %.2e' % (t, pval)

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

# #--
# tissue = samplesheet['TCGA'].str.get_dummies()
#
#
# # d = 'CP466722'
# def regressions(d):
#     df = concat([drug.ix[d], tissue, burden], axis=1).dropna()
#
#     if df.shape[0] > 0:
#         # -- 1st model
#         # Fit model
#         lm = LinearRegression().fit(df.drop([d, 'burden'], axis=1), df[d])
#
#         # Predict
#         y_true, y_pred = df[d].copy(), Series(dict(zip(*(df.index, lm.predict(df.drop([d, 'burden'], axis=1))))))
#
#         # Log likelihood
#         l_lm = log_likelihood(y_true, y_pred)
#
#         # F-statistic
#         f, f_pval = f_statistic(y_true, y_pred, len(y_true), df.drop([d, 'burden'], axis=1).shape[1])
#
#         # R-squared
#         r = r_squared(y_true, y_pred)
#
#         res_1 = {
#             'drug': d, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
#         }
#         print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e, ll: %.2f' % (d, res_1['rsquared'], res_1['f'], res_1['f_pval'], res_1['ll'])
#
#         # -- 2nd model
#         lm = LinearRegression().fit(df.drop([d], axis=1), df[d])
#
#         # Predict
#         y_true, y_pred = df[d].copy(), Series(dict(zip(*(df.index, lm.predict(df.drop([d], axis=1))))))
#
#         # Log likelihood
#         l_lm = log_likelihood(y_true, y_pred)
#
#         # F-statistic
#         f, f_pval = f_statistic(y_true, y_pred, len(y_true), df.drop([d], axis=1).shape[1])
#
#         # R-squared
#         r = r_squared(y_true, y_pred)
#
#         res_2 = {
#             'drug': d, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
#         }
#         print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e, ll: %.2f' % (d, res_2['rsquared'], res_2['f'], res_2['f_pval'], res_2['ll'])
#
#         # -- Likelihood ratio test
#         lr = 2 * (res_2['ll'] - res_1['ll'])
#         p_val = chi2.pdf(np.float(lr), 1)
#
#         res = {'drug': d, 'lr': lr, 'pval': p_val, 'rsquared': r, 'f': f, 'f_pval': f_pval}
#         print res
#
#         return res
#
# ppairs = [regressions(d) for d in drug.index]
# ppairs = DataFrame([i for i in ppairs if i])
# ppairs['fdr'] = multipletests(ppairs['pval'], method='fdr_bh')[1]
# ppairs['targets'] = [';'.join(d_targets.ix[i]) if i in d_targets.index else 'NaN' for i in ppairs['drug']]
# # ppairs.sort('fdr').to_csv('./tables/drug_response.csv', index=False)
# print ppairs[ppairs['fdr'] < .05].sort('fdr')
#
