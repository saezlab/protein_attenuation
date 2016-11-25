#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pickle
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette
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
trans = pivot_table(trans, index='GENE_NAME', columns='SAMPLE_NAME', values='Z_SCORE', fill_value=np.nan, aggfunc=np.mean)
print trans


# --
sig = read_csv('./tables/samples_attenuated_gene_signature.csv', index_col=0)['cor']

burden = trans.ix[sig.index].corrwith(sig).sort_values()
burden.name = 'burden'
print burden.sort_values()


# --
# d = 'CP466722'
def regressions(d):
    df = concat([drug.ix[d], burden], axis=1).dropna()

    # Correlation
    cor, pval = pearsonr(df['burden'], df[d])

    # Fit model
    lm = LinearRegression().fit(df[['burden']], df[d])

    # Predict
    y_true, y_pred = df[d].copy(), Series(dict(zip(*(df.index, lm.predict(df[['burden']])))))

    # Log likelihood
    l_lm = log_likelihood(y_true, y_pred)

    # F-statistic
    f, f_pval = f_statistic(y_true, y_pred, len(y_true), df[['burden']].shape[1])

    # R-squared
    r = r_squared(y_true, y_pred)

    res = {
        'drug': d, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0], 'pearson': cor, 'pearson_pval': pval
    }
    print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e, ll: %.2f' % (d, res['rsquared'], res['f'], res['f_pval'], res['ll'])

    return res

ppairs = DataFrame([regressions(d) for d in drug.index])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs['pearson_fdr'] = multipletests(ppairs['pearson_pval'], method='fdr_bh')[1]
ppairs['targets'] = [';'.join(d_targets.ix[i]) if i in d_targets.index else 'NaN' for i in ppairs['drug']]
ppairs.sort('fdr').to_csv('./tables/drug_response.csv', index=False)
print ppairs[ppairs['fdr'] < .05].sort('pearson', ascending=False)

ds = ['Bortezomib', 'MG-132', 'AUY922', 'SNX-2112', '17-AAG', 'Elesclomol', 'CCT018159']
print ppairs[[i in ds for i in ppairs['drug']]]

# Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

plot_df = ppairs.copy()

plt.scatter(
    x=plot_df['beta'], y=-np.log10(plot_df['fdr']), s=25, linewidths=.3, alpha=.7,
    c=[palette['Overlap'] if f < .05 else sns.light_palette(palette['Overlap']).as_hex()[1] for f in plot_df['fdr']]
)

for fdr, beta, d in plot_df[['fdr', 'beta', 'drug']].values:
    if fdr < .05 and d in ds:
        plt.text(beta, -np.log10(fdr), '%s' % d, fontsize=6)

plt.axhline(-np.log10(0.01), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axhline(-np.log10(0.05), c='#99A3A4', ls='--', lw=.5, alpha=.7)
plt.axvline(0, c='#99A3A4', ls='-', lw=.3, alpha=.7)

plt.ylim(0)

sns.despine()
plt.ylabel('Adj. p-value (-log10)')
plt.xlabel('Beta')
plt.gcf().set_size_inches(4.5, 7)
plt.savefig('./reports/drug_response_associations_volcano.png', bbox_inches='tight', dpi=600)
# plt.savefig('./reports/drug_response_associations_volcano.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
