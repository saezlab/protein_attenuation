#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import chi2
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from cptac.utils import log_likelihood, f_statistic, r_squared
from pymist.utils.corumdb import get_complexes_dict
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pandas import read_csv, pivot_table, concat, Series, DataFrame


# -- CORUM
uniprot = read_uniprot_genename()
corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[g][0] for g in corum_dict[k] if g in uniprot} for k in corum_dict}

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}
print 'corum', len(corum_proteins)


# -- Import data
samplesheet = read_csv('./data/sanger_samplesheet.csv', index_col=0).dropna(subset=['TCGA'])

drug = read_csv('./data/sanger_drug_response_auc.csv', index_col=1).drop('cosmic', axis=1).T

cnv = read_csv('./data/sanger_copy_number.tsv', sep='\t')
cnv['value'] = [1 if i == 'gain' else -1 for i in cnv['MUT_TYPE']]
cnv = pivot_table(cnv, index='gene_name', columns='SAMPLE_NAME', values='value', fill_value=0)

c_samplesheet = read_csv('./data/sanger_drug_samplesheet.csv')
d_targets = c_samplesheet.groupby('Name')['Putative_target'].agg(lambda x: set(x))


# --
c_cnv = cnv.ix[corum_proteins].dropna().replace(-1, 0).sum().sort_values()
c_cnv.name = 'c_cnv'
print c_cnv

tissue = samplesheet['TCGA'].str.get_dummies()
print tissue


# d = 'SNX-2112'
def regressions(d):
    df = concat([drug.ix[[d]].T, tissue, c_cnv], axis=1).dropna()

    if df.shape[0] > 0:
        # -- 1st model
        # Fit model
        lm = LinearRegression().fit(df.drop([d, 'c_cnv'], axis=1), df[d])

        # Predict
        y_true, y_pred = df[d].copy(), Series(dict(zip(*(df.index, lm.predict(df.drop([d, 'c_cnv'], axis=1))))))

        # Log likelihood
        l_lm = log_likelihood(y_true, y_pred)

        # F-statistic
        f, f_pval = f_statistic(y_true, y_pred, len(y_true), df.drop([d, 'c_cnv'], axis=1).shape[1])

        # R-squared
        r = r_squared(y_true, y_pred)

        res_1 = {
            'drug': d, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
        }
        print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e' % (d, res_1['rsquared'], res_1['f'], res_1['f_pval'])

        sm.OLS(df[d], sm.add_constant(df.drop([d, 'c_cnv'], axis=1), has_constant='add')).fit().summary()

        # -- 2nd model
        lm = LinearRegression().fit(df.drop([d], axis=1), df[d])

        # Predict
        y_true, y_pred = df[d].copy(), Series(dict(zip(*(df.index, lm.predict(df.drop([d], axis=1))))))

        # Log likelihood
        l_lm = log_likelihood(y_true, y_pred)

        # F-statistic
        f, f_pval = f_statistic(y_true, y_pred, len(y_true), df.drop([d], axis=1).shape[1])

        # R-squared
        r = r_squared(y_true, y_pred)

        res_2 = {
            'drug': d, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
        }
        print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e' % (d, res_2['rsquared'], res_2['f'], res_2['f_pval'])

        # -- Likelihood ratio test
        lr = 2 * (res_2['ll'] - res_1['ll'])
        p_val = chi2.pdf(np.float(lr), 1)

        res = {'drug': d, 'lr': lr, 'pval': p_val, 'rsquared': r, 'f': f, 'f_pval': f_pval}
        print res

        return res

ppairs = [regressions(d) for d in drug.index]
ppairs = DataFrame([i for i in ppairs if i])
ppairs['fdr'] = multipletests(ppairs['pval'], method='fdr_bh')[1]
ppairs['targets'] = [';'.join(d_targets.ix[i]) if i in d_targets.index else 'NaN' for i in ppairs['drug']]
ppairs.sort('fdr').to_csv('./tables/drug_response.csv', index=False)
print ppairs[ppairs['fdr'] < .05].sort('fdr')

