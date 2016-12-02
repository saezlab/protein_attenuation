#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import matplotlib.pyplot as plt
from cptac.utils import gkn
from pandas import read_csv, DataFrame, concat, Series, pivot_table
from scipy.stats.stats import ttest_ind, spearmanr, pearsonr
from sklearn.linear_model.coordinate_descent import ElasticNetCV
from sklearn.metrics.ranking import roc_auc_score
from sklearn.metrics.regression import r2_score
from sklearn.model_selection import StratifiedShuffleSplit, ShuffleSplit
from sklearn.feature_selection import f_regression
from statsmodels.stats.multitest import multipletests


# -- Import samples attenuation
att_s = read_csv('./tables/samples_correlations.csv', index_col=0)
att_s['cluster'] = np.logical_not(att_s['cluster']).astype(int)
att_s['diff'] = gkn(att_s['diff'])


# -- Import data-sets
# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
genes = list(transcriptomics.index)

scores, fs_m, cv = [], [], StratifiedShuffleSplit(test_size=.3)
for train, test in cv.split(transcriptomics[att_s['diff'].index].T, att_s['cluster']):
    # Samples
    s_train, s_test = list(set(att_s.ix[train].index)), list(set(att_s.ix[test].index))

    # Feature selection
    fs = f_regression(transcriptomics.T.ix[s_train, genes], att_s.ix[s_train, 'diff'], center=False)
    fs = DataFrame({'genes': genes, 'f': fs[0], 'pval': fs[1]}).set_index('genes')
    fs['fdr'] = multipletests(fs['pval'], method='fdr_bh')[1]
    features = list(fs[fs['fdr'] < .05].index)
    print '\nfeatures: ', len(features)

    # Linear regression
    x_train = transcriptomics.T.ix[s_train, features]
    y_train = att_s.ix[s_train, 'diff']

    cv_i = ShuffleSplit(test_size=.2)
    lm = ElasticNetCV(cv=cv_i).fit(x_train, y_train)
    print 'alpha: ', lm.alpha_

    fs_m.append(dict(zip(*(features, lm.coef_))))

    pred = Series(dict(zip(*(s_train, lm.predict(x_train)))))
    p_cor, pval = pearsonr(att_s.ix[s_train, 'diff'], pred[s_train])
    r2 = r2_score(att_s.ix[s_train, 'diff'], pred[s_train])
    print 'train: ', p_cor, r2, pval

    # Evalute
    x_test = transcriptomics.T.ix[s_test, features]
    pred = Series(dict(zip(*(s_test, lm.predict(x_test)))))

    p_cor, pval = pearsonr(att_s.ix[s_test, 'diff'], pred[s_test])
    r2 = r2_score(att_s.ix[s_test, 'diff'], pred[s_test])

    auc = roc_auc_score(att_s.ix[s_test, 'cluster'], pred[s_test])

    scores.extend([{'auc': auc, 'cor': p_cor, 'pval': pval, 'r2': r2}])
    print 'test: ', p_cor, r2, auc, pval

scores = DataFrame(scores)
print scores.median()


# -- Cell lines
# Transcriptomics
c_transcriptomics = read_csv('./data/sanger_gene_experssion_rma.tsv', sep='\t')
c_transcriptomics = pivot_table(c_transcriptomics, index='GENE_NAME', columns='SAMPLE_NAME', values='Z_SCORE', fill_value=np.nan, aggfunc=np.mean)
print c_transcriptomics.shape

# - Model
genes = list(set(transcriptomics.index).intersection(c_transcriptomics.index))

fs = f_regression(transcriptomics.T.ix[att_s.index, genes], att_s.ix[att_s.index, 'diff'], center=False)
fs = DataFrame({'genes': genes, 'f': fs[0], 'pval': fs[1]}).set_index('genes')
fs['fdr'] = multipletests(fs['pval'], method='fdr_bh')[1]
features = list(fs[fs['fdr'] < .05].index)
print '\nfeatures: ', len(features)

# Linear regression
x_train = transcriptomics.T.ix[att_s.index, features]
y_train = att_s.ix[att_s.index, 'diff']

cv_i = ShuffleSplit(test_size=.2)
lm = ElasticNetCV(cv=cv_i).fit(x_train, y_train)
print 'alpha: ', lm.alpha_

# Pred
x_test = c_transcriptomics.T[features]
pred = Series(dict(zip(*(x_test.index, lm.predict(x_test)))))

# Export
pred.to_csv('./tables/cell_lines_predicted_attenuation_score.csv')
print '[INFO] Exported'

