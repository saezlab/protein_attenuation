#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, concat
from sklearn.metrics.ranking import roc_auc_score
from sklearn.model_selection import StratifiedKFold


# -- Import samples attenuation
att_s = read_csv('./tables/samples_correlations.csv', index_col=0)['cluster']
att_s = np.logical_not(att_s).astype(int)


# # -- Samples attenuation gene signature
# p_attenuation_cor = []
# for p in transcriptomics.index:
#     df = concat([cors['diff'], transcriptomics.ix[p]], axis=1).dropna()
#     p_attenuation_cor.append({'gene': p, 'cor': df.corr().ix[0, 1]})
#
# p_attenuation_cor = DataFrame(p_attenuation_cor).set_index('gene')
# p_attenuation_cor.to_csv('./tables/samples_attenuated_gene_signature.csv')
# print p_attenuation_cor.ix[corum_proteins].dropna().sort('cor')


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


scores, cv = [], StratifiedKFold(n_splits=10, random_state=True)
for train, test in cv.split(transcriptomics[att_s.index].T, att_s):
    # Samples
    s_train, s_test = set(att_s[train].index), set(att_s[test].index)

    # Sample attenutaion score
    cors = []
    for s in s_train:
        df = DataFrame({
            'cnv': cnv[s], 'trans': transcriptomics[s], 'prot': proteomics[s]
        }).dropna().corr()

        cors.append({'sample': s, 'cnv_tran': df.ix['cnv', 'trans'], 'cnv_prot': df.ix['cnv', 'prot']})
    cors = DataFrame(cors).dropna().set_index('sample')
    cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']

    # Sig
    sig = transcriptomics.loc[:, s_train].T.corrwith(cors.ix[s_train, 'diff'])

    # Classification
    classif = transcriptomics.loc[:, s_test].corrwith(sig)

    # Evalute
    auc_score = roc_auc_score(att_s[s_test], classif[s_test])
    scores.append(auc_score)
    print auc_score

print np.mean(scores), np.std(scores), np.min(scores), np.max(scores)
