#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
from pandas import DataFrame, Series, read_csv, concat
from sklearn.linear_model import LinearRegression
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import data-sets
# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape


# -- Overlap
proteins, samples = set(proteomics.index).intersection(transcriptomics.index), set(proteomics).intersection(transcriptomics)
proteomics, transcriptomics = proteomics.ix[proteins, samples], transcriptomics.ix[proteins, samples]
print 'proteins', 'samples', len(proteins), len(samples)

# dicts
proteomics_dict = {s: proteomics[s].dropna().to_dict() for s in proteomics}

# -- Complexes proteins
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Activity functions
# s, c, df = 'TCGA-AG-A020', 2174, proteomics_dict
def ztest_complex(s, c, df):
    x1 = [df[s][k] for k in df[s] if k in corum[c]]
    x2 = [df[s][k] for k in df[s] if k not in corum[c]]

    if len(x1) > 1:
        stat = CompareMeans(DescrStatsW(x1), DescrStatsW(x2))

        z_larger, p_larger = stat.ztest_ind(alternative='larger')
        z_smaller, p_smaller = stat.ztest_ind(alternative='smaller')

        z, p = z_larger, p_larger if p_larger < p_smaller else p_smaller

        res = {
            'sample': s, 'complex': c, 'name': corum_n[c],
            'z': z, 'pval': p, 'mean': np.mean(x1), 'targets': len(x1)
        }

        return res


# -- Estimate CORUM activity z-test
c_proteomics = [ztest_complex(s, c, proteomics_dict) for s in proteomics for c in corum]
c_proteomics = DataFrame([i for i in c_proteomics if i])
c_proteomics['fdr'] = multipletests(c_proteomics['pval'],  method='fdr_bh')[1]
c_proteomics.to_csv('./tables/protein_complexes_proteomics_activities.csv')
print c_proteomics.sort('fdr')


# -- Estimate CORUM activity mean
c_proteomics_mean = DataFrame({c: proteomics.ix[corum[c]].mean() for c in corum}).T
c_proteomics_mean.to_csv('./tables/protein_complexes_proteomics_mean_activities.csv')
print c_proteomics_mean.shape
