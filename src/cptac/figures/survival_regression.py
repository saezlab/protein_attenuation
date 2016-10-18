#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac.utils import gkn
from pymist.enrichment.gsea import gsea
from cptac import wd, palette_cnv_number
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from sklearn.linear_model import LinearRegression
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, pivot_table, concat


# -- Imports
# Clinical
clinical = read_csv('./data/tcga_clinical.csv', index_col=0).dropna(subset=['time', 'status'])
clinical = clinical[(clinical['admin.disease_code'] == 'ov') & (clinical['time'] < (365 * 10))]
print 'clinical', clinical.shape

# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# -- Overlap
samples, proteins = set(clinical.index).intersection(proteomics).intersection(transcriptomics).intersection(cnv), set(proteomics.index).intersection(cnv.index)
print 'samples', 'proteins', len(samples), len(proteins)


# -- Proteins correlations
p_correlation = read_csv('./tables/proteins_correlations.csv', index_col=0)
p_correlation['diff'] = p_correlation['CNV_Transcriptomics'] - p_correlation['CNV_Proteomics']
print p_correlation.sort('diff').tail()


# --
cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
cgenes = {
    'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
}


# --
a_proteins = p_correlation[p_correlation['diff'] > (p_correlation['diff'].mean() + 2.5 * p_correlation['diff'].std())]
a_proteins = set(a_proteins.index)

df_proteomics = concat([proteomics.ix[a_proteins, samples].replace(np.nan, 0).T, clinical.ix[samples, ['status', 'time']]], axis=1)

cf_proteomics = CoxPHFitter(normalize=False).fit(df_proteomics, 'time', event_col='status')
print cf_proteomics.print_summary()


df_transcriptomics = concat([transcriptomics.ix[a_proteins, samples].T, clinical.ix[samples, ['status', 'time']]], axis=1).dropna()

cf_transcriptomics = CoxPHFitter(normalize=False).fit(df_transcriptomics, 'time', event_col='status')
print cf_transcriptomics.print_summary()


df_cnv = concat([cnv.ix[a_proteins, samples].T, clinical.ix[samples, ['status', 'time']]], axis=1)

cf_cnv = CoxPHFitter(normalize=False).fit(df_cnv, 'time', event_col='status')
print cf_cnv.print_summary()
