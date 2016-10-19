#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac.utils import gkn
from pymist.enrichment.gsea import gsea
from cptac import wd, palette_cnv_number
from lifelines import CoxPHFitter
from sklearn.decomposition.pca import PCA
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
cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
cgenes = {
    'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
}

p_correlation = read_csv('./tables/proteins_correlations.csv', index_col=0)
p_correlation['diff'] = p_correlation['CNV_Transcriptomics'] - p_correlation['CNV_Proteomics']
p_correlation['type'] = ['suppressor' if i in cgenes['suppressor'] else ('oncogene' if i in cgenes['oncogene'] else 'other') for i in p_correlation.index]
print p_correlation.sort('diff').tail()


# -- Regulatory interactions
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')

associations = {(px, py) for px, py in ppairs_cnv[['px', 'py']].values}
print len(associations)


tms = DataFrame({px: transcriptomics.ix[cgenes['suppressor'], samples].dropna(how='all').T.corrwith(proteomics.ix[px, samples]) for px in set(ppairs_cnv['px']) if px in proteins})
tms = tms.unstack().reset_index()
tms.columns = ['px', 'suppressor', 'corr']
tms = tms[tms['px'] != tms['suppressor']]
print tms.sort('corr')

tms = DataFrame({g: transcriptomics.drop(cgenes['suppressor'], errors='ignore')[list(samples)].T.corrwith(proteomics.ix[g, samples]) for g in cgenes['suppressor'] if g in proteins})



# --
a_proteins = p_correlation[p_correlation['diff'] > (p_correlation['diff'].mean() + 2.5 * p_correlation['diff'].std())]
a_proteins = set(a_proteins.index)

df_proteomics = concat([proteomics.ix[a_proteins, samples].replace(np.nan, 0).T, clinical.ix[samples, ['status', 'time']]], axis=1)

cf_proteomics = CoxPHFitter(normalize=False).fit(df_proteomics, 'time', event_col='status')
print cf_proteomics.print_summary()


df_transcriptomics = concat([transcriptomics.ix[a_proteins, samples].T, clinical.ix[samples, ['status', 'time']]], axis=1).dropna()

cf_transcriptomics = CoxPHFitter(normalize=False).fit(df_transcriptomics, 'time', event_col='status')
print cf_transcriptomics.print_summary()
