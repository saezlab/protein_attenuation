#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
from protein_attenuation.utils import gkn
from pandas import read_csv, DataFrame, pivot_table

# -- Proteomics
proteomics = read_csv('./data/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv', sep='\t', index_col=0)

proteomics = proteomics.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)
proteomics = proteomics.drop(['Mean', 'Median', 'StdDev'], axis=0)

proteomics = proteomics[[c for c in proteomics if c.endswith(' Unshared Log Ratio')]]
proteomics.columns = [c.split(' ')[0] for c in proteomics]

# Match sample ids
brca_samplesheet = read_csv('./data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1.csv')
brca_samplesheet['id'] = ['-'.join(i.split('-')[1:4]) for i in brca_samplesheet['TCGA barcode']]
brca_samplesheet = brca_samplesheet.set_index('id')

proteomics = proteomics.loc[:, [i in brca_samplesheet.index for i in proteomics]]
proteomics.columns = [brca_samplesheet.ix[i, 'TCGA barcode'] for i in proteomics]

# QC status samples
samples_qc = read_csv('./data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv', index_col=2)
samples_qc.index = [i[:12] for i in samples_qc.index]
# samples_qc = samples_qc[samples_qc['QC Status'] == 'pass']

proteomics.columns = [i[:12] for i in proteomics]
proteomics = proteomics.loc[:, [i in samples_qc.index for i in proteomics]]

# Normalise
proteomics = DataFrame({i: gkn(proteomics.ix[i].dropna()).to_dict() for i in proteomics.index}).T

# Export
proteomics.to_csv('./data/cptac_brca_qc_all_proteomics_normalised.csv')
print '[INFO] BRCA tumours protoemics: ', './data/cptac_brca_qc_all_proteomics_normalised.csv'


# -- Transcriptomics
# Import TCGA pancancer rna-seq data
transcriptomics = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)

# Consider only primary tumour samples
transcriptomics = transcriptomics[[i for i in transcriptomics if i[13:16] == '01A']]

# Overlap
transcriptomics = transcriptomics.loc[:, [i[:12] in set(proteomics) for i in transcriptomics]]
transcriptomics.columns = [c[:12] for c in transcriptomics]

# Average replicates
transcriptomics = DataFrame({i: transcriptomics.loc[:, [i]].mean(1) for i in set(transcriptomics)})

# Normalise
transcriptomics = DataFrame({i: gkn(transcriptomics.ix[i].dropna()).to_dict() for i in transcriptomics.index}).T

# Export
transcriptomics.to_csv('./data/cptac_brca_qc_all_transcriptomics_normalised.csv')
print '[INFO] BRCA tumours transcriptomics: ', './data/cptac_brca_qc_all_transcriptomics_normalised.csv'


# -- CNV
# Import whole CNV data-set
cnv = read_csv('/Users/emanuel/Projects/data/cptac/cna_thresholded.tsv', sep='\t')

# Consider only primary tumour samples
cnv = cnv[[i[13:16] == '01A' for i in cnv['barcode']]]

# Sub-set by proteomics samples
cnv = cnv[[i[:12] in set(proteomics) for i in cnv['barcode']]]

# Build matrix - duplicated entries on same sample are discarded if gistic values differ
cnv = cnv.groupby(['barcode', 'hgnc'])['gistic'].agg(lambda x: np.nan if len(set(x)) > 1 else list(x)[0])
cnv = cnv.reset_index()
cnv = cnv.dropna()
cnv = pivot_table(cnv, index='hgnc', columns='barcode', values='gistic', fill_value=0)

# Parse sample ID
cnv.columns = [i[:12] for i in cnv]

# Export
cnv.to_csv('./data/cptac_brca_qc_all_cnv.csv')
print '[INFO] BRCA tumours copy-number: ', './data/cptac_brca_qc_all_cnv.csv'

