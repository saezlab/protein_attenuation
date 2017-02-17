#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import re
import numpy as np
from pandas.stats.misc import zscore
from pandas import read_csv, concat

# -- BRCA
brca = read_csv('./data/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv', sep='\t')

# Map p-site
brca['psite'] = ['%s_%s' % (gene, psite.split(':')[1].upper()) for psite, gene in brca[['Phosphosite', 'Gene']].values]
brca = brca.set_index('psite')

# Selected log ratios
brca = brca[[c for c in brca if c.endswith(' Log Ratio')]]
brca.columns = [c.split(' ')[0] for c in brca]

# Match sample ids
brca_samplesheet = read_csv('./data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1.csv')
brca_samplesheet['id'] = ['-'.join(i.split('-')[1:4]) for i in brca_samplesheet['TCGA barcode']]
brca_samplesheet = brca_samplesheet.set_index('id')

brca = brca.loc[:, [i in brca_samplesheet.index for i in brca]]
brca.columns = [brca_samplesheet.ix[i, 'TCGA barcode'] for i in brca]

# QC status samples
samples_qc = read_csv('./data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv', index_col=2)
samples_qc = samples_qc[samples_qc['QC Status'] == 'pass']

brca = brca.loc[:, [i in samples_qc.index for i in brca]]

# Average p-site log ratios
brca = brca.reset_index().groupby('psite').mean()

# Export
brca.to_csv('./data/brca_phosphoproteomics_processed.csv')
print '[INFO] BRCA phosphoproteomics preprocessed exported' + './data/brca_phosphoproteomics_processed.csv'


# -- HGSC
hgsc = read_csv('./data/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv', sep='\t')

# Map p-site
hgsc['psite'] = ['%s_%s' % (gene, psite.split(':')[1].upper()) for psite, gene in hgsc[['Phosphosite', 'Gene']].values]
hgsc = hgsc.set_index('psite')

# Selected log ratios
hgsc = hgsc[[c for c in hgsc if c.endswith(' Log Ratio')]]
hgsc.columns = [c.split(' ')[0] for c in hgsc]

# Match sample ids
hgsc_samplesheet = read_csv('./data/OV_All_clinical_features_TCGAbiotab_CPTAC_S020.csv')
hgsc_samplesheet['id'] = ['-'.join(i.split('-')[1:4]) for i in hgsc_samplesheet['TCGA barcode']]
hgsc_samplesheet = hgsc_samplesheet.groupby('id')['TCGA barcode'].first()

hgsc = hgsc.loc[:, [i in hgsc_samplesheet.index for i in hgsc]]
hgsc.columns = [hgsc_samplesheet.ix[i] for i in hgsc]

# Average p-site log ratios
hgsc = hgsc.reset_index().groupby('psite').mean()

# Export
hgsc.to_csv('./data/hgsc_phosphoproteomics_processed.csv')
print '[INFO] HGSC phosphoproteomics preprocessed exported' + './data/hgsc_phosphoproteomics_processed.csv'

