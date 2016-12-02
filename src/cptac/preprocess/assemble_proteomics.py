#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from pandas.stats.misc import zscore
from pandas import DataFrame, Series, read_csv, concat


# -- COREAD
coread = read_csv('%s/data/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv' % wd, index_col=0, sep='\t').replace(0, np.nan)
coread = coread.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)

# Select only unambigous mapping peptides
coread = coread[[c for c in coread if c.endswith(' Unshared Area')]]
coread.columns = [c.split(' ')[0] for c in coread]

# Correct for total ion count
coread = coread / coread.sum()

# Log2
coread = np.log2(coread)

# Centre
coread = zscore(coread)

# Match sample ids
coread_samplesheet = read_csv('%s/data/CPTAC_TCGA_ColorectalCancer_select_clinical_data_release1_090413.csv' % wd)
coread_samplesheet['id'] = ['-'.join(i.split('-')[1:5]) for i in coread_samplesheet['TCGA barcode']]
coread_samplesheet = coread_samplesheet.set_index('id')

coread.columns = [coread_samplesheet.ix[i, 'TCGA barcode'] for i in coread]

# Export
coread.to_csv('%s/data/coread_proteomics_processed.csv' % wd)
print '[INFO] COREAD proteomics preprocessed exported'


# -- BRCA
brca = read_csv('%s/data/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv' % wd, sep='\t', index_col=0)

brca = brca.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)
brca = brca.drop(['Mean', 'Median', 'StdDev'], axis=0)

brca = brca[[c for c in brca if c.endswith(' Unshared Log Ratio')]]
brca.columns = [c.split(' ')[0] for c in brca]

# Match sample ids
brca_samplesheet = read_csv('%s/data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1.csv' % wd)
brca_samplesheet['id'] = ['-'.join(i.split('-')[1:4]) for i in brca_samplesheet['TCGA barcode']]
brca_samplesheet = brca_samplesheet.set_index('id')

brca = brca.loc[:, [i in brca_samplesheet.index for i in brca]]
brca.columns = [brca_samplesheet.ix[i, 'TCGA barcode'] for i in brca]

# QC status samples
samples_qc = read_csv('%s/data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv' % wd, index_col=2)
samples_qc = samples_qc[samples_qc['QC Status'] == 'pass']

brca = brca.loc[:, [i in samples_qc.index for i in brca]]

att_s = read_csv('./tables/samples_correlations.csv', index_col=0)
att_s['cluster'] = np.logical_not(att_s['cluster']).astype(int)

plot_df = brca.corr()
plot_df.index = ['Normal' if samples_qc.ix[i, 'QC Status'] == 'pass' else 'Degradation' for i in plot_df.index]
plot_df.columns = ['-' if i[:12] not in att_s.index else ('Attenuated' if att_s.ix[i[:12], 'cluster'] else 'NOT Attenuated') for i in plot_df]

sns.set(font_scale=.5)
sns.clustermap(plot_df, figsize=(15, 15))
plt.savefig('./reports/brca_proteasome_degradation.png', bbox_inches='tight', dpi=300)
plt.close('all')

# Export
brca.to_csv('%s/data/brca_proteomics_processed.csv' % wd)
print '[INFO] BRCA proteomics preprocessed exported'


# -- HGSC
# JHU
hgsc_jhu = read_csv('%s/data/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv' % wd, sep='\t', index_col=0)

hgsc_jhu = hgsc_jhu.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)
hgsc_jhu = hgsc_jhu.drop(['Mean', 'Median', 'StdDev'], axis=0)

hgsc_jhu = hgsc_jhu[[c for c in hgsc_jhu if c.endswith(' Unshared Log Ratio')]]
hgsc_jhu.columns = [c.split(' ')[0] for c in hgsc_jhu]

# PNNL
hgsc_pnnl = read_csv('%s/data/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv' % wd, sep='\t', index_col=0)

hgsc_pnnl = hgsc_pnnl.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)
hgsc_pnnl = hgsc_pnnl.drop(['Mean', 'Median', 'StdDev'], axis=0)

hgsc_pnnl = hgsc_pnnl[[c for c in hgsc_pnnl if c.endswith(' Unshared Log Ratio')]]
hgsc_pnnl.columns = [c.split(' ')[0] for c in hgsc_pnnl]

# Merge two data-sets
hgsc = concat([hgsc_jhu, hgsc_pnnl], axis=1)

# Match sample ids
hgsc_samplesheet = read_csv('%s/data/OV_All_clinical_features_TCGAbiotab_CPTAC_S020.csv' % wd)
hgsc_samplesheet['id'] = ['-'.join(i.split('-')[1:4]) for i in hgsc_samplesheet['TCGA barcode']]
hgsc_samplesheet = hgsc_samplesheet.groupby('id')['TCGA barcode'].first()

hgsc = hgsc.loc[:, [i in hgsc_samplesheet.index for i in hgsc]]
hgsc.columns = [hgsc_samplesheet.ix[i] for i in hgsc]

# Export
hgsc.to_csv('%s/data/hgsc_proteomics_processed.csv' % wd)
print '[INFO] HGSC proteomics preprocessed exported'
