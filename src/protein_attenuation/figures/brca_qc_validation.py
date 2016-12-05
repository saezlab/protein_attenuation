#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture.gaussian_mixture import GaussianMixture
from protein_attenuation import palette_cnv_number
from scipy.stats.stats import ttest_ind
from pandas import read_csv, DataFrame, pivot_table, Series

# -- Imports
# Proteomics
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

# # QC status samples
samples_qc = read_csv('./data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv', index_col=2)
samples_qc.index = [i[:12] for i in samples_qc.index]
# samples_qc = samples_qc[samples_qc['QC Status'] == 'pass']

proteomics.columns = [i[:12] for i in proteomics]
proteomics = proteomics.loc[:, [i in samples_qc.index for i in proteomics]]
print proteomics.shape


# Transcriptomics
# Import TCGA pancancer rna-seq data
transcriptomics = read_csv('/Users/emanuel/Projects/data/cptac/GSE62944_merged_expression_voom.tsv', sep='\t', index_col=0)

# Consider only primary tumour samples
transcriptomics = transcriptomics[[i for i in transcriptomics if i[13:16] == '01A']]

# Overlap
transcriptomics = transcriptomics.loc[:, [i[:12] in set(proteomics) for i in transcriptomics]]
transcriptomics.columns = [c[:12] for c in transcriptomics]

# Average replicates
transcriptomics = DataFrame({i: transcriptomics.loc[:, [i]].mean(1) for i in set(transcriptomics)})
print transcriptomics.shape


# CNV
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
print cnv.shape


# --
samples = set(proteomics).intersection(transcriptomics).intersection(cnv)
proteins = set(proteomics.index).intersection(transcriptomics.index).intersection(cnv.index)
print 'len: ', len(samples)

fail_samples = set(samples_qc[samples_qc['QC Status'] == 'fail'].index).intersection(samples)
pass_samples = set(samples_qc[samples_qc['QC Status'] == 'pass'].index).intersection(samples)
print 'len: ', len(fail_samples)


# -- Attenuated samples
cors = []
for s in samples:
    df = DataFrame({
        'cnv': cnv[s], 'trans': transcriptomics[s], 'prot': proteomics[s]
    }).dropna().corr()

    cors.append({'sample': s, 'cnv_tran': df.ix['cnv', 'trans'], 'cnv_prot': df.ix['cnv', 'prot']})
cors = DataFrame(cors).dropna().set_index('sample')
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
print cors.sort('diff')


# -- GMM attenuation
gmm = GaussianMixture(n_components=2).fit(cors[['diff']])

s_type = Series(dict(zip(*(cors[['diff']].index, gmm.predict(cors[['diff']])))))
clusters = Series(dict(zip(*(range(2), gmm.means_[:, 0]))))

cors['cluster'] = [s_type[i] for i in cors.index]
cors['type'] = ['fail' if i in fail_samples else 'pass' for i in cors.index]
# cors.sort(['cluster', 'diff'], ascending=False).to_csv('./tables/samples_correlations.csv')
# cors = read_csv('./tables/samples_correlations.csv', index_col=0)
print cors.sort(['cluster', 'diff'], ascending=False)

pal = {'pass': palette_cnv_number[2], 'fail': palette_cnv_number[-2]}
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(cors, size=1.25, aspect=2.5)
g = g.map_dataframe(sns.stripplot, 'diff', 'type', orient='h', size=1, jitter=.2, alpha=.4, linewidth=.1, edgecolor='white', palette=pal, split=True)
g = g.map_dataframe(sns.boxplot, 'diff', 'type', orient='h', linewidth=.3, sym='', palette=pal, notch=True)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Attenuation score', 'QC')
plt.title('BRCA tumours attenuation validation')
g.despine(trim=True)
plt.savefig('./reports/protein_attenuation_validation_brca_degradation_samples.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_attenuation_validation_brca_degradation_samples.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Attenuated proteins in tumours
t_cors = read_csv('./tables/proteins_correlations.csv', index_col=0)
p_attenuated = set(t_cors[t_cors['cluster'] == 1].index)

# -- Sample attenutaion score
cors = []
for p in proteins:
    for t, s in [('fail', fail_samples), ('pass', pass_samples)]:
        df = DataFrame({
            'cnv': cnv.ix[p, s], 'trans': transcriptomics.ix[p, s], 'prot': proteomics.ix[p, s]
        }).dropna()

        l = df.shape[0]

        if l > (len(fail_samples) * .5):
            df = df.corr()
            cors.append({'protein': p, 'cnv_tran': df.ix['cnv', 'trans'], 'cnv_prot': df.ix['cnv', 'prot'], 'type': t, 'len': l})

cors = DataFrame(cors).dropna().set_index('protein')
cors['diff'] = cors['cnv_tran'] - cors['cnv_prot']
cors['tumour'] = [
    'Not attenuated' if i not in p_attenuated else ('Attenuated (> %.1f)' % (np.floor((t_cors.ix[i, 'diff'] if t_cors.ix[i, 'diff'] < .5 else .5) * 10) / 10)) for i in cors.index
]
print cors.sort('diff')

# -- Plot
t, pval = ttest_ind(cors[cors['tumour'] != 'Not attenuated']['diff'], cors[cors['tumour'] == 'Not attenuated']['diff'])
print t, pval

order = ['Not attenuated', 'Attenuated (> 0.2)', 'Attenuated (> 0.3)', 'Attenuated (> 0.4)', 'Attenuated (> 0.5)']

pal = {'pass': palette_cnv_number[2], 'fail': palette_cnv_number[-2]}
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(cors, size=1.25, aspect=2.5)
g = g.map_dataframe(sns.stripplot, 'diff', 'tumour', 'type', orient='h', size=1, jitter=.2, alpha=.4, linewidth=.1, edgecolor='white', palette=pal, order=order, split=True)
g = g.map_dataframe(sns.boxplot, 'diff', 'tumour', 'type', orient='h', linewidth=.3, sym='', palette=pal, order=order, notch=True)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Attenuation score', '')
g.add_legend(title='Degradation QC')
plt.title('BRCA tumours protein attenuation validation')
g.despine(trim=True)
plt.savefig('./reports/protein_attenuation_validation_brca_degradation_proteins.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_attenuation_validation_brca_degradation_proteins.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'

