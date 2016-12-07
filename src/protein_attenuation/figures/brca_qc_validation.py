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
proteomics = read_csv('./data/cptac_brca_qc_all_proteomics_normalised.csv', index_col=0)
proteomics = proteomics[proteomics.count(1) > (.5 * proteomics.shape[1])]

transcriptomics = read_csv('./data/cptac_brca_qc_all_transcriptomics_normalised.csv', index_col=0)

cnv = read_csv('./data/cptac_brca_qc_all_cnv.csv', index_col=0)


# -- Overlap
samples = set(proteomics).intersection(transcriptomics).intersection(cnv)
proteins = set(proteomics.index).intersection(transcriptomics.index).intersection(cnv.index)

proteomics, transcriptomics, cnv = proteomics.ix[proteins, samples], transcriptomics.ix[proteins, samples], cnv.ix[proteins, samples]
print 'len: ', len(samples)


# -- Samplesheet
samples_qc = read_csv('./data/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv', index_col=2)
samples_qc.index = [i[:12] for i in samples_qc.index]

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
g = sns.FacetGrid(cors, size=.75, aspect=2)
g = g.map_dataframe(sns.stripplot, 'diff', 'type', orient='h', size=2, jitter=.2, alpha=.4, linewidth=.1, edgecolor='white', palette=pal, split=True)
g = g.map_dataframe(sns.boxplot, 'diff', 'type', orient='h', linewidth=.3, sym='', palette=pal, notch=True)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Attenuation score', 'QC')
plt.title('BRCA tumours protein degradation')
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
plot_df = cors.copy()
plot_df['tumour'] = [i if i == 'Not attenuated' else 'Attenuated' for i in plot_df['tumour']]

order = ['Not attenuated', 'Attenuated']
pal = {'pass': palette_cnv_number[2], 'fail': palette_cnv_number[-2]}
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, size=.75, aspect=2)
g = g.map_dataframe(sns.boxplot, 'diff', 'tumour', 'type', orient='h', linewidth=.3, sym='', palette=pal, order=order, notch=True)
g = g.map(plt.axvline, x=0, ls='-', lw=0.1, c='black', alpha=.5)
g.set_axis_labels('Attenuation score', '')
g.add_legend(title='Degradation QC')
plt.title('BRCA tumours protein attenuation')
g.despine(trim=True)
plt.savefig('./reports/protein_attenuation_validation_brca_degradation_proteins.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/protein_attenuation_validation_brca_degradation_proteins.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
