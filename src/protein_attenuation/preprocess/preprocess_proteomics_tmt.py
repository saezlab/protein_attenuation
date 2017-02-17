#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from protein_attenuation import palette
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from sklearn.linear_model.base import LinearRegression
from pandas import read_csv, DataFrame, concat, Series


# -- Imports
# Proteomics
brca = read_csv('./data/brca_proteomics_processed.csv', index_col=0)
hgsc = read_csv('./data/hgsc_proteomics_processed.csv', index_col=0)

# Samplesheet
samplesheet = Series.from_csv('./data/samplesheet_phosphoproteomics.csv')


# -- Assemble proteomics data-sets
proteomics = concat([brca, hgsc], axis=1)

# Select only primary tumour samples
proteomics = proteomics[[i for i in proteomics if i[15] == 'A']]

# Simplify barcode
proteomics.columns = [i[:12] for i in proteomics]

# Dicard poor corelated
remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)

# Average replicates
proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})

# Overlap with phosphoproteomics
proteomics = proteomics.loc[:, samplesheet.index].dropna(how='all', axis=1)

# Minimum occorence
proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]

# Export
proteomics.to_csv('./data/cptac_proteomics_tmt.csv')
print '[INFO] Proteomics exported: ' + './data/cptac_proteomics_tmt.csv'


# -- Covariates
# Clinical data
clinical = read_csv('./data/tcga_clinical.csv').dropna(subset=['patient.gender', 'patient.days_to_birth'])
clinical['patient.days_to_birth'] *= -1

samples = set(clinical['patient.bcr_patient_barcode']).intersection(samplesheet.index)

clinical_gender = clinical.groupby('patient.bcr_patient_barcode')['patient.gender'].first().to_dict()
clinical_age = clinical.groupby('patient.bcr_patient_barcode')['patient.days_to_birth'].first().to_dict()

# Design matrix
design = samplesheet.ix[samples].str.get_dummies()
design = concat([design, Series({i: clinical_gender[i] for i in samples}).str.get_dummies()], axis=1)
design['days_to_birth'] = [clinical_age[i] for i in samples]


# -- PCA
n_components = 10
pca = PCA(n_components=n_components).fit(proteomics.dropna().T)
pcs = DataFrame(pca.transform(proteomics.dropna().T), index=proteomics.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
pcs = concat([pcs, design], axis=1)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs = plt.figure(figsize=(7, 3)), GridSpec(1, 2, hspace=.3, wspace=.3)

ax = plt.subplot(gs[0])
for t in set(samplesheet):
    sns.regplot(
        pcs.ix[samplesheet[samplesheet == t].index, 'PC1'],
        pcs.ix[samplesheet[samplesheet == t].index, 'PC2'],
    ax=ax, color=palette[t], fit_reg=False, label=t)

ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.set_xlabel('PC1 (%.1f%%)' % (pca.explained_variance_ratio_[0] * 100))
ax.set_ylabel('PC2 (%.1f%%)' % (pca.explained_variance_ratio_[1] * 100))
ax.legend()
sns.despine(trim=True, ax=ax)

ax = plt.subplot(gs[1])
plot_df = DataFrame(zip(['PC%d' % i for i in range(1, n_components + 1)], pca.explained_variance_ratio_), columns=['PC', 'var'])
plot_df['var'] *= 100
sns.barplot('var', 'PC', data=plot_df, color='gray', linewidth=0, ax=ax)
ax.set_xlabel('Explained variance ratio')
ax.set_ylabel('Principal component')
sns.despine(trim=True, ax=ax)
ax.figure.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f%%'))

fig.suptitle('%s samples - PCA before normalisation' % ', '.join(set(samplesheet)))

plt.savefig('./reports/proteomics_tmt_pca_before_correction.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/proteomics_tmt_pca_before_correction.pdf', bbox_inches='tight')
plt.close('all')

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(5, 5))
plt.savefig('./reports/proteomics_tmt_clustermap_before_correction.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/proteomics_tmt_clustermap_before_correction.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation: ', './reports/proteomics_tmt_pca_before_correction.pdf', './reports/proteomics_tmt_clustermap_before_correction.pdf'


# -- Normalise pancan data-set
def rm_batch(x, y, covariates=['BRCA', 'HGSC', 'female', 'male', 'days_to_birth']):
    ys = y.dropna()
    xs = x.ix[ys.index, covariates]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

proteomics = DataFrame({p: rm_batch(design, proteomics.ix[p, design.index]) for p in proteomics.index}).T
proteomics.to_csv('./data/cptac_proteomics_tmt_corrected.csv')
print '[INFO] Proteomics exported (covariates regressed-out): ' + './data/cptac_proteomics_tmt_corrected.csv'


# -- PCA
n_components = 10
pca = PCA(n_components=10).fit(proteomics.round(9).dropna().T)
pcs = DataFrame(pca.transform(proteomics.round(9).dropna().T), index=proteomics.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
pcs = concat([pcs, design], axis=1)

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs = plt.figure(figsize=(7, 3)), GridSpec(1, 2, hspace=.3, wspace=.3)

ax = plt.subplot(gs[0])
for t in set(samplesheet):
    sns.regplot(
        pcs.ix[samplesheet[samplesheet == t].index, 'PC1'],
        pcs.ix[samplesheet[samplesheet == t].index, 'PC2'],
    ax=ax, color=palette[t], fit_reg=False, label=t)

ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.set_xlabel('PC1 (%.1f%%)' % (pca.explained_variance_ratio_[0] * 100))
ax.set_ylabel('PC2 (%.1f%%)' % (pca.explained_variance_ratio_[1] * 100))
ax.legend()
sns.despine(trim=True, ax=ax)

ax = plt.subplot(gs[1])
plot_df = DataFrame(zip(['PC%d' % i for i in range(1, n_components + 1)], pca.explained_variance_ratio_), columns=['PC', 'var'])
plot_df['var'] *= 100
sns.barplot('var', 'PC', data=plot_df, color='gray', linewidth=0, ax=ax)
ax.set_xlabel('Explained variance ratio')
ax.set_ylabel('Principal component')
sns.despine(trim=True, ax=ax)
ax.figure.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f%%'))

fig.suptitle('%s samples - PCA after normalisation' % ', '.join(set(samplesheet)))

plt.savefig('./reports/proteomics_tmt_pca_after_correction.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/proteomics_tmt_pca_after_correction.pdf', bbox_inches='tight')
plt.close('all')

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(5, 5))
plt.savefig('./reports/proteomics_tmt_clustermap_after_correction.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/proteomics_tmt_clustermap_after_correction.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation: ', './reports/proteomics_tmt_pca_after_correction.pdf', './reports/proteomics_tmt_clustermap_after_correction.pdf'
