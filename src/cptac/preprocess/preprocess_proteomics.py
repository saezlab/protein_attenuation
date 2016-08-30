#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from cptac import wd, palette
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from sklearn.linear_model.base import LinearRegression
from pandas import read_csv, DataFrame, concat, Series


# -- Imports
# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()

# Proteomics
brca = read_csv('%s/data/brca_proteomics_processed.csv' % wd, index_col=0)
coread = read_csv('%s/data/coread_proteomics_processed.csv' % wd, index_col=0)
hgsc = read_csv('%s/data/hgsc_proteomics_processed.csv' % wd, index_col=0)


# -- Assemble proteomics data-sets
proteomics = concat([brca, hgsc, coread], axis=1)
proteomics.columns = [i[:15] for i in proteomics]

# Merge replicates and dicard poor corelated
remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)
proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})

# Export
proteomics.to_csv('%s/data/cptac_proteomics.csv' % wd)
print proteomics


# -- Samplesheet
samplesheet = Series({s[:15]: n for n, samples in [('BRCA', brca), ('COREAD', coread), ('HGSC', hgsc)] for s in samples}).drop(remove_samples)
samplesheet.to_csv('%s/data/samplesheet.csv' % wd)
print samplesheet


# -- Covariates
design = samplesheet.str.get_dummies()
design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['AGE'] = [clinical_age[i] for i in design.index]
print design


# -- PCA
n_components = 10
pca = PCA(n_components=n_components).fit(proteomics.dropna().T)
pcs = DataFrame(pca.transform(proteomics.dropna().T), index=proteomics.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
pcs = concat([pcs, design], axis=1)
print pca.explained_variance_ratio_

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

plt.savefig('%s/reports/proteomics_pca_before_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation'

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(5, 5))
plt.savefig('%s/reports/proteomics_clustermap_before_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation'


# -- Normalise pancan data-set
def rm_batch(x, y, covariates=['BRCA', 'COREAD', 'HGSC', 'FEMALE', 'MALE', 'AGE']):
    ys = y.dropna()
    xs = x.ix[ys.index, covariates]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

proteomics = DataFrame({p: rm_batch(design, proteomics.ix[p, design.index]) for p in proteomics.index}).T
proteomics.to_csv('%s/data/cptac_proteomics_corrected.csv' % wd)
print '[INFO] Covariates regressed-out'


# -- PCA
n_components = 10
pca = PCA(n_components=10).fit(proteomics.dropna().T)
pcs = DataFrame(pca.transform(proteomics.dropna().T), index=proteomics.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
pcs = concat([pcs, design], axis=1)
print pca.explained_variance_ratio_

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

plt.savefig('%s/reports/proteomics_pca_after_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation'

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(5, 5))
plt.savefig('%s/reports/proteomics_clustermap_after_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation'
