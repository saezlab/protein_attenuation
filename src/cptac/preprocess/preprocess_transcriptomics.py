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

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq.tsv' % wd, sep='\t', index_col=0)


# -- Covariates
samplesheet = Series.from_csv('%s/data/samplesheet.csv' % wd)

design = Series([samplesheet[i] for i in transcriptomics], index=transcriptomics.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age[i] for i in design.index]
print list(design)


# -- PCA
n_components = 10
pca = PCA(n_components=n_components).fit(transcriptomics.dropna().T)
pcs = DataFrame(pca.transform(transcriptomics.dropna().T), index=transcriptomics.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
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
ax.legend(loc=4)
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

plt.savefig('%s/reports/transcriptomics_pca_before_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation'

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(5, 5))
plt.savefig('%s/reports/transcriptomics_clustermap_before_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation'


# -- Normalise pancan data-set
def rm_batch(x, y, covariates=['BRCA', 'COREAD', 'HGSC', 'FEMALE', 'MALE', 'age']):
    ys = y.dropna()
    xs = x.ix[ys.index, covariates]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

transcriptomics = DataFrame({p: rm_batch(design, transcriptomics.ix[p, design.index]) for p in transcriptomics.index}).T
transcriptomics.to_csv('%s/data/tcga_rnaseq_corrected.tsv' % wd, sep='\t')
print '[INFO] Covariates regressed-out'


# -- PCA
n_components = 10
pca = PCA(n_components=10).fit(transcriptomics.dropna().T)
pcs = DataFrame(pca.transform(transcriptomics.dropna().T), index=transcriptomics.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
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

plt.savefig('%s/reports/transcriptomics_pca_after_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation'

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(5, 5))
plt.savefig('%s/reports/transcriptomics_clustermap_after_correction.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation'
