import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from cptac import wd, palette
from pandas.stats.misc import zscore
from matplotlib.gridspec import GridSpec
from sklearn.decomposition.pca import PCA
from sklearn.linear_model.base import LinearRegression
from pandas import read_csv, DataFrame, concat, Series
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()

# Proteomics
brca = read_csv('%s/data/brca_proteomics_processed.csv' % wd, index_col=0)
coread = read_csv('%s/data/coread_proteomics_processed.csv' % wd, index_col=0)
hgsc = read_csv('%s/data/hgsc_proteomics_processed.csv' % wd, index_col=0)

# Samplesheet
samplesheet = concat([
    DataFrame(zip(*(list(brca), np.repeat('BRCA', len(list(brca))))), columns=['code', 'type']),
    DataFrame(zip(*(list(coread), np.repeat('COREAD', len(list(coread))))), columns=['code', 'type']),
    DataFrame(zip(*(list(hgsc), np.repeat('HGSC', len(list(hgsc))))), columns=['code', 'type'])
]).set_index('code')
samplesheet.to_csv('%s/data/samplesheet.csv' % wd)

# Concatenate all
pancan = concat([brca, hgsc, coread], axis=1)
pancan.to_csv('%s/data/pancan_preprocessed.csv' % wd)
print '[INFO] Exported'


# -- Covariates
design = Series(np.concatenate([
    np.repeat('brca', brca.shape[1]),
    np.repeat('hgsc', hgsc.shape[1]),
    np.repeat('coread', coread.shape[1])
]), index=pancan.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender['-'.join(i.split('-')[:4])[:-1].upper()].lower() for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age['-'.join(i.split('-')[:4])[:-1].upper()] for i in design.index]
design['tmt'] = np.bitwise_or(design['brca'], design['hgsc'])
design['shotgun'] = design['coread']
print list(design)


# -- PCA
n_components = 10
pca = PCA(n_components=n_components).fit(pancan.dropna().T)
pcs = DataFrame(pca.transform(pancan.dropna().T), index=pancan.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
pcs = concat([pcs, design], axis=1)
print pca.explained_variance_ratio_

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs = plt.figure(figsize=(7, 3)), GridSpec(1, 2, hspace=.3, wspace=.3)

ax = plt.subplot(gs[0])
for t in set(samplesheet['type']):
    sns.regplot(
        pcs.ix[samplesheet[samplesheet['type'] == t].index, 'PC1'],
        pcs.ix[samplesheet[samplesheet['type'] == t].index, 'PC2'],
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

fig.suptitle('%s samples - PCA before normalisation' % ', '.join(set(samplesheet['type'])))

plt.savefig('%s/reports/pancan_pca_before_normalisation.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation'

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(3, 3))
plt.savefig('%s/reports/pancan_clustermap_before_normalisation.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA before normalisation'


# -- Normalise pancan data-set
def rm_batch(x, y, covariates=['brca', 'coread', 'hgsc', 'female', 'male', 'age', 'tmt', 'shotgun']):
    ys = y.dropna()
    xs = x.ix[ys.index, covariates]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

pancan = DataFrame({p: rm_batch(design, pancan.ix[p, design.index]) for p in pancan.index}).T
print '[INFO] Covariates regressed-out'


# -- PCA
n_components = 10
pca = PCA(n_components=n_components).fit(pancan.dropna().T)
pcs = DataFrame(pca.transform(pancan.dropna().T), index=pancan.columns, columns=['PC%d' % i for i in range(1, n_components + 1)])
pcs = concat([pcs, design], axis=1)
print pca.explained_variance_ratio_

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs = plt.figure(figsize=(7, 3)), GridSpec(1, 2, hspace=.3, wspace=.3)

ax = plt.subplot(gs[0])
for t in set(samplesheet['type']):
    sns.regplot(
        pcs.ix[samplesheet[samplesheet['type'] == t].index, 'PC1'],
        pcs.ix[samplesheet[samplesheet['type'] == t].index, 'PC2'],
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

fig.suptitle('%s samples - PCA after normalisation' % ', '.join(set(samplesheet['type'])))

plt.savefig('%s/reports/pancan_pca_after_normalisation.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation'

sns.set(style='white', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(pcs.corr(), linewidths=.3, cmap=sns.diverging_palette(220, 10, n=20, as_cmap=True), figsize=(3, 3))
plt.savefig('%s/reports/pancan_clustermap_after_normalisation.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] PCA after normalisation'


# -- Export
pancan.to_csv('%s/data/pancan_preprocessed_normalised.csv' % wd)
print '[INFO] Exported'
