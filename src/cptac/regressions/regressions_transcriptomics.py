import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from scipy.stats.stats import pearsonr, ttest_ind
from matplotlib.gridspec import GridSpec
from statsmodels.stats.weightstats import ztest
from cptac import wd, palette, default_color, palette_cnv_number
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta


# -- Imports
# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

# -- Overlap
genes = set(proteomics.index).intersection(transcriptomics.index)
samples = set(proteomics).intersection(transcriptomics)
print 'genes', 'samples', len(genes), len(samples)


# -- Protein complexes interactions
uniprot = read_uniprot_genename()
uniprot_fasta = read_fasta()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
corum = {(p1, p2) for p1, p2 in corum if p1 in genes and p2 in genes}
print 'corum', len(corum)


# --
# g = 'ZW10'
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in genes}).T
print 'residuals', residuals.shape


# -- Regressions: Py Residuals ~ Px CNV
# px, py = 'ARID1A', 'DPF2'
def regressions(px, py):
    # Protein measurements
    y = residuals.ix[py, samples].dropna()
    x = transcriptomics.ix[[px], y.index].T

    # Fit models
    lm = LinearRegression().fit(x, y)

    # Predict
    y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))

    # Log likelihood
    l_lm = log_likelihood(y_true, y_pred)

    # F-statistic
    f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])

    # R-squared
    r = r_squared(y_true, y_pred)

    res = {
        'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm
    }

    print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
    # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

    return res

ppairs = DataFrame([regressions(px, py) for px, py in corum])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs.to_csv('%s/tables/ppairs_transcriptomics_regulation_all.csv' % wd, index=False)
# ppairs = read_csv('%s/tables/ppairs_transcriptomics_regulation_all.csv' % wd)
print ppairs.sort('fdr')


# -- Protein sequence length
# ppairs = ppairs[(ppairs['cor'] > 0)]
ppairs = ppairs[ppairs['fdr'] < .05]

p_length = DataFrame([{'uniprot': p, 'name': uniprot[p][0], 'length': len(uniprot_fasta[p])} for p in uniprot_fasta if p in uniprot]).groupby('name')['length'].max().to_dict()

plot_df = DataFrame([{'type': t, 'protein': p, 'length': p_length[p]} for px, py in ppairs[['px', 'py']].values for t, p in [('px', px), ('py', py)] if p in p_length])
print plot_df

ttest, pval = ttest_ind(plot_df[plot_df['type'] == 'px']['length'], plot_df[plot_df['type'] == 'py']['length'])
print 'ttest', 'pval', ttest, pval

# Boxplot
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
sns.violinplot(x='type', y='length', data=plot_df, linewidth=.3, cut=0, inner='quartile', split=False, color=palette_cnv_number[0])
sns.stripplot(x='type', y='length', data=plot_df, linewidth=.3, jitter=True, edgecolor='white', split=False, color=default_color)
plt.ylim(0)
sns.despine(trim=True)
plt.ylabel('Protein sequence length (number of AA)')
plt.title('Px (CNV) ~ Py (Residuals)\nT-test: %.2f, p-value: %.2e' % (ttest, pval))
plt.gcf().set_size_inches(2, 4)
plt.savefig('%s/reports/protein_pairs_protein_info_transcriptomics_boxplots.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'

# Histogram
plot_df = DataFrame([{'px': px, 'py': py, 'len_px': p_length[px], 'len_py': p_length[py], 'signif': int(fdr < .05)} for px, py, fdr in ppairs[['px', 'py', 'fdr']].values if px in p_length and py in p_length])
plot_df['diff'] = plot_df['len_px'] - plot_df['len_py']
print plot_df.sort('diff')

z, zpval = ttest_ind(plot_df.loc[plot_df['signif'] == 1, 'diff'].values, plot_df.loc[plot_df['signif'] == 0, 'diff'].values, equal_var=False)

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
for i in [0, 1]:
    values = plot_df[plot_df['signif'] == i]['diff']
    sns.distplot(values, hist=False, kde_kws={'shade': True}, color=palette_cnv_number[i], label='%s (mean = %.2f)' % ('Significant' if i else 'All', np.mean(values)))
plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Protein length difference\nT-test: %.2f, %.2e' % (z, zpval))
plt.xlabel('len(Px) - len(Py)')
sns.despine(trim=True)
plt.legend()
plt.gcf().set_size_inches(4, 2)
plt.savefig('%s/reports/protein_pairs_protein_info_transcriptomics_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
