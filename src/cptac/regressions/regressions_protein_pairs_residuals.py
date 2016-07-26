import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats.stats import pearsonr
from scipy import stats
from matplotlib.gridspec import GridSpec
from cptac import wd, palette, default_color, palette_cnv_number
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0).dropna()
print cnv

# Residuals
p_pairs = read_csv('%s/tables/top_correlated_protein_pairs_proteomics.csv' % wd)
p_pairs = p_pairs[p_pairs['score'] > .5]
p_pairs = {'%s_%s' % (s, t) for p1, p2 in p_pairs[['p1', 'p2']].values for s, t in [(p1, p2), (p2, p1)]}

uniprot = read_uniprot_genename()
corum = get_complexes_pairs()
corum = {'%s_%s' % (uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}

p_pairs = p_pairs.intersection(corum)
print len(p_pairs)

residuals = read_csv('%s/tables/protein_pairs_residuals.csv' % wd, index_col=0).ix[p_pairs]
print residuals

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics

# -- Overlap
samples = set(cnv).intersection(residuals)
print len(samples)


# -- Regressions
# p_pair, gene, df = 'UQCC2_UQCC1', 'UQCC2', residuals
def regressions(p_pair, gene, df):
    print p_pair, gene

    # Protein measurements
    y = df.ix[p_pair, samples].dropna()
    x = cnv.ix[[gene], y.index].T.dropna()

    if len(set(x[gene])) > 0:
        # Fit models
        lm = LinearRegression().fit(x, y)

        # Predict
        y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))

        # Log likelihood
        l_lm = log_likelihood(y_true, y_pred)

        # # Log likelihood ratio
        # lr = 2 * (l_lm1 - l_lm2)
        # lr_pval = stats.chi2.sf(lr, 1)

        # Effect size
        effect = y[x[gene] == 1].mean() - y[x[gene] != 1].mean()

        # F-statistic
        f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])

        # R-squared
        r = r_squared(y_true, y_pred)

        res = {
            'pair': p_pair, 'protein': gene, 'rsquared': r, 'effect_size': effect, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'meas': x[gene].sum()
        }

        if r > .5:
            print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e' % (res['protein'], res['rsquared'], res['f'], res['f_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res

res = [regressions(p_pair, p, residuals) for p_pair in residuals.index for p in p_pair.split('_')]
res = DataFrame([i for i in res if i])
res['f_adjpval'] = multipletests(res['f_pval'], method='fdr_bh')[1]
print res.sort('f_adjpval')


# # -- Export
# res.sort('f_adjpval').to_csv('%s/tables/regressions_overlap_protein_pairs_cnv.csv' % wd, index=False)
# print '[INFO] Tables exported'


# --
p_pair, protein = 'IARS_LARS', 'LARS'
p1, p2 = p_pair.split('_')

plot_df = DataFrame({p1: proteomics.ix[p1], p2: proteomics.ix[p2], 'cnv': cnv.ix[protein], 'res': residuals.ix[p_pair]}).dropna()

# Scatter
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs = GridSpec(1, 2, hspace=.3)

ax = plt.subplot(gs[0])
for c in [0, -1, 1, 2, -2]:
    sns.regplot(
        x=plot_df.ix[plot_df['cnv'] == c, p1], y=plot_df.ix[plot_df['cnv'] == c, p2], ax=ax,
        fit_reg=False, scatter=True, color=palette_cnv_number[c], scatter_kws={'linewidth': .3, 'edgecolor': 'white'}
    )
sns.regplot(x=plot_df[p1], y=plot_df[p2], ax=ax, fit_reg=True, scatter=False, color=default_color)
ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
ax.set_title('pearson: %.2f, p-value: %.2e' % pearsonr(plot_df[p1], plot_df[p2]))
sns.despine(ax=ax)

# Boxplot
ax = plt.subplot(gs[1])
sns.boxplot(plot_df['cnv'], plot_df['res'], palette=palette_cnv_number, ax=ax, linewidth=.3, sym='')
sns.stripplot(plot_df['cnv'], plot_df['res'], palette=palette_cnv_number, ax=ax, linewidth=.3, size=5, jitter=True, edgecolor='white')
ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
sns.despine(ax=ax, trim=True)

plt.gcf().set_size_inches(6, 3)
plt.savefig('%s/reports/regressions_residuals_scatter.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# # -- Plot
# # Significant associations venn diagram
# sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# fig, gs, pos = plt.figure(figsize=(5, 10)), GridSpec(1, 2, hspace=.3), 0
#
# for c in [-2, 2]:
#     ax = plt.subplot(gs[pos])
#
#     associations = {
#         d: {p for p, t, f in res[d][['protein', 'cnv', 'f_adjpval']].values if c == t and f < .05}
#         for d in res}
#
#     venn3(associations.values(), set_labels=associations.keys(), set_colors=[palette[k] for k in associations])
#     venn3_circles(associations.values(), linestyle='solid', color='white')
#
#     ax.set_title('Depletion' if c == -2 else 'Amplification')
#
#     pos += 1
#
# plt.savefig('%s/reports/regressions_overlap_venn.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Done'
#
#
# # -- QQ-plot
# sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
#
# fig, gs, pos = plt.figure(figsize=(3 * len(res), 3)), GridSpec(1, len(res), hspace=.3), 0
# for d in res:
#     plot_df = DataFrame({
#         'x': sorted([-np.log10(np.float(i) / len(res[d])) for i in np.arange(1, len(res[d]) + 1)]),
#         'y': sorted(-np.log10(res[d]['f_pval']))
#     })
#
#     ax = plt.subplot(gs[pos])
#     g = sns.regplot('x', 'y', plot_df, fit_reg=False, ci=False, color=default_color, line_kws={'lw': .3}, ax=ax)
#     g.set_xlim(0)
#     g.set_ylim(0)
#     plt.plot(plt.xlim(), plt.xlim(), 'k--', lw=.3)
#     plt.xlabel('Theoretical -log(P)')
#     plt.ylabel('Observed -log(P)')
#     plt.title('%s ~ Copy Number' % d)
#     sns.despine(trim=True)
#
#     pos += 1
#
# plt.savefig('%s/reports/regressions_overlap_cnv_qqplot.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Plot done'
