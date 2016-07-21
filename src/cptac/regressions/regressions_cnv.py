import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from matplotlib.gridspec import GridSpec
from cptac import wd, palette, default_color
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat


# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv.columns = [i[:15] for i in cnv]
cnv = cnv.applymap(lambda x: 0 if -1 <= x <= 1 else x)
cnv = cnv.loc[:, (cnv != 0).sum() != 0]
print cnv

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected.csv' % wd, index_col=0)
print transcriptomics

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected.csv' % wd, index_col=0)
print proteomics

# Residuals
residuals = read_csv('%s/tables/protein_residuals.csv' % wd, index_col=0)
print residuals

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t').dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# -- Overlap
genes = set(cnv.index).intersection(proteomics.index).intersection(transcriptomics.index).intersection(residuals.index)
samples = set(cnv).intersection(proteomics).intersection(transcriptomics).intersection(residuals).intersection(clinical['SAMPLE_brcID'])
print len(genes), len(samples)


# -- Covariates
samplesheet = Series.from_csv('%s/data/samplesheet.csv' % wd)
samplesheet = {k: v for k, v in samplesheet.to_dict().items()}

design = Series([samplesheet[i] for i in samples], index=samples)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['AGE'] = [clinical_age[i] for i in design.index]
print list(design)


# -- Regressions
# p, c, df, thres = 'ERBB2', 2, transcriptomics, 5
def regressions(p, c, df, thres=5):
    # Protein measurements
    y = df.ix[p, samples].dropna()
    x = concat([cnv.ix[p]], axis=True).ix[y.index].dropna()

    y = y.ix[x.index]

    # Discretise
    x.ix[x[p] != c, p] = 0
    x.ix[x[p] == c, p] = 1

    if x[p].sum() >= thres:
        # Fit models
        lm1 = LinearRegression().fit(x, y)
        lm2 = LinearRegression().fit(np.zeros((len(x), 1)), y)

        # Predict
        y_true_lm1, y_pred_lm1 = y.copy(), Series(dict(zip(*(x.index, lm1.predict(x)))))
        y_true_lm2, y_pred_lm2 = y.copy(), Series(dict(zip(*(x.index, lm2.predict(np.zeros((len(x), 1)))))))

        # np.zeros((len(x), 1))
        # x.drop(p, axis=1)

        # Log likelihood
        l_lm1 = log_likelihood(y_true_lm1, y_pred_lm1)
        l_lm2 = log_likelihood(y_true_lm1, y_pred_lm2)

        # Log likelihood ratio
        lr = 2 * (l_lm1 - l_lm2)
        lr_pval = stats.chi2.sf(lr, 1)

        # Effect size
        effect = y[x[p] == 1].mean() - y[x[p] != 1].mean()

        # F-statistic
        f, f_pval = f_statistic(y_true_lm1, y_pred_lm1, len(y), x.shape[1])

        # R-squared
        r = r_squared(y_true_lm1, y_pred_lm1)

        res = {
            'protein': p, 'cnv': c, 'rsquared': r, 'effect_size': effect, 'f': f, 'f_pval': f_pval, 'll': l_lm1, 'll_small': l_lm2, 'll_ratio': lr, 'lr_pval': lr_pval
        }

        print '%s: Rsquared: %.2f, F: %.2f, F pval: %.2e, ll: %.2f, lr pval: %.2e' % (res['protein'], res['rsquared'], res['f'], res['f_pval'], res['ll'], res['lr_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res


def regressions_dataset(df, adj_pval='bonferroni'):
    res = [regressions(p, c, df) for p in genes for c in [-2, 2]]
    res = DataFrame([i for i in res if i])

    res['f_adjpval'] = multipletests(res['f_pval'], method=adj_pval)[1]
    res['lr_adjpval'] = multipletests(res['lr_pval'], method=adj_pval)[1]

    return res

res = {n: regressions_dataset(df, 'bonferroni') for n, df in [('Transcriptomics', transcriptomics), ('Proteomics', proteomics), ('Residuals', residuals)]}
print res


# -- Export
for d in res:
    res[d].sort('f_adjpval').to_csv('%s/tables/regressions_overlap_%s_cnv.csv' % (wd, d.lower()), index=False)
print '[INFO] Tables exported'


# -- Plot
# Significant associations venn diagram
sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs, pos = plt.figure(figsize=(5, 10)), GridSpec(1, 2, hspace=.3), 0

for c in [-2, 2]:
    ax = plt.subplot(gs[pos])

    associations = {
        d: {p for p, t, f in res[d][['protein', 'cnv', 'f_adjpval']].values if c == t and f < .05}
        for d in res}

    venn3(associations.values(), set_labels=associations.keys(), set_colors=[palette[k] for k in associations])
    venn3_circles(associations.values(), linestyle='solid', color='white')

    ax.set_title('Depletion' if c == -2 else 'Amplification')

    pos += 1

plt.savefig('%s/reports/regressions_overlap_venn.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- QQ-plot
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

fig, gs, pos = plt.figure(figsize=(3 * len(res), 3)), GridSpec(1, len(res), hspace=.3), 0
for d in res:
    plot_df = DataFrame({
        'x': sorted([-np.log10(np.float(i) / len(res[d])) for i in np.arange(1, len(res[d]) + 1)]),
        'y': sorted(-np.log10(res[d]['f_pval']))
    })

    ax = plt.subplot(gs[pos])
    g = sns.regplot('x', 'y', plot_df, fit_reg=False, ci=False, color=default_color, line_kws={'lw': .3}, ax=ax)
    g.set_xlim(0)
    g.set_ylim(0)
    plt.plot(plt.xlim(), plt.xlim(), 'k--', lw=.3)
    plt.xlabel('Theoretical -log(P)')
    plt.ylabel('Observed -log(P)')
    plt.title('%s ~ Copy Number' % d)
    sns.despine(trim=True)

    pos += 1

plt.savefig('%s/reports/regressions_overlap_cnv_qqplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
