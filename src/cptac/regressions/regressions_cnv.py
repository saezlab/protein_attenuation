import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from cptac import wd, palette
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat, merge


# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv.columns = [i[:15] for i in cnv]
cnv = cnv.applymap(lambda x: 0 if -1 <= x <= 1 else x)
cnv = cnv.loc[:, (cnv != 0).sum() != 0]
print cnv

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print transcriptomics

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
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
    x = sm.add_constant(concat([cnv.ix[p], design], axis=True).ix[y.index].dropna(), has_constant='add')

    y = y.ix[x.index]

    # Discretise
    x.ix[x[p] != c, p] = 0
    x.ix[x[p] == c, p] = 1

    if x[p].sum() >= thres:
        # Linear regression
        lm = sm.OLS(y, x).fit()
        print lm.summary()

        # Wald test
        wald = lm.wald_test(p, use_f=False)

        # Effect size
        effect = y[x[p] == 1].mean() - y[x[p] != 1].mean()

        return p, c, x[p].sum(), lm.params[p], effect, wald.pvalue

    else:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


def regressions_dataset(df, adj_pval='bonferroni'):
    res = DataFrame([regressions(p, c, df) for p in genes for c in [-2, 2]], columns=['protein', 'type', 'meas', 'coef', 'effectsize', 'wald']).dropna()
    res['fdr'] = multipletests(res['wald'], method=adj_pval)[1]
    return res

res = {n: regressions_dataset(df, 'fdr_bh') for n, df in [('Transcriptomics', transcriptomics), ('Proteomics', proteomics), ('Residuals', residuals)]}
print res


# -- Export
for d in res:
    res[d].sort('fdr').to_csv('%s/tables/regressions_overlap_%s_cnv.csv' % (wd, d.lower()), index=False)
print '[INFO] Tables exported'


# -- Plot
# Significant associations venn diagram
sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs, pos = plt.figure(figsize=(5, 10)), GridSpec(1, 2, hspace=.3), 0

for c in [-2, 2]:
    ax = plt.subplot(gs[pos])

    associations = {
        d: {p for p, t, f in res[d][['protein', 'type', 'fdr']].values if c == t and f < .05}
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
        'x': sorted(-np.log10(stats.uniform.rvs(size=len(res[d])))),
        'y': sorted(-np.log10(res[d]['wald'].astype(np.float)))
    })

    ax = plt.subplot(gs[pos])
    g = sns.regplot('x', 'y', plot_df, fit_reg=True, ci=False, color='#34495e', line_kws={'lw': .3}, ax=ax)
    g.set_xlim(0)
    g.set_ylim(0)
    plt.xlabel('Theoretical p-value (-log10)')
    plt.ylabel('Observed p-value (-log10)')
    plt.title('%s ~ Copy Number' % d)
    sns.despine(trim=True)

    pos += 1

plt.savefig('%s/reports/regressions_overlap_cnv_qqplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
