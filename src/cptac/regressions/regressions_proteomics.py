import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from cptac.utils import gkn
from cptac import wd, palette_cnv_number, cnv_names
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat, merge


# -- CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
cnv.columns = [i[:15] for i in cnv]
cnv = cnv.applymap(lambda x: 0 if -1 <= x <= 1 else x)
cnv = cnv.loc[:, (cnv != 0).sum() != 0]
print cnv


# -- Proteomics
proteomics = read_csv('%s/data/pancan_proteomics_preprocessed.csv' % wd, index_col=0)
proteomics.columns = [i[:15] for i in proteomics]

remove_samples = {i for i in set(proteomics) if proteomics.loc[:, [i]].shape[1] == 2 and proteomics.loc[:, [i]].corr().ix[0, 1] < .4}
proteomics = proteomics.drop(remove_samples, axis=1)

proteomics = DataFrame({i: proteomics.loc[:, [i]].mean(1) for i in set(proteomics)})
proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]
print proteomics


# -- Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# -- Covariates
samplesheet = Series.from_csv('%s/data/samplesheet.csv' % wd)
samplesheet = {k: v for k, v in samplesheet.to_dict().items()}

design = Series([samplesheet[i] for i in proteomics], index=proteomics.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender[i] for i in design.index}).str.get_dummies()], axis=1)
design['AGE'] = [clinical_age[i] for i in design.index]
print list(design)


# -- Regressions
genes = set(cnv.index).intersection(proteomics.index)
samples = set(cnv).intersection(proteomics)


# p = 'ERBB2'
def regressions(p, c):
    # Protein measurements
    y = proteomics.ix[p, samples].dropna()
    x = sm.add_constant(concat([cnv.ix[p], design], axis=True).ix[y.index].dropna(), has_constant='add')

    y = y.ix[x.index]

    # Discretise
    x.ix[x[p] != c, p] = 0
    x.ix[x[p] == c, p] = 1

    if x[p].sum() >= 5:
        # Linear regression
        lm = sm.OLS(y, x).fit()
        print lm.summary()

        # Wald test
        wald = lm.wald_test(p, use_f=False)

        return p, c, x[p].sum(), lm.params[p], wald.pvalue

    else:
        return np.nan, np.nan, np.nan, np.nan, np.nan

res = DataFrame([regressions(p, c) for p in genes for c in [-2, 2]], columns=['protein', 'type', 'meas', 'coef', 'wald']).dropna()
res['fdr'] = multipletests(res['wald'], method='fdr_bh')[1]
print res.sort_values('fdr')


# -- Export
res.sort('fdr').to_csv('%s/tables/regressions_proteomics_cnv.csv' % wd, index=False)
print '[INFO] Table exported'


# -- QQ-plot
plot_df = DataFrame({
    'x': sorted(-np.log10(stats.uniform.rvs(size=len(res)))),
    'y': sorted(-np.log10(res['wald'].astype(np.float)))
})

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=.75)
g = sns.regplot('x', 'y', plot_df, fit_reg=True, ci=False, color='#34495e', line_kws={'lw': .3})
g.set_xlim(0)
g.set_ylim(0)
plt.xlabel('Theoretical p-value (-log10)')
plt.ylabel('Observed p-value (-log10)')
plt.title('Proteomics ~ Copy Number')
sns.despine(trim=True)
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/regressions_proteomics_cnv_qqplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
