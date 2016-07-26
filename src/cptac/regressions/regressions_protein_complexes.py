import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from matplotlib.gridspec import GridSpec
from cptac import wd, palette, default_color, palette_cnv_number
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
print cnv

# Protein complexes activities
c_activity = read_csv('%s/tables/protein_complexes_activities.tsv' % wd, sep='\t', index_col=0)
c_activity = pivot_table(c_activity, index='complex', columns='sample', values='mean')

uniprot = read_uniprot_genename()
complex_name = get_complexes_name()
complex_proteins = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(cnv.index) for k, v in get_complexes_dict().items()}
print c_activity

# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t').dropna(subset=['VITAL_STATUS', 'DAYS_TO_LAST_FOLLOWUP'])
clinical['VITAL_STATUS'] = [0 if i == 'Alive' else 1 for i in clinical['VITAL_STATUS']]
clinical = clinical[clinical['DAYS_TO_LAST_FOLLOWUP'] > 1]
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# -- Overlap
genes = set(cnv.index)
complexes = set(c_activity.index)
samples = set(cnv).intersection(c_activity).intersection(clinical['SAMPLE_brcID'])
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
# protein, p_complex, c, df, thres = 'ERBB2', 100, 2, c_activity, 5
def regressions(protein, p_complex, c, df, thres=5):
    # Protein measurements
    y = df.ix[p_complex, samples].dropna()
    x = concat([cnv.ix[protein]], axis=True).ix[y.index].dropna()

    y = y.ix[x.index]

    # Discretise
    x.ix[x[protein] != c, protein] = 0
    x.ix[x[protein] == c, protein] = 1

    if x[protein].sum() >= thres:
        # Fit models
        lm = LinearRegression().fit(x, y)

        # Predict
        y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))

        # # Log likelihood
        # l_lm = log_likelihood(y_true, y_pred)

        # # Log likelihood ratio
        # lr = 2 * (l_lm1 - l_lm2)
        # lr_pval = stats.chi2.sf(lr, 1)

        # Effect size
        effect = y[x[protein] == 1].mean() - y[x[protein] != 1].mean()

        # F-statistic
        f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])

        # R-squared
        r = r_squared(y_true, y_pred)

        res = {
            'protein': protein, 'complex': p_complex, 'cnv': c, 'rsquared': r, 'effect_size': effect, 'f': f, 'f_pval': f_pval
        }

        print '%s, %s: Rsquared: %.2f, F: %.2f, F pval: %.2e' % (res['protein'], complex_name[p_complex], res['rsquared'], res['f'], res['f_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res

res = [regressions(protein, p_complex, c, c_activity) for p_complex in complexes for protein in complex_proteins[p_complex] for c in [-2, 2]]
res = DataFrame([i for i in res if i])
res['f_adjpval'] = multipletests(res['f_pval'], method='fdr_bh')[1]
print res.sort('f_adjpval')


# -- QQ-plot
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})

plot_df = DataFrame({
    'x': sorted([-np.log10(np.float(i) / len(res)) for i in np.arange(1, len(res) + 1)]),
    'y': sorted(-np.log10(res['f_pval']))
})

g = sns.regplot('x', 'y', plot_df, fit_reg=False, ci=False, color=default_color, line_kws={'lw': .3})
g.set_xlim(0)
g.set_ylim(0)
plt.plot(plt.xlim(), plt.xlim(), 'k--', lw=.3)
plt.xlabel('Theoretical -log(P)')
plt.ylabel('Observed -log(P)')
plt.title('Protein Complex ~ Copy Number')
sns.despine(trim=True)

plt.gcf().set_size_inches(5, 5)
plt.savefig('%s/reports/regressions_overlap_protein_complexes_qqplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


