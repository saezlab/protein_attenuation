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
c_activity = pivot_table(c_activity, index='complex', columns='sample', values='z')

uniprot = read_uniprot_genename()
complex_name = get_complexes_name()
complex_proteins = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(cnv.index) for k, v in get_complexes_dict().items()}
print c_activity

# Proteomics
samplesheet = Series.from_csv('%s/data/samplesheet.csv' % wd)
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics


# -- Overlap
genes = set(cnv.index)
complexes = set(c_activity.index)
samples = set(cnv).intersection(c_activity)
print len(genes), len(samples)


# -- Regressions
# protein, p_complex, df = 'SMAD2', 3038, c_activity
def regressions(protein, p_complex, df):
    # Protein measurements
    y = df.ix[p_complex, samples].dropna()
    x = cnv.ix[[protein], y.index].T.dropna()

    if len(set(x[protein])) > 0:
        # Fit models
        lm = LinearRegression().fit(x, y)

        # Predict
        y_pred = Series(dict(zip(*(x.index, lm.predict(x)))))

        # # Log likelihood
        # l_lm = log_likelihood(y_true, y_pred)

        # # Log likelihood ratio
        # lr = 2 * (l_lm1 - l_lm2)
        # lr_pval = stats.chi2.sf(lr, 1)

        # Effect size
        effect = y[x[protein] == 1].mean() - y[x[protein] != 1].mean()

        # F-statistic
        f, f_pval = f_statistic(y, y_pred, len(y), x.shape[1])

        # R-squared
        r = r_squared(y, y_pred)

        res = {
            'protein': protein, 'complex': p_complex, 'rsquared': r, 'effect_size': effect, 'f': f, 'f_pval': f_pval
        }

        print '%s, %s: Rsquared: %.2f, F: %.2f, F pval: %.2e' % (res['protein'], complex_name[p_complex], res['rsquared'], res['f'], res['f_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res

res = [regressions(protein, p_complex, c_activity) for p_complex in complexes for protein in complex_proteins[p_complex]]
res = DataFrame([i for i in res if i])
res['f_adjpval'] = multipletests(res['f_pval'], method='fdr_bh')[1]
res['n_proteins'] = [len(complex_proteins[i].intersection(proteomics.index)) for i in res['complex']]
print res.sort('f_adjpval').head(40)


# --
p_complex, protein = 2217, 'MRE11A'

plot_df = DataFrame({'complex': c_activity.ix[p_complex], 'cnv': cnv.ix[protein], 'type': samplesheet})
plot_df = concat([plot_df, proteomics.ix[complex_proteins[p_complex]].dropna(how='all').T], axis=1).dropna()

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs = GridSpec(1, 3, hspace=.3)

# Scatter
ax = plt.subplot(gs[0])
plot_df_ = plot_df.drop(['cnv', 'complex', 'type'], axis=1).unstack().reset_index()
plot_df_.columns = ['protein', 'sample', 'proteomics']
plot_df_['type'] = [samplesheet.ix[i] for i in plot_df_['sample']]
g = sns.stripplot('sample', 'proteomics', hue='type', data=plot_df_, order=list(plot_df.sort('complex').index), ax=ax, palette=palette, size=3, jitter=True, edgecolor='white', linewidth=.3)
ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
sns.despine(trim=True, bottom=True)
ax.set_xticklabels([])
plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
ax.legend(loc=2)
ax.set_ylabel('Protein abundance')
ax.set_xlabel('Samples (ordered by z-test)')

# Boxplot
ax = plt.subplot(gs[1])
sns.boxplot(plot_df['cnv'], plot_df['complex'], palette=palette_cnv_number, ax=ax, linewidth=.3, sym='')
sns.stripplot(plot_df['cnv'], plot_df['complex'], palette=palette_cnv_number, ax=ax, linewidth=.3, size=5, jitter=True, edgecolor='white')
ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
sns.despine(ax=ax, trim=True)
ax.set_xlabel('%s - Copy number variation' % protein)
ax.set_ylabel('Protein complex (z-test)')

# Heatmap
ax = plt.subplot(gs[2])
pal = sns.diverging_palette(220, 20, n=7, as_cmap=True)
sns.heatmap(plot_df.drop(['cnv', 'complex', 'type'], axis=1).corr(), ax=ax, cmap=pal, center=0, annot=True, fmt='.2f')

plt.suptitle(complex_name[p_complex])
plt.gcf().set_size_inches(12, 4)
plt.savefig('%s/reports/regressions_complex_cnv.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
