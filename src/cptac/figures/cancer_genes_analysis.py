#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

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
from cptac import palette, default_color, palette_cnv_number
from cptac.utils import log_likelihood, f_statistic, r_squared
from sklearn.linear_model import LinearRegression
from matplotlib_venn import venn3, venn3_circles
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.stringdb import get_stringdb
from pymist.utils.map_peptide_sequence import read_uniprot_genename, read_fasta


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- Tumour suppressors and oncogenes
cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
cgenes = {
    'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
}
print 'cgenes', len(cgenes)


# -- Protein correlations
p_correlations = read_csv('./tables/proteins_correlations.csv', index_col=0)
print 'p_correlations', len(p_correlations)


# -- Plot scatter of correlations
ax_min = np.min([p_correlations['CNV_Transcriptomics'].min() * 1.10, p_correlations['Transcriptomics_Proteomics'].min() * 1.10])
ax_max = np.min([p_correlations['CNV_Transcriptomics'].max() * 1.10, p_correlations['Transcriptomics_Proteomics'].max() * 1.10])

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'lines.linewidth': .75})
g = sns.jointplot(
    'CNV_Transcriptomics', 'CNV_Proteomics', p_correlations, 'scatter', color='#808080', xlim=[ax_min, ax_max], ylim=[ax_min, ax_max],
    space=0, s=15, edgecolor='w', linewidth=.1, marginal_kws={'hist': False, 'rug': False}, stat_func=None, alpha=.3
)
g.plot_marginals(sns.kdeplot, shade=True, color='#595959', lw=.3)

g.ax_joint.axhline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.axvline(0, ls='-', lw=0.1, c='black', alpha=.3)
g.ax_joint.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', lw=.3)

g.x = p_correlations['CNV_Transcriptomics']
g.y = p_correlations['CNV_Proteomics']
g.plot_joint(sns.kdeplot, cmap=sns.light_palette('#595959', as_cmap=True), legend=False, shade=False, shade_lowest=False, n_levels=9, alpha=.8, lw=.1)

g.x = p_correlations[[i in cgenes['suppressor'] for i in p_correlations.index]]['CNV_Transcriptomics']
g.y = p_correlations[[i in cgenes['suppressor'] for i in p_correlations.index]]['CNV_Proteomics']
g.plot_joint(sns.regplot, fit_reg=False, color=palette_cnv_number[-2])
g.plot_marginals(sns.kdeplot, color=palette_cnv_number[-2], shade=True, legend=False)

g.x = p_correlations[[i in cgenes['oncogene'] for i in p_correlations.index]]['CNV_Transcriptomics']
g.y = p_correlations[[i in cgenes['oncogene'] for i in p_correlations.index]]['CNV_Proteomics']
g.plot_joint(sns.regplot, fit_reg=False, color=palette_cnv_number[2])
g.plot_marginals(sns.kdeplot, color=palette_cnv_number[2], shade=True, legend=False)

plt.gcf().set_size_inches(3, 3)

g.set_axis_labels('CNV ~ Transcriptomics', 'CNV ~ Proteomics')
plt.savefig('./reports/correlation_difference_lmplot_corr_cancer_genes.pdf', bbox_inches='tight')
plt.savefig('./reports/correlation_difference_lmplot_corr_cancer_genes.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Plot done'


# -- Protein complexes interactions
uniprot = read_uniprot_genename()

# string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in get_stringdb(900) for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
# print 'string', len(string)

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print 'corum', len(corum)


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
proteins = {p for (px, py) in corum for p in [px, py]}.intersection(proteomics.index)
print 'samples', len(samples)


# --
# g = 'ZW10'
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in set(transcriptomics.index).intersection(proteomics.index)}).T
print 'residuals', residuals.shape


# -- Tumour suppressors and oncogenes
cgenes = read_csv('./tables/Cancer5000.oncodriveROLE.0.3-0.7.txt', sep='\t').append(read_csv('./tables/HCD.oncodriveROLE.0.3-0.7.txt', sep='\t'))
cgenes = {
    'suppressor': set(cgenes[cgenes['oncodriveROLE'] == 'Loss of function']['SYM']),
    'oncogene': set(cgenes[cgenes['oncodriveROLE'] == 'Activating']['SYM'])
}
print 'cgenes', len(cgenes)

# interactions = {(px, py) for px, py in read_csv('./tables/network-112015_symbol.txt', sep='\t').dropna()[['V4', 'V5']].values if px in cgenes['suppressor'] or px in cgenes['oncogene']}
interactions = {(px, py) for px, py in read_csv('./tables/network-112015_symbol.txt', sep='\t').dropna()[['V4', 'V5']].values if px in cgenes['suppressor']}
print 'interactions', len(interactions)

sources, targets = {px for px, py in interactions}, {py for px, py in interactions}
print 'sources', 'targets', len(sources), len(targets)


# -- Regressions: Py Residuals ~ Px CNV
# px, py = 'SMARCA4', 'ZEB1'
def regressions(px, py):
    if px in cnv.index and py in residuals.index and px != py:
        # Protein measurements
        y = residuals.ix[py, samples].dropna()
        x = cnv.ix[[px], y.index].T.dropna()

        y = y.ix[x.index]

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
            'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
        }

        print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res

ppairs = [regressions(px, py) for px in sources for py in targets]
ppairs = DataFrame([i for i in ppairs if i])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs['type'] = ['suppressor' if i in cgenes['suppressor'] else 'oncogene' for i in ppairs['px']]
print ppairs.sort(['fdr', 'f_pval'])


ppairs[ppairs['fdr'] < .05].groupby('px')['py'].agg(lambda x: set(x)).to_dict()

sns.boxplot(x='type', y='beta', data=ppairs)
sns.boxplot(x='type', y='beta', data=ppairs[ppairs['fdr'] < .05])
