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

ppairs = [regressions(px, py) for px, py in corum if px in cgenes['suppressor'] or px in cgenes['oncogene']]
ppairs = DataFrame([i for i in ppairs if i])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs['type'] = ['suppressor' if i in cgenes['suppressor'] else 'oncogene' for i in ppairs['px']]
print ppairs.sort(['fdr', 'f_pval'])


ppairs[ppairs['fdr'] < .05].groupby('px')['py'].agg(lambda x: set(x)).to_dict()

sns.boxplot(x='type', y='beta', data=ppairs)
sns.boxplot(x='type', y='beta', data=ppairs[ppairs['fdr'] < .05])
