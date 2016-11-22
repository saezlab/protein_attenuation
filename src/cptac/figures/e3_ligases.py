#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from igraph import plot
from cptac import palette, palette_cnv_number, default_color
from matplotlib.gridspec import GridSpec
from scipy.stats.stats import pearsonr
from matplotlib_venn import venn2, venn2_circles
from scipy.stats import chi2
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from cptac.utils import log_likelihood, f_statistic, r_squared
from pandas import DataFrame, Series, read_csv, concat

# -- Import data-sets
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape

# Residuals
residuals = read_csv('./data/residuals_protein_transcript.csv', index_col=0)
print 'residuals', residuals.shape


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
print 'samples', len(samples)


# -- Improt regression results
ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = {(px, py) for px, py in ppairs_trans[ppairs_trans['fdr'] < .05][['px', 'py']].values}
print len(ppairs_trans)

ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
ppairs_cnv = ppairs_cnv[[(px, py) in ppairs_trans for px, py in ppairs_cnv[['px', 'py']].values]]
print ppairs_cnv.sort('fdr')


# -- Residuals
# g = 'ZW10'
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in set(transcriptomics.index).intersection(proteomics.index)}).T
print 'residuals', residuals.shape


# -- Regressions: Py Residuals ~ Px CNV + Pe3
e3_ligases = ['HUWE1', 'UBR1', 'UBR2', 'UBR3', 'UBR5', 'RNF126', 'CNOT4', 'LTN1', 'ZNF598', 'UBE4B', 'UBE4A', 'SYVN1', 'STUB1']


# px, py, pe3 = 'RPA2', 'RPA3', 'HUWE1'
def regressions(px, py, pe3):
    if py in residuals.index and px in cnv.index and pe3 in transcriptomics.index:
        # -- 1st model
        # Protein measurements
        y = residuals.ix[py, samples].dropna()
        x = cnv.ix[[px], y.index].T

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

        res_1 = {
            'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta': lm.coef_[0]
        }
        print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res_1['rsquared'], res_1['f'], res_1['f_pval'])

        # -- 2nd model
        # Protein measurements
        y = residuals.ix[py, samples].dropna()
        x = concat([cnv.ix[[px], y.index].T, transcriptomics.ix[[pe3], y.index].T], axis=1).ix[y.index]

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

        res_2 = {
            'px': px, 'py': py, 'pe3': pe3, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm, 'beta_px': lm.coef_[0], 'beta_pe3': lm.coef_[1]
        }
        print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res_2['rsquared'], res_2['f'], res_2['f_pval'])

        # -- Likelihood ratio test
        lr = 2 * (res_2['ll'] - res_1['ll'])
        p_val = chi2.pdf(lr, 1)

        res = {'px': px, 'py': py, 'pe3': pe3, 'beta_px': lm.coef_[0], 'beta_pe3': lm.coef_[1], 'lr': lr, 'pval': p_val}

        return res

ppairs = [regressions(px, py, pe3) for px, py in ppairs_cnv[['px', 'py']].values for pe3 in e3_ligases]
ppairs = DataFrame([i for i in ppairs if i])
ppairs['fdr'] = multipletests(ppairs['pval'], method='fdr_bh')[1]
print ppairs[ppairs['beta_pe3'] < 0].sort('fdr')


# -- Scatter
def ppair_correlation(px, py):
    x, y = zip(*proteomics.ix[[px, py]].T.dropna().values)
    return pearsonr(x, y)

plot_df = ppairs[(ppairs['beta_pe3'] < 0) & (ppairs['fdr'] < .05)]
plot_df['cor'] = [ppair_correlation(px, py)[0] for px, py in plot_df[['pe3', 'py']].values]

# px, py = 'COG3', 'COG2'
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
gs, pos = GridSpec(len(plot_df), 2, hspace=.5, wspace=.5), 0
for px, py in plot_df.sort('cor', ascending=False)[['pe3', 'py']].values:
    #
    ax = plt.subplot(gs[pos])

    df = DataFrame({
        'y': residuals.ix[py, samples],
        'x': transcriptomics.ix[px, samples]
    }).dropna()

    sns.regplot(df['x'], df['y'], ax=ax, color=default_color, fit_reg=True, scatter=True, truncate=True, line_kws={'linewidth': .3})

    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('E3: %s (proteomics)' % px)
    ax.set_ylabel('Py: %s (residuals)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(df['x'], df['y']))
    ax.set_ylim(df['y'].min() * 1.05, df['y'].max() * 1.05)

    #
    ax = plt.subplot(gs[pos + 1])

    df = DataFrame({
        'y': proteomics.ix[py, samples],
        'x': transcriptomics.ix[px, samples]
    }).dropna()

    sns.regplot(df['x'], df['y'], ax=ax, color=default_color, fit_reg=True, scatter=True, truncate=True, line_kws={'linewidth': .3})

    sns.despine(ax=ax)
    ax.axhline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
    ax.set_xlabel('Py: %s (proteomics)' % px)
    ax.set_ylabel('E3: %s (proteomics)' % py)
    ax.set_title('Pearson\'s r: %.2f, p-value: %.2e' % pearsonr(df['x'], df['y']))
    ax.set_ylim(df['y'].min() * 1.05, df['y'].max() * 1.05)

    pos += 2

plt.gcf().set_size_inches(6, 3 * len(plot_df))
plt.savefig('./reports/regressions_associations_scatter_e3ligases.png', bbox_inches='tight', dpi=300)
plt.savefig('./reports/regressions_associations_scatter_e3ligases.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
