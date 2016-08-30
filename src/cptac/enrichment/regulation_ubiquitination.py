#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests
from cptac import wd, palette, default_color
from sklearn.metrics.ranking import roc_curve, auc
from sklearn.linear_model import LinearRegression
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# -- Overlap
proteins, samples = set(proteomics.index).intersection(transcriptomics.index), set(proteomics).intersection(transcriptomics)
print 'proteins', 'samples', len(proteins), len(samples)

# -- Gene-sets
uniprot = read_uniprot_genename()

# Uniprot PTMs lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('%s/files/uniprot_ptms_proteins/%s' % (wd, f), sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('%s/files/uniprot_ptms_proteins/' % wd) if f.startswith('uniprot_')
}
print 'ptms', len(ptms)


# -- Residuals
# p = 'LAMB1'
def protein_residual(p):
    y = proteomics.ix[p, samples].dropna()
    x = transcriptomics.ix[[p], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[p] - lm.intercept_

    return y_
residuals = DataFrame({p: protein_residual(p) for p in proteins}).T
print 'residuals', residuals.shape


# --
e3_ligases = set(read_csv('%s/files/putative_e3_ligases.csv' % wd).dropna()['gene']).intersection(proteins)
print 'e3_ligases', len(e3_ligases)


# ligase, protein = 'PPIL2', 'LAMB3'
def cor_protein(ligase, protein):
    x, y = zip(*DataFrame({'x': proteomics.ix[ligase, samples], 'y': residuals.ix[protein, samples]}).dropna().values)

    cor, pval = pearsonr(x, y)
    res = {'protein': protein, 'ligase': ligase, 'cor': cor, 'pval': pval}

    return res

cor_df = DataFrame([cor_protein(ligase, protein) for protein in proteins.difference(e3_ligases) for ligase in e3_ligases])
cor_df['fdr'] = multipletests(cor_df['pval'], method='fdr_bh')[1]
cor_df['ubq'] = [int(i in ptms['ubl_conjugation']) for i in cor_df['protein']]
print cor_df.sort('cor')


# -- Enrichment
(float(cor_df[(cor_df['cor'] < 0) & (cor_df['fdr'] < .05)]['ubq'].sum()) / cor_df['ubq'].sum())

(float(cor_df[(cor_df['cor'] > 0) & (cor_df['fdr'] < .05)]['ubq'].sum()) / cor_df['ubq'].sum())

curve_fpr, curve_tpr, _ = roc_curve(cor_df['ubq'], cor_df['cor'])

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
plt.plot(curve_fpr, curve_tpr, label='AUC %0.2f' % auc(curve_fpr, curve_tpr), c=default_color)
plt.plot([0, 1], [0, 1], 'k--', lw=.3)
sns.despine(trim=True)
plt.legend(loc='lower right')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.gcf().set_size_inches(4, 4)
plt.savefig('%s/reports/protein_ubiquitination_aroc.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
