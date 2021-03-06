#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from sklearn.linear_model import LinearRegression
from pandas import DataFrame, Series, read_csv
from statsmodels.stats.multitest import multipletests
from protein_attenuation.utils import log_likelihood, f_statistic, r_squared, read_uniprot_genename, read_fasta, get_complexes_pairs


# -- Imports
# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)

# Residuals
residuals = read_csv('./data/residuals_protein_transcript.csv', index_col=0)


# -- Overlap
genes = set(proteomics.index).intersection(transcriptomics.index)
samples = set(proteomics).intersection(transcriptomics)
print 'genes', 'samples', len(genes), len(samples)


# -- Protein complexes interactions
uniprot = read_uniprot_genename()
uniprot_fasta = read_fasta()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
corum = {(p1, p2) for p1, p2 in corum if p1 in genes and p2 in genes}
print 'corum', len(corum)


# -- Regressions: Py Residuals ~ Px CNV
# px, py = 'ARID1A', 'DPF2'
def regressions(px, py):
    # Protein measurements
    y = residuals.ix[py, samples].dropna()
    x = proteomics.ix[[px], y.index].T.dropna()

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
        'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'll': l_lm
    }

    print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
    # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

    return res

ppairs = DataFrame([regressions(px, py) for px, py in corum])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs.sort('fdr').to_csv('./tables/ppairs_proteomics_regulation_all.csv', index=False)
print '[INFO] Protein complex associations (Residuals ~ Proteomics): ', './tables/ppairs_proteomics_regulation_all.csv'
