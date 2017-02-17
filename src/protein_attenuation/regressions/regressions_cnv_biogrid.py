#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from pandas import DataFrame, Series, read_csv
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from protein_attenuation.utils import f_statistic, r_squared, read_uniprot_genename, read_fasta, get_complexes_pairs


# -- Imports
# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq_corrected_normalised.csv', index_col=0)

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)

# Residuals
residuals = read_csv('./data/residuals_protein_transcript.csv', index_col=0)


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)


# -- CORUM
uniprot = read_uniprot_genename()
uniprot_fasta = read_fasta()

corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}


# -- BioGRID
biogrid = read_csv('./files/BIOGRID-ALL-3.4.145.tab2.txt', sep='\t')

biogrid = biogrid[(biogrid['Organism Interactor A'] == 9606) & (biogrid['Organism Interactor B'] == 9606) & (biogrid['Experimental System Type'] == 'physical')]

biogrid_exp = {'Affinity Capture-MS', 'Affinity Capture-Western', 'Co-crystal Structure', 'Co-purification', 'Reconstituted Complex'}
biogrid = biogrid[[i in biogrid_exp for i in biogrid['Experimental System']]]

biogrid = {(s, t) for p1, p2 in biogrid[['Official Symbol Interactor A', 'Official Symbol Interactor B']].values for s, t in [(p1, p2), (p2, p1)] if s != t}


# -- Concat interactions set
interactions = corum.union(biogrid)
print len(interactions)


# -- Regressions: Py Residuals ~ Px CNV
# px, py = 'ARID1A', 'DPF2'
def regressions(px, py):
    if py in residuals.index and px in cnv.index:
        # Protein measurements
        y = residuals.ix[py, samples].dropna()
        x = cnv.ix[[px], y.index].T

        # Fit models
        lm = LinearRegression().fit(x, y)

        # Predict
        y_true, y_pred = y.copy(), Series(dict(zip(*(x.index, lm.predict(x)))))

        # F-statistic
        f, f_pval = f_statistic(y_true, y_pred, len(y), x.shape[1])

        # R-squared
        r = r_squared(y_true, y_pred)

        res = {
            'px': px, 'py': py, 'rsquared': r, 'f': f, 'f_pval': f_pval, 'beta': lm.coef_[0]
        }

        # print 'Px (%s), Py (%s): Rsquared: %.2f, F: %.2f, F pval: %.2e' % (px, py, res['rsquared'], res['f'], res['f_pval'])
        # print sm.OLS(y, sm.add_constant(x, has_constant='add')).fit().summary()

        return res

ppairs = [regressions(px, py) for px, py in interactions]
ppairs = DataFrame([i for i in ppairs if i])
ppairs['fdr'] = multipletests(ppairs['f_pval'], method='fdr_bh')[1]
ppairs.sort('fdr').to_csv('./tables/ppairs_cnv_regulation_all_tmt.csv', index=False)
print '[INFO] Protein complex associations (Residuals ~ Copy-number): ', './tables/ppairs_cnv_regulation_all.csv'
