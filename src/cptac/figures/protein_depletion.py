#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from cptac import palette, palette_dbs
from cptac.utils import read_gmt
from scipy.stats.stats import ttest_ind
from pandas import read_csv, DataFrame, Series
from pymist.enrichment.gsea import gsea
from scipy.stats.stats import spearmanr
from pymist.utils.stringdb import get_stringdb
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


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


# -- Overlap
samples = set(cnv).intersection(proteomics).intersection(transcriptomics)
genes = set(cnv.index).intersection(transcriptomics.index).intersection(proteomics.index)
print 'samples', 'genes', len(samples), len(genes)


# -- CORUM
uniprot = read_uniprot_genename()

corum = get_complexes_pairs()
corum = {uniprot[g][0] for p in corum for g in p if g in uniprot}.intersection(genes)

corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[g][0] for g in corum_dict[k] if g in uniprot and uniprot[g][0] in genes} for k in corum_dict}

corum_proteins = {p for k in corum_dict for p in corum_dict[k]}

corum_name = get_complexes_name()
print 'corum', len(corum)


# -- GO terms
msigdb_go_bp = read_gmt('./files/c5.bp.v5.1.symbols.gmt')
msigdb_go_cc = read_gmt('./files/c5.cc.v5.1.symbols.gmt')
msigdb_go_mf = read_gmt('./files/c5.mf.v5.1.symbols.gmt')
print 'msigdb_go_mf', 'msigdb_go_cc', 'msigdb_go_bp', len(msigdb_go_mf), len(msigdb_go_cc), len(msigdb_go_bp)


# -- Uniprot PTMs lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('./files/uniprot_ptms_proteins/%s' % f, sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('./files/uniprot_ptms_proteins/') if f.startswith('uniprot_')
}
print 'ptms', len(ptms)


# -- Correlations
res = {}
for g in genes:
    df = DataFrame({'CNV': cnv.ix[g, samples], 'Transcriptomics': transcriptomics.ix[g, samples], 'Proteomics': proteomics.ix[g, samples]}).dropna().corr()

    res[g] = {
        'CNV_Transcriptomics': df.ix['CNV', 'Transcriptomics'],
        'CNV_Proteomics': df.ix['CNV', 'Proteomics'],
        'Transcriptomics_Proteomics': df.ix['Transcriptomics', 'Proteomics']
    }

    print g

res = DataFrame(res).T
res['diff'] = res['CNV_Transcriptomics'] - res['CNV_Proteomics']
print res


# --
dataset, permutations = res['diff'].to_dict(), 10

# dataset, signature, permutations=1000, plot_name=None, plot_title='', y1_label='Enrichment score', y2_label='Data value'
ptm_gsea = DataFrame({k: gsea(dataset, ptms[k], permutations, plot_name='./reports/correlation_difference_ptms_gsea_%s.pdf' % k) for k in ptms}, index=['escore', 'pval']).T

