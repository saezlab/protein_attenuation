import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, default_color, palette_cnv_number
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest
from pandas import DataFrame, Series, read_csv, concat
from sklearn.linear_model import LinearRegression
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import data-sets
# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# -- Overlap
genes = set(proteomics.index).intersection(transcriptomics.index)
samples = set(proteomics).intersection(transcriptomics)
print 'genes', 'samples', len(genes), len(samples)


# Regress-out transcriptomics from proteomics
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in genes}).T
print 'residuals', residuals.shape


# -- Complexes proteins
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(genes) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 0}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Estimate CORUM
# s, c = 'TCGA-AG-A020-01', 5615
def ztest_complex(s, c):
    x1 = residuals[s].ix[corum[c]].dropna()
    x2 = residuals[s].drop(corum[c], errors='ignore').dropna()

    if len(x1) > 1:
        z, pvalue = ztest(x1, x2)
        return s, c, corum_n[c], z, pvalue, x1.mean(), len(x1)

c_activity = [ztest_complex(s, c) for s in samples for c in corum]
c_activity = DataFrame([i for i in c_activity if i], columns=['sample', 'complex', 'name', 'z', 'pval', 'mean', 'targets'])
c_activity['fdr'] = multipletests(c_activity['pval'],  method='fdr_bh')[1]
c_activity.to_csv('%s/tables/protein_complexes_activities.tsv' % wd, sep='\t')
print c_activity.sort('mean')
