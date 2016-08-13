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
proteins, samples = set(proteomics.index).intersection(transcriptomics.index), set(proteomics).intersection(transcriptomics)
proteomics, transcriptomics = proteomics.ix[proteins, samples], transcriptomics.ix[proteins, samples]
print 'proteins', 'samples', len(proteins), len(samples)


# -- Complexes proteins
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Activity functions
# s, c = 'TCGA-AG-A020-01', 5615
def ztest_complex(s, c, df):
    x1 = df[s].ix[corum[c]].dropna()
    x2 = df[s].drop(corum[c], errors='ignore').dropna()

    if len(x1) > 1:
        z, pvalue = ztest(x1, x2)
        return s, c, corum_n[c], z, pvalue, x1.mean(), len(x1)


# -- Estimate CORUM activity z-test
# Proteomics
c_proteomics = [ztest_complex(s, c, proteomics) for s in proteomics for c in corum]
c_proteomics = DataFrame([i for i in c_proteomics if i], columns=['sample', 'complex', 'name', 'z', 'pval', 'mean', 'targets'])
c_proteomics['fdr'] = multipletests(c_proteomics['pval'],  method='fdr_bh')[1]
c_proteomics.to_csv('%s/tables/protein_complexes_proteomics_activities.tsv' % wd, sep='\t')
print c_proteomics.sort('fdr')

# Transcriptomics
c_transcriptomics = [ztest_complex(s, c, transcriptomics) for s in transcriptomics for c in corum]
c_transcriptomics = DataFrame([i for i in c_transcriptomics if i], columns=['sample', 'complex', 'name', 'z', 'pval', 'mean', 'targets'])
c_transcriptomics['fdr'] = multipletests(c_transcriptomics['pval'],  method='fdr_bh')[1]
c_transcriptomics.to_csv('%s/tables/protein_complexes_transcriptomics_activities.tsv' % wd, sep='\t')
print c_transcriptomics.sort('fdr')


# -- Estimate CORUM activity mean
# Proteomics
c_proteomics_mean = DataFrame({c: proteomics.ix[corum[c]].mean() for c in corum}).T
c_proteomics_mean.to_csv('%s/tables/protein_complexes_proteomics_mean_activities.tsv' % wd, sep='\t')
print c_proteomics_mean.shape

# Transcriptomics
c_transcriptomics_mean = DataFrame({c: transcriptomics.ix[corum[c]].mean() for c in corum}).T
c_transcriptomics_mean.to_csv('%s/tables/protein_complexes_transcriptomics_mean_activities.tsv' % wd, sep='\t')
print c_transcriptomics_mean.shape
