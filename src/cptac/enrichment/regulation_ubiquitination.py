import igraph
import numpy as np
from pandas.stats.misc import zscore
import seaborn as sns
from sklearn.cross_validation import ShuffleSplit
from sklearn.preprocessing.imputation import Imputer
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from yeast_phospho.utilities import pearson, spearman
from scipy.stats.stats import spearmanr
from cptac.utils import gkn
from cptac import wd, palette
from matplotlib.gridspec import GridSpec
from pymist.utils.stringdb import get_stringdb
from pymist.utils.biogriddb import get_biogriddb
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from sklearn.linear_model import LinearRegression, Lasso, LassoCV, Ridge, RidgeCV, ElasticNet, ElasticNetCV
from sklearn.decomposition import FactorAnalysis, PCA
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

# -- Cor
t_p_cor = DataFrame({s: pearson(transcriptomics.ix[proteins, s], proteomics.ix[proteins, s]) for s in samples}, index=['cor', 'pval', 'meas']).T
t_p_cor = t_p_cor[t_p_cor['cor'] > .1]

samples = samples.intersection(t_p_cor.index)

proteomics.ix[proteins, samples].T.corrwith(transcriptomics.ix[proteins, samples].T).hist()


# -- Residuals
# p = 'LAMB1'
def protein_residual(p):
    y = proteomics.ix[p, samples].dropna()
    x = transcriptomics.ix[[p], y.index].dropna().T

    lm = LinearRegression().fit(x, y)

    return y - lm.predict(x)
residuals = DataFrame({p: protein_residual(p) for p in proteins}).T
print 'residuals', residuals.shape


# --
e3_ligases = set(read_csv('%s/files/putative_e3_ligases.csv' % wd).dropna()['gene']).intersection(proteins)

x = proteomics.ix[e3_ligases, samples].T.copy()
y = residuals.drop(e3_ligases, axis=0, errors='ignore').ix[:, x.index].copy().T

x = DataFrame(Imputer().fit_transform(x), index=x.index, columns=x.columns)
y = DataFrame(Imputer().fit_transform(y), index=y.index, columns=y.columns)

lm = LinearRegression().fit(x, y.ix[x.index])
coefs = DataFrame(lm.coef_, index=y.columns, columns=x.columns)
print coefs

df = coefs.unstack().reset_index()
df.columns = ['feature', 'protein', 'coef']
print df.sort('coef')

# print df[df['e3_ligase'] == 'FBXW11'].sort('coef')
# print df.sort('coef')

# df = DataFrame({'coef': coefs.mean(1)})
# df['e3_ligases'] = [int(i in e3_ligases) for i in df.index]

# print pearson(proteomics.ix['TP53', samples], proteomics.ix['UBE3A', samples])
# print pearson(residuals.ix['TP53', samples], proteomics.ix['UBE3A', samples])