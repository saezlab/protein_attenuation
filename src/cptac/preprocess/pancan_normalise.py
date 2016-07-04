import itertools as it
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from gdsc import wd as gdsc_wd
from mtkirc.utils import gkn
from pandas.stats.misc import zscore
from sklearn.linear_model.base import LinearRegression
from pandas import read_csv, DataFrame, concat, Series
from scipy.stats.stats import spearmanr
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from statsmodels.stats.multitest import multipletests
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()


# -- Import CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}


# -- Imports
# Clinical data
clinical = read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')
clinical_gender = clinical.groupby('SAMPLE_brcID')['GENDER'].first().to_dict()
clinical_age = clinical.groupby('SAMPLE_brcID')['AGE'].first().to_dict()


# Proteomics
brca = read_csv('%s/tables/brca_proteomics_processed.csv' % wd, index_col=0)
coread = read_csv('%s/tables/coread_proteomics_processed.csv' % wd, index_col=0)
hgsc = read_csv('%s/tables/hgsc_proteomics_processed.csv' % wd, index_col=0)


# Concatenate all
pancan = concat([brca, coread, hgsc], axis=1)


# -- Normalise pancan data-set
design = Series(np.concatenate([
    np.repeat('brca', brca.shape[1]),
    np.repeat('coread', coread.shape[1]),
    np.repeat('hgsc', hgsc.shape[1])
]), index=pancan.columns)
design = design.str.get_dummies()

design = concat([design, Series({i: clinical_gender['-'.join(i.split('-')[:4])[:-1].upper()] for i in design.index}).str.get_dummies()], axis=1)
design['age'] = [clinical_age['-'.join(i.split('-')[:4])[:-1].upper()] for i in design.index]


def rm_batch(x, y):
    ys = y.dropna()
    xs = x.ix[ys.index]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

pancan_n = DataFrame({p: rm_batch(design, pancan.ix[p, design.index]) for p in pancan.index}).T
pancan_n = zscore(pancan_n)
print pancan_n


# -- Export
pancan_n.to_csv('%s/table/pancan_normalised.csv' % wd)
print '[INFO] Exported'
