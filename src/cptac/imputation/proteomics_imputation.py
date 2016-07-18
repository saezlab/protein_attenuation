import igraph
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
from cptac import wd, palette
from matplotlib.gridspec import GridSpec
from pymist.utils.stringdb import get_stringdb
from pymist.utils.biogriddb import get_biogriddb
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Correlated protein pairs
p_pairs = read_csv('%s/tables/top_correlated_protein_pairs.csv' % wd)
p_pairs = p_pairs[p_pairs['score'] > .5]
p_pairs = {(p1, p2) for p1, p2 in p_pairs[['p1', 'p2']].values}
print len(p_pairs)


# -- Create network
network_i = igraph.Graph(directed=False)

# Add nodes
network_i.add_vertices(list({i for p in p_pairs for i in p}))

# Add edges
network_i.add_edges([(p1, p2) for p1, p2 in p_pairs])

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
print DataFrame(Series(dict(zip(*(network_i.vs['name'], network_i.degree())))).sort_values())
print network_i.summary()


# -- Proteomics
proteomics = read_csv('%s/data/pancan_proteomics_preprocessed_normalised.csv' % wd, index_col=0)
annot = Series.from_csv('%s/data/samplesheet.csv' % wd)
print proteomics


# -- Impute
# p = 'AKT1'
imputed = {}
for p in proteomics.index:
    if p in network_i.vs['name']:
        c_proteins = set(network_i.vs[network_i.neighborhood(p)]['name']).difference({p}).intersection(proteomics.index)

        y = proteomics.ix[p].dropna()
        x = proteomics.ix[c_proteins, y.index].T.replace(np.nan, 0)

        y_ = set(proteomics.loc[p, proteomics.ix[p].isnull()].index)
        x_ = proteomics.ix[c_proteins, y_].T.replace(np.nan, 0)

        if len(y_) > 0:
            lm = sm.OLS(y, sm.add_constant(x)).fit()

            y_ = Series(dict(zip(*(y_, lm.predict(sm.add_constant(x_, has_constant='add'))))))

            imputed[p] = concat([y, y_]).to_dict()

        else:
            imputed[p] = proteomics.ix[p].to_dict()

    else:
        imputed[p] = proteomics.ix[p].to_dict()

imputed = DataFrame(imputed).T
print '%.2f' % (float(proteomics.count().sum()) / np.prod(proteomics.shape) * 100)
print '%.2f' % (float(imputed.count().sum()) / np.prod(imputed.shape) * 100)


imputed.to_csv('%s/data/pancan_proteomics_preprocessed_normalised_imputed.csv' % wd)
print '[INFO] Exported'
