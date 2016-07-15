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
proteomics = read_csv('%s/data/pancan_preprocessed_normalised.csv' % wd, index_col=0)
annot = read_csv('%s/data/samplesheet.csv' % wd, index_col=0)
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

imputed.to_csv('%s/data/pancan_preprocessed_normalised_imputed.csv' % wd)
print '[INFO] Exported'


# --
uniprot = read_uniprot_genename()

imputed = read_csv('%s/data/pancan_preprocessed_normalised_imputed.csv' % wd, index_col=0).dropna()
print imputed.shape

# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print len(corum)

# String
string = get_stringdb(900)
string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}
print len(string)

# Biogrid
biogrid = get_biogriddb()
biogrid = {(s, t) for p1, p2 in biogrid for s, t in [(p1, p2), (p2, p1)]}
print len(biogrid)

# Omnipath
omnipath = read_csv('%s/files/omnipathdb.txt' % wd, sep='\t')
omnipath = {(uniprot[s][0], uniprot[t][0]) for s, t in omnipath[['source', 'target']].values if s in uniprot and t in uniprot}
print len(omnipath)


# -- Protein-protein correlation
cor_res = {}
for name, d_df in [('Pancancer', imputed)]:
    df = d_df.T.corr(method='pearson')
    df.values[np.tril_indices(df.shape[0], 0)] = np.nan
    df.index.name = None
    df = df.unstack().reset_index().dropna()
    df.columns = ['p1', 'p2', 'cor']

    df['CORUM'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['STRING'] = [1 if ((p1, p2) in string) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['BioGRID'] = [1 if ((p1, p2) in biogrid) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['OmniPath'] = [1 if ((p1, p2) in omnipath) else 0 for p1, p2 in df[['p1', 'p2']].values]

    df['score'] = df['cor'].abs()

    cor_res[name] = {}
    for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
        curve_fpr, curve_tpr, _ = roc_curve(df[db], df['score'])
        cor_res[name][db] = (curve_fpr, curve_tpr)

print '[INFO] Done'


# -- Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs, pos = plt.figure(figsize=(13, 3)), GridSpec(1, 4, hspace=.3, wspace=.3), 0

for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
    ax = plt.subplot(gs[pos])

    for name in ['Pancancer']:
        curve_fpr, curve_tpr = cor_res[name][db]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title(db)

    pos += 1

plt.savefig('%s/reports/pancan_aroc_imputed.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
