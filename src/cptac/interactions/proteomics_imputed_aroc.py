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

# -- Get overlapping proteins
brca = read_csv('%s/data/brca_proteomics_processed.csv' % wd, index_col=0)
brca = brca[brca.count(1) > (brca.shape[1] * .5)]

hgsc = read_csv('%s/data/hgsc_proteomics_processed.csv' % wd, index_col=0)
hgsc = hgsc[hgsc.count(1) > (hgsc.shape[1] * .5)]

coread = read_csv('%s/data/coread_proteomics_processed.csv' % wd, index_col=0)
coread = coread[coread.count(1) > (coread.shape[1] * .5)]

ov_prot = set(brca.index).intersection(hgsc.index).intersection(coread.index)

# -- Import imputed proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised_imputed.csv' % wd, index_col=0).ix[ov_prot]
proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]
print proteomics.shape

transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0).ix[ov_prot].dropna()
print transcriptomics.shape


# -- Import protein-protein interactions
uniprot = read_uniprot_genename()

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
for name, d_df in [('Proteomics', proteomics), ('Transcriptomics', transcriptomics)]:
    df = d_df.T.corr(method='pearson')
    df.values[np.tril_indices(df.shape[0], 0)] = np.nan
    df.index.name = None
    df = df.unstack().reset_index().dropna()
    df.columns = ['p1', 'p2', 'cor']

    df['CORUM'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['STRING'] = [1 if ((p1, p2) in string) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['BioGRID'] = [1 if ((p1, p2) in biogrid) else 0 for p1, p2 in df[['p1', 'p2']].values]
    df['OmniPath'] = [1 if ((p1, p2) in omnipath) else 0 for p1, p2 in df[['p1', 'p2']].values]

    cor_res[name] = {}
    for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
        curve_fpr, curve_tpr, _ = roc_curve(df[db], df['cor'])
        cor_res[name][db] = (curve_fpr, curve_tpr)

print '[INFO] Done'


# -- Plot
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.set_style({'xtick.direction': 'in', 'ytick.direction': 'in'})

fig, gs, pos = plt.figure(figsize=(13, 3)), GridSpec(1, 4, hspace=.3, wspace=.3), 0

for db in ['CORUM', 'STRING', 'BioGRID', 'OmniPath']:
    ax = plt.subplot(gs[pos])

    for name in ['Proteomics', 'Transcriptomics']:
        curve_fpr, curve_tpr = cor_res[name][db]
        plt.plot(curve_fpr, curve_tpr, label='%s (AUC %0.2f)' % (name, auc(curve_fpr, curve_tpr)), c=palette[name])

    ax.plot([0, 1], [0, 1], 'k--', lw=.3)
    sns.despine(trim=True)
    ax.legend(loc='lower right')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title(db)

    pos += 1

plt.savefig('%s/reports/proteomics_ppi_aroc_imputed.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
