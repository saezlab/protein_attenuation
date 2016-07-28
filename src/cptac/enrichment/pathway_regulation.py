import igraph
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
from cptac import wd, palette, default_color
from cptac.utils import read_gmt
from scipy.stats.stats import pearsonr
from matplotlib.gridspec import GridSpec
from pymist.utils.stringdb import get_stringdb
from pymist.utils.biogriddb import get_biogriddb
from sklearn.metrics.ranking import roc_curve, auc
from pymist.utils.corumdb import get_complexes_pairs
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import data-sets
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0).dropna()
print 'proteomics', proteomics.shape

transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0).ix[proteomics.index].dropna()
print 'transcriptomics', transcriptomics.shape

# -- Import protein-protein interactions
uniprot = read_uniprot_genename()

# CORUM
corum = set()
for p1, p2 in get_complexes_pairs():
    if (p2, p1) not in corum:
        corum.add((p1, p2))
corum = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in corum if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}

# String
string = set()
for p1, p2 in get_stringdb(900):
    if (p2, p1) not in string:
        string.add((p1, p2))
string = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in string if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}

p_pairs = corum.union(string)
print 'p_pairs', len(p_pairs)

# -- Import KEGG pathways
msigdb_kegg = read_gmt('%s/files/c2.cp.kegg.v5.1.symbols.gmt' % wd)
# msigdb_kegg = read_gmt('%s/files/h.all.v5.1.symbols.gmt' % wd)
print 'msigdb_kegg', len(msigdb_kegg)

# -- Overlap
samples = set(proteomics).intersection(transcriptomics)
proteins = set(proteomics.index).intersection(transcriptomics.index)
print 'samples', 'proteins', len(samples), len(proteins)


# --
# p = 'KEGG_GLYCOLYSIS_GLUCONEOGENESIS'
# p1, p2 = 'LDHB', 'LDHA'
cor_df = []
for p in msigdb_kegg:
    for p1, p2 in p_pairs:
        if p1 in proteins and p2 in proteins:
            if p1 in msigdb_kegg[p] and p2 in msigdb_kegg[p]:
                samples = list(concat([proteomics.ix[[p1, p2]].T, transcriptomics.ix[[p1, p2]].T], axis=1).dropna().index)

                for t, df in [('Proteomics', proteomics), ('Transcriptomics', transcriptomics)]:
                    cor, pval = pearsonr(df.ix[p1, samples], df.ix[p2, samples])
                    cor_df.append({'pathway': p, 'p1': p1, 'p2': p2, 'cor': cor, 'pval': pval, 'type': t})

cor_df = DataFrame(cor_df)
cor_df['name'] = [i[5:].replace('_', ' ').lower() for i in cor_df['pathway']]
print cor_df


# --
pathway_pairs = Series(dict(zip(*(np.unique(cor_df['name'], return_counts=True)))))
pathway_pairs = set(pathway_pairs[pathway_pairs > 6].index)

plot_df = cor_df[[i in pathway_pairs for i in cor_df['name']]]

order = Series({p: cor_df.loc[(cor_df['name'] == p) & (cor_df['type'] == 'Proteomics'), 'cor'].mean() - cor_df.loc[(cor_df['name'] == p) & (cor_df['type'] == 'Transcriptomics'), 'cor'].mean() for p in set(plot_df['name'])}).dropna().sort_values()

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
g = sns.boxplot('cor', 'name', 'type', palette=palette, data=cor_df, linewidth=.3, order=order.index, orient='h', sym='')
sns.stripplot('cor', 'name', 'type', palette=palette, data=cor_df, linewidth=.3, order=order.index, orient='h', jitter=True, edgecolor='white', size=3, split=True)
plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
sns.despine(trim=True)
plt.gcf().set_size_inches(5, 15)
plt.savefig('%s/reports/pathway_regulation.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

