import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, palette_dbs
from pymist.utils.stringdb import get_stringdb
from pymist.utils.biogriddb import get_biogriddb
from pymist.utils.corumdb import get_complexes_pairs
from pandas import DataFrame, Series, read_csv, pivot_table
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()


# -- Import p-p interactions
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


# -- Import normalised pancan proteomics
pancan = read_csv('%s/data/pancan_proteomics_preprocessed_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .50)]
print pancan


# -- Correlation
df = pancan.T.corr(method='pearson')
df.values[np.tril_indices(df.shape[0], 0)] = np.nan
df.index.name = None
df = df.unstack().reset_index().dropna()
df.columns = ['p1', 'p2', 'cor']

df['corum'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in df[['p1', 'p2']].values]
df['string'] = [1 if ((p1, p2) in string) else 0 for p1, p2 in df[['p1', 'p2']].values]
df['biogrid'] = [1 if ((p1, p2) in biogrid) else 0 for p1, p2 in df[['p1', 'p2']].values]
df['omnipath'] = [1 if ((p1, p2) in omnipath) else 0 for p1, p2 in df[['p1', 'p2']].values]

df['score'] = df['cor'].abs()
print df


# -- Histograms
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=.75)

g = sns.distplot(df['cor'], color=palette_dbs['All'], hist=False, kde_kws={'shade': True}, label='All')
sns.distplot(df.ix[df['corum'] == 1, 'cor'], color=palette_dbs['CORUM'], hist=False, kde_kws={'shade': True}, label='CORUM')
sns.distplot(df.ix[df['string'] == 1, 'cor'], color=palette_dbs['STRING'], hist=False, kde_kws={'shade': True}, label='STRING')
sns.distplot(df.ix[df['biogrid'] == 1, 'cor'], color=palette_dbs['BioGRID'], hist=False, kde_kws={'shade': True}, label='BioGRID')
sns.distplot(df.ix[df['omnipath'] == 1, 'cor'], color=palette_dbs['OmniPath'], hist=False, kde_kws={'shade': True}, label='OmniPath')

plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Proteomics pancancer correlation histograms')
plt.xlabel('Pearson\'s r')
plt.ylabel('Density')
g.set_xlim(-1, 1)
sns.despine(trim=True)
plt.gcf().set_size_inches(4, 2)
plt.savefig('%s/reports/pancan_correlation_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Export top correlated protein pairs
df[df['score'] > .4].to_csv('%s/tables/top_correlated_protein_pairs.csv' % wd, index=False)
print '[INFO] Exported'
