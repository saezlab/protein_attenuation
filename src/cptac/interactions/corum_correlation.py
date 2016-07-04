import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from scipy.stats import pearsonr, spearmanr
from pymist.utils.stringdb import get_stringdb
from pymist.utils.corumdb import get_complexes_pairs
from pandas import DataFrame, Series, read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()


# -- Import p-p interactions
# CORUM
corum = get_complexes_pairs()
corum = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in corum for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}

# String
string = get_stringdb(900)
string = {(uniprot[s][0], uniprot[t][0]) for p1, p2 in string for s, t in [(p1, p2), (p2, p1)] if s in uniprot and t in uniprot}


# -- Import normalised pancan proteomics
pancan = read_csv('%s/tables/pancan_normalised.csv' % wd, index_col=0)
pancan = pancan[pancan.count(1) > (pancan.shape[1] * .5)]
print pancan


# -- Correlation
df = pancan.T.corr(method='pearson')
df.values[np.tril_indices(df.shape[0], 0)] = np.nan
df.index.name = None
df = df.unstack().reset_index().dropna()
df.columns = ['p1', 'p2', 'cor']

df['corum'] = [1 if ((p1, p2) in corum) else 0 for p1, p2 in df[['p1', 'p2']].values]
df['string'] = [1 if ((p1, p2) in string) else 0 for p1, p2 in df[['p1', 'p2']].values]
df['score'] = df['cor'].abs()
print df


# -- Histograms
palette = {'All': '#808080', 'CORUM': '#e74c3c', 'STRING': '#34495e'}

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=.75)

g = sns.distplot(df['cor'], color=palette['All'], hist=False, kde_kws={'shade': True}, label='All')
sns.distplot(df.ix[df['corum'] == 1, 'cor'], color=palette['CORUM'], hist=False, kde_kws={'shade': True}, label='CORUM')
sns.distplot(df.ix[df['string'] == 1, 'cor'], color=palette['STRING'], hist=False, kde_kws={'shade': True}, label='STRING')

plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Protein-protein pearson\'s correlation')
plt.xlabel('Pearson\'s r')
plt.ylabel('Density')
g.set_xlim(-1, 1)
sns.despine(trim=True)
plt.gcf().set_size_inches(3, 2)
plt.savefig('%s/reports/protein_correlation_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Export top correlated protein pairs
df[df['score'] > .5].to_csv('%s/tables/top_correlated_protein_pairs.csv' % wd, index=False)
print '[INFO] Exported'
