import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from pandas import DataFrame, Series, read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict


# -- Import associations
ppairs = read_csv('%s/tables/ppairs_cnv_regulation.csv' % wd)
ppairs = ppairs[ppairs['cor'] > 0][['px', 'py', 'fdr', 'cor']]
proteins = set(ppairs['px']).union(ppairs['py'])
print ppairs.sort('fdr')


# -- Complexes proteins
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- ppairs clusters
clusters = read_csv('%s/tables/ppairs_signif_complexes_cluster.csv' % wd, index_col=0)
cmap = dict(zip(*(set(clusters['cluster']), sns.color_palette('Paired', 10).as_hex())))
clusters['group'] = [cmap[i] for i in clusters['cluster']]


# --
ppairs['corum'] = [[c for c in corum if px in corum[c] and py in corum[c]] for px, py in ppairs[['px', 'py']].values]
pairs_complexes = DataFrame({'%s_%s' % (p1x, p1y): {'%s_%s' % (p2x, p2y): int(len(set(c1).intersection(c2)) > 0) for p2x, p2y, c2 in ppairs[['px', 'py', 'corum']].values} for p1x, p1y, c1 in ppairs[['px', 'py', 'corum']].values})

cmap = sns.light_palette('#e74c3c', n_colors=2, as_cmap=True)
sns.set(style='white', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3}, font_scale=0.75)
sns.clustermap(pairs_complexes, cmap=cmap, figsize=(8, 8), linewidths=.3, col_colors=clusters['group'])
plt.savefig('%s/reports/ppairs_signif_complexes_clustermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# --
ppairs['log10(fdr)'] = -np.log10(ppairs['fdr'])
ppairs['pair'] = ['%s_%s' % (px, py) for px, py in ppairs[['px', 'py']].values]
ppairs['group'] = [clusters.ix[i, 'group'] for i in ppairs['pair']]
ppairs['cluster'] = [clusters.ix[i, 'cluster'] for i in ppairs['pair']]

clusters_legend = {c: '; '.join([corum_n[i] for i in set.intersection(*map(set, ppairs.ix[ppairs['cluster'] == c, 'corum'].values))]) for c in set(ppairs['cluster'])}
ppairs['cluster complexes'] = [clusters_legend[i] for i in ppairs['cluster']]

sns.set(style='white', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .0, 'xtick.direction': 'out', 'ytick.direction': 'out'}, font_scale=0.75)
g = sns.PairGrid(ppairs, x_vars=['log10(fdr)', 'cor'], y_vars=['pair'])
g.map(sns.barplot, linewidth=0, palette=ppairs['group'])

for ax, title in zip(g.axes.flat, ['Association FDR (-log10)\nPx (CNV) ~ Py (Residuals)', 'Pearson\nPx (Proteomics) ~ Py (Proteomics)']):

    # Set a different title for each axes
    ax.set(title=title)

    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

sns.despine(left=True, bottom=True)

plt.gcf().set_size_inches(4, 8)
plt.savefig('%s/reports/ppairs_signif_complexes_barplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

# --
palette = ppairs.groupby('cluster complexes')['group'].first().to_dict()
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'}, font_scale=0.75)
sns.boxplot(x='log10(fdr)', y='cluster complexes', data=ppairs, orient='h', sym='', linewidth=.3, palette=palette)
sns.swarmplot(x='log10(fdr)', y='cluster complexes', data=ppairs, orient='h', linewidth=.3, edgecolor='white', palette=palette, size=3)
plt.axvline(-np.log10(0.01), ls='-', lw=0.3, c='gray', alpha=.3, label='FDR < 1%')
plt.title('Protein-pairs significantly associated (FDR < 5%)\nPy (residuals) ~ Px (CNV)')
plt.legend(loc=4)
sns.despine(trim=True)
plt.gcf().set_size_inches(2, 3)
plt.savefig('%s/reports/ppairs_signif_complexes_boxplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# --
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=0.75)
sns.lmplot(x='log10(fdr)', y='cor', hue='cluster complexes', data=ppairs, palette=palette, legend=False, fit_reg=False)
sns.despine()
plt.ylabel('Pearson\nPx (Proteomics) ~ Py (Proteomics)')
plt.xlabel('Association FDR (-log10)\nPx (CNV) ~ Py (Residuals)')
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/ppairs_signif_complexes_scatter.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

# --
ppairs[['px', 'py', 'fdr', 'cor', 'cluster', 'cluster complexes']].to_csv('%s/tables/ppairs_signif_associations_annotated.txt' % wd, sep='\t', index=False)
print '[INFO] Exported'
