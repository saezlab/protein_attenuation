import os
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pandas.stats.misc import zscore
from scipy.stats.stats import pearsonr
from matplotlib_venn import venn2, venn2_circles
from cptac.utils import hypergeom_test, read_gmt
from cptac import wd, palette
from statsmodels.stats.multitest import multipletests
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.corumdb import get_complexes_pairs
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import protein-pairs
uniprot = read_uniprot_genename()
uniprot_proteins = {v[0] for v in uniprot.values()}

# CORUM
corum = set()
for p1, p2 in get_complexes_pairs():
    if (p2, p1) not in corum:
        corum.add((p1, p2))
corum = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in corum if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}
print 'corum', len(corum)

# Overlap
p_pairs = corum
p_pairs_proteins = {p for p1, p2 in p_pairs for p in [p1, p2]}
print 'p_pairs', len(p_pairs)


# -- Import data-sets
# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

# Overlap
proteins = set(transcriptomics.index).intersection(proteomics.index)
samples = set(transcriptomics).intersection(proteomics)
print len(proteins), len(samples)

# -- Turnover rates
turnover = read_csv('%s/files/proteins_turnovers_preprocessed.csv' % wd).dropna(subset=['Uniprot IDs human'])
turnover = DataFrame([{'protein': i, 'p_halflife': r, 't_halflife': t} for p, r, t in turnover[['Uniprot IDs human', 'Protein half-life average [h]', 'mRNA half-life average [h]']].values for i in p.split(';')])
turnover['protein'] = [uniprot[i][0] for i in turnover['protein']]
turnover = turnover.groupby('protein').max().dropna()

# -- Trans vs Prot
# p = 'SPR'
def pt_correlation(p):
    x = proteomics.ix[p, samples].dropna()
    y = transcriptomics.ix[p, x.index]
    cor, pval = pearsonr(x, y)
    return {'cor': cor, 'pval': pval}
pt_cor = DataFrame({p: pt_correlation(p) for p in proteins}).T
pt_cor['fdr'] = multipletests(pt_cor['pval'], method='fdr_bh')[1]
print pt_cor


# -- Import gene-sets
# Uniprot PTMs lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('%s/files/uniprot_ptms_proteins/%s' % (wd, f), sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('%s/files/uniprot_ptms_proteins/' % wd) if f.startswith('uniprot_')
}
print 'ptms', len(ptms)

# GO terms
msigdb_go_bp = read_gmt('%s/files/c5.bp.v5.1.symbols.gmt' % wd)
msigdb_go_cc = read_gmt('%s/files/c5.cc.v5.1.symbols.gmt' % wd)
msigdb_go_mf = read_gmt('%s/files/c5.mf.v5.1.symbols.gmt' % wd)
print 'msigdb_go_mf', 'msigdb_go_cc', 'msigdb_go_bp', len(msigdb_go_mf), len(msigdb_go_cc), len(msigdb_go_bp)

# Pathways
msigdb_cp = read_gmt('%s/files/c2.cp.v5.1.symbols.gmt' % wd)
msigdb_kegg = read_gmt('%s/files/c2.cp.kegg.v5.1.symbols.gmt' % wd)
msigdb_cgp = read_gmt('%s/files/c2.cgp.v5.1.symbols.gmt' % wd)
print 'msigdb_cp', 'msigdb_kegg', 'msigdb_cgp', len(msigdb_cp), len(msigdb_kegg), len(msigdb_cgp)


# -- Protein pairs correlation
# px, py = 'SDHA', 'SDHB'
def protein_correlation(px, py):
    if px in proteins and py in proteins:
        samples = list(concat([proteomics.ix[[px, py]].T, transcriptomics.ix[[px, py]].T], axis=1).dropna().index)

        p_cor, p_pval = pearsonr(proteomics.ix[px, samples], proteomics.ix[py, samples])
        t_cor, t_pval = pearsonr(transcriptomics.ix[px, samples], transcriptomics.ix[py, samples])

        return {
            'px': px, 'py': py,
            'p_cor': p_cor, 'p_pval': p_pval,
            't_cor': t_cor, 't_pval': t_pval,
            'len': len(samples)
        }

cor_df = {'%s_%s' % (px, py): protein_correlation(px, py) for px, py in p_pairs}
cor_df = DataFrame({i: cor_df[i] for i in cor_df if cor_df[i]}).T

cor_df['diff'] = [p_cor - t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]
cor_df['sum'] = [p_cor + t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]

cor_df['t_fdr'] = multipletests(cor_df['t_pval'], method='fdr_bh')[1]
cor_df['p_fdr'] = multipletests(cor_df['p_pval'], method='fdr_bh')[1]
print cor_df.sort('p_fdr')


# --
bkg = set(pt_cor[pt_cor['fdr'] > .05].index).intersection(p_pairs_proteins)

p_pairs = {p for pp in cor_df[(cor_df['p_fdr'] < .05) & (cor_df['p_cor'] > .0) & (cor_df['t_fdr'] > .05)].index for p in pp.split('_')}.intersection(bkg)
t_pairs = {p for pp in cor_df[(cor_df['t_fdr'] < .05) & (cor_df['t_cor'] > .0) & (cor_df['p_fdr'] > .05)].index for p in pp.split('_')}.intersection(bkg)
print 'bkg', 'p_pairs', 't_pairs', len(bkg), len(p_pairs), len(t_pairs)


def enrichment_hypergeom(signature_type, signature):
    # Hypergeometric
    p_pval, p_len = hypergeom_test(signature.intersection(bkg), bkg, p_pairs)
    t_pval, t_len = hypergeom_test(signature.intersection(bkg), bkg, t_pairs)

    s_len = len(signature.intersection(bkg))

    return {
        'type': signature_type,
        'p_pval': p_pval, 'p_len': p_len, 'p_perc': float(p_len) / s_len * 100,
        't_pval': t_pval, 't_len': t_len, 't_perc': float(t_len) / s_len * 100,
        's_len': s_len
    }

hyper = DataFrame({
    signature_name: enrichment_hypergeom(signature_type, signatures[signature_name]) for signature_type, signatures in [('cc', msigdb_go_cc), ('bp', msigdb_go_bp)]
    for signature_name in signatures if len(signatures[signature_name].intersection(bkg)) > 0
}).T.dropna()
hyper['p_fdr'] = multipletests(hyper['p_pval'], method='fdr_bh')[1]
hyper['t_fdr'] = multipletests(hyper['t_pval'], method='fdr_bh')[1]
print hyper.sort(['p_len', 'p_perc'], ascending=False)

plot_df = hyper.loc[(hyper[['t_fdr', 'p_fdr']] < .05).sum(1) == 1, ['p_perc', 't_perc']].astype(float).T
plot_df.columns = [i.lower().replace('_', ' ') for i in plot_df]
plot_df.index = ['Proteomics' if i.startswith('p_') else 'Transcriptomics' for i in plot_df.index]

sns.clustermap(plot_df, cmap=sns.diverging_palette(220, 20, n=7, as_cmap=True), row_colors=[palette[i] for i in plot_df.index], center=0)
plt.gcf().set_size_inches(35, 3)
plt.savefig('%s/reports/protein_pairs_goterms_clustermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


plot_df = {
    'Transcriptomics': t_pairs,
    'Proteomics': p_pairs
}

order = ['Transcriptomics', 'Proteomics']
venn2([plot_df[i] for i in order], set_labels=order, set_colors=[palette[k] for k in order])
venn2_circles([plot_df[i] for i in order], linestyle='solid', color='white')
plt.savefig('%s/reports/protein_pairs_correlation_venn.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Half-life
f_halflife = turnover.dropna().copy()

p_lb, p_ub = np.percentile(f_halflife['p_halflife'], 25), np.percentile(f_halflife['p_halflife'], 75)
t_lb, t_ub = np.percentile(f_halflife['t_halflife'], 25), np.percentile(f_halflife['t_halflife'], 75)
f_halflife['p_halflife'] = ['Unstable' if i < p_lb else ('Stable' if i > p_ub else np.nan) for i in f_halflife['p_halflife']]
f_halflife['t_halflife'] = ['Unstable' if i < t_lb else ('Stable' if i > t_ub else np.nan) for i in f_halflife['t_halflife']]
f_halflife = f_halflife.dropna()


plot_df = cor_df[['t_cor', 'p_cor']].unstack().reset_index()
plot_df.columns = ['cor_type', 'pair', 'cor']
plot_df = plot_df[[p.split('_')[0] in f_halflife.index and p.split('_')[1] in f_halflife.index for p in plot_df['pair']]]
plot_df = DataFrame([{
     'pair': pp, 'c_type': c_type, 'cor': cor, 'l_type': l_type, 'halflife': ' - '.join(f_halflife.ix[pp.split('_'), '%s_halflife' % l_type])}
for c_type, pp, cor in plot_df[['cor_type', 'pair', 'cor']].values for l_type in ['p', 't']])

plot_df['c_type'] = ['Transcriptomics' if i.startswith('t_') else 'Proteomics' for i in plot_df['c_type']]
plot_df['l_type'] = ['mRNA' if i == 't' else 'Protein' for i in plot_df['l_type']]
print plot_df


order = ['Stable - Stable', 'Stable - Unstable', 'Unstable - Stable', 'Unstable - Unstable']
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.FacetGrid(plot_df, row='l_type', xlim=[-1, 1], size=1.5, aspect=.8, legend_out=False)
# g.map(sns.violinplot, 'cor', 'halflife', 'c_type', palette=palette, order=order, cut=0, orient='h', inner='quartile', linewidth=.3)
g.map(sns.boxplot, 'cor', 'halflife', 'c_type', palette=palette, order=order, orient='h', sym='', notch=True, linewidth=.3)
g.map(sns.stripplot, 'cor', 'halflife', 'c_type', palette=palette, order=order, split=True, orient='h', size=2, edgecolor='white', linewidth=.3, jitter=True, alpha=.2)
g.map(plt.axvline, x=0, ls='--', lw=0.3, c='black', alpha=.5)
g.set_axis_labels('Pearson\'s r', '')
g.add_legend()
g.despine(trim=True)
g.set_titles(row_template='{row_name}')
plt.gcf().set_size_inches(3, 4)
plt.savefig('%s/reports/protein_pairs_correlation_halflife_violin.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
