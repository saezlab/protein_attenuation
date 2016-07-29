import os
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
from pymist.enrichment.gsea import gsea
from scipy.stats.stats import pearsonr, ttest_ind
from mtkirc.slapenrich import slapenrich
from sklearn.metrics.ranking import roc_curve, auc
from cptac.utils import hypergeom_test, read_gmt
from cptac import wd, default_color, palette_cnv_number
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest
from pymist.utils.stringdb import get_stringdb
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import protein-pairs
uniprot = read_uniprot_genename()

# CORUM
corum = set()
for p1, p2 in get_complexes_pairs():
    if (p2, p1) not in corum:
        corum.add((p1, p2))
corum = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in corum if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}

corum_proteins = {p for p1, p2 in corum for p in [p1, p2]}
print 'corum', len(corum)

# String
string = set()
for p1, p2 in get_stringdb(900):
    if (p2, p1) not in string:
        string.add((p1, p2))
string = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in string if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}

string_proteins = {p for p1, p2 in string for p in [p1, p2]}
print 'string', len(string)

# Overlap
p_pairs = corum.union(string)
p_pairs_proteins = corum_proteins.union(string_proteins)
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
# p1, p2 = 'SDHA', 'SDHB'
def protein_correlation(p1, p2):
    if p1 in proteins and p2 in proteins:
        samples = list(concat([proteomics.ix[[p1, p2]].T, transcriptomics.ix[[p1, p2]].T], axis=1).dropna().index)

        p_cor, p_pval = pearsonr(proteomics.ix[p1, samples], proteomics.ix[p2, samples])
        t_cor, t_pval = pearsonr(transcriptomics.ix[p1, samples], transcriptomics.ix[p2, samples])

        return {'p1': p1, 'p2': p2, 'p_cor': p_cor, 'p_pval': p_pval, 't_cor': t_cor, 't_pval': t_pval, 'len': len(samples)}

cor_df = {'%s_%s' % (p1, p2): protein_correlation(p1, p2) for p1, p2 in p_pairs}
cor_df = DataFrame({i: cor_df[i] for i in cor_df if cor_df[i]}).T

cor_df['diff'] = [p_cor - t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]
cor_df['sum'] = [p_cor + t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]

cor_df['t_fdr'] = multipletests(cor_df['t_pval'], method='fdr_bh')[1]
cor_df['p_fdr'] = multipletests(cor_df['p_pval'], method='fdr_bh')[1]
print cor_df.sort('diff')


# -- Gene set enrichment
p_pairs_trans = cor_df[(cor_df['t_fdr'] < .05) & (cor_df['t_cor'] > 0) & (cor_df['p_fdr'] > .05)]
p_pairs_prot = cor_df[(cor_df['p_fdr'] < .05) & (cor_df['p_cor'] > 0) & (cor_df['t_fdr'] > .05)]


# signature_type, signature_name, signature, background = 'mf', 'PROTEIN_PHOSPHATASE_BINDING', msigdb_go_mf['PROTEIN_PHOSPHATASE_BINDING'], p_pairs_proteins
def enrichment_hypergeom(signature_type, signature_name, signature, background):
    # Hypergeometric
    t_pval, t_len = hypergeom_test(signature, background, {p for i in p_pairs_trans.index for p in i.split('_')})
    p_pval, p_len = hypergeom_test(signature, background, {p for i in p_pairs_prot.index for p in i.split('_')})

    # Effect size
    idx = Series({i: int(np.any([p in signature for p in i.split('_')])) for i in cor_df.index})

    t_mean = cor_df.ix[idx[idx == 0].index, 't_cor'].mean() - cor_df.ix[idx[idx == 1].index, 't_cor'].mean()
    p_mean = cor_df.ix[idx[idx == 0].index, 'p_cor'].mean() - cor_df.ix[idx[idx == 1].index, 'p_cor'].mean()

    return {
        'type': signature_type,
        'name': signature_name,
        't_pval': t_pval,
        't_len': t_len,
        't_mean': t_mean,
        'p_pval': p_pval,
        'p_len': p_len,
        'p_mean': p_mean,
        's_len': len(signature.intersection(background))
    }

hyper = DataFrame([enrichment_hypergeom(signature_type, signature_name, signatures[signature_name], p_pairs_proteins) for signature_type, signatures in [
    ('mf', msigdb_go_mf), ('cc', msigdb_go_cc), ('bp', msigdb_go_bp), ('ptms', ptms), ('kegg', msigdb_kegg)
] for signature_name in signatures]).dropna()

hyper['t_fdr'] = multipletests(hyper['t_pval'], method='fdr_bh')[1]
hyper['p_fdr'] = multipletests(hyper['p_pval'], method='fdr_bh')[1]
print hyper


# -- Plot enrichment
plot_df = hyper.loc[(hyper['t_fdr'] < .05) | (hyper['p_fdr'] < .05), ['name', 't_mean', 'p_mean']].set_index('name')
plot_df = plot_df[(plot_df.abs() > .1).sum(1) != 0]

cmap = sns.diverging_palette(220, 20, n=7, as_cmap=True)
sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
sns.clustermap(plot_df, cmap=cmap, center=0)
plt.gcf().set_size_inches(3, 20)
plt.savefig('%s/reports/protein_pairs_geneset_enrichment_clustermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
