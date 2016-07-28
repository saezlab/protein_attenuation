import os
import pypath
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


# -- Imports
uniprot = read_uniprot_genename()

# # CORUM
# corum = set()
# for p1, p2 in get_complexes_pairs():
#     if (p2, p1) not in corum:
#         corum.add((p1, p2))
# corum = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in corum if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}
#
# corum_dict = get_complexes_dict()
# corum_dict = {k: {uniprot[p][0] for p in corum_dict[k] if p in uniprot} for k in corum_dict}
#
# corum_name = get_complexes_name()
#
# corum_proteins = {p for p1, p2 in corum for p in [p1, p2]}
# print len(corum)

# String
string_thres = {'low': 150, 'medium': 400, 'high': 700, 'highest': 900}

string = set()
for p1, p2 in get_stringdb(string_thres['highest']):
    if (p2, p1) not in string:
        string.add((p1, p2))
string = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in string if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}

string_proteins = {p for p1, p2 in string for p in [p1, p2]}
print 'string', len(string)


# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

# uniprot Protein lists
ptms = {'_'.join(f[:-4].split('_')[1:]):
            {uniprot[i][0] for i in read_csv('%s/files/uniprot_ptms_proteins/%s' % (wd, f), sep='\t')['Entry'] if i in uniprot}
    for f in os.listdir('%s/files/uniprot_ptms_proteins/' % wd) if f.startswith('uniprot_')
}
print 'ptms', len(ptms)

# GO terms
go_terms_bp = read_gmt('%s/files/c5.bp.v5.1.symbols.gmt' % wd)
go_terms_cc = read_gmt('%s/files/c5.cc.v5.1.symbols.gmt' % wd)
go_terms_mf = read_gmt('%s/files/c5.mf.v5.1.symbols.gmt' % wd)
print 'go_terms_mf', 'go_terms_cc', 'go_terms_bp', len(go_terms_mf), len(go_terms_cc), len(go_terms_bp)

msigdb_cp = read_gmt('%s/files/c2.cp.v5.1.symbols.gmt' % wd)
msigdb_kegg = read_gmt('%s/files/c2.cp.kegg.v5.1.symbols.gmt' % wd)
msigdb_cgp = read_gmt('%s/files/c2.cgp.v5.1.symbols.gmt' % wd)
print 'msigdb_cp', 'msigdb_kegg', 'msigdb_cgp', len(msigdb_cp), len(msigdb_kegg), len(msigdb_cgp)


# -- Overlap
proteins = set(transcriptomics.index).intersection(proteomics.index)
samples = set(transcriptomics).intersection(proteomics)
print len(proteins), len(samples)


# -- Protein pairs correlation
# p1, p2 = 'SDHA', 'SDHB'
def protein_correlation(p1, p2):
    if p1 in proteins and p2 in proteins:
        samples = list(concat([proteomics.ix[[p1, p2]].T, transcriptomics.ix[[p1, p2]].T], axis=1).dropna().index)

        p_cor, p_pval = pearsonr(proteomics.ix[p1, samples], proteomics.ix[p2, samples])
        t_cor, t_pval = pearsonr(transcriptomics.ix[p1, samples], transcriptomics.ix[p2, samples])

        return {'p1': p1, 'p2': p2, 'p_cor': p_cor, 'p_pval': p_pval, 't_cor': t_cor, 't_pval': t_pval, 'len': len(samples)}

cor_df = {'%s_%s' % (p1, p2): protein_correlation(p1, p2) for p1, p2 in string}
cor_df = DataFrame({i: cor_df[i] for i in cor_df if cor_df[i]}).T

cor_df['diff'] = [p_cor - t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]
cor_df['sum'] = [p_cor + t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]

cor_df['t_fdr'] = multipletests(cor_df['t_pval'], method='fdr_bh')[1]
cor_df['p_fdr'] = multipletests(cor_df['p_pval'], method='fdr_bh')[1]

print cor_df.sort('diff')


# --
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=.75)
g = sns.distplot(cor_df['diff'], color=default_color, hist=False, kde_kws={'shade': True})
plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Correlation difference (median: %.2f)' % cor_df['diff'].median())
plt.xlabel('pearson(proteomics) - pearson(transcriptomics)')
plt.ylabel('Density')
sns.despine(trim=True)
plt.gcf().set_size_inches(6, 3)
plt.savefig('%s/reports/corum_pairs_correlation_difference_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
p_pairs_trans = cor_df[(cor_df['t_fdr'] < .05) & (cor_df['t_cor'] > 0) & (cor_df['p_fdr'] > .05)]
p_pairs_prot = cor_df[(cor_df['p_fdr'] < .05) & (cor_df['p_cor'] > 0) & (cor_df['t_fdr'] > .05)]


def enrichment_hypergeom(set_type, set_name, subset_type, pair_set, background, subset):
    h_pval, h_len = hypergeom_test(pair_set, background, {p for i in subset.index for p in i.split('_')})
    return {'type': set_type, 'set_name': set_name, 'subset_type': subset_type, 'pval': h_pval, 'len': h_len}

hyper = DataFrame([enrichment_hypergeom(set_type, set_name, subset_type, pair_set, string_proteins, subset) for set_type, set_paths in [
    ('mf', go_terms_mf), ('cc', go_terms_cc), ('bp', go_terms_bp),
    ('ptms', ptms),
    ('cp', msigdb_cp)
    # ('kegg', msigdb_kegg)
] for set_name, pair_set in set_paths.items() for subset_type, subset in [('transcriptomics', p_pairs_trans), ('proteomics', p_pairs_prot)]])
# print go_term_enrch[go_term_enrch['len'] > 20].sort('aroc')
hyper['fdr'] = multipletests(hyper['pval'], method='fdr_bh')[1]
print hyper.sort('fdr')


hyper_trans = hyper[(hyper['subset_type'] == 'transcriptomics') & (hyper['fdr'] < .05)]
hyper_prot = hyper[(hyper['subset_type'] == 'proteomics') & (hyper['fdr'] < .05)]

hyper_ov = set(hyper_trans['set_name']).intersection(hyper_prot['set_name'])

print hyper_trans[[i not in hyper_ov for i in hyper_trans['set_name']]].sort('fdr')
print hyper_prot[[i not in hyper_ov for i in hyper_prot['set_name']]].sort('fdr')


# --
# type_set, pair_set = ('isopeptide_bond', ptms['isopeptide_bond'])
def enrichment_ttest(type_set, pair_set):
    indx = Series({i: int(np.all([p in pair_set for p in i.split('_')])) for i in cor_df.index})
    p_bck, p_set = cor_df.ix[indx[indx == 0].index, 'diff'], cor_df.ix[indx[indx == 1].index, 'diff']

    t_stat, t_pval = ttest_ind(p_bck, p_set, equal_var=False)
    return {'type': type_set, 't_stat': t_stat, 't_pval': t_pval, 'size': len(p_set), 'diff': p_bck.mean() - p_set.mean()}

go_term_enrch = DataFrame({k.lower(): enrichment_ttest(t, v) for t, df in [
    # ('mf', go_terms_mf),
    # ('cc', go_terms_cc),
    # ('bp', go_terms_bp)
    ('ptms', ptms)
    # ('cp', msigdb_cp)
    # ('kegg', msigdb_kegg)
] for k, v in df.items()}).T.dropna()
print go_term_enrch.sort('diff')


# -- Enrichment analysis
# pair_set = go_terms_bp['UBIQUITIN_CYCLE']
# type_set, pair_set = ('oxidation', ptms['oxidation'])
def enrichment_aroc(type_set, pair_set):
    df = cor_df[['diff']].copy()
    df['TP'] = [int(np.all([p in pair_set for p in i.split('_')])) for i in df.index]

    curve_fpr, curve_tpr, _ = roc_curve(df['TP'], df['diff'])
    curve_auc = auc(curve_fpr, curve_tpr)

    return {'aroc': curve_auc, 'type': type_set, 'len': df['TP'].sum()}

go_term_enrch = DataFrame({k.lower(): enrichment_aroc(t, v) for t, df in [
    # ('mf', go_terms_mf),
    # ('cc', go_terms_cc),
    # ('bp', go_terms_bp)
    ('ptms', ptms)
    # ('cp', msigdb_cp)
    # ('kegg', msigdb_kegg)
] for k, v in df.items()}).T.dropna()
print go_term_enrch.sort('aroc')


# type_set, pair_set = 'kegg', msigdb_kegg['kegg_n_glycan_biosynthesis'.upper()]
def enrichment_hypergeom(type_set, pair_set, background, thres=.0):
    tran_proteins = {p for i in cor_df[(cor_df['t_fdr'] < .05) & (cor_df['p_fdr'] > .05) & (cor_df['diff'].abs() > thres)].index for p in i.split('_')}
    prot_proteins = {p for i in cor_df[(cor_df['t_fdr'] > .05) & (cor_df['p_fdr'] < .05) & (cor_df['diff'].abs() > thres)].index for p in i.split('_')}

    # tran_proteins = {p for i in cor_df[cor_df['diff'] < -thres].index for p in i.split('_')}
    # prot_proteins = {p for i in cor_df[cor_df['diff'] > thres].index for p in i.split('_')}

    tran_pval, tran_len = hypergeom_test(pair_set, background, tran_proteins)
    prot_pval, prot_len = hypergeom_test(pair_set, background, prot_proteins)

    if tran_len <= 1 and prot_len <= 1:
        return {'type': type_set, 'trans_pval': np.nan, 'prot_pval': np.nan, 'trans_len': tran_len, 'prot_len': prot_len}
    else:
        return {'type': type_set, 'trans_pval': tran_pval, 'prot_pval': prot_pval, 'trans_len': tran_len, 'prot_len': prot_len}

hyper = DataFrame({k.lower(): enrichment_hypergeom(t, v, string_proteins) for t, df in [
    # ('mf', go_terms_mf), ('cc', go_terms_cc), ('bp', go_terms_bp)
    ('ptms', ptms)
    # ('cp', msigdb_cp)
    # ('kegg', msigdb_kegg)
] for k, v in df.items()}).T.dropna()
# print go_term_enrch[go_term_enrch['len'] > 20].sort('aroc')
print hyper.sort('trans_pval')
print hyper.sort('prot_pval')
