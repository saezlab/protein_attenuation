import os
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
from pymist.enrichment.gsea import gsea
from scipy.stats.stats import pearsonr, ttest_ind, spearmanr
from matplotlib.gridspec import GridSpec
from mtkirc.slapenrich import slapenrich
from sklearn.metrics.ranking import roc_curve, auc
from cptac.utils import hypergeom_test, read_gmt
from cptac import wd, default_color, palette_cnv_number, palette
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest
from pymist.utils.stringdb import get_stringdb
from pandas import DataFrame, Series, read_csv, concat, pivot_table
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name
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

# String
string = set()
for p1, p2 in get_stringdb(900):
    if (p2, p1) not in string:
        string.add((p1, p2))
string = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in string if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}
print 'string', len(string)

# Overlap
p_pairs = corum.union(string)
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

        p1_cor, p1_pval = pearsonr(proteomics.ix[p1, samples], transcriptomics.ix[p1, samples])
        p2_cor, p2_pval = pearsonr(proteomics.ix[p2, samples], transcriptomics.ix[p2, samples])

        return {
            'p1': p1, 'p2': p2,
            'p_cor': p_cor, 'p_pval': p_pval,
            't_cor': t_cor, 't_pval': t_pval,
            'p1_cor': p1_cor, 'p1_pval': p1_pval,
            'p2_cor': p2_cor, 'p2_pval': p2_pval,
            'len': len(samples)
        }

cor_df = {'%s_%s' % (p1, p2): protein_correlation(p1, p2) for p1, p2 in p_pairs}
cor_df = DataFrame({i: cor_df[i] for i in cor_df if cor_df[i]}).T

cor_df['diff'] = [p_cor - t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]
cor_df['sum'] = [p_cor + t_cor for p_cor, t_cor in cor_df[['p_cor', 't_cor']].values]

cor_df['t_fdr'] = multipletests(cor_df['t_pval'], method='fdr_bh')[1]
cor_df['p_fdr'] = multipletests(cor_df['p_pval'], method='fdr_bh')[1]

cor_df['p1_fdr'] = multipletests(cor_df['p1_pval'], method='fdr_bh')[1]
cor_df['p2_fdr'] = multipletests(cor_df['p2_pval'], method='fdr_bh')[1]

cor_df['p_diff'] = [p1_cor - p2_cor for p1_cor, p2_cor in cor_df[['p1_cor', 'p2_cor']].values]
print cor_df.sort('diff')


# --
p_pairs_prot = cor_df[(cor_df['p_fdr'] < .05) & (cor_df['p_cor'] > .0) & (cor_df['t_fdr'] > .05)]
protein_prot = {p for p1, p2, p1_fdr, p2_fdr in p_pairs_prot[['p1', 'p2', 'p1_fdr', 'p2_fdr']].values for p, fdr in [(p1, p1_fdr), (p2, p2_fdr)] if fdr > .05}

p_pairs_trans = cor_df[(cor_df['t_fdr'] < .05) & (cor_df['t_cor'] > .0) & (cor_df['p_fdr'] > .05)]
protein_trans = {p for p1, p2, p1_fdr, p2_fdr in p_pairs_trans[['p1', 'p2', 'p1_fdr', 'p2_fdr']].values for p, fdr in [(p1, p1_fdr), (p2, p2_fdr)] if fdr > .05}
print 'protein_prot', 'protein_trans', len(protein_prot), len(protein_trans)

protein_prot_xor = protein_prot.difference(protein_trans)
protein_trans_xor = protein_trans.difference(protein_prot)
print 'protein_prot_xor', 'protein_trans_xor', len(protein_prot_xor), len(protein_trans_xor)


# signature_type, signature_name, signature, background = 'mf', 'PROTEIN_PHOSPHATASE_BINDING', msigdb_go_mf['PROTEIN_PHOSPHATASE_BINDING'], p_pairs_proteins
def enrichment_hypergeom(signature_type, signature_name, signature):
    # Hypergeometric
    p_pval, p_len = hypergeom_test(signature, p_pairs_proteins, protein_prot)
    t_pval, t_len = hypergeom_test(signature, p_pairs_proteins, protein_trans)

    # XOR
    p_xor_pval, p_xor_len = hypergeom_test(signature, p_pairs_proteins, protein_prot_xor)
    t_xor_pval, t_xor_len = hypergeom_test(signature, p_pairs_proteins, protein_trans_xor)

    return {
        'type': signature_type,
        'name': signature_name,
        'p_pval': p_pval, 'p_len': p_len,
        't_pval': t_pval, 't_len': t_len,
        'p_xor_pval': p_xor_pval, 'p_xor_len': p_xor_len,
        't_xor_pval': t_xor_pval, 't_xor_len': t_xor_len,
        's_len': len(signature.intersection(uniprot_proteins))
    }

hyper = DataFrame([enrichment_hypergeom(signature_type, signature_name, signatures[signature_name]) for signature_type, signatures in [
    # ('mf', msigdb_go_mf), ('cc', msigdb_go_cc), ('bp', msigdb_go_bp), ('ptms', ptms), ('kegg', msigdb_kegg)
    ('cc', msigdb_go_cc), ('bp', msigdb_go_bp)
] for signature_name in signatures]).dropna()
hyper['p_fdr'] = multipletests(hyper['p_pval'], method='fdr_bh')[1]
hyper['t_fdr'] = multipletests(hyper['t_pval'], method='fdr_bh')[1]
hyper['p_xor_fdr'] = multipletests(hyper['p_xor_pval'], method='fdr_bh')[1]
hyper['t_xor_fdr'] = multipletests(hyper['t_xor_pval'], method='fdr_bh')[1]
print hyper.sort('t_fdr')


# --
prot_go_terms = {n: db[n].intersection(protein_prot_xor) for n in hyper[hyper['p_xor_fdr'] < .05].sort('p_xor_fdr')['name'] for db in [msigdb_go_cc, msigdb_go_bp] if n in db}
trans_go_terms = {n: db[n].intersection(protein_trans_xor) for n in hyper[hyper['t_xor_fdr'] < .05].sort('t_xor_fdr')['name'] for db in [msigdb_go_cc, msigdb_go_bp] if n in db}

prot_go_terms_proteins = {p for n in prot_go_terms for p in prot_go_terms[n]}
trans_go_terms_proteins = {p for n in trans_go_terms for p in trans_go_terms[n]}

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
pos, gs = 0, GridSpec(1, 4, hspace=.2, wspace=.75)
for db_n, db in [('proteomics', prot_go_terms), ('transcriptomics', trans_go_terms)]:
    ax = plt.subplot(gs[pos])
    plot_df = DataFrame([{'geneset': n, 'overlap': len(db[n])} for n in db]).sort('overlap', ascending=False)
    plot_df['geneset'] = [i.lower().replace('_', ' ') for i in plot_df['geneset']]
    sns.barplot('overlap', 'geneset', data=plot_df, color=default_color, ax=ax, orient='h', ci=None, linewidth=0)
    sns.despine(ax=ax)
    ax.set_xlabel('Overlap')
    ax.set_title('Protein-pairs correlating at %s level\nSignificant gene-sets overlap' % db_n)
    pos += 1

    ax = plt.subplot(gs[pos])
    plot_df = DataFrame(zip(*(np.unique([p for n in db for p in db[n]], return_counts=True))), columns=['gene', 'counts']).sort('counts', ascending=False)
    sns.barplot('counts', 'gene', data=plot_df, color=default_color, ax=ax, orient='h', ci=None, linewidth=0)
    sns.despine(ax=ax)
    ax.set_xlabel('Counts')
    ax.set_title('Protein-pairs correlating at %s level\nSignificant gene-sets gene counts' % db_n)
    pos += 1

plt.gcf().set_size_inches(20, 7)
plt.savefig('%s/reports/protein_pairs_go_terms_barplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plot protein pair
[(p1, p2) for p1, p2 in p_pairs if p1 in trans_go_terms_proteins and p2 in trans_go_terms_proteins and (p1 == 'NDUFS8' or p2 == 'NDUFS8')]

p1, p2 = 'RPL29', 'RPSA'

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
gs = GridSpec(2, 2, hspace=.2, wspace=.2)

samples = list(concat([proteomics.ix[[p1, p2]].T, transcriptomics.ix[[p1, p2]].T], axis=1).dropna().index)

ax = plt.subplot(gs[0])
sns.regplot(proteomics.ix[p1, samples], transcriptomics.ix[p1, samples], ax=ax, color=default_color)
ax.set_title(p1)
ax.set_xlabel('Proteomics')
ax.set_ylabel('Transcriptomics')
ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)

ax = plt.subplot(gs[1])
sns.regplot(proteomics.ix[p2, samples], transcriptomics.ix[p2, samples], ax=ax, color=default_color)
ax.set_title(p2)
ax.set_xlabel('Proteomics')
ax.set_ylabel('Transcriptomics')
ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)

ax = plt.subplot(gs[2])
sns.regplot(proteomics.ix[p1, samples], proteomics.ix[p2, samples], ax=ax, color=palette['Proteomics'])
ax.set_title('Proteomics')
ax.set_xlabel(p1)
ax.set_ylabel(p2)
ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)

ax = plt.subplot(gs[3])
sns.regplot(transcriptomics.ix[p1, samples], transcriptomics.ix[p2, samples], ax=ax, color=palette['Transcriptomics'])
ax.set_title('Transcriptomics')
ax.set_xlabel(p1)
ax.set_ylabel(p2)
ax.axvline(0, ls='--', lw=0.1, c='black', alpha=.3)
ax.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)

plt.gcf().set_size_inches(8, 8)
plt.savefig('%s/reports/protein_pairs_cor_scatter.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
