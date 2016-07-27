import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
from pymist.utils.read_gmt import read_gmt
from pymist.enrichment.gsea import gsea
from scipy.stats.stats import pearsonr
from mtkirc.slapenrich import slapenrich
from sklearn.metrics.ranking import roc_curve, auc
from cptac import wd, default_color, palette_cnv_number
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_pairs, get_complexes_dict, get_complexes_name
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# CORUM
uniprot = read_uniprot_genename()
corum = set()
for p1, p2 in get_complexes_pairs():
    if (p2, p1) not in corum:
        corum.add((p1, p2))
corum = {(uniprot[p1][0], uniprot[p2][0]) for p1, p2 in corum if p1 in uniprot and p2 in uniprot and uniprot[p1][0] != uniprot[p2][0]}

corum_dict = get_complexes_dict()
corum_dict = {k: {uniprot[p][0] for p in corum_dict[k] if p in uniprot} for k in corum_dict}
corum_dict = {k: {'%s_%s' % (p1, p2) for p1, p2 in corum if p1 in corum_dict[k] and p2 in corum_dict[k]} for k in corum_dict}

corum_name = get_complexes_name()
print len(corum)

#
go_terms_bp = read_gmt('%s/files/c5.bp.v5.1.symbols.gmt' % wd)
go_terms_bp = {k: {'%s_%s' % (p1, p2) for p1, p2 in corum if p1 in go_terms_bp[k] and p2 in go_terms_bp[k]} for k in go_terms_bp}
print len(go_terms_bp)

go_terms_cc = read_gmt('%s/files/c5.cc.v5.1.symbols.gmt' % wd)
go_terms_cc = {k: {'%s_%s' % (p1, p2) for p1, p2 in corum if p1 in go_terms_cc[k] and p2 in go_terms_cc[k]} for k in go_terms_cc}
print len(go_terms_cc)

go_terms_mf = read_gmt('%s/files/c5.mf.v5.1.symbols.gmt' % wd)
go_terms_mf = {k: {'%s_%s' % (p1, p2) for p1, p2 in corum if p1 in go_terms_mf[k] and p2 in go_terms_mf[k]} for k in go_terms_mf}
print len(go_terms_mf)


# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print transcriptomics

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print proteomics

# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
print cnv

# PAX DB
protein_ppm_map = Series.from_csv('%s/files/9606-paxdb_uniprot.txt' % wd, sep='\t')
protein_ppm_map = protein_ppm_map.reset_index()
protein_ppm_map = protein_ppm_map.groupby('index')[0].agg(lambda x: set(x)).to_dict()
protein_ppm_map = {k: protein_ppm_map[k].intersection(uniprot) for k in protein_ppm_map}

protein_ppm = read_csv('%s/files/9606-WHOLE_ORGANISM-integrated.txt' % wd, sep='\t')
protein_ppm = protein_ppm[[i in protein_ppm_map.index for i in protein_ppm['string_external_id']]]
protein_ppm['uniprot'] = [protein_ppm_map[i] for i in protein_ppm['string_external_id']]
protein_ppm = DataFrame([{'gene': uniprot[p][0], 'ppm': v} for v, k in protein_ppm[['abundance', 'uniprot']].values for p in k])

protein_ppm_dup = Series(dict(zip(*(np.unique(protein_ppm['gene'], return_counts=True)))))
protein_ppm_dup = set(protein_ppm_dup[protein_ppm_dup > 1].index)

protein_ppm = protein_ppm[[i not in protein_ppm_dup for i in protein_ppm['gene']]].set_index('gene')
print protein_ppm

# -- Overlap
proteins = set(transcriptomics.index).intersection(proteomics.index)
samples = set(transcriptomics).intersection(proteomics)
print len(proteins), len(samples)


# --
# p1, p2 = 'SDHA', 'SDHB'
def protein_correlation(p1, p2):
    if p1 in proteins and p2 in proteins:
        samples = list(concat([proteomics.ix[[p1, p2]].T, transcriptomics.ix[[p1, p2]].T], axis=1).dropna().index)

        p_cor, p_pval = pearsonr(proteomics.ix[p1, samples], proteomics.ix[p2, samples])
        t_cor, t_pval = pearsonr(transcriptomics.ix[p1, samples], transcriptomics.ix[p2, samples])

        return {'p1': p1, 'p2': p2, 'p_cor': p_cor, 'p_pval': p_pval, 't_cor': t_cor, 't_pval': t_pval, 'len': len(samples)}

res = [protein_correlation(p1, p2) for p1, p2 in corum]
res = DataFrame([i for i in res if i])
res['diff'] = [p_cor - t_cor for p_cor, t_cor in res[['p_cor', 't_cor']].values]
res['sum'] = [p_cor + t_cor for p_cor, t_cor in res[['p_cor', 't_cor']].values]
res['pair'] = ['%s_%s' % (p1, p2) for p1, p2 in res[['p1', 'p2']].values]
res = res.set_index('pair')
res = res[(res['p_cor'] > 0) & (res['t_cor'] > 0)]
res['ppm'] = [protein_ppm.ix[i.split('_'), 'ppm'].mean() for i in res.index]
print res.sort('diff')


# --
sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=.75)
g = sns.distplot(res['diff'], color=default_color, hist=False, kde_kws={'shade': True})
plt.axvline(0, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Correlation difference (median: %.2f)' % res['diff'].median())
plt.xlabel('pearson(proteomics) - pearson(transcriptomics)')
plt.ylabel('Density')
sns.despine(trim=True)
plt.gcf().set_size_inches(6, 3)
plt.savefig('%s/reports/corum_pairs_correlation_difference_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
corum_aroc = {}
for k in corum_dict:
    df = res[['diff']]
    df['TP'] = [int(i in corum_dict[k]) for i in df.index]

    curve_fpr, curve_tpr, _ = roc_curve(df['TP'], df['diff'])
    curve_auc = auc(curve_fpr, curve_tpr)

    corum_aroc[k] = curve_auc
corum_aroc = DataFrame(Series(corum_aroc).dropna().sort_values(), columns=['aroc'])
corum_aroc['name'] = [corum_name[i] for i in corum_aroc.index]
print corum_aroc

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=.75)
g = sns.distplot(corum_aroc['aroc'], color=default_color, hist=False, kde_kws={'shade': True})
g.set(xlim=(corum_aroc['aroc'].min(), corum_aroc['aroc'].max()))
plt.axvline(0.5, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Protein complex (median: %.2f)' % corum_aroc['aroc'].median())
plt.xlabel('AROC')
plt.ylabel('Density')
sns.despine()
plt.gcf().set_size_inches(6, 3)
plt.savefig('%s/reports/corum_pairs_aroc_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
corum_aroc_sum = {}
for k in corum_dict:
    df = res[['sum']]
    df['TP'] = [int(i in corum_dict[k]) for i in df.index]

    curve_fpr, curve_tpr, _ = roc_curve(df['TP'], df['sum'])
    curve_auc = auc(curve_fpr, curve_tpr)

    corum_aroc_sum[k] = curve_auc
corum_aroc_sum = DataFrame(Series(corum_aroc_sum).dropna().sort_values(), columns=['aroc'])
corum_aroc_sum['name'] = [corum_name[i] for i in corum_aroc_sum.index]
print corum_aroc_sum

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=.75)
g = sns.distplot(corum_aroc_sum['aroc'], color=default_color, hist=False, kde_kws={'shade': True})
g.set(xlim=(corum_aroc_sum['aroc'].min(), corum_aroc_sum['aroc'].max()))
plt.axvline(0.5, ls='--', lw=0.3, c='black', alpha=.5)
plt.title('Protein complex (median: %.2f)' % corum_aroc_sum['aroc'].median())
plt.xlabel('AROC')
plt.ylabel('Density')
sns.despine()
plt.gcf().set_size_inches(6, 3)
plt.savefig('%s/reports/corum_pairs_aroc_sum_histogram.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# --
go_terms = go_terms_mf.copy()

goterms_diff_aroc = {}
for k in go_terms:
    if len(go_terms[k]) > 0:
        df = res[['diff']]
        df['TP'] = [int(i in go_terms[k]) for i in df.index]

        curve_fpr, curve_tpr, _ = roc_curve(df['TP'], df['diff'])
        curve_auc = auc(curve_fpr, curve_tpr)

        goterms_diff_aroc[k] = curve_auc
goterms_diff_aroc = DataFrame(Series(goterms_diff_aroc).dropna().sort_values(), columns=['aroc'])
print goterms_diff_aroc
