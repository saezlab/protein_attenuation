import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, default_color, palette_cnv_number
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ztest
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pymist.utils.map_peptide_sequence import read_uniprot_genename

# -- Proteomics
prot = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print prot

# -- Complexes proteins
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
print len(corum)


# -- Estimate CORUM
def ztest_complex(s, c):
    x1 = prot[s].ix[corum[c]].dropna()
    x2 = prot[s].drop(corum[c], errors='ignore').dropna()

    if len(x1) > 1:
        z, pvalue = ztest(x1, x2)
        return s, c, z, pvalue, x1.mean(), len(x1)

c_activity = [ztest_complex(s, c) for s in prot for c in corum]
c_activity = DataFrame([i for i in c_activity if i], columns=['sample', 'complex', 'z', 'pval', 'mean', 'targets'])
c_activity['FDR'] = multipletests(c_activity['pval'],  method='fdr_bh')[1]
c_activity.to_csv('%s/tables/protein_complexes_activities.csv' % wd)
print c_activity.sort('FDR')

# -- CNV validation
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)

plot_df = [(s, c, p, z, cnv.ix[p, s]) for s, c, z, m in c_activity[['sample', 'complex', 'z', 'mean']].values for p in corum[c] if p in cnv.index and s in cnv.columns]
plot_df = DataFrame(plot_df, columns=['sample', 'complex', 'protein', 'z', 'cnv'])

sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
plt.axhline(0, ls='--', lw=0.1, c='black', alpha=.3)
sns.violinplot('cnv', 'z', data=plot_df, palette=palette_cnv_number, linewidth=.3, cut=0, inner='quartiles')
sns.despine(trim=True)
plt.gcf().set_size_inches(3, 3)
plt.savefig('%s/reports/c_activities_boxplots.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
