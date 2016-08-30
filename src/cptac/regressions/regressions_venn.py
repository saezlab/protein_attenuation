import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette
from matplotlib_venn import venn3, venn3_circles
from pandas import DataFrame, Series, read_csv


# -- Import associations
ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all.csv')
ppairs_cnv = ppairs_cnv[ppairs_cnv['fdr'] < .05]
print 'ppairs_cnv', ppairs_cnv.shape

ppairs_trans = read_csv('./tables/ppairs_transcriptomics_regulation_all.csv')
ppairs_trans = ppairs_trans[ppairs_trans['fdr'] < .05]
print 'ppairs_trans', ppairs_trans.shape

ppairs_prot = read_csv('./tables/ppairs_proteomics_regulation_all.csv')
ppairs_prot = ppairs_prot[ppairs_prot['fdr'] < .05]
print 'ppairs_prot', ppairs_prot.shape


# -- Overlap
sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

associations = {
    'CNV': {(px, py) for px, py in ppairs_cnv[['px', 'py']].values},
    'Transcriptomics': {(px, py) for px, py in ppairs_trans[['px', 'py']].values},
    'Proteomics': {(px, py) for px, py in ppairs_prot[['px', 'py']].values},
}

venn3(associations.values(), set_labels=associations.keys(), set_colors=[palette[k] for k in associations])
venn3_circles(associations.values(), linestyle='solid', color='white')

plt.savefig('./reports/regressions_associations_venn.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
