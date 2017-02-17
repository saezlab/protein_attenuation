#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv
from matplotlib_venn import venn2, venn2_circles
from protein_attenuation import palette


# -- Improt regression results
ppairs_cnv = read_csv('./tables/ppairs_cnv_regulation_all_tmt.csv')
ppairs_phospho = read_csv('./tables/ppairs_phospho_regulation_all_tmt.csv')


# -- Venn: overlap between transcriptomics and CNV
sns.set(style='white', font_scale=1., rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

associations = {
    'Copy-number variation': {(px, py) for px, py, fdr in ppairs_cnv[['px', 'py', 'fdr']].values if fdr < .05},
    'Phosphoproteomics': {(px.split('_')[0], py) for px, py, fdr in ppairs_phospho[['px', 'py', 'fdr']].values if fdr < .05}
}

venn2(associations.values(), set_labels=associations.keys(), set_colors=[palette[k] for k in associations])
venn2_circles(associations.values(), linestyle='solid', color='white')

plt.savefig('./reports/regressions_associations_venn_tmt.pdf', bbox_inches='tight')
plt.savefig('./reports/regressions_associations_venn_tmt.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Venn overlap between transcriptomics and CNV: ', './reports/regressions_associations_venn.pdf'


associations['Copy-number variation'].intersection(associations['Phosphoproteomics'])

px, py = 'MRE11A', 'RAD50'

ppairs_phospho[[x.startswith(px) and (y == py) for x, y in ppairs_phospho[['px', 'py']].values]]
