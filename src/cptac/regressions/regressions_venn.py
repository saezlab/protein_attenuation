import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, palette
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn3, venn3_circles
from pandas import DataFrame, Series, read_csv


# -- Import associations
transcriptomics = read_csv('%s/tables/regressions_transcriptomics_cnv.csv' % wd)
proteomics = read_csv('%s/tables/regressions_proteomics_cnv.csv' % wd)
residuals = read_csv('%s/tables/regressions_residuals_cnv.csv' % wd)


# -- Overlap
sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
fig, gs, pos = plt.figure(figsize=(5, 10)), GridSpec(1, 2, hspace=.3), 0

for c in [-2, 2]:
    ax = plt.subplot(gs[pos])

    associations = {
        n: {p for p, t, f in df[['protein', 'type', 'fdr']].values if c == t and f < .05}
        for n, df in [('transcriptomics', transcriptomics), ('proteomics', proteomics), ('residuals', residuals)]}

    venn3(associations.values(), set_labels=associations.keys(), set_colors=[palette_datasets[k] for k in associations])
    venn3_circles(associations.values(), linestyle='solid', color='white')

    ax.set_title('Depletion' if c == -2 else 'Amplification')

    pos += 1

plt.savefig('%s/reports/regressions_venn.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
