#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from protein_attenuation import palette
from protein_attenuation.utils import jaccard
from pandas import DataFrame, Series, read_csv


# -- Imports
# Samplesheet
samplesheet = Series.from_csv('./data/samplesheet.csv')

# Proteomics
proteomics = read_csv('./data/cptac_proteomics.csv', index_col=0)

# Transcriptomics
transcriptomics = read_csv('./data/tcga_rnaseq.csv', index_col=0)

# CNV
cnv = read_csv('./data/tcga_cnv.csv', index_col=0)

# Clinical data
clinical = read_csv('./data/tcga_clinical.csv').dropna(subset=['patient.gender', 'patient.days_to_birth'])


# -- Overlap table
overlap = DataFrame({
    'Proteomics': {s: int(s in proteomics) for s in proteomics},
    'Transcriptomics': {s: int(s in transcriptomics) for s in proteomics},
    'Copy-number variation': {s: int(s in cnv) for s in proteomics},
    'Clinical': {s: int(s in clinical['patient.bcr_patient_barcode'].values) for s in proteomics},
    'Tumour': {s: samplesheet[s] for s in proteomics}
})

overlap['Overlap'] = (overlap[['Transcriptomics', 'Copy-number variation', 'Clinical', 'Proteomics']].sum(1) == 4).astype(int)
overlap.to_csv('./tables/overlap_table.csv')
print '[INFO] Samples data-sets overlap table: ', './tables/overlap_table.csv'


# -- Samples count
plot_df = DataFrame([{'Tumour': t, 'Data': d, 'Counts': overlap[overlap['Tumour'] == t][d].sum()} for t in set(samplesheet) for d in ['Transcriptomics', 'Copy-number variation', 'Clinical', 'Overlap', 'Proteomics']])

hue_order = ['Copy-number variation', 'Transcriptomics', 'Proteomics', 'Clinical', 'Overlap']
sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
g = sns.factorplot(x='Tumour', y='Counts', hue='Data', data=plot_df, kind='bar', palette=palette, ci=None, hue_order=hue_order, legend_out=False, lw=.3)
g.despine()
g.add_legend(title='', label_order=hue_order)
g.set_ylabels('Number of samples')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/overlap_samples.pdf', bbox_inches='tight')
plt.savefig('./reports/overlap_samples.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Samples data-sets overlap: ', './reports/overlap_samples.pdf'


# -- Proteins overlap
plot_df = []
for s in set(proteomics).intersection(transcriptomics):
    jt = jaccard(set(transcriptomics.index), set(proteomics[s].dropna().index))
    jp = jaccard(set(proteomics[s].dropna().index), set(transcriptomics[s].index))

    res = {'Sample': s, 'Tumour': samplesheet[s], 'Transcript': jt, 'Protein': jp}

    plot_df.append(res)

plot_df = DataFrame(plot_df)

order = ['COREAD', 'HGSC', 'BRCA']
sns.set(style='ticks', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})

sns.boxplot(x='Tumour', y='Transcript', data=plot_df, palette=palette, linewidth=.5, notch=True, sym='', order=order)
sns.stripplot(x='Tumour', y='Transcript', data=plot_df, palette=palette, linewidth=.3, jitter=True, size=2, alpha=.5, order=order)
for x, t in zip(*([0, 1, 2], order)):
    plt.text(x, plot_df[plot_df['Tumour'] == t]['Transcript'].max() * 1.05, 'N = %d' % (plot_df['Tumour'] == t).sum(), ha='center', va='bottom', fontsize=5)

plt.axhline(0.5, ls='--', lw=0.3, c='black', alpha=.5)
plt.ylim(0, .7)
sns.despine(trim=True)
plt.ylabel('Proteomic coverage of the expressed transcriptome\n(Jaccard index)')
plt.gcf().set_size_inches(1.5, 3)
plt.savefig('./reports/overlap_proteins.pdf', bbox_inches='tight')
plt.savefig('./reports/overlap_proteins.png', bbox_inches='tight', dpi=300)
plt.close('all')
print '[INFO] Proteins coverage: ', './reports/overlap_proteins.pdf'
