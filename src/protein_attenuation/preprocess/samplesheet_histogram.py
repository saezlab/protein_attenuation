#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from protein_attenuation import wd, palette
from matplotlib_venn import venn2, venn2_circles
from pandas import DataFrame, Series, read_csv

# -- Import
samplesheet = Series.from_csv('%s/data/samplesheet.csv' % wd)

# Clinical data
clinical = set(read_csv('%s/data/clinical_data.tsv' % wd, sep='\t')['SAMPLE_brcID'])
print len(clinical)

# CNV
cnv = set(read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)).intersection(clinical)
print len(cnv)

# Transcriptomics
transcriptomics = set(read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)).intersection(clinical)
print len(transcriptomics)

# Proteomics
proteomics = set(read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)).intersection(clinical)
print len(proteomics)


# -- Histogram
plot_df = Series(dict(zip(*(np.unique(samplesheet, return_counts=True))))).reset_index()
plot_df.columns = ['tissue', 'samples']

sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.set_style({'xtick.direction': 'in', 'ytick.direction': 'in'})
sns.barplot('tissue', 'samples', palette=palette, data=plot_df, linewidth=0, ci=None)
sns.despine(bottom=True)
plt.xlabel('Tissue of origin')
plt.ylabel('Number of samples')
plt.title('Proteomics samples')
plt.gcf().set_size_inches(4, 6)
plt.savefig('%s/reports/samplesheet_histogram_tissue.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'


# -- Overlap
groups = {'CNV': cnv, 'Transcriptomics': transcriptomics}

sns.set(style='white', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
venn2(groups.values(), set_labels=groups.keys(), set_colors=[palette[k] for k in groups])
venn2_circles(groups.values(), linestyle='solid', color='white')

plt.savefig('%s/reports/samplesheet_venn_datasets.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
