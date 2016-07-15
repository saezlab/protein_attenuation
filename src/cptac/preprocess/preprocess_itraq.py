import os
import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, data
from cptac.utils import gkn
from pandas.stats.misc import zscore
from pandas import DataFrame, Series, read_csv, concat
from pymist.utils.map_peptide_sequence import read_fasta


columns = ['Protein', 'iTRAQ114', 'iTRAQ115', 'iTRAQ116', 'iTRAQ117']

t_dir = 'TCGA Ovarian Cancer'

ss = read_csv('%s/TCGA Ovarian Cancer/CPTAC_TCGA_Ovarian_Cancer_iTRAQ_Sample_Mapping.txt' % data, sep='\t', index_col=2)

psms = []
# s_dir = ss.index[0]
for s_dir in ss.index:
    # Assemble PSM fraction files
    psm = concat([read_csv('%s/%s/%s/PSMs/TSV/%s' % (data, t_dir, s_dir, f), sep='\t', index_col=False)[columns] for f in os.listdir('%s/%s/%s/PSMs/TSV/' % (data, t_dir, s_dir)) if f.endswith('.raw.cap.psm')])

    # Remove ambigously mapping peptides
    psm = psm[[len(i.split(';')) == 1 for i in psm[columns[0]]]]
    psm[columns[0]] = [i.split('|')[3] for i in psm[columns[0]]]

    # Remove intensity QC scores
    psm[columns[1:]] = psm[columns[1:]].applymap(lambda x: x.split('/')[0]).astype(np.float).replace(0, np.nan)

    # Match
    psm.columns = [columns[0]] + [ss.ix[s_dir, '%s-Biospecimen' % c[5:]] for c in columns[1:]]

    # Log 2 transform
    psm[psm.columns[1:]] = np.log2(psm[psm.columns[1:]])

    # Calculate ratios
    psm[psm.columns[1:]] = psm[psm.columns[1:]].subtract(psm['Internal Reference'], axis=0)

    # Drop interal reference
    psm = psm.drop('Internal Reference', axis=1)

    # Groupby protein
    psm = psm.groupby(columns[0]).max()

    # Append
    psms.append(psm)

    print s_dir

psms = concat(psms, axis=1)
psms.to_csv('%s/tables/pancan_proteomics_preprocessed_intensities.csv' % wd)
print psms
