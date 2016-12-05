#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import numpy as np
from protein_attenuation.utils import gkn
from pandas import read_csv, DataFrame
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Uniprot
uniprot = read_uniprot_genename()

# -- Cell lines: proteomics
proteomics = read_csv('./data/brca_cell_lines_proteomics.csv')

# Drop tumour samples
proteomics = proteomics[[c for c in proteomics if 'Tumor' not in c]]

# Parse uniprot IDs
proteomics['Uniprot Entry'] = [i.split('|')[1] for i in proteomics['Uniprot Entry']]
proteomics = proteomics[[i in uniprot for i in proteomics['Uniprot Entry']]]
proteomics['Uniprot Entry'] = [uniprot[i][0] for i in proteomics['Uniprot Entry']]
proteomics = proteomics.groupby('Uniprot Entry').mean()

# Log2
proteomics = np.log2(proteomics)

# Average replicates
proteomics.columns = [i.split(' ')[0] for i in proteomics]
proteomics = DataFrame({c: proteomics[c].mean(1) for c in set(proteomics)})

# Proteins measured in at least 50% of the samples
proteomics = proteomics[proteomics.count(1) > proteomics.shape[1] * .5]

# Normalise
proteomics = DataFrame({p: gkn(proteomics.ix[p].dropna()) for p in proteomics.index}).T

# Export
proteomics.to_csv('./data/brca_cell_lines_proteomics_preprocessed.csv')
print 'proteomics', proteomics
