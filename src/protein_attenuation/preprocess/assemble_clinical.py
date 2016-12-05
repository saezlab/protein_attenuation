#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

import os
import numpy as np
from pandas import read_csv, concat

# Source data folder
dwd = '/Users/emanuel/Projects/data/tcga/clinical/'

# Read directories
dirs = [i for i in os.listdir(dwd) if i.startswith('gdac.broadinstitute.')]
print len(dirs)

# Clinical data files
files = ['%s/%s/%s.merged_only_clinical_clin_format.txt' % (dwd, d, d.split('_')[1].split('.')[0]) for d in dirs]

# Read tables
dfs = []
for f in files:
    df = read_csv(f, sep='\t', index_col=0)
    df.ix['days_to_last_followup'] = df.ix[[i for i in df.index if 'days_to_last_followup' in i]].astype(np.float).max()
    df.ix['days_to_death'] = df.ix[[i for i in df.index if 'days_to_death' in i]].astype(np.float).max()
    df.ix['patient.bcr_patient_barcode'] = [i.upper() for i in df.ix['patient.bcr_patient_barcode']]
    dfs.append(df)

print len(dfs)

# Concatenate
df = concat(dfs, axis=1, join='inner').T.set_index('patient.bcr_patient_barcode')
print df.shape

# Annotate
df = df.dropna(subset=['patient.vital_status'])
df['time'] = [dlf if np.isnan(dd) else dd for dd, dlf in df[['days_to_death', 'days_to_last_followup']].values]
df['status'] = [0 if i.lower() == 'alive' else 1 for i in df['patient.vital_status']]
print df.shape

# Export
df.to_csv('%s/tcga_clinical.csv' % dwd)
df.to_csv('./data/tcga_clinical.csv')
print '[INFO] Exported'
