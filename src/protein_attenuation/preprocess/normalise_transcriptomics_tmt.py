#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from protein_attenuation.utils import gkn
from pandas import DataFrame, read_csv


# -- Normalise transcriptomics
for f in ['tcga_rnaseq_tmt_corrected']:
    transcriptomics = read_csv('./data/%s.csv' % f, index_col=0)
    transcriptomics = DataFrame({i: gkn(transcriptomics.ix[i].dropna()).to_dict() for i in transcriptomics.index}).T
    transcriptomics.to_csv('./data/%s_normalised.csv' % f)
    print '[INFO] Done: %s' % f
