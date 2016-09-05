#!/usr/bin/env python
# Copyright (C) 2016  Emanuel Goncalves

from cptac.utils import gkn
from pandas import DataFrame, read_csv


# -- Normalise proteomics
for f in ['cptac_proteomics', 'cptac_proteomics_corrected']:
    proteomics = read_csv('./data/%s.csv' % f, index_col=0)
    proteomics = proteomics[proteomics.count(1) > (proteomics.shape[1] * .5)]
    proteomics = DataFrame({i: gkn(proteomics.ix[i].dropna()).to_dict() for i in proteomics.index}).T
    proteomics.to_csv('./data/%s_normalised.csv' % f)
    print '[INFO] Done: %s' % f
