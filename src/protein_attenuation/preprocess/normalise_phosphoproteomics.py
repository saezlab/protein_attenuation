#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

from protein_attenuation.utils import gkn
from pandas import DataFrame, read_csv


# -- Normalise proteomics
for f in ['cptac_phosphoproteomics', 'cptac_phosphoproteomics_corrected']:
    pp = read_csv('./data/%s.csv' % f, index_col=0)
    pp = pp[pp.count(1) > (pp.shape[1] * .5)]
    pp = DataFrame({i: gkn(pp.ix[i].dropna()).to_dict() for i in pp.index}).T
    pp.to_csv('./data/%s_normalised.csv' % f)
    print '[INFO] Done: %s' % f
