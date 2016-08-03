import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd
from pandas import DataFrame, Series, read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- ID conversion maps
uniprot_human = read_uniprot_genename()
uniprot_mouse = read_uniprot_genename(os='Mus musculus')
print 'uniprot_human', 'uniprot_mouse', len(uniprot_human), len(uniprot_mouse)

mouse_human = {m: {h for h in uniprot_human if uniprot_mouse[m][1].split('_')[0] == uniprot_human[h][1].split('_')[0]} for m in uniprot_mouse}
mouse_human = {m: mouse_human[m] for m in mouse_human if len(mouse_human[m]) == 1}
print 'mouse_human', len(mouse_human)

# -- Turnover rates
turnover = read_csv('%s/files/proteins_turnovers.csv' % wd).dropna(subset=['Uniprot IDs'])
turnover['Uniprot IDs human'] = [';'.join([list(mouse_human[p])[0] for p in i.split(';') if p in mouse_human]) for i in turnover['Uniprot IDs']]
turnover.to_csv('%s/files/proteins_turnovers_preprocessed.csv' % wd)
print turnover
