import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import wd, palette_cnv_number
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from sklearn.linear_model import LinearRegression
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict
from pandas import DataFrame, Series, read_csv, pivot_table, concat

# -- Imports
# CNV
cnv = read_csv('%s/data/tcga_cnv.tsv' % wd, sep='\t', index_col=0)
print 'cnv', cnv.shape

# Transcriptomics
transcriptomics = read_csv('%s/data/tcga_rnaseq_corrected_normalised.csv' % wd, index_col=0)
print 'transcriptomics', transcriptomics.shape

# Proteomics
proteomics = read_csv('%s/data/cptac_proteomics_corrected_normalised.csv' % wd, index_col=0)
print 'proteomics', proteomics.shape

# Clinical data
clinical = read_csv('%s/data/tcga_clinical.csv' % wd, index_col=0).dropna(subset=['time', 'status'])
clinical.index = ['%s-01' % i for i in clinical.index]
print 'clinical', clinical.shape

# Proteomics complex activities
c_prot = read_csv('%s/tables/protein_complexes_proteomics_activities.tsv' % wd, sep='\t')
c_prot = pivot_table(c_prot, index='complex', columns='sample', values='z')
print 'c_prot', c_prot.shape


# -- Complexes proteins
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot} for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Overlap
samples = set(cnv).intersection(clinical.index).intersection(transcriptomics).intersection(proteomics)
genes = set(cnv.index).intersection(transcriptomics.index).intersection(proteomics.index)
print 'samples', 'genes', len(samples), len(genes)


# Residuals
def protein_residual(g):
    y = proteomics.ix[g, samples].dropna()
    x = transcriptomics.ix[[g], y.index].T

    lm = LinearRegression().fit(x, y)
    y_ = y - lm.coef_[0] * x[g] - lm.intercept_

    return y_
residuals = DataFrame({g: protein_residual(g) for g in genes}).T
print 'residuals', residuals.shape


# -- Signif Px -> Py associations
ppairs = read_csv('./tables/ppairs_cnv_regulation.csv')
print ppairs


# px, py = 'COG3', 'COG2'
def survival_regression(px, py):
    df = {}
    for p in {p for c in corum if px in corum[c] and py in corum[c] for p in corum[c]}:
        if p in proteomics.index:
            x = proteomics.ix[p, samples].dropna()
            p25, p75 = np.percentile(x, 25), np.percentile(x, 75)
            df[p] = x.map(lambda i: 1 if i > p75 else (-1 if i < p25 else 0))

    df = DataFrame(df).mode(axis=1).apply(lambda x: min(x.min(), x.max(), key=abs), axis=1)
    df = concat([df, clinical.ix[samples, ['time', 'status']]], axis=1).dropna()

    cf = CoxPHFitter(normalize=False)
    cf.fit(df, 'time', event_col='status')
    print cf.print_summary()

    cindex = concordance_index(df['time'], -cf.predict_partial_hazard(df[[0]]).values.ravel(), df['status'])

    res = {'px': px, 'py': py, 'cindex': cindex}
    return res

c_survival = DataFrame([survival_regression(px, py) for px, py in ppairs[['px', 'py']].values])
print c_survival.sort('cindex')


px, py = c_survival.sort('cindex').tail(1)[['px', 'py']].values[0]


df = {}
for p in {p for c in corum if px in corum[c] and py in corum[c] for p in corum[c]}:
    if p in proteomics.index:
        x = proteomics.ix[p, samples].dropna()
        p25, p75 = np.percentile(x, 25), np.percentile(x, 75)
        df[p] = x.map(lambda i: 1 if i > p75 else (-1 if i < p25 else 0))

df = DataFrame(df).mode(axis=1).apply(lambda x: min(x.min(), x.max(), key=abs), axis=1)
df = concat([df, clinical.ix[samples, ['time', 'status']]], axis=1).dropna()

samples_up = set(df[df[0] == 1].index)
samples_down = set(df[df[0] == -1].index)

logrank = logrank_test(
    df.ix[samples_up, 'time'], df.ix[samples_down, 'time'],
    df.ix[samples_up, 'status'], df.ix[samples_down, 'status']
)


sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .0, 'xtick.direction': 'out', 'ytick.direction': 'out'}, font_scale=0.75)
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(df.ix[samples_up, 'time'], event_observed=df.ix[samples_up, 'status'], label='Up')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[1], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(df.ix[samples_down, 'time'], event_observed=df.ix[samples_down, 'status'], label='Down')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[-1], ax=ax, show_censors=True, censor_styles={'ms': 5})

sns.despine(ax=ax)

ax.set_title('Logrank test: %.2e' % (logrank.p_value))
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(4, 3)
plt.savefig('%s/reports/survival_protein_complex.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
