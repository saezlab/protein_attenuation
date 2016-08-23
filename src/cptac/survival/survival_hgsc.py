import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette_cnv_number, palette
from pandas import DataFrame, Series, read_csv, concat
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict


# -- Imports
# Clinical data
cohorts = read_csv('./data/hgsc_cohort_sample.csv').groupby('bcr_patient_barcode').agg(lambda x: '_'.join(set(x)))

clinical = read_csv('./data/hgsc_clinical.csv', index_col=0)
clinical['cohort'] = [cohorts.ix[i, 'PCC'] for i in clinical.index]
clinical.index = ['%s-01' % i for i in clinical.index]
clinical['time'] = [i for i in clinical['daystodeath_or_LFU']]
clinical['status'] = [0 if i == 'LIVING' else 1 for i in clinical['vital_status']]

# Proteomics
proteomics = read_csv('./data/cptac_proteomics_corrected_normalised.csv', index_col=0)
print 'proteomics', proteomics.shape


# -- Overlap
proteins, samples = set(proteomics.index), set(clinical.index).intersection(proteomics)
print 'proteins', 'samples', len(proteins), len(samples)


# -- Protein complexes
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteomics.index) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Signif Px -> Py associations
ppairs = read_csv('./tables/ppairs_cnv_regulation.csv')
p_complexes = {c for px, py in ppairs[['px', 'py']].values for c in corum if px in corum[c]}
print ppairs

# --
# c = 351

c_survival = {}
for c in p_complexes:
    df = proteomics.ix[corum[c], samples].dropna().T
    c_proteins = list(df)

    if df.shape[1] > 0:
        df = concat([df, clinical.ix[samples, ['time', 'status']]], axis=1).dropna()

        phm = PHReg(df['time'], df[c_proteins], df['status'], ties='efron').fit_regularized(alpha=.01)
        phm_cindex = concordance_index(df['time'], -phm.predict(df[c_proteins], pred_type='hr').predicted_values, df['status'])
        print phm.summary()
        print phm_cindex

        # cf = CoxPHFitter(normalize=False).fit(df, 'time', event_col='status')
        # cf_cindex = concordance_index(df['time'], -cf.predict_partial_hazard(df[c_proteins]), df['status'])
        # print cf.print_summary()

        c_survival[c] = {'complex': c, 'c_index': phm_cindex, 'len': len(c_proteins)}
c_survival = DataFrame(c_survival).T
c_survival['name'] = [corum_n[i] for i in c_survival.index]
print c_survival.sort('c_index', ascending=False)


# -- Plot
c = 846

df = proteomics.ix[corum[c], samples].dropna().T
c_proteins = list(df)

df = concat([df, clinical.ix[samples, ['time', 'status']]], axis=1).dropna()

phm = PHReg(df['time'], df[c_proteins], df['status'], ties='efron').fit_regularized(alpha=.01)
phm_cindex = concordance_index(df['time'], -phm.predict(df[c_proteins], pred_type='hr').predicted_values, df['status'])
print phm.summary()
print phm_cindex


# df = c_survival.ix[c, 'model'].predict_partial_hazard(df[c_proteins])[0]
df = Series(dict(zip(*(df.index, phm.predict(df[c_proteins], pred_type='hr').predicted_values))))
df = df.map(lambda i: 1 if i > np.percentile(df, 50) else -1)
df = concat([df, clinical.ix[samples, ['time', 'status']]], axis=1).dropna()

samples_up = set(df[df[0] == 1].index)
samples_down = set(df[df[0] == -1].index)

logrank = logrank_test(
    df.ix[samples_up, 'time'], df.ix[samples_down, 'time'],
    df.ix[samples_up, 'status'], df.ix[samples_down, 'status']
)

sns.set(style='ticks', context='paper', rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'in', 'ytick.direction': 'in'}, font_scale=0.75)
ax = plt.subplot(111)

kmf = KaplanMeierFitter()

kmf.fit(df.ix[samples_up, 'time'], event_observed=df.ix[samples_up, 'status'], label='Up')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[1], ax=ax, show_censors=True, censor_styles={'ms': 5})

kmf.fit(df.ix[samples_down, 'time'], event_observed=df.ix[samples_down, 'status'], label='Down')
kmf.plot(ci_force_lines=False, color=palette_cnv_number[-1], ax=ax, show_censors=True, censor_styles={'ms': 5})

sns.despine(ax=ax)

ax.set_title('%s\nLogrank test: %.2e' % (corum_n[c], logrank.p_value))
ax.set_xlabel('Timeline (days)')
ax.set_ylabel('Survival fraction')

plt.gcf().set_size_inches(5, 3)
plt.savefig('./reports/survival_protein_complex.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Done'
