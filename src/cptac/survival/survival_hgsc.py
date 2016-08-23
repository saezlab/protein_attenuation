import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from cptac import palette_cnv_number, palette
from pandas import DataFrame, Series, read_csv, concat
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index, k_fold_cross_validation
from lifelines.statistics import logrank_test
from statsmodels.duration.hazard_regression import PHReg
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.utils.corumdb import get_complexes_name, get_complexes_dict


# -- Clinical
clinical = read_csv('./data/hgsc_clinical.csv', index_col=0)
clinical.index = ['%s-01' % i for i in clinical.index]

clinical['time'] = [i for i in clinical['daystodeath_or_LFU']]
clinical['status'] = [0 if i == 'LIVING' else 1 for i in clinical['vital_status']]

clinical = clinical[clinical['time'] < (10 * 365)]
print clinical.shape


# -- Proteomics
hgsc_samplesheet = read_csv('./data/OV_All_clinical_features_TCGAbiotab_CPTAC_S020.csv')
hgsc_samplesheet['id'] = ['-'.join(i.split('-')[1:4]) for i in hgsc_samplesheet['TCGA barcode']]
hgsc_samplesheet = hgsc_samplesheet.groupby('id')['TCGA barcode'].first()

# PNNL
hgsc_pnnl = read_csv('./data/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv', sep='\t', index_col=0)

hgsc_pnnl = hgsc_pnnl.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)
hgsc_pnnl = hgsc_pnnl.drop(['Mean', 'Median', 'StdDev'], axis=0)

hgsc_pnnl = hgsc_pnnl[[c for c in hgsc_pnnl if c.endswith(' Unshared Log Ratio')]]
hgsc_pnnl.columns = [c.split(' ')[0] for c in hgsc_pnnl]

hgsc_pnnl = hgsc_pnnl.loc[:, [i in hgsc_samplesheet.index for i in hgsc_pnnl]]
hgsc_pnnl.columns = [hgsc_samplesheet.ix[i][:15] for i in hgsc_pnnl]

hgsc_pnnl = hgsc_pnnl.dropna()
print hgsc_pnnl.shape

# JHU
hgsc_jhu = read_csv('./data/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv', sep='\t', index_col=0)

hgsc_jhu = hgsc_jhu.drop(['Description', 'Organism', 'Chromosome', 'Locus'], axis=1)
hgsc_jhu = hgsc_jhu.drop(['Mean', 'Median', 'StdDev'], axis=0)

hgsc_jhu = hgsc_jhu[[c for c in hgsc_jhu if c.endswith(' Unshared Log Ratio')]]
hgsc_jhu.columns = [c.split(' ')[0] for c in hgsc_jhu]

hgsc_jhu = hgsc_jhu.loc[:, [i in hgsc_samplesheet.index for i in hgsc_jhu]]
hgsc_jhu.columns = [hgsc_samplesheet.ix[i][:15] for i in hgsc_jhu]

hgsc_jhu = hgsc_jhu.dropna()
print hgsc_jhu.shape


# -- Overlap
proteins = list(set(hgsc_pnnl.dropna().index).intersection(hgsc_jhu.dropna().index))

hgsc_pnnl = hgsc_pnnl.ix[proteins, list(set(hgsc_pnnl).intersection(clinical.index))]
print hgsc_pnnl.shape

hgsc_jhu = hgsc_jhu.ix[proteins, list(set(hgsc_jhu).difference(hgsc_pnnl).intersection(clinical.index))]
print hgsc_jhu.shape


# -- Protein complexes
uniprot = read_uniprot_genename()
corum = {k: {uniprot[p][0] for p in v if p in uniprot}.intersection(proteins) for k, v in get_complexes_dict().items()}
corum = {k: corum[k] for k in corum if len(corum[k]) > 1}
corum_n = get_complexes_name()
print 'corum', len(corum)


# -- Signif Px -> Py associations
ppairs = read_csv('./tables/ppairs_cnv_regulation.csv')
p_complexes = {(px, c) for px, py in ppairs[['px', 'py']].values for c in corum if px in corum[c] and py in corum[c]}
print ppairs


# --
# c = 1181
c_survival = {}
for px, c in p_complexes:
    c_proteins = list(corum[c].intersection(proteins))

    # Train
    x_train = hgsc_pnnl.ix[c_proteins].T
    x_train = concat([x_train, clinical[['time', 'status']]], axis=1).dropna()

    model = PHReg(x_train['time'], x_train[c_proteins], x_train['status'], ties='efron').fit_regularized(alpha=0.01)

    c_proteins = Series(dict(zip(*(x_train.columns, model.params))))
    c_proteins = c_proteins[c_proteins != 0]

    if len(c_proteins) > 0:
        c_proteins = list(set(c_proteins.index))

        # Train
        x_train = hgsc_pnnl.ix[c_proteins].T
        x_train = concat([x_train, clinical[['time', 'status']]], axis=1).dropna()

        x_test = hgsc_jhu.ix[c_proteins].T
        x_test = concat([x_test, clinical[['time', 'status']]], axis=1).dropna()

        model = CoxPHFitter(normalize=False, penalizer=.0, tie_method='Efron').fit(x_train, 'time', event_col='status')

        hazards = model.predict_partial_hazard(x_test[c_proteins])
        cindex = concordance_index(x_test['time'], -hazards, x_test['status'])

        print model.print_summary()
        print corum_n[c], cindex

        c_survival[c] = {'complex': c, 'px': px, 'cindex': cindex, 'len': len(c_proteins), 'hazards': hazards}

c_survival = DataFrame(c_survival).T
c_survival['name'] = [corum_n[i] for i in c_survival.index]
print c_survival.sort('cindex', ascending=False)[['cindex', 'name', 'len', 'px']]


# -- Plot
c = 1181

# df = c_survival.ix[c, 'model'].predict_partial_hazard(df[c_proteins])[0]
df = c_survival.ix[c, 'hazards'][0]
df = df.map(lambda i: 1 if i > np.percentile(df, 50) else -1)
df = concat([df, clinical[['time', 'status']]], axis=1).dropna()

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


# --
px = 'RPA2'

df = DataFrame([hz[0].map(lambda i: 1 if i > np.percentile(hz[0], 50) else -1) for hz in c_survival[c_survival['px'] == px]['hazards']])
df = df.mode().apply(lambda x: min(x.min(), x.max(), key=abs))
df = concat([df, clinical[['time', 'status']]], axis=1).dropna()

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
