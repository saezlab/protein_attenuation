# -- Paths
wd = '/Users/emanuel/Projects/projects/cptac/'
data = '/Volumes/EBI BACKUP/CPTAC'

# -- Colors
default_color = '#808080'

palette = {'BRCA': '#f39c12', 'COREAD': '#bdc3c7', 'HGSC': '#f198a4', 'Pancancer': '#34495e'}
palette_dbs = {'All': '#d8d8d8', 'CORUM': '#e39e54', 'STRING': '#d64d4d', 'BioGRID': '#4d7358'}
palette_survival = {'high': '#d62d20', 'low': '#d8d8d8'}
palette_binary = {0: '#767676', 1: '#e74c3c'}
palette_cnv = {'neutral': '#767676', 'depletion': '#e74c3c', 'amplification': '#2ecc71'}
palette_cnv_number = {-2: '#e74c3c', -1: '#F0938A', 0: '#d5d5d5', 1: '#7bc399', 2: '#239c56'}

# -- Default filters
genomic_mod = {
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Indel',
    'Missense',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Nonstop_Mutation',
    'Splice_Site_Del',
    'Splice_Site_Ins',
    'Splice_Site_SNP'
}