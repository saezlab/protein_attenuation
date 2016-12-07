# Copyright (C) 2016  Emanuel Goncalves

####	Analysis
# Samples overlap
overlap:
	ipython src/protein_attenuation/overview/overlap.py



####	Regulatory associations
# Calculate the protein residuals having regressed-out transcriptomics measurements
residuals:
	ipython src/protein_attenuation/regressions/residuals_protein_transcript.py

# Regression associations within protein complexes (using different omics as input feature)
regressions:
	ipython src/protein_attenuation/regressions/regressions_cnv.py
	ipython src/protein_attenuation/regressions/regressions_transcriptomics.py
	ipython src/protein_attenuation/regressions/regressions_proteomics.py	



####	Data preprocessing
# Assemble TCGA and CPTAC data-sets for the overlapping samples. 
# COREAD (pmid 25043054), HGSC (pmid 27372738), BRCA (pmid 27251275)
assemble_datasets:
	ipython src/protein_attenuation/preprocess/assemble_clinical.py
	ipython src/protein_attenuation/preprocess/assemble_proteomics.py
	ipython src/protein_attenuation/preprocess/assemble_transcriptomics.py
	ipython src/protein_attenuation/preprocess/assemble_copy_number.py

# Merge, preprocess and normalise TCGA and CPTAC data-sets
preprocess_datasets:
	ipython src/protein_attenuation/preprocess/preprocess_proteomics.py
	ipython src/protein_attenuation/preprocess/preprocess_transcriptomics.py

	ipython src/protein_attenuation/preprocess/normalise_proteomics.py
	ipython src/protein_attenuation/preprocess/normalise_transcriptomics.py

# Preprocess cell lines protoemics. HGSC (pmid 27561551), BRCA (pmid 25892236)
preprocess_proteomics_celllines:
	ipython src/protein_attenuation/preprocess/preprocess_proteomics_hgsc_cell_lines.py
	ipython src/protein_attenuation/preprocess/preprocess_proteomics_brca_cell_lines.py

# Preprocess, merge and normalise proteomics, transcriptomics and copy-number for all the BRCA tumours
# inclusively samples failing in QC due to high protein degradation
# BRCA (pmid 27251275)
preprocess_brca_tumours:
	ipython src/protein_attenuation/preprocess/preprocess_brca_tumours_qc.py



####	Miscellaneous
# Merge highly similar CORUM complexes (overlaping with the measured proteins)
corum_redundancy:
	ipython src/protein_attenuation/corum_complexes_redundancy.py



#### General
clean:
	rm -rf '*.pyc'
	rm -rf 'gurobi.log'
	rm -rf '.DS_Store'

help:
	@echo "\t sfd"
