# Copyright (C) 2016  Emanuel Goncalves

####	Data preprocessing
assemble_datasets:
	# Assemble TCGA and CPTAC data-sets for the overlapping samples. 
	# COREAD (pmid 25043054), HGSC (pmid 27372738), BRCA (pmid 27251275)
	ipython src/protein_attenuation/preprocess/assemble_clinical.py
	ipython src/protein_attenuation/preprocess/assemble_proteomics.py
	ipython src/protein_attenuation/preprocess/assemble_transcriptomics.py
	ipython src/protein_attenuation/preprocess/assemble_copy_number.py

preprocess_datasets:
	# Merge, preprocess and normalise TCGA and CPTAC data-sets
	ipython src/protein_attenuation/preprocess/preprocess_proteomics.py
	ipython src/protein_attenuation/preprocess/preprocess_transcriptomics.py

	ipython src/protein_attenuation/preprocess/normalise_proteomics.py
	ipython src/protein_attenuation/preprocess/normalise_transcriptomics.py

preprocess_proteomics_celllines:
	# Preprocess cell lines protoemics. HGSC (pmid 27561551), BRCA (pmid 25892236)
	ipython src/protein_attenuation/preprocess/preprocess_proteomics_hgsc_cell_lines.py
	ipython src/protein_attenuation/preprocess/preprocess_proteomics_brca_cell_lines.py

preprocess_brca_tumours:
	# Preprocess, merge and normalise proteomics, transcriptomics and copy-number for all the BRCA tumours
	# inclusively samples failing in QC due to high protein degradation
	# BRCA (pmid 27251275)
	ipython src/protein_attenuation/preprocess/preprocess_brca_tumours_qc.py

corum_redundancy:
	# Merge highly similar CORUM complexes (overlaping with the measured proteins)
	ipython src/protein_attenuation/corum_complexes_redundancy.py


####	Analysis
overlap:
	# Samples overlap
	ipython src/protein_attenuation/overview/overlap.py

regressions:
	# Regulatory associations: Regression associations within protein complexes (using different omics as input feature)
	ipython src/protein_attenuation/regressions/regressions_cnv.py
	ipython src/protein_attenuation/regressions/regressions_transcriptomics.py
	ipython src/protein_attenuation/regressions/regressions_proteomics.py
	ipython src/protein_attenuation/regressions/regression_associations.py

residuals:
	# Regulatory associations: Calculate the protein residuals having regressed-out transcriptomics measurements
	ipython src/protein_attenuation/regressions/residuals_protein_transcript.py

correlation:
	# Protein pairwise correlation
	ipython src/protein_attenuation/correlation/protein_pairs_correlation.py
	ipython src/protein_attenuation/correlation/protein_clustering.py

attenuation:
	# Protein attenuation potential
	ipython src/protein_attenuation/attenuation/protein_attenuation.py
	ipython src/protein_attenuation/attenuation/samples_protein_attenuation.py
	ipython src/protein_attenuation/attenuation/protein_attenuation_cell_lines.py
	ipython src/protein_attenuation/attenuation/samples_protein_attenuation_brca_qc_failed.py
	ipython src/protein_attenuation/attenuation/proteasome_inhibition.py
	ipython src/protein_attenuation/attenuation/drug_response.py


#### General
clean:
	rm -rf '*.pyc'
	rm -rf 'gurobi.log'
	rm -rf '.DS_Store'
