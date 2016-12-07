# Copyright (C) 2016  Emanuel Goncalves

WD=./

# Scripts to assemble TCGA and CPTAC data-sets for the overlapping samples
assemble_datasets:
	ipython src/protein_attenuation/preprocess/assemble_clinical.py
	ipython src/protein_attenuation/preprocess/assemble_proteomics.py
	ipython src/protein_attenuation/preprocess/assemble_transcriptomics.py
	ipython src/protein_attenuation/preprocess/assemble_copy_number.py

# Scripts to merge, preprocess and normalise TCGA and CPTAC data-sets
preprocess_datasets:
	ipython src/protein_attenuation/preprocess/preprocess_proteomics.py
	ipython src/protein_attenuation/preprocess/preprocess_transcriptomics.py

	ipython src/protein_attenuation/preprocess/normalise_proteomics.py
	ipython src/protein_attenuation/preprocess/normalise_transcriptomics.py

# Script to merge highly similar CORUM complexes (overlaping with the measured proteins)
corum_redundancy:
	ipython src/protein_attenuation/corum_complexes_redundancy.py

clean:
	rm -rf '*.pyc'

help:
	@echo "\t sfd"
