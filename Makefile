# Copyright (C) 2016  Emanuel Goncalves

WD=./

overlap:
	python src/protein_attenuation/figures/overlap.py

clean:
	rm -rf gurobi.log
	rm -rf .DS_Store

help:
	@echo "\t sfd"