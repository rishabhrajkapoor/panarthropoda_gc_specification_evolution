# panarthropoda_gc_specification_evolution
Code and associated data for ancestral state reconstruction of germ plasm specification mode across Panarthropoda

## Dependencies 
R version 4.2.3
R packages: phytools version 2.4.21; ape version 5.7.1; phylolm version 2.6.5;

Python 3.9.7
Python packages: ete3 version 3.1.2; pandas version 2.2.3; numpy version 1.23.5; scipy 1.7.3

Software accessed via biocontainers singularity containers
Transdecoder version 5.7.1;
blast version 2.9.0;
hmmer version 3.4;
busco version 5.8.3;
IQ-Tree version 2.4.0;

Additional software downloads:
IQ-Tree version 3.0.1;
PAML version 4.10.9

## Steps
phylogram_constraints.ipynb: construct a full species-level phylogram concatenating literature-derived phylogenies in the "constraint_tree" folder. Outputs 'constraint_trees/all_combined_species.nwk' used for Figure 9 and supplementary figure 1 of the manuscript

time_calibrated_busco_phylogenomics.ipynb: Downloads genomic/transcriptomic datasets from phylogenetic_data_with_substitutions.tsv. Passes predicted proteins through the BUSCOphylogenomics pipeline, then to 

