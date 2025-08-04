# panarthropoda_gc_specification_evolution
This repository contains code and associated data for ancestral state reconstruction of germ plasm specification mode across Panarthropoda.
## supplementary_tables
Supplementary Tables S3, S4 and S9 (remaining tables are available in the Supplementary Materials document).
### Supplementary Table S3. Hypothesized transitions in germ plasm specification in Panarthropoda.
### Supplementary Table S4. Genomic and transcriptomic datasets used for phylogenetic inference.
### Supplementary Table S9. Results of search for oskar orthologs.
## Dependencies
R version 4.2.3
 
R packages: phytools version 2.4.21; ape version 5.7.1; phylolm version 2.6.5

Python 3.9.7 

Python packages: ete3 version 3.1.2; pandas version 2.2.3; numpy version 1.23.5; scipy 1.7.3

Software accessed via BioContainers Singularity containers: Transdecoder version 5.7.1; BLAST version 2.9.0; hmmer version 3.4; BUSCO version 5.8.3; IQ-Tree version 2.4.0. 

Additional software downloads: IQ-Tree version 3.0.1; PAML version 4.10.9

## Steps
### phylogram_constraints.ipynb:
Constructs a full species-level phylogram concatenating literature-derived phylogenies in the "constraint_tree" folder. Outputs 'constraint_trees/all_combined_species.nwk' used for Figure 5 and Supplementary Figure S1 of the manuscript.
### time_calibrated_busco_phylogenomics.ipynb:
Downloads genomic/transcriptomic datasets from phylogenetic_data_with_substitutions.tsv. Passes predicted proteins through the BUSCOphylogenomics pipeline (BUSCO_py/supermatrix/SUPERMATRIX.phylip), then to IQtree for constrained maximum likelihood topology search using a literature-derived constraint tree from phylogram_constraints.ipynb (outputs in BUSCO_py/supermatrix/ with "constrained_optimized" prefix), then to IQtree for the IQ2MC pipeline for hessian matrix computation (outputs in "BUSCO_py/supermatrix/" with mcmc_prep_no_partition prefix ), and finally PAML for divergence time estimation with MCMCtree (outputs in "BUSCO_py/supermatrix/run1R_mcmc_prep_no_partition.mcmctree.ctl").
### oskar_search.ipynb:
Searches for oskar homologs in genomic/transcriptomic datasets with hmmsearch and outputs oskar_data.tsv and oskar_proteins.fa with Oskar protein outputs. hmmsearch outputs in hmmsearch_results.
### maximum_likelihood_ancestral_reconstruction.R:
Fits two-state substitution models in phytools and performs maximum likelihood ancestry inference.
### germ_band_correlations.R:
Performs Pagel's test and phylogenetic logistic regression to test the association between germ cell specification mode and germ band type.
### oskar_fitphyloglm.R:
Performs phylogenetic logistic regression to test the association between oskar presence in a genome and germ band type.


