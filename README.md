# panarthropoda_gc_specification_evolution
Code and associated data for ancestral state reconstruction of germ plasm specification mode across Panarthropoda
## supplementary_tables
Supplementary tables S3, S4 and S9 (remaining tables are are available in the paper supplement).

### Supplementary Table S9. Results of search for oskar orthologs. 
BUSCO scores may differ from Table S4 for the 7 unannotated genomes (Clogmia albipunctata, Mayetiola destructor, 
Macrocentus cingulum, Epiphyas postvittana, Nymphalis Antopa, Pyrrhocoris apterus, Tribolium 
confusum), as the BUSCO scores and oskar presence are determined from AUGUSTUS genome 
annotations. “Blondel_oskar_count” column indicates the total number of oskar orthologs 
previously identified across all datasets from this species (Blondel et al., 2021) Empty cells 
indicate that the species was not analyzed in Blondel et al., 2021 on oskar ortholog evolution. “oskar_present” column indicates that oskar was detected in the current analysis and/or 
in our previous work (Blondel et al., 2021). "data_type" indicates whether the proteins predictions were obtained directly from protein fastas on NCBI ("NCBI-annotated genome"), from TransDecoder run on transcriptome shotgun assembly nucleotide sequences ("TSA"), or from ab initio annotations generated with AUGUSTUS ("AUGUSTUS-annotated genome"). We were unable to find sequence data for species in the "sub_species" columns, and instead substituted the genomes/transcriptomes for species in the same taxonomic family. "oskar_accession" columns provide either NCBI protein accessions for NCBI-annotated genomes ("https://www.ncbi.nlm.nih.gov/protein/"), TSA nucleotide accessions for TSAs (https://www.ncbi.nlm.nih.gov/genbank/tsa/), or AUGUSTUS-generated accessions for AUGUSTUS-annotated genomes. All oskar protien sequences are availabel on github ("oskar_proteins.fa"). 


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
### phylogram_constraints.ipynb: 
construct a full species-level phylogram concatenating literature-derived phylogenies in the "constraint_tree" folder. Outputs 'constraint_trees/all_combined_species.nwk' used for Figure 9 and supplementary figure 1 of the manuscript. 

### time_calibrated_busco_phylogenomics.ipynb: 
Downloads genomic/transcriptomic datasets from phylogenetic_data_with_substitutions.tsv. Passes predicted proteins through the BUSCOphylogenomics pipeline (BUSCO_py/supermatrix/SUPERMATRIX.phylip) , then to IQtree for constrained maximum likelihood topology search using a literature-derived constraint tree from phylogram_constraints.ipynb (outputs in BUSCO_py/supermatrix/ with "constrained_optimized" prefix), then to IQtree for the IQ2MC pipeline for hessian matrix computation (outputs in "BUSCO_py/supermatrix/" with mcmc_prep_no_partition prefix ), and finally PAML for divergence time estimation with MCMCtree (outputs in "BUSCO_py/supermatrix/run1R_mcmc_prep_no_partition.mcmctree.ctl").

### oskar_search.ipynb
Search for oskar homologs in genomic/transcriptomic datasets with hmmsearch and outputs oskar_data.tsv and oskar_proteins.fa with oskar outputs. hmmsearch outputs in hmmsearch_results

### maximum_likelihood_ancestral_reconstruction.R
Fits two-state substitution models in phytools and performs maximum likelihood ancestry inference.

### germ_band_correlations.R
Perform's Pagel's test and phylogenetic logistic regression to test the association between germ cell specification mode and germ band type.

### oskar_fitphyloglm.r
Perform's  phylogenetic logistic regression to test the association between oskar presence and germ band type.
