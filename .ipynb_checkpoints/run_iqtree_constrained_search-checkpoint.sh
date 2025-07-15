#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J BUSCO # Job name
#SBATCH -n 64 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem 100G
#SBATCH -t 01-00:00
#SBATCH --mail-user=rkapoor@g.harvard.edu
#SBATCH --mail-type=END 
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

cd BUSCO_py/supermatrix

singularity exec /cvmfs/singularity.galaxyproject.org/i/q/iqtree:2.4.0--h503566f_0 iqtree -m LG+G4 -s 'SUPERMATRIX.phylip' -p 'SUPERMATRIX.partitions.nex' -g '../constraint_tree_iqtree_species.tree' -pre 'constrained_optimized' -nt AUTO