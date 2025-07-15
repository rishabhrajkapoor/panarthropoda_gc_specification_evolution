#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J BUSCO # Job name
#SBATCH -n 50 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem 200G
#SBATCH -t 00-00:30
#SBATCH --mail-user=rkapoor@g.harvard.edu
#SBATCH --mail-type=END 
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

n=$1

dir="ncbi_dataset/data/$n/*.fna"

singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.3--pyhdfd78af_1 busco -i $dir -l arthropoda_odb12 -c 50 -m genome -f -o ./BUSCO_outputs/$n --offline 