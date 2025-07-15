#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J BUSCO # Job name
#SBATCH -n 50 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem 50G
#SBATCH -t 00-01:00
#SBATCH --mail-user=rkapoor@g.harvard.edu
#SBATCH --mail-type=END 
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL

n=$1
dir="/net/bos-nfsisilon/ifs/rc_labs/extavour_lab/rkapoor/IF_project_iteration2/TSA_transdec/$n/$n.fasta.transdecoder.pep"
singularity exec /cvmfs/singularity.galaxyproject.org/b/u/busco:5.8.3--pyhdfd78af_1 busco -i $dir -l arthropoda_odb12 -c 50 -m prot -f -o ./BUSCO_outputs/$n --offline 