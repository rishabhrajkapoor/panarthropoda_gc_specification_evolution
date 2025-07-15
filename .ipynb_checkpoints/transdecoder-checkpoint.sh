#!/bin/bash

n=$1
##directory to store transdecoder results
dir="/net/bos-nfsisilon/ifs/rc_labs/extavour_lab/rkapoor/IF_project_iteration2/TSA_transdec/"$n
mkdir dir
singularity exec /cvmfs/singularity.galaxyproject.org/t/r/transdecoder:5.7.1--pl5321hdfd78af_0 TransDecoder.LongOrfs -t "TSA/"$n".fasta" --output_dir $dir --complete_orfs_only
singularity exec /cvmfs/singularity.galaxyproject.org/t/r/transdecoder:5.7.1--pl5321hdfd78af_0  TransDecoder.Predict -t "TSA/"$n".fasta" --output_dir $dir