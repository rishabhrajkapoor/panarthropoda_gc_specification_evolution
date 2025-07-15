#!/bin/bash

curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$1/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED&filename=$1.zip" -H "Accept: application/zip"
unzip -o $1.zip
rm $1.zip