#!/bin/sh

# Requires Wget

type wget >/dev/null 2>&1 || { echo >&2 "Wget is required"; exit 1; }

wget -O datasets/nrdb90.gz  "ftp://ftp.ebi.ac.uk/pub/databases/nrdb90/nrdb90.gz"
gzip -d datasets/nrdb90.gz 

wget -O datasets/rna.fa.gz  "ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/RNA/rna.fa.gz"
gzip -d datasets/rna.fa.gz


