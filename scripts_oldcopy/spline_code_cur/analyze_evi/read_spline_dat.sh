#!/bin/bash
#$ -pe omp 1
#$ -l mem_total=98G
#$ -N read_bin
#$ -j y


Rscript ./read_bin.R $1 $2 $3 $4 $5



