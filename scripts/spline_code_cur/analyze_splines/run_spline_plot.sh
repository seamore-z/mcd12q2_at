#!/bin/bash
#$ -pe omp 1
#$ -N plot_splines
#$ -V
#$ -j y
#$ -l h_rt=48:00:00


in_dir="../MCD12I5"
in_stem="MCD12I5"
start_year=2001
num_years=16

for t in `cat tiles.txt`;
do

Rscript ./read_extracts.R $in_dir $in_stem $t $start_year $num_years
done
