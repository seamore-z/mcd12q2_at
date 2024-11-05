#!/bin/bash
#$ -pe omp 1
#$ -N plot_splines
#$ -V
#$ -j y


in_dir="/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs_test2"
in_stem="c6_str5"
start_year=2000
num_years=6

for t in `cat tiles.txt`;
do

Rscript ./read_extracts.R $in_dir $in_stem $t $start_year $num_years
done
