#!/bin/bash
#$ -pe omp 1
#$ -N plot_splines
#$ -V
#$ -j y


in_dir="../outputs_viirs"
in_stem="vi_str8"
start_year=2011
num_years=4

for t in `cat tiles.txt`;
do

Rscript ./read_extracts_comp.R $in_dir $in_stem $t $start_year $num_years
done
