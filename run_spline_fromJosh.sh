#!/bin/bash

module load hdf

#$ -N smooth_spline
#$ -V
#$ -m ae
#$ -pe omp 16
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_total=98G

if [ -z "$1" ]
then
	echo Usage \"./run_spline_fromJosh.sh in_tile
	exit 1
else
	in_tile=$1
fi

echo "Writing to local disk on $TMPDIR"

## where the ouputs will be produced
# out_dir="../outputs_global"
out_year=$(($2+1))
out_dir="/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/$in_tile"
## name for the outputs will start with this string
out_stem="sz_v5"

## where the step extract is held - optional
# extract_path="/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/run_dir/all_step.csv"
extract_path="NONE"
## we don't read in every day of inputs instead this number is used to determine every nth observation we will read in 
stride=5
## the number of rows read in at one time - this will greatly impact the memory usage of the process 
## a chunk of 800 will use about 75GB of memory, a chunk of 600 will use about 55GB of memory, etc
chunk_size=600
## whether to output the quantile summaries of each band - used for the land cover classification inputs
quant_flag=1
## the number of years to read in at once to be splined - this should stay constant
## will produce num_yrs-2 outputs
num_yrs=5
## the starting year we will read in, the outputs will start at start_yr+1
start_yr=$2
first_yr=$(( start_yr + 1 ))

## location of input nbar and qa data
# nbar_dir="/projectnb/modislc/data/mcd12_in/c6/mcd43a4"
nbar_dir="/projectnb2/modislc/projects/sat/data/mcd43a4_link"
# qa_dir="/projectnb/modislc/data/mcd12_in/c6/mcd43a2"
qa_dir="/projectnb2/modislc/projects/sat/data/mcd43a2_link"




echo "Starting processing for $in_tile for output year $first_yr."
## to only process the file if it doesnt exist
#if [ ! -f $out_dir/${out_stem}.${in_tile}.$first_yr.evi2.bip ] 
#then
#	echo "No outputs for tile $in_tile. Will process"

	# echo ../modis_smooth_t/bin/modis_smooth_t --tile_id $in_tile --chunk $chunk_size --smooth 1,2,4,5,6,7,8 --prefilter 8 --residual 8 --quantile $quant_flag --start_year $start_yr --num_yrs $num_yrs --stride $stride --out_stem $out_stem --out_dir $out_dir --nbar_dir $nbar_dir --qa_dir $qa_dir --extract_path $extract_path
	echo /projectnb/modislc/users/seamorez/MCD12Q2/scripts/spline_code_cur/bin/modis_smooth_t --tile_id $in_tile --chunk $chunk_size --smooth 1,2,4,5,6,7,8 --prefilter 8 --residual 8 --quantile $quant_flag --start_year $start_yr --num_yrs $num_yrs --stride $stride --out_stem $out_stem --out_dir $out_dir --nbar_dir $nbar_dir --qa_dir $qa_dir --extract_path $extract_path
	# ../modis_smooth_t/bin/modis_smooth_t --tile_id $in_tile --chunk $chunk_size --smooth 1,2,4,5,6,7,8 --prefilter 8 --residual 8 --quantile $quant_flag --start_year $start_yr --num_yrs $num_yrs --stride $stride --out_stem $out_stem --out_dir $out_dir --nbar_dir $nbar_dir  --qa_dir $qa_dir --extract_path $extract_path
	/projectnb/modislc/users/seamorez/MCD12Q2/scripts/spline_code_cur/bin/modis_smooth_t --tile_id $in_tile --chunk $chunk_size --smooth 1,2,4,5,6,7,8 --prefilter 8 --residual 8 --quantile $quant_flag --start_year $start_yr --num_yrs $num_yrs --stride $stride --out_stem $out_stem --out_dir $out_dir --nbar_dir $nbar_dir  --qa_dir $qa_dir --extract_path $extract_path
#fi 

echo "Now finished processing tile $in_tile."
    


