# Get NBAR data
tiles <- c('h10v02')  # Test tile 1
# tiles <- c('h07v03','h08v03','h09v02','h09v03','h10v02','h10v03','h11v02',
#            'h12v01','h12v02','h13v01','h13v02')  # AK Tundra tiles

setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/e4ftl01.cr.usgs.gov/')
for(i in seq(length(tiles))){
  for(year in 2016:2023){
    system(paste('qsub -N NB2',tiles[i],'_',year,' -V -pe omp 2 -l h_rt=12:00:00 /projectnb/modislc/users/seamorez/MCD12Q2/run_nbara2.sh ', tiles[i],year,sep=''))
    system(paste('qsub -N NB4',tiles[i],'_',year,' -V -pe omp 2 -l h_rt=12:00:00 /projectnb/modislc/users/seamorez/MCD12Q2/run_nbara4.sh ', tiles[i],year,sep=''))
  }  
}


# Get NBAR omitted
setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/e4ftl01.cr.usgs.gov/')
for(i in seq(length(tiles))){
  for(year in 2016:2023){
    system(paste('qsub -N om2_',tiles[i],'_',year,' -V -pe omp 2 -l h_rt=12:00:00 /projectnb/modislc/users/seamorez/MCD12Q2/run_nbara2_omit.sh ', tiles[i],year,sep=''))
    system(paste('qsub -N om4_',tiles[i],'_',year,' -V -pe omp 2 -l h_rt=12:00:00 /projectnb/modislc/users/seamorez/MCD12Q2/run_nbara4_omit.sh ', tiles[i],year,sep=''))
  }
}

# Run "change_folder_name.R" code
setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/')
for(i in seq(length(tiles))){
  for(year in 2016:2023){
    system(paste('qsub -N cfld',tiles[i],'_',year,' -V -pe omp 2 -l h_rt=12:00:00 /projectnb/modislc/users/seamorez/MCD12Q2/run_change_folder.sh ', tiles[i], year,sep=''))
  }
}

# Run spline code 
# tiles <- c('h10v02')

setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/')
for(i in seq(length(tiles))){
  #Create output directory if it doesn't exist
  for(start_yr in 2015:2022){
    outDir <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/',tiles[i],'/',(start_yr+1))
    if (!dir.exists(outDir)) {dir.create(outDir)}
    system(paste('qsub /projectnb/modislc/users/seamorez/MCD12Q2/run_spline_fromJosh.sh ', tiles[i],' ',start_yr,sep=''))    
  }
}


# To see splined results
setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/log/')
for(i in seq(length(tiles))){
  for(year in 2016:2023){
    for(j in 0:0){
      system(paste('qsub -N ck_sp -V -pe omp 4 -l h_rt=03:00:00 -l mem_total=32G /projectnb/modislc/users/seamorez/MCD12Q2/run_ck_spline.sh ',tiles[i],year,j,sep=''))  
    }
  }  
}


# Sub AnnualPhenology
setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/pheno/')
for(i in seq(length(tiles))){
  #Create output directory if it doesn't exist
  # outDir <- paste0('/projectnb/modislc/users/mkmoon/LCD_C6/v7/pheno/bu/',tiles[i],'/')
  # if (!dir.exists(outDir)) {dir.create(outDir)}
  for(year in 2016:2023){
    system(paste('qsub -V -pe omp 16 -l h_rt=06:00:00 -l mem_total=98G /projectnb/modislc/users/seamorez/MCD12Q2/MCD12Q2C6_subAnnualPhenology.sh ',tiles[i],' ',year,sep=''))  
  }
}


###################
# Check NBAR
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/nbar/')
for(i in 1:3){
  system(paste('qsub -N NbCk',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_check_nbar.sh ',tiles[i],sep=''))  
}


###################
# Get MCD12Q1 2020
setwd('/projectnb/modislc/projects/sat/data/mcd12q/q1')
for(year in 2002:2020){
  system(paste('qsub -N Q1',year,' -V -pe omp 8 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_12q1.sh ',year,sep=''))  
}


###################
# Compare with C5 - quality check
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/C5_C6_1/')
for(i in 1:6){
    for(j in 1:10){
      system(paste('qsub -V -l h_rt=01:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_compare.sh ',i,j,sep=''))    
    }
}

# Compare with C5 - hist
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/hist/')
for(j in 1:10){
  system(paste('qsub -V -l h_rt=02:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_compare_hist.sh ',j,sep=''))    
}

# Compare with C5 - 1to1
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/1_to_1/')
for(j in 1:10){
  system(paste('qsub -V -l h_rt=02:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_compare_1to1.sh ',j,sep=''))    
}


#### Load C6 for phenocam comparison
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/load_c6/')
for(i in 1:2){
  system(paste('qsub -V -pe omp 4 -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_load_c6.sh ',i,sep=''))    
}


#### MSLSP comparison
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/hls/')
for(i in 1:6){
  for(j in 1:4){
    system(paste('qsub -V -pe omp 2 -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_hls.sh ',i,j,sep=''))      
  }
}




####################################################
# V7
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v7/')
system('qsub -V -pe omp 2 -j y -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/download_wget.sh')    


####################################################
# V8
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v8/figure/')
for(tt in 1:3){
  for(yy in 2003:2019){
    system(paste('qsub -V -pe omp 2 -j y -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_v8.sh ',tt,yy,sep=''))        
  }
}

