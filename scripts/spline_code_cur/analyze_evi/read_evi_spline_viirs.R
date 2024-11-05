rm(list=ls())

args = commandArgs(trailingOnly=T)
debug_flag = 0
if(debug_flag==0) {
  in_stem = args[1]
  in_dir = args[2]
  in_year = as.numeric(args[3])
  in_tile = args[4]
  cur_pix = as.numeric(args[5]) 
} else {
  in_stem = "c6_str4_0.3"
  in_dir = "/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs_test2"
  in_year = 2005
  in_tile = "h12v04"
  cur_pix = 2501
}

##############################
## start function definitions #
##############################

## this function reads an integer image stack
read_int_file <- function(in_file, dsize, nbands, bsq_flag) {
  
  ## get the total dimensions of file = nrow*ncol*dsize*nbands
  f_size <- file.info(in_file)$size
  ## open the file
  f <- file(in_file,"rb")
  ## dsize should be an int() 2 bytes
  tot_size <- (f_size/dsize)
  ## can deal with signed or unsigned integers
  if(dsize < 3) {
  	temp <- readBin(f,integer(),n=tot_size,endian= "little",size=dsize,signed=T)
	} else {
	temp <- readBin(f,integer(),n=tot_size,endian= "little",size=dsize)
   }
  close(f)
  
  ## if bsq_flag is 1 then we read by row
  ## else we read by column
   ## re-order the temp array into a matrix
   if(bsq_flag==1) { byr_flag = FALSE } else { byr_flag = TRUE }
   temp <- matrix(temp, ncol=nbands, byrow=byr_flag)
 
  return(temp)
}

plot_pix <- function(evi2,res,flag,cur_pix, flag_on) {
  
  evi2_pix = evi2[cur_pix,]
  res_pix = res[cur_pix,]
  flag_pix = flag[cur_pix,]
  
  fill_val = 32767
  evi2_pix[evi2_pix == fill_val] = NA
  res_pix[res_pix == fill_val] = NA
  
  evi2_pix = evi2_pix/10000
  res_pix = res_pix/10000
  
  orig_dat = evi2_pix - res_pix
  
  ndays = 365
  qa_col = array(NA,ndays)
  snow_pch = array(NA,ndays)
  
  ## make a dates array for the x-axis of our plots
  dates = seq(1,ndays)
  
  pix_y = floor((cur_pix-1) / 2400)
  pix_x = cur_pix - (pix_y*2400) - 1
  
  if(flag_on == 1) {
    qa_col[flag_pix==0] = "blue"
    qa_col[flag_pix==1] = "cyan"
    qa_col[flag_pix==2] = "orange"
    qa_col[flag_pix==3] = "red"
    qa_col[flag_pix==4] = "purple"
    qa_col[flag_pix>4] = NA
    
    snow_pch[flag_pix==0] = 15
    snow_pch[flag_pix==1] = 0
    snow_pch[flag_pix==2] = 15
    snow_pch[flag_pix==3] = 0
    snow_pch[flag_pix==4] = 0
  } else {
    snow_pch = 15
    qa_col = "blue"
  }

  temp_tit = paste("Pixel ",as.integer(cur_pix)," (x,y) (",pix_x,",",pix_y,")",sep="")
  ## output the nbar points color coded by the qa values across our data range
  plot( dates, orig_dat, col = qa_col,pch=snow_pch,cex=0.7,main=temp_tit,
        ylim=c(0.1,0.5),axes=T,ylab="EVI2",xlab="Time")
  lines( dates, evi2_pix)
  
  return()
}  ## end plot pix


##############################
## end function definitions #
##############################

nbands = 365
dsize = 2

ref_file = paste(in_dir,"/",in_stem,".",in_tile,".",in_year,".evi2.bip",sep="")
res_file = paste(in_dir,"/",in_stem,".",in_tile,".",in_year,".evi2.residual.bip",sep="")
flag_file = paste(in_dir,"/",in_stem,".",in_tile,".",in_year,".evi2.flag.bip",sep="")
if(file.exists(ref_file) && file.exists(res_file)) {
  evi2 <- read_int_file(in_file = ref_file, dsize = dsize, nbands = nbands, bsq_flag = 0)
  res <- read_int_file(in_file = res_file, dsize = dsize, nbands = nbands, bsq_flag = 0)
} else {
  print(paste("No input files found ",ref_file,sep=""))
}
flag_on = 1
if(file.exists(flag_file)) {
  flag <- read_int_file(in_file = flag_file, dsize = dsize, nbands = nbands, bsq_flag = 0)
} else {
  flag_on = 0
  flag = NULL
}

out_name = paste("./",in_stem,"_",cur_pix,".jpg",sep="")
jpeg(out_name,height=4,width=6,res=400,units="in")
plot_pix(evi2,res,flag,cur_pix,flag_on)
dev.off()

