rm(list=ls())

# read in the arguments
args = commandArgs(trailingOnly=T)

debug = 0
if(debug==0) {
  in_dir = args[1]
  in_stem = args[2]
  tile = args[3]
  start_year = as.numeric(args[4])
  num_years = as.numeric(args[5])
} else {
  in_dir = "../MCD12I5"
  in_stem = "MCD12I5"
  tile= "h12v04"
  ## should be actual year 
  start_year = 2001
  num_years = 5
}


tile_h = as.numeric(substring(tile,2,3))
tile_v = as.numeric(substring(tile,5,6))

#in_files = vector("list",length(in_stem))

in_files=array(NA,num_years)
for(y in 1:num_years) {
    search_str = paste(in_stem,".A",(start_year+y-1),"001.",tile,"*.txt",sep="")
  
    temp_list = list.files(path=in_dir,pattern=glob2rx(search_str),full.names=T,include.dirs=F)
    if(length(temp_list==1)) {
        in_files[y] = temp_list
    }
 #     print(in_files[[s]][y])
}



in_file = paste("./step_out.csv",sep="")
if(file.exists(in_file)) {
  step_dat = read.table(in_file,sep=",",na.strings="?")
}

tile_ind = step_dat[,1] == tile_h & step_dat[,2] == tile_v
tile_dat = step_dat[tile_ind,]


## now we reformat the arrays so they contain a continuous time series
## break up the arrays into band_dat and flags and pre and post splining
## each array has 3 columns to be ignored then we want to extract out the middle year of 3
## then the data are organized in row chunks of processing 
#  with each pixel containing 8 bands including bands(1,2,4,5,6,7) and EVI2 + snow_flag
## there is one set of that pixel pre-processing and one post splining

## first we select out the columns we want
start_col = 3 + 365 + 1
end_col = 3 + 365*2
nbands = 8

all_dat = vector("list",num_years)
for(y in 1:num_years) {
    if(!is.na(in_files[y])) {
      all_dat[[y]] = read.table(in_files[y],sep=",")
    }
    else {
	    print(paste("File ",in_files[y]," does not exist.",sep=""))
      quit(save="no")
    }
}
## the first 3 are site id info
all_ids = all_dat[[1]][,1]
un_ids = unique(all_ids)
band_nums = all_dat[[1]][,3]
npix = length(all_ids)
num_ids = length(un_ids)
ncols = dim(all_dat[[1]])[[2]]

## now we make a set of arrays to hold the outputs
tot_out_dates = 365*num_years
all_post = vector("list",num_ids)
all_pre = vector("list",num_ids)
pre_out = array(NA,dim=c(nbands,tot_out_dates))
post_out = array(NA,dim=c(nbands,tot_out_dates))
for(i in 1:num_ids) {
  temp_ind = all_ids == un_ids[i]

  for(y in 1:num_years) {
    temp_dat = all_dat[[y]][temp_ind,start_col:end_col]
    y_start = (y-1)*365 + 1
    y_end = y*365
    for(b in 1:nbands) {
      cur_dat = as.numeric(temp_dat[b,])
      na_ind = cur_dat==32767
      cur_dat[na_ind] = NA
      pre_out[b,y_start:y_end] = cur_dat
    }
    for(b in (nbands+1):(nbands*2)) {
      cur_dat = as.numeric(temp_dat[b,])
      na_ind = cur_dat==32767
      cur_dat[na_ind] = NA
      post_out[(b-nbands),y_start:y_end] = cur_dat
    }
  }
  all_pre[[i]] = pre_out
  all_post[[i]] = post_out
}

y_arr = as.integer(un_ids/2400)
x_arr = un_ids-(y_arr*2400)

b_names = c("Red","NIR","Green","SWIR1","SWIR2","SWIR3","EVI2","flag")
mycols <- c("#636363", "#31A354","#3182BD", "#DE2D26","#756BB1")

date_array = rep(seq(1,365),num_years)
year_array = date_array
for(y in 1:num_years) {
  start = ((y-1)*365)+1
  end = y*365
  year_array[start:end] = start_year+y
}
all_dates = paste(year_array,date_array,sep="-")
flag_labs = c("Spline-smoothed", "Non-Snow Observed", "Non-Snow Interpolated", "Non-Snow Min Filled", "Snow Observed", "Snow Interpolated")

plot_lims = array(NA,dim=c((nbands-1),2))
for(b in 1:(nbands-1)) {
  plot_lims[b,1] = 0
  plot_lims[b,2] = 1
}
plot_lims[nbands-1,1] = 0.1
plot_lims[nbands-1,2]=0.8

date_out_int = seq(75,tot_out_dates,75)
all_dates_out = all_dates[date_out_int] 
num_dates = length(all_dates)
plot_legend=1

full_plot = 0
mean_plot = 1
if(full_plot == 1) {
  ## i loop goes here
  out_f = paste("./",tile,".pdf",sep="")
  pdf(out_f,width=12,height=10)
  par(mar=c(4,3,2,1)+0.1,mgp=c(2,1,0),cex.axis=0.8,cex.lab=0.8,cex.main=0.9,las=0,mfrow=c(4,2))
  #par(mar=c(4,3,2,1)+0.1,mgp=c(2,1,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,las=0)
  for(i in 1:num_ids) {
  
    snow_flags = all_post[[i]][8,]
    cols <- mycols[snow_flags + 1]
    pchs <- snow_flags + 1
    
   # b = 7
    for(b in 1:(nbands-1)) {
      main_tit = paste("Tile, ",tile,", band ",b_names[b],", pix (",x_arr[i],",",y_arr[i],"), site id ",tile_dat[i,6],", igbp = ",tile_dat[i,16],sep="")
      plot(all_post[[i]][b,]/10000,type="l",main=main_tit,xlab="Date",ylab="Reflectance",ylim=plot_lims[b,],xlim=c(1,tot_out_dates),axes=F)
      points(all_pre[[i]][b,]/10000,cex=0.5,col=cols,pch=pchs)
      axis(1,at=date_out_int,labels=all_dates_out,)
      axis(2,at=seq(plot_lims[b,1],plot_lims[b,2],0.1),labels=seq(plot_lims[b,1],plot_lims[b,2],0.1))
      box()
  
      if(plot_legend) {
        legend("topleft", legend=flag_labs, pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white",cex=0.8)
      }
   }
   ## makes an empty plot to fill the last spot
   plot.new()
   legend("top", legend=flag_labs, pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white",cex=0.7)
  
  }
  dev.off()
}

site_ids = unique(tile_dat[,6])
num_site_ids = length(site_ids)
## do it differently to get one set of plots for whole site
if(mean_plot == 1) {
  ## i loop goes here
  out_f = paste("./mean.",tile,".",in_stem,".pdf",sep="")
  pdf(out_f,width=12,height=10)
  par(mar=c(4,3,2,1)+0.1,mgp=c(2,1,0),cex.axis=0.8,cex.lab=0.8,cex.main=0.9,las=0,mfrow=c(4,2))
  #par(mar=c(4,3,2,1)+0.1,mgp=c(2,1,0),cex.axis=0.7,cex.lab=0.7,cex.main=0.7,las=0)
  for(i in 1:num_site_ids) {
    
    temp_ind = tile_dat[,6]==site_ids[i]
    snow_flags = array(NA,num_dates)
    for(n in 1:num_dates) {
      for(j in 1:length(temp_ind)) {
        if(temp_ind[j]) {
          snow_flags[n]=max(snow_flags[n],all_post[[j]][8,n],na.rm=T)
        }
      }
    }
    cols <- mycols[snow_flags + 1]
    pchs <- snow_flags + 1
    cur_x = floor(mean(x_arr[temp_ind]))
    cur_y = floor(mean(y_arr[temp_ind]))
    cur_igbp = mean(tile_dat[temp_ind,16])
    # b = 7
    for(b in 1:(nbands-1)) {
      
      ## need to massage to get in right format
      count_pre = array(0,num_dates)
      sum_pre = array(0,num_dates)
      count_post = array(0,num_dates)
      sum_post = array(0,num_dates)
      for(n in 1:num_dates) {
        for(j in 1:length(temp_ind)) {
          if(temp_ind[j]) {
            if(!is.na(all_pre[[j]][b,n])) {
              sum_pre[n]= sum_pre[n] + all_pre[[j]][b,n]
              count_pre[n]= count_pre[n] + 1
            }
            if(!is.na(all_post[[j]][b,n])) {
              sum_post[n]= sum_post[n] + all_post[[j]][b,n]
              count_post[n]= count_post[n] + 1
            }
          }  
        }
      }
      pre_dat = sum_pre/count_pre
      post_dat = sum_post/count_post
      
      main_tit = paste("Tile, ",tile,", band ",b_names[b],", pix (",cur_x,",",cur_y,"), site id ",site_ids[i],", igbp = ",cur_igbp,sep="")
      plot(post_dat/10000,type="l",main=main_tit,xlab="Date",ylab="Reflectance",ylim=plot_lims[b,],xlim=c(1,tot_out_dates),axes=F)
      points(pre_dat/10000,cex=0.5,col=cols,pch=pchs)
      axis(1,at=date_out_int,labels=all_dates_out,)
      axis(2,at=seq(plot_lims[b,1],plot_lims[b,2],0.1),labels=seq(plot_lims[b,1],plot_lims[b,2],0.1))
      box()
      

    }
    ## makes an empty plot to fill the last spot
    plot.new()
     if(plot_legend) {
        legend("top", legend=flag_labs, pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white",cex=1)
     }
    
    
  }
  dev.off()
}

