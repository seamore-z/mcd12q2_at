library(stringr)

make_plot <- function(in_name,out_name) {

## now we read in the file with the values - premade fun to read table with sep=","
my.data <- read.csv( file = in_name, 
                     colClasses = "numeric",
                     header=FALSE
)
# each row is a pixel - first 5 columns go pixid, y, x, pixnum, date (0-365*num_years)
## next 5 go orig_evi2, orig_ndsi, orig_qa, snow, out_evi2

# separate headers from values:
pix.headers <- my.data[ ,c(1:5)]
pix.data <- my.data[ , c(6:10) ] 

nyears = 5
ndays = 365*nyears

## npix are the number of pixels in our sample
npix = length(pix.headers[,1])/ndays

## setup the output arrays that will be used to plot the points - each will have npix rows and 365*nyears cols
evi2 = array(NA,dim=c(npix,ndays))
ndsi = array(NA,dim=c(npix,ndays))
qa = array(NA,dim=c(npix,ndays))
snow = array(NA,dim=c(npix,ndays))
smooth.data = array(NA,dim=c(npix,ndays))
    
y <- rep(1:nyears,each=365)
d<-1:365
  
## we make colors for the qa values
qa_col <- vector (length = ndays)
snow_pch <- vector (length = ndays)
## make a dates array for the x-axis of our plots
dates = seq(1,(ndays))
## setup the x axis as one label per year at the midpoint of each year
years <- seq(2004,2008)
date_int <- seq(182,ndays,365)

## these are graphing parameters that can be edited
## there is a half year buffer so the buffer at the start and end will be removed
temp_xlim = c(date_int[1],date_int[nyears])
## will be the ylimit - could be up to (0,1)
temp_ylim = c(0,1)

## data range contains all the data we will plot
data_range = seq(temp_xlim[1],temp_xlim[2])

## fill missing data
pix.data[pix.data[,1]==32767,1] = NA
pix.data[pix.data[,2]==32767,2] = NA

## now we setup all inputs and we open the visual device
## in this case a pdf creation function

#open output file
pdf(out_name)

#plot graph for each pixel and band
for(i in 1:npix){
  pix_start = (i-1)*ndays +1
  pix_end = i*ndays 

  smooth.data <- pix.data[pix_start:pix_end,5]/10000
  evi2 <- pix.data[pix_start:pix_end ,1]/10000
  ndsi <- pix.data[pix_start:pix_end ,2]/10000
  qa <- pix.data[pix_start:pix_end,3]
  snow <- pix.data[pix_start:pix_end,4]
  qa_col[qa==0] = "blue"
  qa_col[qa==1] = "cyan"
  qa_col[qa==2] = "forestgreen"
  qa_col[qa==3] = "orange"
  qa_col[qa==4] = "red"
  qa_col[qa>4] = NA
  
  snow_pch[snow==0] = 19
  snow_pch[snow==1] = 1
  snow_pch[snow==2] = 2
  snow_pch[snow > 2] = NA
  
  temp_tit = paste("PixelID: ", pix.headers[pix_start,1],
                   "Row: ", pix.headers[pix_start,2], 
                   "Col: ", pix.headers[pix_start,3])
  
  ## output the nbar points color coded by the qa values across our data range
  plot( dates[data_range], 
        evi2[data_range], 
        col = qa_col[data_range],
        pch=snow_pch[data_range],cex=0.6,
        main=temp_tit,
        ylim=temp_ylim,
        xlim=temp_xlim,
        axes=F,
        ylab="Reflectance",
        xlab="Time")
  
  axis(1,at=date_int,labels=years,pos=temp_ylim)
  ## adds the smoothed splines as a line over the NBAR points
  lines(dates[data_range], 
        smooth.data[data_range], 
        lwd=1.5, col="gray40")
  ## increment tick marks on y axis by 0.05
  axis(2,at=seq(temp_ylim[1],temp_ylim[2],0.1),labels=seq(temp_ylim[1],temp_ylim[2],0.1),pos=temp_xlim[1])
  points(dates[data_range],ndsi[data_range],col="red",pch=17,cex=0.6)
  #  abline(v=365)
  
} ## end i loop
## turn off the pdf device
dev.off()
}  ## end function

make_plot2 <- function(c6_in,c5_in,out_name) {

## now we read in the file with the values - premade fun to read table with sep=","
c6_dat <- read.csv( file = c6_in, 
                     colClasses = "numeric",
                     header=FALSE)

## now we read in the file with the values - premade fun to read table with sep=","
c5_dat <- read.csv( file = c5_in, 
                     colClasses = "numeric",
                     header=FALSE)
# all headers for pixels have ID > 0 in the first column 
# and all rows with values start with 0

# separate headers from values:
head_c6 <- c6_dat[ c6_dat[,1] > 0, (1:5) ]
pix_c6 <- c6_dat[ c6_dat[,1] == 0, (2:5) ] 

head_c5 <- c5_dat[ c5_dat[,1] > 0, (1:5) ]
pix_c5 <- c5_dat[ c5_dat[,1] == 0, (2:5) ] 

nbands = 7
## npix are the number of pixels - rows
npix = length(head_c6[,1])/nbands

if(npix != length(head_c5[,1])/nbands) {
	print("Files do not have same number of lines! Abort!")
	return(NULL)
}

nyears = 3
## the band that corresponds to MODIS BAND 1
## can look this up if we are looping through multiple bands
#band = switch(pix.band + 1, "red","green","blue")

## now we break up each row/pixel by all the features
## first two columns are pix_num and a dummy value (0)
## then order of nbar in, qa in, output for each (365) day for 14 years
## setup the output arrays that will be used to plot the points - each will have npix rows and 365*nyears cols

ndays = 365*nyears
    
y <- rep(1:nyears,each=365)
d<-1:365
  
## we make colors for the qa values
qa_col_c6 <- vector (length = ndays)
snow_pch_c6 <- vector (length = ndays)
qa_col_c5 <- vector (length = ndays)
snow_pch_c5 <- vector (length = ndays)
## make a dates array for the x-axis of our plots
dates = seq(1,(ndays))
## setup the x axis as one label per year at the midpoint of each year
years <- seq(2002,2005)
date_int <- seq(182,ndays,365)

## these are graphing parameters that can be edited
## there is a half year buffer so the buffer at the start and end will be removed
temp_xlim = c(date_int[1],date_int[nyears])
temp_xlim = c(1,ndays)
## will be the ylimit - could be up to (0,1)
temp_ylim = array(0,dim=c(nbands,2))
temp_ylim[1,]=c(0,0.2)
temp_ylim[2,]=c(0,0.5)
temp_ylim[3,]=c(0,0.2)
temp_ylim[4,]=c(0,0.2)
temp_ylim[5,]=c(0,0.5)
temp_ylim[6,]=c(0,0.5)
temp_ylim[7,]=c(0,0.3)

## data range contains all the data we will plot
data_range = seq(temp_xlim[1],temp_xlim[2])
## now we setup all inputs and we open the visual device
## in this case a pdf creation function

#open output file
pdf(out_name,width=9,height=9)

par(mar=c(4,3,2,1)+0.1,mgp=c(2,1,0),cex.axis=0.8,cex.lab=0.8,cex.main=0.9,las=0,mfrow=c(2,2))
#plot graph for each pixel and band
for(i in 1:npix){
for(b in 1:nbands) {

  pix_start = ndays*npix*(b-1) + (i-1)*ndays +1
  pix_end = ndays*npix*(b-1) + i*ndays 
  smooth_c5 <- pix_c5[pix_start:pix_end,4]/10000
  nbar_c5 <- pix_c5[pix_start:pix_end ,1]/10000
  qa_c5 <- pix_c5[pix_start:pix_end,2]
  snow_c5 <- pix_c5[pix_start:pix_end,3]

  smooth_c6 <- pix_c6[pix_start:pix_end,4]/10000
  nbar_c6 <- pix_c6[pix_start:pix_end ,1]/10000
  qa_c6 <- pix_c6[pix_start:pix_end,2]
  snow_c6 <- pix_c6[pix_start:pix_end,3]

  qa_col_c5[qa_c5==0] = "blue"
  qa_col_c5[qa_c5==1] = "cyan"
  qa_col_c5[qa_c5==2] = "forestgreen"
  qa_col_c5[qa_c5==3] = "orange"
  qa_col_c5[qa_c5==4] = "red"
  qa_col_c5[qa_c5>4] = NA
  
  snow_pch_c5[snow_c5==0] = 19
  snow_pch_c5[snow_c5==1] = 1
  snow_pch_c5[snow_c5 > 1] = NA

  qa_col_c6[qa_c6==0] = "blue"
  qa_col_c6[qa_c6==1] = "cyan"
  qa_col_c6[qa_c6==2] = "forestgreen"
  qa_col_c6[qa_c6==3] = "orange"
  qa_col_c6[qa_c6==4] = "red"
  qa_col_c6[qa_c6>4] = NA
  
  snow_pch_c6[snow_c6==0] = 15
  snow_pch_c6[snow_c6==1] = 0
  snow_pch_c6[snow_c6 > 1] = NA

  if(b == 1) {
	smooth_r_c5 = smooth_c5
	nbar_r_c5 = nbar_c5
	qa_r_c5 = qa_c5

	smooth_r_c6 = smooth_c6
	nbar_r_c6 = nbar_c6
	qa_r_c6 = qa_c6
} else if(b==2) {
	smooth_n_c5 = smooth_c5
	nbar_n_c5 = nbar_c5
	qa_n_c5 = qa_c5

	smooth_n_c6 = smooth_c6
	nbar_n_c6 = nbar_c6
	qa_n_c6 = qa_c6
}
  
  head_loc = ((b-1)*npix) + i
  temp_tit = paste("PixelID: ", head_c5[head_loc,1],
                   "Row: ", head_c5[head_loc,2], 
                   "Col: ", head_c5[head_loc,3], 
                   "Band: ",head_c5[head_loc,4]+1)
  
  ## output the nbar points color coded by the qa values across our data range
  plot( dates[data_range], 
        nbar_c5[data_range], 
        col = qa_col_c5[data_range],
        pch=snow_pch_c5[data_range],cex=0.5,
        main=temp_tit,
        ylim=temp_ylim[b,],
        xlim=temp_xlim,
        axes=F,
        ylab="Reflectance",
        xlab="Time")
  points( dates[data_range], 
        nbar_c6[data_range], 
        col = qa_col_c6[data_range],
        pch=snow_pch_c6[data_range],cex=0.5)
  axis(1,at=date_int,labels=years,pos=temp_ylim[b,1])
  ## adds the smoothed splines as a line over the NBAR points
  lines(dates[data_range], 
        smooth_c5[data_range], 
        lwd=1.5, col="gray20")
  lines(dates[data_range], 
        smooth_c6[data_range], 
        lwd=1.5, col="gray50",lty=2)
  ## increment tick marks on y axis by 0.05
  axis(2,at=seq(temp_ylim[b,1],temp_ylim[b,2],0.05),labels=seq(temp_ylim[b,1],temp_ylim[b,2],0.05),pos=temp_xlim[1])
  #  abline(v=365)
  if(b==1 || b==5) {
  legend(x=temp_xlim[1]+100,y=temp_ylim[b,2],legend=c("no snow c5","snow c5","no snow c6","snow c6","c5_fit","c6_fit"),col=c("black","black","black","black","gray20","gray50"),pch=c(19,1,15,0,NA,NA),lty=c(NA,NA,NA,NA,1,2),cex=0.5)
  legend(x=temp_xlim[1]+400,y=temp_ylim[b,2],legend=c("qa = 0","qa = 1","qa = 2","qa =3"),fill=c("blue","cyan","forestgreen","orange"),cex=0.5)
}
}  ## end bands
temp_tit = paste("PixelID: ", head_c5[head_loc,1],
                   "Row: ", head_c5[head_loc,2], 
                   "Col: ", head_c5[head_loc,3], 
                   "EVI2")

   comb_qa_c5 = qa_r_c5 < 4 & qa_n_c5 < 4
   comb_qa_c6 = qa_r_c6 < 4 & qa_n_c6 < 4
   temp_dat = cbind(smooth_r_c5,smooth_n_c5)
  smooth_evi_c5 = calc_evi2(temp_dat)
   temp_dat = cbind(smooth_r_c6,smooth_n_c6)
  smooth_evi_c6 = calc_evi2(temp_dat)
temp_dat = cbind(nbar_r_c5,nbar_n_c5)
  evi_c5 = calc_evi2(temp_dat)
temp_dat = cbind(nbar_r_c6,nbar_n_c6)
  evi_c6 = calc_evi2(temp_dat)

  qa_col_c5[comb_qa_c5] = "red"
  qa_col_c5[!comb_qa_c5] = NA

  qa_col_c6[comb_qa_c6] = "blue"
  qa_col_c6[!comb_qa_c6] = NA

  ## output the EVI points color coded by the qa values across our data range
  plot( dates[data_range], 
        evi_c5[data_range], 
        col = qa_col_c5[data_range],
        pch=19,cex=0.5,
        main=temp_tit,
        ylim=c(0,1),
        xlim=temp_xlim,
        axes=F,
        ylab="EVI2",
        xlab="Time")
  points( dates[data_range], 
        evi_c6[data_range], 
        col = qa_col_c6[data_range],
        pch=19,cex=0.5)
  axis(1,at=date_int,labels=years,pos=0)
  ## adds the smoothed splines as a line over the NBAR points
  lines(dates[data_range], 
        smooth_evi_c5[data_range], 
        lwd=1.5, col="gray20")
  lines(dates[data_range], 
        smooth_evi_c6[data_range], 
        lwd=1.5, col="gray50",lty=2)
  ## increment tick marks on y axis by 0.05
  axis(2,at=seq(0,1,0.1),labels=seq(0,1,0.1),pos=temp_xlim[1])
  legend(x=temp_xlim[1],y=1,legend=c("c5 points","c6 points","c5_fit","c6_fit"),col=c("red","blue","gray20","gray50"),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),cex=0.5)
## for a blank plot - placeholder
 # plot( dates[data_range],nbar_c5[data_range], type="n",ylim=temp_ylim[b,],xlim=temp_xlim, axes=F, ylab=NA, xlab=NA)
}  ## end pix
## turn off the pdf device
dev.off()

return(1)
}

## expects 2 bands (red,nir) and returns 1
calc_evi2 <- function(in_dat) {
	npix = length(in_dat[,1])
	temp_ind = in_dat[,1] < 1 & in_dat[,2] < 1
	out = array(NA,npix)
	out[temp_ind] = 2.5*((in_dat[temp_ind,2]-in_dat[temp_ind,1]) / (in_dat[temp_ind,2] + (2.4*in_dat[temp_ind,1]) + 1))

	return(out)

}


## setup input/output files
in_tile = "h19v08"
in_dir = "../inputs"
in_name = paste(in_tile,"_pixels.out",sep="")
out_dir = "."
out_name = paste(in_tile,"evi2_c5.pdf",sep="_")

full_in_name = paste(in_dir,in_name,sep="/")
full_out_name = paste(out_dir,out_name,sep="/")

make_plot(full_in_name,full_out_name)




