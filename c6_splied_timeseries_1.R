library(sp)
library(raster)
#library(gdalUtils)
library(terra)
library(rgdal)
library(doMC)
library(doParallel)


args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- c(as.numeric(substr(args[3],7,10)))
xx <- as.numeric(substr(args[3],11,11))


#tile = 'h10v02'
#year = c(2017, 2018, 2019, 2020, 2021)
mcd12a4_path <- paste('/projectnb/modislc/projects/sat/data/e4ftl01.cr.usgs.gov/MOTA/MCD43A4.061')
search_str <- paste('*43A4.A',year[1],'*.',tile,'.061*',sep='')
files <- list.files(path=mcd12a4_path, pattern=glob2rx(search_str),full.names=T,include.dirs=F,recursive=T)

#test <- raster(get_subdatasets(files[1])[[1]])
rast_obj <- rast(files[1])
subdatasets <- names(rast_obj)
test <- raster(rast_obj[subdatasets[1]])

# x <- -163.25640; y_orig <- 61.30227 64.919067,-166.376788
x <- -166.376788; y_orig <- 64.919067
# x <- -162.655197; y_orig <- 60.923207
# x <- -166.375628; y_orig <- 64.981045
# x <- -166.358837; y_orig <- 64.893833
# x <- -165.414587; y_orig <- 64.743895

points <- cbind(x,y_orig)
v <- SpatialPoints(points, proj4string = CRS("+proj=longlat +datum=WGS84"))
y <- spTransform(v, CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"))

numPix <- length(test)
imgNum <- setValues(test, 1:numPix)
z <- extract(imgNum, y)
sam <- z #make lat/lon coord sample((2400*2400),30)


## Extract splined EVI 
ncol <- 2400
nrow <- 2400
nbands <- 365
cnt <- ncol*nrow*nbands

data <- readBin(paste('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/',tile,'/sz_v5.',tile,'.',as.character(year),'.evi2.bip',sep=''),
                what="integer",n=cnt,size=2,endian="little")
data <- array(data,c(nbands, ncol, nrow))
data <- aperm(data, c(3,2,1)) #for transposing
data <- brick(data)

splined <- matrix(NA,365,1)
for(i in 1:365){
  rast <- data[[i]]
  
  splined[i,] <- rast[sam]
}

registerDoMC()
evi_mm <- matrix(NA,365,1)
ndsi_mm <- matrix(NA,365,1)
for(i in 1:365){
# evi_mm <- foreach(i=1:365,.combine=rbind) %dopar% {
  
  #sds <- get_subdatasets(files[i])
  sds <- rast(files[i])
  # sds <- names(rast_obj)
  # sds <- raster(rast_obj[subdatasets[i]])
  
  b1 <- raster(sds[[8]])[sam]
  b2 <- raster(sds[[9]])[sam]
  b4 <- raster(sds[[11]])[sam]
  b6 <- raster(sds[[13]])[sam]
  
  evi_mm[i,1] <- 2.5*((b2-b1)/(b2 + (2.4*b1) +1)) # calculate EVI2 
  ndsi_mm[i,1] <- (b4-b6)/(b4+b6)                 # calculate NDSI 

  # print(i)
}


setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/spline_test/')
png(file=paste('sz_v5_',as.character(year),'_',as.character(x),'_',as.character(y_orig),'.png',sep=''),res=300,units='in',width=13.5,height=6.5)
par(oma=c(1,1,1,1),mar=c(4,4,0,0))
plot(splined[,1]/10000,ylim=c(-0.2,0.8),
     xlab='Date',ylab='EVI2',lwd=1.5,
     cex.lab=1.3,type='l')
points(evi_mm,col='red',cex=.5)
dev.off()

setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/spline_test/')
png(file=paste('sz_v5_',as.character(year),'_',as.character(x),'_',as.character(y_orig),'_ndsi','.png',sep=''),res=300,units='in',width=13.5,height=6.5)
par(oma=c(1,1,1,1),mar=c(4,4,0,0))
plot(splined[,1]/10000,ylim=c(-0.8,1.1),
     xlab='Date',ylab='NDSI',lwd=1.5,
     cex.lab=1.3,type='l')
points(ndsi_mm,col='blue',cex=.5)
dev.off()


## Testing spline smoothing param spar ##
## Note: no weights on this spline     ##
evi_dd <- matrix(NA,365,1)
for(i in 1:365){
  if (!is.na(ndsi_mm[i, 1]) && ndsi_mm[i,1] < -0.3) {
    evi_dd[i,1] <- evi_mm[i,1]
  } else {
    #if (!is.na(evi_mm[i,1]) && evi_mm[i,1] < 0.3) {
      #evi_dd[i,1] <- 0.15
    #}
  }
}
evi_dd_qa <- evi_dd[evi_dd!=0.1500000]; evi_dd_qa <- evi_dd_qa[!is.na(evi_dd_qa)]
evi_min <- quantile(evi_dd_qa, probs = 0.10)
if (evi_min>0.15) {
  evi_dd[evi_dd==0.15] <- evi_min
}
  
dates <- c(1:365)
dates[is.na(evi_dd)] <- NA
dates <- na.omit(dates)
evi_dd <- na.omit(evi_dd)

spl <- smooth.spline(dates, evi_dd, spar=.55)

pred_dates <- seq(1,365,length.out = 365)
xSmooth <- predict(spl, as.numeric(pred_dates))$y

setwd('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/MCD43A4/data/spline/spline_test/')
png(file=paste('sz_v5_',as.character(year),'_',as.character(x),'_',as.character(y_orig),'_tspine','.png',sep=''),res=300,units='in',width=13.5,height=6.5)
par(oma=c(1,1,1,1),mar=c(4,4,0,0))
plot(pred_dates, xSmooth,ylim=c(0,0.8),
     xlab='Date',ylab='EVI2',lwd=1.5,
     cex.lab=1.3,type='l')
points(dates, evi_dd,col='green',cex=.5)
dev.off()
