library(sp)
library(raster)
library(gdalUtils)
library(rgdal)

args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- as.numeric(substr(args[3],7,10))
# tile <- 'h11v09'
# year <- 2019

sam <- sample((2400*2400),30)
num <- 1

## Extract splined EVI 
ncol <- 2400
nrow <- 2400
nbands <- 365
cnt <- ncol*nrow*nbands


splined <- matrix(NA,(365*3),2)

# BU
for(i in 1:(365*3)){
  if(i==1){
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/spline/c61_v7.',tile,'.',(year-1),'.evi2.bip',sep=''),
                    what="integer",n=cnt,size=2,endian="little")
    # data <- readBin(files[1],what="integer",n=cnt,size=2,endian="little")
    data <- array(data,c(nbands, ncol, nrow))
    data <- aperm(data, c(3,2,1)) #for transposing
    data <- brick(data)
    rast <- data[[i]]
  }else if(i>1 & i<(365*1+1)){
    rast <- data[[i]]
  }else if(i==(365*1+1)){
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/spline/c61_v7.',tile,'.',(year-0),'.evi2.bip',sep=''),
                    what="integer",n=cnt,size=2,endian="little")
    # data <- readBin(files[2],what="integer",n=cnt,size=2,endian="little")
    data <- array(data,c(nbands, ncol, nrow))
    data <- aperm(data, c(3,2,1)) #for transposing
    data <- brick(data)
    rast <- data[[i-(365*1)]]
  }else if(i>(365*1+1) & i<(365*2+1)){
    rast <- data[[i-(365*1)]]
  }else if(i==(365*2+1)){
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/spline/c61_v7.',tile,'.',(year+1),'.evi2.bip',sep=''),
                    what="integer",n=cnt,size=2,endian="little")
    # data <- readBin(files[3],what="integer",n=cnt,size=2,endian="little")
    data <- array(data,c(nbands, ncol, nrow))
    data <- aperm(data, c(3,2,1)) #for transposing
    data <- brick(data)
    rast <- data[[i-(365*2)]]
  }else{
    rast <- data[[i-(365*2)]]
  }

  splined[i,1] <- as.matrix(rast)[sam[num]]
  if(i%%100==0) print(i)
}

rm(data)

# BU_1
for(i in 1:(365*3)){
  if(i==1){
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v7/spline/bu/c61_v7.',tile,'.',(year-1),'.evi2.bip',sep=''),
                    what="integer",n=cnt,size=2,endian="little")
    # data <- readBin(files[1],what="integer",n=cnt,size=2,endian="little")
    data <- array(data,c(nbands, ncol, nrow))
    data <- aperm(data, c(3,2,1)) #for transposing
    data <- brick(data)
    rast <- data[[i]]
  }else if(i>1 & i<(365*1+1)){
    rast <- data[[i]]
  }else if(i==(365*1+1)){
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v7/spline/bu/c61_v7.',tile,'.',(year-0),'.evi2.bip',sep=''),
                    what="integer",n=cnt,size=2,endian="little")
    # data <- readBin(files[2],what="integer",n=cnt,size=2,endian="little")
    data <- array(data,c(nbands, ncol, nrow))
    data <- aperm(data, c(3,2,1)) #for transposing
    data <- brick(data)
    rast <- data[[i-(365*1)]]
  }else if(i>(365*1+1) & i<(365*2+1)){
    rast <- data[[i-(365*1)]]
  }else if(i==(365*2+1)){
    data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/v7/spline/bu/c61_v7.',tile,'.',(year+1),'.evi2.bip',sep=''),
                    what="integer",n=cnt,size=2,endian="little")
    # data <- readBin(files[3],what="integer",n=cnt,size=2,endian="little")
    data <- array(data,c(nbands, ncol, nrow))
    data <- aperm(data, c(3,2,1)) #for transposing
    data <- brick(data)
    rast <- data[[i-(365*2)]]
  }else{
    rast <- data[[i-(365*2)]]
  }
  
  splined[i,2] <- as.matrix(rast)[sam[num]]
  if(i%%100==0) print(i)
}

rm(data)

##########################
# Phenology
mcd12q1_path <- paste('/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01',sep='')
file <- list.files(path=mcd12q1_path,pattern=glob2rx(paste('*',tile,'*',sep='')),full.names=T)
sds <- get_subdatasets(file)
lct <- raster(sds[1])


phemat<- matrix(NA,7,2)

phe <- c('gup','mid-gup','mat','peak','gdn','mid-gdn','dor')
phe_na <- c(2,3,5,4,6:11)
phe_bu <- c(5:11,4,3,2)

doy_offset <- as.integer(as.Date(paste((year-1),'-12-31',sep='')) - as.Date("1970-1-1"))

for(vv in c(1,3,5,7)){
  ## 
  path <- '/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/pheno'
  files <- list.files(path=path,pattern=glob2rx(paste('*',tile,'*2019',sep='')),full.names=T)
  rast <- raster(files,band=phe_bu[vv])
  values(rast)[values(rast)>32700] <- NA
  values(rast)[values(lct)==17] <- NA
  rast_na <- rast - doy_offset
  
  ## 
  path <- '/projectnb/modislc/users/mkmoon/LCD_C6/v7/pheno/bu'
  files <- list.files(path=path,pattern=glob2rx(paste('*',tile,'*2019',sep='')),full.names=T)
  rast <- raster(files,band=phe_bu[vv])
  values(rast)[values(rast)>32700] <- NA
  values(rast)[values(lct)==17] <- NA
  rast_bu <- rast - doy_offset
  
  phemat[vv,1] <- as.matrix(rast_na)[sam[num]]
  phemat[vv,2] <- as.matrix(rast_bu)[sam[num]]
 
  print(vv) 
}




# # NASA
# path <- '/projectnb/modislc/users/mkmoon/LCD_C6/v7/spline/nasa'
# for(i in 1:(365*3)){
#   if(i==1){
#     sstr <- paste('*A',(year-1),'001.',tile,'*.bip',sep='')
#     files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
#     data <- readBin(files,what="integer",n=cnt,size=2,endian="little")
#     # data <- readBin(files[1],what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i]]
#   }else if(i>1 & i<(365*1+1)){
#     rast <- data[[i]]
#   }else if(i==(365*1+1)){
#     sstr <- paste('*A',(year-0),'001.',tile,'*.bip',sep='')
#     files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
#     data <- readBin(files,what="integer",n=cnt,size=2,endian="little")
#     # data <- readBin(files[2],what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*1)]]
#   }else if(i>(365*1+1) & i<(365*2+1)){
#     rast <- data[[i-(365*1)]]
#   }else if(i==(365*2+1)){
#     sstr <- paste('*A',(year+1),'001.',tile,'*.bip',sep='')
#     files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
#     data <- readBin(files,what="integer",n=cnt,size=2,endian="little")
#     # data <- readBin(files[3],what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*2)]]
#   }else{
#     rast <- data[[i-(365*2)]]
#   }
#   
#   splined[i,2] <- as.matrix(rast)[sam[num]]
#   if(i%%10==0) print(i)
# }



sstr <- sprintf("%010d",sam[num])
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/spline/ts/')
png(file=paste('c6_v7_',tile,'_',year,'_',sstr,'.png',sep=''),
    res=300,units='in',width=13.5,height=6.5)
par(oma=c(2,2,2,2),mar=c(3,3,1,0))
  plot(splined[,1]/10000,ylim=c(0,0.85),
     type='o',xlab='Date',ylab='EVI2',lwd=1.5,
     cex.lab=1.3,axe=F,ann=F,col='blue')
  points(splined[,2]/10000,col='red',cex=0.5)
axis(1,at=c(1,(365*1+1),(365*2+1),(365*3+1),
            c(as.character(year-1),as.character(year),as.character(year+1)),cex.axis=1.3))
axis(2,at=seq(0,1,0.2),cex.axis=1.3)
mtext('EVI2',2,cex=1.3,line=2.5)
abline(v=c(1,(365*1+1),(365*2+1),(365*3+1)),lty=1)

abline(v=(phemat[c(1,3,5,7),1]+365),lty=2,lwd=2,col='blue')
abline(v=(phemat[c(1,3,5,7),2]+365),lty=3,lwd=2,col='red')


box(lty=1)
dev.off()


# splined <- matrix(NA,(365*10),1)
# for(i in 1:(365*10)){
#   if(i==1){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2004.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i]]
#   }else if(i>1 & i<(365*1+1)){
#     rast <- data[[i]]
#   }else if(i==(365*1+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2005.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*1)]]
#   }else if(i>(365*1+1) & i<(365*2+1)){
#     rast <- data[[i-(365*1)]]
#   }else if(i==(365*2+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2006.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*2)]]
#   }else if(i>(365*2+1) & i<(365*3+1)){
#     rast <- data[[i-(365*2)]]
#   }else if(i==(365*3+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2007.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*3)]]
#   }else if(i>(365*3+1) & i<(365*4+1)){
#     rast <- data[[i-(365*3)]]
#   }else if(i==(365*4+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2008.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*4)]]
#   }else if(i>(365*4+1) & i<(365*5+1)){
#     rast <- data[[i-(365*4)]]
#   }else if(i==(365*5+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2009.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*5)]]
#   }else if(i>(365*5+1) & i<(365*6+1)){
#     rast <- data[[i-(365*5)]]
#   }else if(i==(365*6+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2010.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*6)]]
#   }else if(i>(365*6+1) & i<(365*7+1)){
#     rast <- data[[i-(365*6)]]
#   }else if(i==(365*7+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2011.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*7)]]
#   }else if(i>(365*7+1) & i<(365*8+1)){
#     rast <- data[[i-(365*7)]]
#   }else if(i==(365*8+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2012.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*8)]]
#   }else if(i>(365*8+1) & i<(365*9+1)){
#     rast <- data[[i-(365*8)]]
#   }else if(i==(365*9+1)){
#     data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/c6_03to16.',tile,'.2013.evi2.bip',sep=''),
#                     what="integer",n=cnt,size=2,endian="little")
#     data <- array(data,c(nbands, ncol, nrow))
#     data <- aperm(data, c(3,2,1)) #for transposing
#     data <- brick(data)
#     rast <- data[[i-(365*9)]]
#   }else{
#     rast <- data[[i-(365*9)]]
#   }
# 
#   splined[i,1] <- as.matrix(rast)[sam[num]]
#   if(i%%10==0) print(i)
# }
# 
# sstr <- sprintf("%010d",sam[num])
# setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/splined/')
# png(file=paste('c6_1_',tile,'_',sstr,'.png',sep=''),
#     res=300,units='in',width=13.5,height=6.5)
# par(oma=c(2,2,2,2),mar=c(3,3,1,0))
#   plot(splined[,1]/10000,ylim=c(0,0.85),
#      type='o',xlab='Date',ylab='EVI2',lwd=1.5,
#      cex.lab=1.3,axe=F,ann=F)
# axis(1,at=c(1,(365*1+1),(365*2+1),(365*3+1),(365*4+1),(365*5+1),(365*6+1),(365*7+1),(365*8+1),(365*9+1),(365*10+1)),
#      c('2004','2005','2006','2007','2008','2009','2010','2011','2012'),
#      cex.axis=1.3)
# axis(2,at=seq(0,1,0.2),cex.axis=1.3)
# mtext('EVI2',2,cex=1.3,line=2.5)
# abline(v=c(1,(365*1+1),(365*2+1),(365*3+1),(365*4+1),(365*5+1),(365*6+1),(365*7+1),(365*8+1),(365*9+1),(365*10+1)),lty=5)
# box(lty=1)
# dev.off()


