library(raster)
library(sp)

#NEED TO CHANGE YEAR
ncol <- 2400
nrow <- 2400
nbands <- 365
cnt <- ncol*nrow*nbands

data <- readBin('/projectnb/modislc/projects/sat/data/spline/h10v03/tg_v1.h10v03.2019.evi2.bip',
                what="integer",n=cnt,size=2,endian="little")
data <- array(data,c(nbands, ncol, nrow))
data <- aperm(data, c(3,2,1)) #for transposing
data <- brick(data)


saveRDS(data, file='/projectnb/modislc/projects/sat/data/spline/h10v03/splineEVIRData/data2019.rds')
#save(data, file='/projectnb/modislc/projects/sat/data/spline/h10v03/data2019.RData')


path <- paste('/projectnb/modislc/projects/sat/data/spline')
search_str <- paste('*.Rds',year,tile,sep='') #Search for files with specific year and tile
files <- list.files(path=path, pattern=glob2rx(search_str),full.names=T,include.dirs=F,recursive=T) #List those files
