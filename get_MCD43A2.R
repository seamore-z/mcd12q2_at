# Get NBAR (MCD43A4 C6) data from USGS
args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- as.numeric(substr(args[3],7,10))

# MCD43A2
setwd('/projectnb/modislc/projects/sat/data/')
if(year==2000){
  mod_dates <- seq(1,365)
  for(i in 55:365){
    d = mod_dates[i]
    d_str = sprintf("%03d",d)
    dates <- as.Date(d,format='%Y.%m.%d',origin=paste((year-1),'.12.31',sep=''))
    mm <- substr(dates,6,7)
    dd <- substr(dates,9,10)
    
    url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD43A2.061/',year,'.',mm,'.',dd,'/',sep='')
    
    system(paste('wget --user=seamore --password=stxt105*T -l1 -r --no-parent -A "MCD43A2.A',year,d_str,'.',tile,'.061.*.hdf" ',url,sep=''))
    
    print(i)
  }
}else{
  mod_dates <- seq(1,365)
  for(i in 1:365){
    d = mod_dates[i]
    d_str = sprintf("%03d",d)
    dates <- as.Date(d,format='%Y.%m.%d',origin=paste((year-1),'.12.31',sep=''))
    mm <- substr(dates,6,7)
    dd <- substr(dates,9,10)
    
    url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD43A2.061/',year,'.',mm,'.',dd,'/',sep='')
    
    system(paste('wget --user=seamore --password=stxt105*T -l1 -r --no-parent -A "MCD43A2.A',year,d_str,'.',tile,'.061.*.hdf" ',url,sep=''))
    
    print(i)
  }
}




# #####
# args <- commandArgs()
# print(args)
# 
# chunk <- as.numeric(args[3])
# 
# 
# flist <- as.matrix(read.table('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/download2.txt'))
# if(chunk==1){
#   for(i in 1:500){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==2){
#   for(i in 501:1000){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==3){
#   for(i in 1001:1500){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==4){
#   for(i in 1501:2000){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==5){
#   for(i in 2001:2500){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==6){
#   for(i in 2501:3000){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==7){
#   for(i in 3001:3500){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==8){
#   for(i in 3501:4000){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==9){
#   for(i in 4001:4500){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==10){
#   for(i in 4501:5000){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else if(chunk==11){
#   for(i in 5001:5500){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }else{
#   for(i in 5501:5840){
#     setwd(paste('/projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/',substr(flist[i],79,82),'/',substr(flist[i],83,85),sep=''))
#     system(paste('wget --user=mkmoon --password=M159k258! ',flist[i],sep=''))   
#   }
# }
    
  

