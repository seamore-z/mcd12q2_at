args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- as.numeric(substr(args[3],7,10))
# tile <- 'h11v04'

# Check number of files: 1280
sstr <- paste('*',tile,'*hdf',sep='')

print(sstr)

lnga2 <- length(list.files('/projectnb/modislc/projects/sat/data/e4ftl01.cr.usgs.gov/MOTA/MCD43A2.061/',pattern=glob2rx(sstr),recursive=T,full.names=T))
lnga4 <- length(list.files('/projectnb/modislc/projects/sat/data/e4ftl01.cr.usgs.gov/MOTA/MCD43A4.061/',pattern=glob2rx(sstr),recursive=T,full.names=T))

print(paste(lnga2,'; ',lnga4))


# if(lnga2==1095 & lnga4==1095){
  
  # system('rm -r /projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a2_link/*')
  # system('rm -r /projectnb/modislc/users/mkmoon/MuSLI/c6_2018/mcd43a4_link/*')
  # 
  # # Make directory
  # for(yy in 2001:2016){
  #   fld_yy <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a2_link/',yy,sep='')
  #   system(paste('mkdir ',fld_yy,sep=''))
  #   fld_yy <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a4_link/',yy,sep='')
  #   system(paste('mkdir ',fld_yy,sep=''))
  #   for(dd in 1:365){
  #     doy <- sprintf("%03d",dd)
  #     fld_dd <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a2_link/',yy,'/',doy,sep='')
  #     system(paste('mkdir ',fld_dd,sep=''))
  #     fld_dd <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a4_link/',yy,'/',doy,sep='')
  #     system(paste('mkdir ',fld_dd,sep=''))
  #   }
  #   print(yy)
  # }
  
  # Copy and change folder names
#  if(year==2000|year==2004|year==2008|year==2012|year==2016|year==2020){
    nday <- 365
 # }else{
  #  nday <- 365  
  #}
  
  fld_new <- paste('/projectnb/modislc/projects/sat/data/mcd43a4_link/',year,sep='')
  if (!dir.exists(fld_new)) {dir.create(fld_new)}

  fld_new <- paste('/projectnb/modislc/projects/sat/data/mcd43a2_link/',year,sep='')
  if (!dir.exists(fld_new)) {dir.create(fld_new)}
  
  mod_dates <- seq(1,365)
  for(i in 1:nday){
    d = mod_dates[i]
    d_str = sprintf("%03d",d)
    dates <- as.Date(d,format='%Y.%m.%d',origin=paste((year-1),'.12.31',sep=''))
    mm <- substr(dates,6,7)
    dd <- substr(dates,9,10)
    doy <- sprintf("%03d",i)
    
    sstr <- paste('*',tile,'*hdf',sep='')
    
    fld_old <- list.files(paste('/projectnb/modislc/projects/sat/data/e4ftl01.cr.usgs.gov/MOTA/MCD43A4.061/',year,'.',mm,'.',dd,sep=''),pattern=glob2rx(sstr),full.names=T)
    fld_new <- paste('/projectnb/modislc/projects/sat/data/mcd43a4_link/',year,'/',doy,sep='')
    if (!dir.exists(fld_new)) {dir.create(fld_new)}
    system(paste('cp ',fld_old,' ',fld_new,sep=''))    
    
    fld_old <- list.files(paste('/projectnb/modislc/projects/sat/data/e4ftl01.cr.usgs.gov/MOTA/MCD43A2.061/',year,'.',mm,'.',dd,sep=''),pattern=glob2rx(sstr),full.names=T)
    fld_new <- paste('/projectnb/modislc/projects/sat/data/mcd43a2_link/',year,'/',doy,sep='')
    if (!dir.exists(fld_new)) {dir.create(fld_new)}
    system(paste('cp ',fld_old,' ',fld_new,sep=''))    
    
    if(i%%10==0) print(paste(year,'; DOY ',i))
  }
  
# }

# from 2017 to 2020 (365x3+185=1280)
length(list.files('/projectnb/modislc/projects/sat/data/mcd43a2_link/',pattern=glob2rx(sstr),recursive=T))
length(list.files('/projectnb/modislc/projects/sat/data/mcd43a4_link/',pattern=glob2rx(sstr),recursive=T))



# ###############################
# sstr <- '*A2019*'
# aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a2_link/2020',pattern=glob2rx(sstr),recursive=T,full.names=T)
# for(i in 1:length(aa)){
#   file.rename(aa[i],paste0(substr(aa[i],1,75),'2020',substr(aa[i],80,121)))  
# }
# aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a4_link/2020',pattern=glob2rx(sstr),recursive=T,full.names=T)
# for(i in 1:length(aa)){
#   file.rename(aa[i],paste0(substr(aa[i],1,75),'2020',substr(aa[i],80,121)))  
# }
# 
# sstr <- '*A2020*'
# aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a2_link/2021',pattern=glob2rx(sstr),recursive=T,full.names=T)
# for(i in 1:length(aa)){
#   file.rename(aa[i],paste0(substr(aa[i],1,75),'2021',substr(aa[i],80,121)))  
# }
# aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/mcd43a4_link/2021',pattern=glob2rx(sstr),recursive=T,full.names=T)
# for(i in 1:length(aa)){
#   file.rename(aa[i],paste0(substr(aa[i],1,75),'2021',substr(aa[i],80,121)))  
# }

# sstr <- '*A2018*'
# aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/input/mcd43a2/i2018/2017',pattern=glob2rx(sstr),recursive=T,full.names=T)
# for(i in 1:length(aa)){
#   file.rename(aa[i],paste0(substr(aa[i],1,87),'2017',substr(aa[i],92,150)))
# }
#sstr <- '*A2020*'
#aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/input/mcd43a2/h12v04/i2020/2021',pattern=glob2rx(sstr),recursive=T,full.names=T)

#for(i in 1:length(aa)){
# file.rename(aa[i],paste0(substr(aa[i],1,94),'2021',substr(aa[i],99,150)))
#}



# sstr <- '*A2018*'
# aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/input/mcd43a4/i2018/2017',pattern=glob2rx(sstr),recursive=T,full.names=T)
# for(i in 1:length(aa)){
#   file.rename(aa[i],paste0(substr(aa[i],1,87),'2017',substr(aa[i],92,150)))
# }
#sstr <- '*A2020*'
#aa <- list.files('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/input/mcd43a4/h12v04/i2020/2021',pattern=glob2rx(sstr),recursive=T,full.names=T)
#for(i in 1:length(aa)){
#  file.rename(aa[i],paste0(substr(aa[i],1,94),'2021',substr(aa[i],99,150)))
#}

