# -----------------------------------------------------------------------
##
## This script: downloads CHELSA climate data (PR and TAS)
##
# -----------------------------------------------------------------------
# Download CHELSA climate data 
#     TAS (mean daily air temp)
# FAQ: https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf
# Download: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F

library(tidyverse)  #basic use
library(sf)         #spatial manipulations
library(sp)         #spatial manipulations
library(raster)     #spatial manipulations
library(terra)      #spatial manipulations
library(utils)      #for the read.table and for loop read-in 



# CHELSA_monthly_tas_1980-2019 --------------------------------------------
##  Monthly Historic = 1979-2019; one .tiff per month
##  NOT the file downloaded from the "climatologies" section. 

# 1. List of urls of files to download (filelist)
filelist <- utils::read.table("data/climate_databases/CHELSA_1980-2018_global_MeanDailyAirTemp_MONTHLY/envidatS3paths_CHELSA_1980-2018_V2_Global_Monthly_MeanDailyAirTemp.txt") 

# Run these lines to see what part of file path is being read
ii=1 
url <- filelist[ii,]
strsplit(url,"/")[[1]][10]

# 2. Directory to download the files to (pathhead)
pathhead <- paste(getwd(), "data/climate_databases/CHELSA_1980-2018_global_MeanDailyAirTemp_MONTHLY/", sep="/") 
pathhead

# 3. Loop over filelist and do download (requires stable internet connection)
for(ii in 1:nrow(filelist)){ 
  url <- filelist[ii,]
  path <- paste(pathhead, strsplit(url, "/")[[1]][10], sep="/")
  try(download.file(url, path, method="curl", quiet=T))
}




# CHELSA_climatologies_ssp370_tas__2011-2040__2041-2070__2071-2100 --------

# 1. List of urls of files to download (filelist)
filelist <- utils::read.table("data/climate_databases/CHELSA_climatologies_ssp370_tas__2011-2040__2041-2070__2071-2100/CHELSA_climatologies_ssp370_tas__2011-2040__2041-2070__2071-2100__envidatS3paths.txt") 

# Run these lines to see what part of file path is being read
ii=1 
url <- filelist[ii,]
strsplit(url,"/")[[1]][13]

# 2. Directory to download the files to (pathhead)
pathhead <- paste(getwd(), "data/climate_databases/CHELSA_climatologies_ssp370_tas__2011-2040__2041-2070__2071-2100/", sep="/") 
pathhead

# 3. Loop over filelist and do download (requires stable internet connection)
for(ii in 1:nrow(filelist)){ 
  url <- filelist[ii,]
  path <- paste(pathhead, strsplit(url, "/")[[1]][13], sep="/")
  try(download.file(url, path, method="curl", quiet=T))
}



