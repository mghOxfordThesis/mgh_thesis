#################################################################################################
#
#  This script curates the climate variables and collates the raster layers into a list. 
#  The temporary .tif files are deleted. 
#  
#  This script is run after:
#  1. 01_data2rasters.R
#
#################################################################################################

#################################################################################################
print("<MVSE> curating climate variables on geo-pixels ...")

#################################################################################################
# use one raster file to determine dimensions and save these for later
filepath<- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_temp_", ym_df[1, 1], ym_df[1, 2],".tif")
rast <- raster(x=filepath)
filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_dimRaster.Rdata")
rdim <- dim(rast)
save(rdim, file=filenameOut)

#################################################################################################
print("<MVSE> temperature ...")

curate_temp <- function(yr, mo) {
  filepath<- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_temp_",yr,mo,".tif")
  rast <- raster(x=filepath)
  rcl <- matrix(c(c(-Inf, 0, 0),        ## (-Inf, 0) <- 0
                  c(45, +Inf, 45)),     ## (45, Inf) <- 45
                ncol=3, byrow=TRUE)
  rast <- reclassify(rast, rcl, right=NA)
  return(rast)
}
temp_rast_list <- purrr::map(1:nrow(ym_df), function(x) curate_temp(yr=ym_df[x, 1], mo=ym_df[x, 2]))
if (!overwrite && length(list.files(path=RAW_CLIMATE_RASTER_FOLDER, pattern=paste0(region, tag, "_temp_list.Rdata")))==1) {
  tmp <- temp_rast_list
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")) #temp_rast_list
  temp_rast_list <- c(temp_rast_list, tmp)
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, "_date_ordering.Rdata")) #date_ordering
  temp_rast_list <- temp_rast_list[date_ordering]
}
fileout <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")
save(temp_rast_list, file=fileout)

#################################################################################################
print("<MVSE> dew point ...")

curate_hum <- function(yr, mo) {
  filepath<- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_hum_",yr,mo,".tif")
  rast <- raster(x=filepath)
  rcl <- matrix(c(-Inf, 0, 0),      ## (-Inf, 0) <- 0
                ncol=3, byrow=TRUE)
  rast <- reclassify(rast, rcl, right=FALSE)
  return(rast)
}
humrel_rast_list <- purrr::map(1:nrow(ym_df), function(x) curate_hum(yr=ym_df[x, 1], mo=ym_df[x, 2]))
if (!overwrite && length(list.files(path=RAW_CLIMATE_RASTER_FOLDER, pattern=paste0(region, tag, "_humrel_list.Rdata")))==1) {
  tmp <- humrel_rast_list
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_humrel_list.Rdata")) #humrel_rast_list
  humrel_rast_list <- c(humrel_rast_list, tmp)
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, "_date_ordering.Rdata")) #date_ordering
  humrel_rast_list <- humrel_rast_list[date_ordering]
}
fileout <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_humrel_list.Rdata")
save(humrel_rast_list, file=fileout)

#################################################################################################
print("<MVSE> precipitation ...")

curate_pr <- function(yr, mo) {
  filepath<- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_prec_",yr,mo,".tif")
  rast <- raster(x=filepath)
  rcl <- matrix(c(-Inf, 0, 0),      ## (-Inf, 0) <- 0
                ncol=3, byrow=TRUE)
  rast <- reclassify(rast, rcl, right=FALSE)
  return(rast)
}
prec_rast_list <- purrr::map(1:nrow(ym_df), function(x) curate_pr(yr=ym_df[x, 1], mo=ym_df[x, 2]))
if (!overwrite && length(list.files(path=RAW_CLIMATE_RASTER_FOLDER, pattern=paste0(region, tag, "_prec_list.Rdata")))==1) {
  tmp <- prec_rast_list
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_prec_list.Rdata")) #prec_rast_list
  prec_rast_list <- c(prec_rast_list, tmp)
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, "_date_ordering.Rdata")) #date_ordering
  prec_rast_list <- prec_rast_list[date_ordering]
}
fileout <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_prec_list.Rdata")
save(prec_rast_list, file=fileout)

#################################################################################################
print("<MVSE> deleting temporary raster files ...")

remove_files <- function(yr, mo) {
  ## temperature
  filepath<- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_temp_",yr,mo,".tif")
  file.remove(filepath)
  ##vapor/humidity
  filepath<-paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_hum_",yr,mo,".tif")
  file.remove(filepath)
  ##precipitation
  filepath<- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_prec_",yr,mo,".tif")
  file.remove(filepath)
}
junk <- purrr::map(1:nrow(ym_df), function(x) remove_files(yr=ym_df[x, 1], mo=ym_df[x, 2]))

#file.remove(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_date_ordering.Rdata"))
