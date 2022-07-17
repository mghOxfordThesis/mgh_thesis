#################################################################################################
#
#  This script reads and crops the climate data according to the specified region boundary 
#  and temporarily saves the climate data for each month as .tif files. 
#
#################################################################################################

#################################################################################################
print("<MVSE> reading maps into raster files  ...")

#################################################################################################
print("<< temperature >>")

new_dates <- c() ##will store the dates, assumes all sources of data have the same dates and order

for (ii in seq_along(ym_ind)) {
  pre1 <- raster(paste0(COPERNICUS_DATA_FOLDER, "temperature.grib.nc"), band=ym_ind[ii])
  print(paste(ii, " -- ", as.character(pre1@z)))
  new_dates <- c(new_dates, as.character(pre1@z))
  
  pre1 <- brick(pre1) ##easier to manage
  rast <- raster::rotate(pre1) ##transform to -180 180
  
  # crop_extent <- spTransform(crop_extent, CRS(proj4string(rast)))
  cr <- crop(rast, crop_extent, snap="out") 
  rast <- rasterize(crop_extent, cr, getCover=TRUE)
  rast[as.vector(rast<0.1)] <- NA # only consider pixels which the polygon covers at least 10% of area (arbitrary)
  rast <- mask(x=cr, mask=rast)
  rast <- rast-273.15 # adjust to degrees Celcius
  filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, 
                        "_temp_", ym_df[ii, 1], ym_df[ii, 2], ".tif")
  writeRaster(rast, filename=filenameOut,  overwrite=TRUE)
}

#################################################################################################
# save information on dates
if (!overwrite && length(list.files(path=RAW_CLIMATE_RASTER_FOLDER, pattern=paste0(region, tag, "_dates.Rdata")))==1) {
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_dates.Rdata")) #dates
  dates <- dates$date
  dates <- c(dates, new_dates)
  date_ordering <- match(sort(dates), dates)
  dates <- sort(dates)
  filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, "_date_ordering.Rdata")
  save(date_ordering, file=filenameOut)
  new_date_positions <- match(new_dates, dates)
} else {
  dates <- new_dates
  new_date_positions <- match(new_dates, dates)
}
dates <- data.frame(date=dates)
filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag,"_", "dates.Rdata")
save(dates, file=filenameOut)
filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, "_new_date_positions.Rdata")
save(new_date_positions, file=filenameOut)

#################################################################################################
print("<< humidity >>")

for (ii in seq_along(ym_ind)) {
  pre1 <- raster(paste0(COPERNICUS_DATA_FOLDER,"dewpoint.grib.nc"), band=ym_ind[ii])
  print(paste(ii, " -- ", as.character(pre1@z)))
  
  pre1 <- brick(pre1) ##easier to manage
  rast <- raster::rotate(pre1) ##transform to -180 180
  
  # crop_extent <- sp::spTransform(crop_extent, CRS(proj4string(rast)))
  cr <- crop(rast, crop_extent, snap="out") 
  rast <- rasterize(crop_extent, cr, getCover=TRUE)
  rast[as.vector(rast<0.1)] <- NA # only consider pixels which the polygon covers at least 10% of area (arbitrary)
  rast <- mask(x=cr, mask=rast)
  filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, 
                        "_hum_", ym_df[ii, 1], ym_df[ii, 2], ".tif")
  writeRaster(rast, filename=filenameOut,  overwrite=TRUE)
}

#################################################################################################
print("<< precipitation >>")

for (ii in seq_along(ym_ind)) {
  pre1 <- raster(paste0(COPERNICUS_DATA_FOLDER,"precipitation.grib.nc"), band=ym_ind[ii])
  print(paste(ii, " -- ", as.character(pre1@z)))
  
  pre1 <- brick(pre1) ##easier to manage
  rast <- raster::rotate(pre1) ##transform to -180 180
  
  # crop_extent <- sp::spTransform(crop_extent, CRS(proj4string(rast)))
  cr <- crop(rast, crop_extent, snap="out") 
  rast <- rasterize(crop_extent, cr, getCover=TRUE)
  rast[as.vector(rast<0.1)] <- NA # only consider pixels which the polygon covers at least 10% of area (arbitrary)
  rast <- mask(x=cr, mask=rast)
  filenameOut <- paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region,
                        "_prec_", ym_df[ii, 1], ym_df[ii, 2], ".tif")
  writeRaster(rast, filename=filenameOut,  overwrite=TRUE)
}
