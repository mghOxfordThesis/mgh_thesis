#################################################################################################
#
#  This script converts dewpoint to relative humidity. 
#
#  This script is run after:
#  1. 01_data2rasters.R
#  2. 02_curate_calculate.R
# 
#################################################################################################

#################################################################################################
print("<MVSE> estimating humidity (%) ...")

load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata"))    # temp_rast_list
load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_humrel_list.Rdata"))  # humrel_rast_list

# convert dewpoint to relative humidity
compute_rel_hum <- function(x) {
  temp <- temp_rast_list[[x]]
  dew <- humrel_rast_list[[x]]
  
  dew <- dew-273.15 ## convert dew to C (comes in K)
  lahr <- 100*(exp((17.625*dew)/(243.04+dew))/exp((17.625*temp)/(243.04+temp))) 
  lahr <- reclassify(lahr, cbind(-Inf, 0, 0), right=NA)        ##from [-Inf to 0]  inclusive <- 0
  lahr <- reclassify(lahr, cbind(100, +Inf, 100), right=NA)    ##from [100 to Inf] inclusive <- 100
  return(lahr)
}
load(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag, "_new_date_positions.Rdata")) #new_date_positions
humrel_rast_list[new_date_positions] <- purrr::map(new_date_positions, compute_rel_hum)
fileout<- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_humrel_list.Rdata")
save(humrel_rast_list, file=fileout)
file.remove(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/", region, tag,"_new_date_positions.Rdata"))
#unlink(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/"), recursive = TRUE)

