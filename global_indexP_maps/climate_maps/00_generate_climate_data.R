#################################################################################################
#
#  This script is the "master script" that generates climate data (temperature, humidity and 
#  precipitation) at the pixel level for a region of interest (e.g. a country). For our analyses
#  we use the climate data obtained from Copernicus (i.e. temperature.grib.nc, 
#  dewpoint.grib.nc, precipitation.grib.nc) and spatial polygons obtained from GADM. 
#
#  After setting up the required directories, this script runs five other scripts in succession:
#  1. 01_data2rasters.R
#  2. 02_curate_calculate.R
#  3. 03_do_relative_humidity_range.R
#  4. 04_plotting_climate_maps.R
#  5. 05_organize_for_MVSE.R
#
#################################################################################################

#################################################################################################
# load libraries
library(pacman)
pacman::p_load(sf, pbapply, rgdal, tidyverse, raster, ncdf4, 
               RColorBrewer, rasterVis, lubridate)

# directories
EXTERNAL_HARD_DRIVE_PATH <- ""
SHP_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "maps_data_SHP/")
COPERNICUS_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "maps_data_Copernicus2/")
RAW_CLIMATE_RASTER_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "climate_data/raw_raster_data/")
RAW_CLIMATE_PLOT_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "climate_data/raster_plots/")
MVSE_INPUT_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "climate_data/mvse_input_data/")


#################################################################################################
# get_climate_data takes a region name (region) (e.g. Tanzania), a SpatialPolygonsDataFrame 
# from a shapefile (crop_extent), vector of focal years (years), a naming tag to be used as 
# an identifier on the filename along with the region (tag), and a logical flag to indicate
# whether existing climate data for the selected region should be overwritten if it exists 
# (overwrite). It retrieves the climate data (temperature,relative humidity and rainfall) 
# for the region of interest along with some intermediate output necessary to perform
# Index P estimation at the pixel level for the region. 
get_climate_data <- function(region, crop_extent, years, tag="", overwrite=FALSE) {
  # verify that data is available for the requested years 
  available_years <- 1981:2019
  if (length(intersect(years, available_years))!=length(years)) {
    unavailable_years <- setdiff(years, intersect(years, available_years))
    stop(cat("Climate data for the following years are not available: ", 
             unavailable_years))
  }
  
  # set up directories
  RAW_CLIMATE_RASTER_FOLDER <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, "/")
  RAW_CLIMATE_PLOT_FOLDER <- paste0(RAW_CLIMATE_PLOT_FOLDER, region, "/")
  MVSE_INPUT_DATA_FOLDER <- paste0(MVSE_INPUT_DATA_FOLDER, region, "/")
  if (!file.exists(RAW_CLIMATE_RASTER_FOLDER)) dir.create(RAW_CLIMATE_RASTER_FOLDER)
  if (!file.exists(RAW_CLIMATE_PLOT_FOLDER)) dir.create(RAW_CLIMATE_PLOT_FOLDER)
  if (!file.exists(MVSE_INPUT_DATA_FOLDER)) dir.create(MVSE_INPUT_DATA_FOLDER)
  if (!file.exists(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/"))) dir.create(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary/"))
  if (tag!="") tag <- paste0("_", tag)
  
  # check the years for which climate data has already been generated
  if (!overwrite) {
    available_files <- list.files(path=MVSE_INPUT_DATA_FOLDER)
    if (length(available_files)==0) {required_years <- years}
    
    else if (paste0(region, "_list_cell_THP.Rdata") %in% available_files) {
      load(paste0(MVSE_INPUT_DATA_FOLDER, region, "_list_cell_THP.Rdata")) #list_cell_THP
      dates <- list_cell_THP[[1]]$date
      completed_years <- unique(lubridate::year(dates))
      required_years <- setdiff(years, intersect(years, completed_years))
      if (length(required_years)==0) {
        cat("## climate data for the specified years has already been generated.\n")
        cat("## set overwrite=TRUE to regenerate the climate data \n")
        return()
      }
    } else {
      cat("## climate data already generated for", intersect(years, completed_years), "\n")
      required_years <- years
    }
  } else {
    required_years <- years
  }
  
  # indices for the multi-layer raster files
  ym_ind <- as.vector(sapply(required_years, 
                             function(x) seq(12*(which(available_years==x)-1)+1, 
                                             12*(which(available_years==x)))))
  ym_df <- data.frame(year=required_years, month=rep(1:12, each=length(required_years))) %>% 
    arrange(desc(-year), desc(-month))
  
  # run the required code
  source('01_data2rasters.R', local=TRUE)
  source("02_curate_calculate.R", local=TRUE)
  source("03_do_relative_humidity_range.R", local=TRUE)
  source('04_plotting_climate_maps.R', local=TRUE)
  source('05_organize_for_MVSE.R', local=TRUE)
  return()
}

#################################################################################################
# Here is an example of how it might be used 

# read in the shapefile for the country of interest, which in this case is Haiti
crop_extent <- readOGR(file.path(SHP_DATA_FOLDER, "Haiti", "gadm40_HTI_0.shp"))

# obtain the climate data
get_climate_data(region="Haiti", crop_extent=crop_extent, years=1981:2019, overwrite=TRUE)
