#################################################################################################
#
#  This script is the "master script" that estimates Index P for each pixel within an area of 
#  interest (typically a country). A collection of summary statistics of Index P are calculated. 
#  The climate data must be generated via 00_generate_climate_data.R before this script can be 
#  run. 
#
#  After setting up the required directories, this script runs five other scripts in succession:
#  1. 01_run_MVSE_per_cell.R
#  2. 01_run_MVSE_region.R
#  3. 02_run_MVSE_per_cell.R
#  4. 03_MVSE2raster.R
#  5. 04_draw_maps.R
#
#################################################################################################


#################################################################################################
# load libraries
library(pacman)
pacman::p_load(MVSE)
pacman::p_load(parallel, raster, tidyverse, RColorBrewer, rasterVis,
               cowplot, ggpubr, lubridate, rgdal, sf, lme4, Rcpp, terra)

# folder path
EXTERNAL_HARD_DRIVE_PATH <- ""
RAW_CLIMATE_RASTER_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "../climate_maps/climate_data/raw_raster_data/")
MVSE_INPUT_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "../climate_maps/climate_data/mvse_input_data/")
MVSE_OUTPUT_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "mvse_output/")
SHP_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "../climate_maps/maps_data_SHP/")

#################################################################################################
# gen_MVSE_maps takes a region name (region), an associated tag for MVSE output files (tag),  
# the tag for the requested climate files (climate_tag), the requested years (years), the 
# arguments to construct a mvse_model object (mvse_model_args), arguments for the fitting procedure, 
# (fitting_args), the number of threads/cores (no_clus), a spatial polygons object for the region 
# of interest (crop_extent). It generates a collection of summary statistics of the estimated 
# Index P for the region of interest. overwrite=FALSE simply checks whether MVSE has already been 
# simulated for years of interest and if so, it will throw an error. `all_plots=TRUE` results in 
# only a subset of the diagnostic plots being generated. `delete_raw_output=TRUE` ensures that 
# the raw Index P distribution data is deleted once the summary statistics have been generated. 
gen_MVSE_maps <- function(region, tag="", climate_tag="", years, mvse_model_args, 
                          fitting_args, no_clus, overwrite=FALSE, all_plots=TRUE, delete_raw_output=FALSE, crop_extent) {
  # check if climate data for the region of interest exists
  RAW_CLIMATE_RASTER_FOLDER <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, "/")
  MVSE_INPUT_DATA_FOLDER <- paste0(MVSE_INPUT_DATA_FOLDER, region, "/")
  if (!dir.exists(MVSE_INPUT_DATA_FOLDER))
    stop(cat("## Input climate data not available for ", region, "\n"))
  if (!dir.exists(RAW_CLIMATE_RASTER_FOLDER))
    stop(cat("## Raw raster data not available for ", region, "\n"))
  
  # ensure years are provided
  if (is.null(years)) {
    stop(cat("## Must request some years. \n"))
  }
  
  # ensure that provided years are continuous
  if (!all(diff(years)==1)) {
    stop(cat("## Provided years must be continuous over a given interval. \n"))
  }
  
  # check if climate data is available for years of interest
  if (climate_tag!="") climate_tag <- paste0("_", climate_tag)
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, climate_tag, "_dates.Rdata")) #dates
  dates <- as.character(dates$date)
  available_years <- unique(lubridate::year(dates))
  available_years <- sort(available_years)
  if (length(base::intersect(available_years, years))!=length(years))
    stop(cat("## Climate data on requested years not available.", "Data for", available_years, 
             "are available. \n"))
  
  # create required directories
  if (tag!="") tag <- paste0("_", tag)
  MVSE_OUTPUT_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, region, "/")
  MVSE_OUTPUT_DATA_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, "data/")
  MVSE_OUTPUT_PLOTS_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, "plots/")
  if (!dir.exists(MVSE_OUTPUT_FOLDER)) dir.create(MVSE_OUTPUT_FOLDER)
  if (!dir.exists(MVSE_OUTPUT_DATA_FOLDER)) dir.create(MVSE_OUTPUT_DATA_FOLDER)
  if (!dir.exists(MVSE_OUTPUT_PLOTS_FOLDER)) dir.create(MVSE_OUTPUT_PLOTS_FOLDER)
  if (!dir.exists(paste0(MVSE_OUTPUT_DATA_FOLDER, "temporary/"))) dir.create(paste0(MVSE_OUTPUT_DATA_FOLDER, "temporary/"))
  if (!dir.exists(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"))) dir.create(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"))
  if (!dir.exists(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/"))) dir.create(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/"))
  if (!dir.exists(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/pixel_summary/"))) dir.create(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/pixel_summary/"))
  
  # check if simulations have already been run for the years of interest
  if (!overwrite) {
    if (file.exists(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_indexP_month_region.Rdata"))) {
      load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_indexP_month_region.Rdata")) #region_indexP_sims
      sim_dates <- region_indexP_sims$date
      simulated_years <- unique(lubridate::year(sim_dates))
      if (length(intersect(years, simulated_years))==length(years)) {
        stop(cat("## MVSE already performed for the years of interest \n## set overwrite=TRUE to overwrite the previously simulated data. \n"))
      }
    }
  }
  
  # curtail data to the years of interest
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, climate_tag, "_dates.Rdata")) #dates
  dates <- as.character(dates$date)
  dates <- dates[which(lubridate::year(dates) %in% years)]
  load(paste0(MVSE_INPUT_DATA_FOLDER, region, climate_tag, "_list_cell_THP.Rdata")) #list_cell_THP
  for (ii in seq_along(list_cell_THP)) {
    list_cell_THP[[ii]] <- filter(list_cell_THP[[ii]], year %in% years)
  }
  
  # create required directories
  overall_time <- Sys.time()
  source('01_run_MVSE_per_cell.R', local=TRUE)
  source('01_run_MVSE_region.R', local=TRUE)
  source('02_run_MVSE_per_cell.R', local=TRUE)
  source('03_MVSE2raster.R', local=TRUE)
  source('04_draw_maps.R', local=TRUE)
  overall_time <- Sys.time() - overall_time
  cat("## overall time:", overall_time, attributes(overall_time)$units, "\n")
  return()
}

#################################################################################################]
# here is an example of how it might be run. We take Haiti as an example

# define the arguments for the MVSE package
mvse_model_args <- list(model_category="aegypti")
fitting_args <- list(iter=0.2*10^5, warmup=0.3, init=c(rho=1, eta=10, alpha=3), 
                     gauJump=c(rho=2, eta=5), samples=1000, seed=123)

# load the simplified polygons object for Haiti
load(file=file.path(SHP_DATA_FOLDER, "Haiti", "simplified_gadm40_HTI_0.Rdata"))
crop_extent <- st_as_sf(crop_extent)
st_crs(crop_extent) <- 4326

# generate summary statistics of Index P for the region of interest
gen_MVSE_maps(region=country, years=1981:2019, mvse_model_args=mvse_model_args,
              fitting_args=fitting_args, no_clus=1, overwrite=TRUE, all_plots=FALSE, delete_raw_output=TRUE,
              crop_extent=crop_extent)
