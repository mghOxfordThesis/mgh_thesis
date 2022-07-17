#################################################################################################
#
# We generate time series of the climate data for a municipality in Brazil. 
# It takes about 30 mins to obtain the climate data for each municipality. Given that there
# are 5570 municipalities (~116 days), it is best to organize climate data for several 
# municipalities in parallel. 
#
# Must run the following scripts prior to this script: 
#  1. 01_organise_population_data.R
#  2. 01_organise_adm2_sf_data.R
#  3. transform_grib_2_nc.sh
#
# Outputs: A .csv file with the monthly temperature, humidity and rainfall data from 2000 to 
# 2014 for each municipality. (5570 csv files in the CLIMATE_DATA_FOLDER)
# 
#################################################################################################

#################################################################################################
# Definition of required packages, directories and local data

# libraries
library(sf)
library(tidyverse)
library(raster)
library(ncdf4)

# directories
CASE_DATA_FOLDER <- "bra2_cases/"
CLIMATE_DATA_FOLDER <- "bra2_climate/"
SPATIAL_POLYGONS_DATA_FOLDER <- "bra2_sf/"
COPERNICUS_DATA_FOLDER <- "../global_indexP_maps/climate_maps/maps_data_Copernicus2/"
RAW_CLIMATE_RASTER_FOLDER <- "temporary/"
CLIMATE_MAPS_FOLDER <- "../global_indexP_maps/climate_maps/"
MVSE_INPUT_DATA_FOLDER <- "temporary/"

# local data
spatial_polygons_data <- readRDS(file.path(SPATIAL_POLYGONS_DATA_FOLDER, "simplified_bra2_sf_data.rds"))
case_data <- readRDS(file.path(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds"))

#################################################################################################
# Identify municipalities that do not have climate data

all_municipalities <- spatial_polygons_data %>% 
  st_drop_geometry() %>% 
  dplyr::select(state, city) %>% 
  left_join(unique(dplyr::select(case_data, city_id, state, city)), by=c("state", "city"))
has_climate_data <- function(this_state, this_city) {
  city_id <- filter(all_municipalities, state==this_state, city==this_city) %>%
    dplyr::select(city_id) %>% 
    unique()
  climate_data_file <- paste0("BR_mun_", city_id, "_climate_entire_region.csv")
  file_path <- file.path(CLIMATE_DATA_FOLDER, climate_data_file)
  if (file.size(file_path)<8000 || !file.exists(file_path)) { # if all entries NA then <8000 bytes
    return(FALSE)
  }
  return(TRUE)
}
municipalities_without_climate_data <- all_municipalities %>% 
  apply(1, function(x) has_climate_data(this_state=x[1], this_city=x[2])) %>%
  cbind(all_municipalities, has_climate_data=.) %>%
  filter(!has_climate_data)

#################################################################################################
# Generate climate data (temperature, humidity and precipitation) for each municipality

# `generate_climate_data` takes as arguments a state name (this_state), a municipality name (this_city), a municipality ID (this_city_id), 
# the years of interest (required_years) and a logical flag for whether any existing data should be overwritten (overwrite). 
generate_climate_data <- function(this_state, this_city, this_city_id, required_years=2000:2014, overwrite=TRUE) {
  # check if climate data is already available
  has_climate_data <- (municipalities_without_climate_data %>%
    filter(state==this_state, city==this_city) %>%
    nrow())==0
  if (!overwrite && has_climate_data) {
    stop(paste0(this_state, this_city, this_city_id, " already has climate data."))
  }

  # generate the climate data
  print(paste0("## starting ", this_state, this_city, this_city_id, "..."))
  if (!file.exists(RAW_CLIMATE_RASTER_FOLDER)) {dir.create(RAW_CLIMATE_RASTER_FOLDER)}
  if (!file.exists(paste0(RAW_CLIMATE_RASTER_FOLDER, "temporary"))) {dir.create(file.path(RAW_CLIMATE_RASTER_FOLDER, "temporary"))}
  crop_extent <- spatial_polygons_data %>% filter(state==this_state, city==this_city) %>% sf::as_Spatial()
  available_years <- 1981:2019
  ym_ind <- as.vector(sapply(required_years, 
                             function(x) seq(12*(which(available_years==x)-1)+1, 
                                             12*(which(available_years==x)))))
  ym_df <- data.frame(year=required_years, month=rep(1:12, each=length(required_years))) %>% 
    arrange(desc(-year, -month))
  region <- this_city_id
  tag <- ""
  source(file.path(CLIMATE_MAPS_FOLDER, '01_data2rasters.R'), local=TRUE)
  source(file.path(CLIMATE_MAPS_FOLDER, "02_curate_calculate.R"), local=TRUE)
  source(file.path(CLIMATE_MAPS_FOLDER, "03_do_relative_humidity_range.R"), local=TRUE)
  source(file.path(CLIMATE_MAPS_FOLDER, "05_organize_for_MVSE.R"), local=TRUE)
  
  # compute the average climate data of region
  load(paste0(MVSE_INPUT_DATA_FOLDER, region, tag, "_list_cell_THP.Rdata")) #list_cell_THP
  dates <- list_cell_THP[[1]]$date
  years <- list_cell_THP[[1]]$year
  months <- list_cell_THP[[1]]$month
  days <- list_cell_THP[[1]]$day
  T_list <- list()
  H_list <- list()
  R_list <- list()
  for (ii in seq_along(list_cell_THP)) {
    T_list[[ii]] <- list_cell_THP[[ii]][["T"]]
    H_list[[ii]] <- list_cell_THP[[ii]][["H"]]
    R_list[[ii]] <- list_cell_THP[[ii]][["R"]]
  }
  avg_T <- rowMeans(do.call(cbind, T_list), na.rm=TRUE)
  avg_H <- rowMeans(do.call(cbind, H_list), na.rm=TRUE)
  avg_R <- rowMeans(do.call(cbind, R_list), na.rm=TRUE)
  climate_data <- data.frame(date=dates, T=avg_T, H=avg_H, R=avg_R, year=years, month=months, day=days)
  outfile_name <- paste0("BR_mun_", this_city_id, "_climate_entire_region.csv")
  write.csv(climate_data, file.path(CLIMATE_DATA_FOLDER, outfile_name), row.names=FALSE)
  
  # delete temporary climate files
  unlink(file.path(RAW_CLIMATE_RASTER_FOLDER), recursive=TRUE)
  
  return()
}

# generate climate data for required regions
input_mat <- municipalities_without_climate_data %>% 
  dplyr::select(state, city, city_id) %>%
  as.matrix()
parallel::mcmapply(FUN=generate_climate_data, this_state=input_mat[, 1], this_city=input_mat[, 2], 
                   this_city_id=input_mat[, 3], mc.cores=2)
