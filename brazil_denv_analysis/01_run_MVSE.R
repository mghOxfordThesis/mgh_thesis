#################################################################################################
#
# We perform the mvse fitting for all the municipalities in Brazil for which climate data is 
# available. A .rds object containing the sampling output for Index P is generated and saved 
# for each municipality. 
# 
# The following scripts must be run and their output must be generated before this script 
# can be run:
#  1. 01_organise_population_data.R
#  2. 01_organise_adm2_sf_data.R
#  3. transform_grib_2_nc.sh
#  4. 02_generate_climate_data.R
#
#################################################################################################

#################################################################################################
# Definition of required packages, directories and local data

# libraries
library(pacman)
p_load(MVSE, tidyverse)

# input directories
CASE_DATA_FOLDER <- "data/bra2_cases/"
CLIMATE_DATA_FOLDER <- "data/bra2_climate/"

# output directories
RAW_OUTPUT_DATA_FOLDER <- "output/data/raw_output/"

#################################################################################################
# `run_MVSE` takes a city name (this_city), a state name (this_state), a city ID (this_city_id), 
# a list of arguments for mvse model (mvse_model_args) and a list of arguments for fitting the 
# mvse model (fitting_args) as arguments and runs MVSE on that region. `overwrite is a logical 
# which determines whether existing (if available) is overwritten or not. The sampling output 
# of Index P are saved as .rds files in RAW_OUTPUT_DATA_FOLDER. `
run_MVSE <- function(this_city, this_state, this_city_id, mvse_model_args, fitting_args, overwrite=TRUE) {
  # retrieve case data for region
  if (!file.exists(file.path(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")))
    stop(cat("## Input DENV case data not available. \n"))
  
  # check if the index P has already been generated 
  if (!overwrite && file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata"))) {
    cat("## Index P already computed for ", this_city_id, ", ", this_state, ", ", this_city, "\n")
    return()
  }
  
  # retrieve climate data for region 
  climate_data_file <- paste0("BR_mun_", this_city_id, "_climate_entire_region.csv")
  if (file.size(paste0(CLIMATE_DATA_FOLDER, climate_data_file))<8000) {
    cat("## No climate data for ", this_city_id, ", ", this_state, ", ", this_city, "\n")
    return()
  }
  climate_data <- read.csv(file=paste0(CLIMATE_DATA_FOLDER, climate_data_file))
  
  # construct mvse model
  if (is.null(mvse_model_args$model_category) || (mvse_model_args$model_category!="aegypti")) {
    this_mvse_model <- mvse_model(priors=mvse_model_args$priors, climate_data=climate_data, warning=FALSE)
  } else {
    this_mvse_model <- mvse_model(model_category="aegypti", climate_data=climate_data, warning=FALSE)
  }
  
  # run sampling procedure
  this_mvse_fit <- mvse::fitting(this_mvse_model, iter=fitting_args$iter, warmup=fitting_args$warmup, 
                                 seed=fitting_args$seed, init=fitting_args$init, 
                                 gauJump=fitting_args$gauJump, samples=fitting_args$samples, verbose=FALSE)
  sim_p <- MVSE::extract(this_mvse_fit, pars="indexP")$indexP
  fileout <- paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata")
  save(sim_p, file=fileout)
  
}

#################################################################################################
# define mvse model arguments and fitting arguments. 
mvse_model_args <- list(model_category="aegypti")
fitting_args <- list(iter=10^5, warmup=0.3, init=c(rho=1, eta=10, alpha=3), 
                     gauJump=c(rho=2, eta=5), samples=1000, seed=123)

# run this for all municipalities. 
regions <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")) %>% 
  select(state, city, city_id) %>% 
  filter(!grepl(pattern="Município ignorado", x=city, fixed=TRUE)) %>% 
  unique()
{ apply(regions, 1, function(x) run_MVSE(this_city=x[2], this_state=x[1], this_city_id=x[3], 
                mvse_model_args=mvse_model_args, fitting_args=fitting_args)) }
