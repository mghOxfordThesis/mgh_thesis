#################################################################################################
#
#  This script is estimates the posterior distribution of the Index P time series for the 
#  the entire region by aggregating the climate data across all the pixels. 
# 
#  This script is run after:
#  1. 00_simplify_SHP_files.R (optional)
#  2. 00_generate_MVSE_maps.R
#  3. 01_run_MVSE_per_cell.R
#
#################################################################################################

########################################################################################
cat("## Step 1b: run MCMC for entire region ... \n")
step_1b_time <- Sys.time()

########################################################################################
# just to get raster matrix dims
fileInDim <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_dimRaster.Rdata")
load(fileInDim) # rdim
ori_l<- rdim[1]
ori_c<- rdim[2]
cat(paste("## data with",ori_l,"lines","and", ori_c, "columns \n"))

########################################################################################
cat("## loading MVSE input data \n")

# compute mean for each climate variable across the region at each time point
temp_data <- matrix(NA, nrow=length(dates), ncol=length(list_cell_THP))
humrel_data <- matrix(NA, nrow=length(dates), ncol=length(list_cell_THP))
prec_data <- matrix(NA, nrow=length(dates), ncol=length(list_cell_THP))
for (ii in seq_along(list_cell_THP)) {
  temp_data[, ii] <- list_cell_THP[[ii]][["T"]]
  humrel_data[, ii] <- list_cell_THP[[ii]][["H"]]
  prec_data[, ii] <- list_cell_THP[[ii]][["R"]]
}
mean_temp <- rowMeans(temp_data, na.rm=TRUE)
mean_humrel <- rowMeans(humrel_data, na.rm=TRUE)
mean_prec <- rowMeans(prec_data, na.rm=TRUE)
region_climate_data <- data.frame(date=list_cell_THP[[1]]$date, 
                                  year=list_cell_THP[[1]]$year, 
                                  month=list_cell_THP[[1]]$month, 
                                  day=list_cell_THP[[1]]$day, 
                                  T=mean_temp,  
                                  H=mean_humrel, 
                                  R=mean_prec)

########################################################################################
cat('## run MVSE for entire region...\n')
if (is.null(mvse_model_args$model_category) || (mvse_model_args$model_category!="aegypti")) {
  this_mvse_model <- mvse_model(priors=mvse_model_args$priors, climate_data=region_climate_data)
} else {
  this_mvse_model <- mvse_model(model_category="aegypti", climate_data=region_climate_data)
}
this_mvse_fit <- fitting(this_mvse_model, iter=fitting_args$iter, warmup=fitting_args$warmup, 
                         seed=fitting_args$seed, init=fitting_args$init, 
                         gauJump=fitting_args$gauJump, samples=fitting_args$samples, 
                         verbose=TRUE)
region_indexP_sims <- MVSE::extract(this_mvse_fit, pars="indexP")$indexP
# region_Q_sims <- MVSE::extract(this_mvse_fit, pars="Q")$Q
# region_V0_sims <- MVSE::extract(this_mvse_fit, pars="V0")$V0
fileout <- paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_", "indexP_month_region.Rdata")
save(region_indexP_sims, file=fileout)
# fileout <- paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_", "Q_month_region.Rdata")
# save(region_Q_sims, file=fileout)
# fileout <- paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_", "V0_month_region.Rdata")
# save(region_V0_sims, file=fileout)

step_1b_time <- Sys.time() - step_1b_time
cat("## step 1b time:", step_1b_time, attributes(step_1b_time)$units, "\n")
cat("\n")
