#################################################################################################
#
#  This script is calculates summary statistics of the posterior distribution of Index P. 
#  A full list of the summary statistics is provided with the description of the data set. 
# 
#  This script is run after:
#  1. 00_simplify_SHP_files.R (optional)
#  2. 00_generate_MVSE_maps.R
#  3. 01_run_MVSE_per_cell.R
#  4. 01_run_MVSE_region.R
#  5. 02_run_MVSE_per_cell.R
#
#################################################################################################

##################################################################
cat("## Step 3: making MVSE rasters and organizing MVSE output ...\n")
step_3_time <- Sys.time()

##################################################################
## just to get raster matrix dims
fileInDim <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_dimRaster.Rdata")
load(fileInDim) # rdim
ori_l <- rdim[1]
ori_c <- rdim[2]
N <- ori_l*ori_c

##################################################################
cat('## organize data per pixel...\n')

# compute summary statistics for each pixel for each month 
# (mean, median, sd, 95% credible interval)
# the output is saved as chunks to overcome any memory issues (these chucks are then rejoined).
summarise_index <- function(index) {
  raw_output_files <- list.files(path=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"),
                                 pattern=paste0(region, tag, "_", index, "_month_per_cell_parallel*"))
  num_files <- length(raw_output_files)
  for (ii in 1:num_files) {
    raw_output_files[[ii]] <- paste0(region, tag, "_", index, "_month_per_cell_parallel", ii, ".Rdata")
  }
  for (ii in seq_along(raw_output_files)) {
    load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", raw_output_files[[ii]])) #allPs, allQs or allV0s
    if (index=="indexP") data <- allPs
    else if (index=="Q") data <- allQs
    else data <- allV0s
    
    mean_month_list <- list()
    median_month_list <- list()
    sd_month_list <- list()
    l95_month_list <- list()
    u95_month_list <- list()
    
    for (kk in seq_along(data)) {
      if (length(data[[kk]])==1 && is.na(data[[kk]])) {
        mean_month_list[[kk]] <- rep(NA, length(dates))
        median_month_list[[kk]] <- rep(NA, length(dates))
        sd_month_list[[kk]] <- rep(NA, length(dates))
        l95_month_list[[kk]] <- rep(NA, length(dates))
        u95_month_list[[kk]] <- rep(NA, length(dates))
      }
      else {
        this_data <- as.matrix(data[[kk]][, -1])
        mean_month_list[[kk]] <- round(rowMeans(this_data), 3)
        median_month_list[[kk]] <- round(apply(this_data, 1, median), 3)
        sd_month_list[[kk]] <- round(apply(this_data, 1, sd), 3)
        l95_month_list[[kk]] <- round(apply(this_data, 1, function(x) quantile(x, probs=0.025, names=FALSE)), 3)
        u95_month_list[[kk]] <- round(apply(this_data, 1, function(x) quantile(x, probs=0.975, names=FALSE)), 3)
      }
    }
    monthly_mean <- do.call(cbind, mean_month_list)
    monthly_median <- do.call(cbind, median_month_list)
    monthly_sd <- do.call(cbind, sd_month_list)
    monthly_l95 <- do.call(cbind, l95_month_list)
    monthly_u95 <- do.call(cbind, u95_month_list)
    output <- list(monthly_mean=monthly_mean, monthly_median=monthly_median,
                   monthly_sd=monthly_sd, monthly_l95=monthly_l95,
                   monthly_u95=monthly_u95)
    filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/pixel_summary/", region, tag, "_", index, "_month_summary_per_pixel", ii, ".Rdata")
    eval(parse(text=paste0(index, "_month_summary_per_pixel", "<-", "output")))
    eval(parse(text=paste0("save(", index, "_month_summary_per_pixel", ", ", "file=filename)")))
    if (index=="indexP") rm(allPs)
    else if (index=="Q") rm(allQs)
    else rm(allV0s)
    rm(data);
  }
  return()
}
{ summarise_index(index="indexP") }
# { summarise_index(index="Q") }
# { summarise_index(index="V0") }

# collate monthly summary statistics for each pixel into a list of length 5 (monthly_mean, monthly_medinan, etc.)
# with matrices of dimension (length(dates), number of pixels)
collate_summary_statistics <- function(index) {
  # get monthly summary data for each pixel
  pixel_summary_files <- list.files(path=paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/pixel_summary/"),
                                    pattern=paste0(region, tag, "_", index, "_month_summary_per_pixel*"))
  pixel_summary_files <- pixel_summary_files[grepl(".Rdata$", pixel_summary_files)]
  num_files <- length(pixel_summary_files)
  for (ii in 1:num_files) {
    pixel_summary_files[[ii]] <- paste0(region, tag, "_", index, "_month_summary_per_pixel", ii, ".Rdata")
  }
  
  if (num_files==1) {
    load(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/pixel_summary/", pixel_summary_files[[1]])) #indexP_month_summary_per_pixel
    eval(parse(text=paste0("index_monthly_summary", "<-", index, "_month_summary_per_pixel")))
  } else {
    index_monthly_summary_list <- list()
    for (ii in seq_along(pixel_summary_files)) {
      load(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/pixel_summary/", pixel_summary_files[[ii]])) #indexP_month_summary_per_pixel
      eval(parse(text=paste0("tmp", "<-", index, "_month_summary_per_pixel")))
      index_monthly_summary_list[[ii]] <- tmp
    }
    index_monthly_summary <- list()
    for (ii in seq_along(index_monthly_summary_list[[1]])) {
      tmp <- purrr::map(seq_along(index_monthly_summary_list), function(x) return(index_monthly_summary_list[[x]][[ii]]))
      index_monthly_summary[[ii]] <- do.call(cbind, tmp)
    }
    names(index_monthly_summary) <- names(index_monthly_summary_list[[1]])
  }
  return(index_monthly_summary)
}
indexP_monthly_summary <- collate_summary_statistics(index="indexP")
#Q_monthly_summary <- collate_summary_statistics(index="Q")
#V0_monthly_summary <- collate_summary_statistics(index="V0")

# generate the raster bricks of the summary statistics for each month
create_monthly_rasters <- function(index_monthly_summary, index) {
  # load raster template
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")) #temp_rast_list
  rastEx <- temp_rast_list[[1]] # use temperature raster list as a template
  
  for (ii in seq_along(index_monthly_summary)) {
    data <- index_monthly_summary[[ii]]
    data_name <- names(index_monthly_summary)[ii]
    output <- list()
    for (jj in 1:nrow(data)) {
      rastEx[] <- matrix(data[jj, ], ncol=ori_c, nrow=ori_l, byrow=TRUE)
      output[[jj]] <- rastEx
    }
    output <- raster::brick(output)
    names(output) <- as.character(dates)
    output <- terra::rast(output)
    filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", data_name, "_rasters.tif")
    eval(parse(text=paste0(index, "_", data_name, "_rasters", "<-", "output")))
    eval(parse(text=paste0("terra::writeRaster(", index, "_", data_name, "_rasters", ", ", "file=filename, overwrite=TRUE)")))
  }
  return()
}
{
  create_monthly_rasters(indexP_monthly_summary, "indexP")
  #create_monthly_rasters(Q_monthly_summary, "Q")
  #create_monthly_rasters(V0_monthly_summary, "V0")
}

# create a dataframe with all monthly summary data (date, long, lat, ...)
get_full_data <- function(index_monthly_summary, index) {
  all_data_output_df <- data.frame()
  for (ii in seq_along(index_monthly_summary)) {
    data_name <- names(index_monthly_summary)[ii]
    output <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", data_name, "_rasters.tif"))

    data_output_list <- list()
    for (t in 1:length(dates)) {
      spts <- rasterToPoints(output[[t]], spatial=TRUE)
      spts <- as.data.frame(spts)
      spts <- cbind(spts, dates[t])
      colnames(spts) <- c(data_name, "lon", "lat", "date")
      data_output_list[[t]] <- spts
    }
    if (ii==1) all_data_output_df <- as.data.frame(do.call(rbind, data_output_list))
    else all_data_output_df <- left_join(all_data_output_df,
                                         as.data.frame(do.call(rbind, data_output_list)),
                                         by=c("lon", "lat", "date"))
  }
  all_data_output_df <- dplyr::select(all_data_output_df, c("date", "lon", "lat", dplyr::everything()))
  filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_monthly_summary_dataframe.Rdata")
  eval(parse(text=paste0(index, "_monthly_summary_dataframe", "<-", "all_data_output_df")))
  eval(parse(text=paste0("save(", index, "_monthly_summary_dataframe", ", ", "file=filename)")))
}
{
  get_full_data(indexP_monthly_summary, "indexP")
  #get_full_data(Q_monthly_summary, "Q")
  #get_full_data(V0_monthly_summary, "V0")
}

# compute summary statistics for each pixel for each year
# (mean, median, sd, 95% credible interval)
# saved as raster bricks
create_yearly_rasters  <- function(index) {
  raw_output_files <- list.files(path=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"),
                                 pattern=paste0(region, tag, "_", index, "_month_per_cell_parallel*"))
  num_files <- length(raw_output_files)
  for (ii in 1:num_files) {
    raw_output_files[[ii]] <- paste0(region, tag, "_", index, "_month_per_cell_parallel", ii, ".Rdata")
  }
  
  mean_year_list <- list()
  median_year_list <- list()
  sd_year_list <- list()
  l95_year_list <- list()
  u95_year_list <- list()
  
  mean_months_above_one_list <- list()
  median_months_above_one_list <- list()
  sd_months_above_one_list <- list()
  l95_months_above_one_list <- list()
  u95_months_above_one_list <- list()
  
  years <- lubridate::year(dates) %>% unique() %>% sort()
  jj <- 1
  for (ii in seq_along(raw_output_files)) {
    load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", raw_output_files[[ii]])) #allPs, allQs or allV0s
    if (index=="indexP") data <- allPs
    else if (index=="Q") data <- allQs
    else data <- allV0s
    
    for (kk in seq_along(data)) {
      if (length(data[[kk]])==1 && is.na(data[[kk]])) {
        mean_year_list[[jj]] <- rep(NA, length(years))
        median_year_list[[jj]] <- rep(NA, length(years))
        sd_year_list[[jj]] <- rep(NA, length(years))
        l95_year_list[[jj]] <- rep(NA, length(years))
        u95_year_list[[jj]] <- rep(NA, length(years))
        
        mean_months_above_one_list[[jj]] <- rep(NA, length(years))
        median_months_above_one_list[[jj]] <- rep(NA, length(years))
        sd_months_above_one_list[[jj]] <- rep(NA, length(years))
        l95_months_above_one_list[[jj]] <- rep(NA, length(years))
        u95_months_above_one_list[[jj]] <- rep(NA, length(years))
      } else {
        this_data <- list()
        months_above_one <- list()
        for (t in seq_along(years)) {
          this_data[[t]] <- data[[kk]] %>% mutate(year=lubridate::year(date)) %>% filter(year==years[t])
          this_data[[t]] <- as.matrix(this_data[[t]] %>% dplyr::select(-date, -year))
          months_above_one[[t]] <- apply(this_data[[t]], 2, function(x) length(which(x>1)))
        }
        mean_year_list[[jj]] <- round(sapply(this_data, mean), 3)
        median_year_list[[jj]] <- round(sapply(this_data, median), 3)
        sd_year_list[[jj]] <- round(sapply(this_data, sd), 3)
        l95_year_list[[jj]] <- round(sapply(this_data, function(x) quantile(x, probs=0.025, names=FALSE)), 3)
        u95_year_list[[jj]] <- round(sapply(this_data, function(x) quantile(x, probs=0.975, names=FALSE)), 3)
        
        
        mean_months_above_one_list[[jj]] <- round(sapply(months_above_one, mean), 3)
        median_months_above_one_list[[jj]] <- round(sapply(months_above_one, median), 3)
        sd_months_above_one_list[[jj]] <- round(sapply(months_above_one, sd), 3)
        l95_months_above_one_list[[jj]] <- round(sapply(months_above_one, function(x) quantile(x, probs=0.025, names=FALSE)), 3)
        u95_months_above_one_list[[jj]] <- round(sapply(months_above_one, function(x) quantile(x, probs=0.975, names=FALSE)), 3)
      }
      jj <- jj + 1
    }
    if (index=="indexP") rm(allPs)
    else if (index=="Q") rm(allQs)
    else rm(allV0s)
    rm(data);
  }
  yearly_mean <- do.call(cbind, mean_year_list)
  yearly_median <- do.call(cbind, median_year_list)
  yearly_sd <- do.call(cbind, sd_year_list)
  yearly_l95 <- do.call(cbind, l95_year_list)
  yearly_u95 <- do.call(cbind, u95_year_list)
  yearly_months_above_one_mean <- do.call(cbind, mean_months_above_one_list)
  yearly_months_above_one_median <- do.call(cbind, median_months_above_one_list)
  yearly_months_above_one_sd <- do.call(cbind, sd_months_above_one_list)
  yearly_months_above_one_l95 <- do.call(cbind, l95_months_above_one_list)
  yearly_months_above_one_u95 <- do.call(cbind, u95_months_above_one_list)
  index_yearly_summary <- list(yearly_mean=yearly_mean, yearly_median=yearly_median, yearly_sd=yearly_sd, 
                               yearly_l95=yearly_l95, yearly_u95=yearly_u95, yearly_months_above_one_mean=yearly_months_above_one_mean, 
                               yearly_months_above_one_median=yearly_months_above_one_median, yearly_months_above_one_sd=yearly_months_above_one_sd, 
                               yearly_months_above_one_l95=yearly_months_above_one_l95, yearly_months_above_one_u95=yearly_months_above_one_u95)
  
  # load raster template
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")) #temp_rast_list
  rastEx <- temp_rast_list[[1]] # use temperature raster list as a template
  
  for (ii in seq_along(index_yearly_summary)) {
    data <- index_yearly_summary[[ii]]
    data_name <- names(index_yearly_summary)[ii]
    output <- list()
    for (jj in 1:nrow(data)) {
      rastEx[] <- matrix(data[jj, ], ncol=ori_c, nrow=ori_l, byrow=TRUE)
      output[[jj]] <- rastEx
    }
    output <- raster::brick(output)
    names(output) <- as.character(years)
    output <- terra::rast(output)
    filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", data_name, "_rasters.tif")
    eval(parse(text=paste0(index, "_", data_name, "_rasters", "<-", "output")))
    eval(parse(text=paste0("terra::writeRaster(", index, "_", data_name, "_rasters", ", ", "file=filename, overwrite=TRUE)")))
  }
  return()
}
{create_yearly_rasters("indexP")}
#{ create_yearly_rasters("Q")}
#{ create_yearly_rasters("V0")}

# compute summary statistics for each pixel for the entire period
# (mean, median, sd, 95% credible interval)
# saved as Raster objects 
create_entire_period_raster  <- function(index) {
  raw_output_files <- list.files(path=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"),
                                 pattern=paste0(region, tag, "_", index, "_month_per_cell_parallel*"))
  num_files <- length(raw_output_files)
  for (ii in 1:num_files) {
    raw_output_files[[ii]] <- paste0(region, tag, "_", index, "_month_per_cell_parallel", ii, ".Rdata")
  }
  
  mean_list <- list()
  median_list <- list()
  sd_list <- list()
  l95_list <- list()
  u95_list <- list()
  
  jj <- 1
  for (ii in seq_along(raw_output_files)) {
    load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", raw_output_files[[ii]])) #allPs, allQs or allV0s
    if (index=="indexP") data <- allPs
    else if (index=="Q") data <- allQs
    else data <- allV0s
    
    for (kk in seq_along(data)) {
      if (length(data[[kk]])==1 && is.na(data[[kk]])) {
        mean_list[[jj]] <- NA
        median_list[[jj]] <- NA
        sd_list[[jj]] <- NA
        l95_list[[jj]] <- NA
        u95_list[[jj]] <- NA
      } else {
        this_data <- as.matrix(data[[kk]][, -1])
        mean_list[[jj]] <- round(mean(this_data), 3)
        median_list[[jj]] <- round(median(this_data), 3)
        sd_list[[jj]] <- round(sd(this_data), 3)
        l95_list[[jj]] <- round(quantile(this_data, probs=0.025, names=FALSE), 3)
        u95_list[[jj]] <- round(quantile(this_data, probs=0.975, names=FALSE), 3)
      }
      jj <- jj + 1
    }
    if (index=="indexP") rm(allPs)
    else if (index=="Q") rm(allQs)
    else rm(allV0s)
    rm(data);
  }
  entire_period_mean <- unlist(mean_list)
  entire_period_median <- unlist(median_list)
  entire_period_sd <- unlist(sd_list)
  entire_period_l95 <- unlist(l95_list)
  entire_period_u95 <- unlist(u95_list)
  index_entire_period_summary <- list(entire_period_mean=entire_period_mean, entire_period_median=entire_period_median, 
                                      entire_period_sd=entire_period_sd, 
                                      entire_period_l95=entire_period_l95, entire_period_u95=entire_period_u95)
  
  # load raster template
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")) #temp_rast_list
  rastEx <- temp_rast_list[[1]] # use temperature raster list as a template
  
  for (ii in seq_along(index_entire_period_summary)) {
    data <- index_entire_period_summary[[ii]]
    data_name <- names(index_entire_period_summary)[ii]
    rastEx[] <- matrix(data, ncol=ori_c, nrow=ori_l, byrow=TRUE)
    output <- rastEx
    output <- terra::rast(output)
    filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", data_name, "_raster.tif")
    eval(parse(text=paste0(index, "_", data_name, "_raster", "<-", "output")))
    eval(parse(text=paste0("terra::writeRaster(", index, "_", data_name, "_raster", ", ", "file=filename, overwrite=TRUE)")))
  }
  return()
}
{ create_entire_period_raster("indexP")}
#{ create_entire_period_raster("Q")}
#{ create_entire_period_raster("V0")}

# generate summary statistic rasters for a typical year per pixel
# (i.e. mean typical year, median typical year, sd typical year, 95% credible interval typical year)
# (i.e. mean peak timing typical year, etc.)
# (i.e. mean trough timing typical year, etc.)
# (i.e. mean number of months above one typical year, etc.)
# saved as Raster bricks
gen_typical_year_rasters <- function(index) {
  raw_output_files <- list.files(path=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"),
                                 pattern=paste0(region, tag, "_", index, "_month_per_cell_parallel*"))
  num_files <- length(raw_output_files)
  for (ii in 1:num_files) {
    raw_output_files[[ii]] <- paste0(region, tag, "_", index, "_month_per_cell_parallel", ii, ".Rdata")
  }
  
  months_in_year <- 12
  typical_year_data_per_cell <- list()
  typical_year_peak_month_data_per_cell <- list()
  typical_year_trough_month_data_per_cell <- list()
  typical_year_months_above_one_data_per_cell <- list()
  jj <- 1
  for (ii in raw_output_files) {
    load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", ii)) #allPs, allQs or allV0s
    if (index=="indexP") data <- allPs
    else if (index=="Q") data <- allQs
    else data <- allV0s
    
    for (kk in seq_along(data)) {
      if (length(data[[kk]])==1 && is.na(data[[kk]])) {
        typical_year_data_per_cell[[jj]] <- NA
        typical_year_peak_month_data_per_cell[[jj]] <- rep(NA, 5)
        typical_year_trough_month_data_per_cell[[jj]] <- rep(NA, 5)
        typical_year_months_above_one_data_per_cell[[jj]] <- rep(NA, 5)
        jj <- jj+1
        next
      }
      
      this_data <- as.matrix(data[[kk]][, -1])
      num_years <- nrow(this_data)/months_in_year
      typical_years <- matrix(NA, nrow=months_in_year, ncol=ncol(this_data))
      for (tt in 1:months_in_year) {
        typical_years[tt, ] <- colMeans(this_data[seq(tt, nrow(this_data), months_in_year), ])
      }
      typical_year <- apply(typical_years, 1, median)
      typical_year_mean <- apply(typical_years, 1, mean)
      typical_year_sd <- apply(typical_years, 1, sd)
      typical_year_l95 <- apply(typical_years, 1, function(x) quantile(x, probs=0.025, names=FALSE))
      typical_year_u95 <- apply(typical_years, 1, function(x) quantile(x, probs=0.975, names=FALSE))
      
      typical_year_peak_months <- apply(typical_years, 2, function(x) which(x==max(x))[1])
      typical_year_peak_month_mean <- mean(typical_year_peak_months)
      typical_year_peak_month_median <- median(typical_year_peak_months)
      typical_year_peak_month_sd <- sd(typical_year_peak_months)
      typical_year_peak_month_lb <- min(typical_year_peak_months)
      typical_year_peak_month_ub <- max(typical_year_peak_months)
      
      typical_year_trough_months <- apply(typical_years, 2, function(x) which(x==min(x))[1])
      typical_year_trough_month_mean <- mean(typical_year_trough_months)
      typical_year_trough_month_median <- median(typical_year_trough_months)
      typical_year_trough_month_sd <- sd(typical_year_trough_months)
      typical_year_trough_month_lb <- min(typical_year_trough_months)
      typical_year_trough_month_ub <- max(typical_year_trough_months)
      
      typical_year_number_months_above_one <- apply(typical_years, 2, function(x) length(which(x>1)))
      typical_year_number_months_above_one_mean <- mean(typical_year_number_months_above_one)
      typical_year_number_months_above_one_median <- median(typical_year_number_months_above_one)
      typical_year_number_months_above_one_sd <- sd(typical_year_number_months_above_one)
      typical_year_number_months_above_one_l95 <- quantile(typical_year_number_months_above_one, probs=0.025, names=FALSE)
      typical_year_number_months_above_one_u95 <- quantile(typical_year_number_months_above_one, probs=0.975, names=FALSE)
      
      
      df <- data.frame(month=1:months_in_year, mean=typical_year, median=typical_year, sd=typical_year_sd, l95=typical_year_l95, u95=typical_year_u95)
      typical_year_data_per_cell[[jj]] <- df
      peak_month_data <- c(peak_month_mean=typical_year_peak_month_mean, peak_month_median=typical_year_peak_month_median, 
                           peak_month_sd=typical_year_peak_month_sd, peak_month_lb=typical_year_peak_month_lb, 
                           peak_month_ub=typical_year_peak_month_ub)
      typical_year_peak_month_data_per_cell[[jj]] <- peak_month_data
      trough_month_data <- c(trough_month_mean=typical_year_trough_month_mean, trough_month_median=typical_year_trough_month_median, 
                             trough_month_sd=typical_year_trough_month_sd, trough_month_lb=typical_year_trough_month_lb, 
                             trough_month_ub=typical_year_trough_month_ub)
      typical_year_trough_month_data_per_cell[[jj]] <- trough_month_data
      months_above_one_data <- c(months_above_one_mean=typical_year_number_months_above_one_mean,
                                 months_above_one_median=typical_year_number_months_above_one_median, 
                                 months_above_one_sd=typical_year_number_months_above_one_sd, 
                                 months_above_one_l95=typical_year_number_months_above_one_l95, 
                                 months_above_one_u95=typical_year_number_months_above_one_u95)
      typical_year_months_above_one_data_per_cell[[jj]] <- months_above_one_data
      
      jj <- jj + 1
    }
    if (index=="indexP") rm(allPs)
    else if (index=="Q") rm(allQs)
    else rm(allV0s)
    rm(data);
  }
  typical_year_peak_month_data_per_cell <- do.call(rbind, typical_year_peak_month_data_per_cell)
  typical_year_trough_month_data_per_cell <- do.call(rbind, typical_year_trough_month_data_per_cell)
  typical_year_months_above_one_data_per_cell <- do.call(rbind, typical_year_months_above_one_data_per_cell)
  
  # load raster template
  load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")) #temp_rast_list
  rastEx <- temp_rast_list[[1]] # use temperature raster list as a template
  
  output_names <- c("mean", "median", "sd", "l95", "u95", 
                    "peak_month_mean", "peak_month_median", "peak_month_sd", "peak_month_lb", 
                    "peak_month_ub", "trough_month_mean", "trough_month_median", "trough_month_sd", 
                    "trough_month_lb", "trough_month_ub", "months_above_one_mean", "months_above_one_median", 
                    "months_above_one_sd", "months_above_one_l95", "months_above_one_u95")
  for (ii in output_names) {
    if (ii %in% c("mean", "median", "sd", "l95", "u95")) {
      output <- list()
      for (t in 1:months_in_year) {
        cell_data <- rep(NA, ori_c*ori_l)
        for (jj in seq_along(typical_year_data_per_cell)) {
          if (all(is.na(typical_year_data_per_cell[[jj]]))) cell_data[jj] <- NA
          else cell_data[jj] <- typical_year_data_per_cell[[jj]][[ii]][t]
        }
        rastEx[] <- matrix(cell_data, ncol=ori_c, nrow=ori_l, byrow=TRUE)
        output[[t]] <- rastEx
      }
      output <- raster::brick(output)
      names(output) <- as.character(1:12)
    } else if (ii %in% c("peak_month_mean", "peak_month_median", "peak_month_sd", "peak_month_lb", "peak_month_ub")) {
      cell_data <- rep(NA, ori_c*ori_l)
      for (jj in 1:nrow(typical_year_peak_month_data_per_cell)) {
        cell_data[jj] <- typical_year_peak_month_data_per_cell[jj, ii]
      }
      rastEx[] <- matrix(cell_data, ncol=ori_c, nrow=ori_l, byrow=TRUE)
      output <- rastEx
    } else if (ii %in% c("trough_month_mean", "trough_month_median", "trough_month_sd", "trough_month_lb", "trough_month_ub")) {
      cell_data <- rep(NA, ori_c*ori_l)
      for (jj in 1:nrow(typical_year_trough_month_data_per_cell)) {
        cell_data[jj] <- typical_year_trough_month_data_per_cell[jj, ii]
      }
      rastEx[] <- matrix(cell_data, ncol=ori_c, nrow=ori_l, byrow=TRUE)
      output <- rastEx
    } else {
      cell_data <- rep(NA, ori_c*ori_l)
      for (jj in 1:nrow(typical_year_months_above_one_data_per_cell)) {
        cell_data[jj] <- typical_year_months_above_one_data_per_cell[jj, ii]
      }
      rastEx[] <- matrix(cell_data, ncol=ori_c, nrow=ori_l, byrow=TRUE)
      output <- rastEx
    }
    output <- terra::rast(output)
    filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", "typical_year_", ii, "_rasters.tif")
    eval(parse(text=paste0(index, "_", "typical_year_", ii, "_rasters", "<-", "output")))
    eval(parse(text=paste0("terra::writeRaster(", index, "_", "typical_year_", ii, "_rasters", ", ", "file=filename, overwrite=TRUE)")))
  }
}
{
  gen_typical_year_rasters(index="indexP")
  #gen_typical_year_rasters(index="Q")
  #gen_typical_year_rasters(index="V0")
}

# generate the summary rasters for entire period
# (i.e. mean, median, sd, l95, u95)
# gen_summary_rasters <- function(index) {
#   raw_output_files <- list.files(path=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/"),
#                                  pattern=paste0(region, tag, "_", index, "_month_per_cell_parallel*"))
#   num_files <- length(raw_output_files)
#   for (ii in 1:num_files) {
#     raw_output_files[[ii]] <- paste0(region, tag, "_", index, "_month_per_cell_parallel", ii, ".Rdata")
#   }
#   
#   summary_index_per_cell <- list()
#   jj <- 1
#   for (ii in raw_output_files) {
#     load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", ii)) #allPs, allQs or allV0s
#     if (index=="indexP") data <- allPs
#     else if (index=="Q") data <- allQs
#     else data <- allV0s
#     
#     for (kk in seq_along(data)) {
#       if (length(data[[kk]])==1 && is.na(data[[kk]])) {summary_index_per_cell[[jj]] <- NA; jj<- jj+1; next;}
#       
#       this_data <- as.matrix(data[[kk]][, -1])
#       summary_index <- colMeans(this_data)
#       
#       summary_index_median <- median(summary_index)
#       summary_index_mean <- mean(summary_index)
#       summary_index_sd <- sd(summary_index)
#       summary_index_l95 <- quantile(summary_index, probs=0.025, names=FALSE)
#       summary_index_u95 <- quantile(summary_index, probs=0.975, names=FALSE)
#       
#       summary <- c(mean=summary_index_mean, median=summary_index_median, sd=summary_index_sd, 
#                    l95=summary_index_l95, u95=summary_index_u95)
#       summary_index_per_cell[[jj]] <- summary
#       jj <- jj + 1
#     }
#     if (index=="indexP") rm(allPs)
#     else if (index=="Q") rm(allQs)
#     else rm(allV0s)
#     rm(data); 
#   }
#   
#   # load raster template
#   load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")) #temp_rast_list
#   rastEx <- temp_rast_list[[1]] # use temperature raster list as a template
#   
#   output_names <- c("mean", "median", "sd", "l95", "u95")
#   for (ii in output_names) {
#     output <- list()
#     cell_data <- rep(NA, ori_c*ori_l)
#     for (jj in seq_along(summary_index_per_cell)) {
#       if (all(is.na(summary_index_per_cell[[jj]]))) cell_data[jj] <- NA
#       else cell_data[jj] <- summary_index_per_cell[[jj]][ii]
#     }
#     rastEx[] <- matrix(cell_data, ncol=ori_c, nrow=ori_l, byrow=TRUE)
#     output <- rastEx
#     output <- terra::rast(output)
#     filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", "summary_", ii, "_raster.tif")
#     eval(parse(text=paste0(index, "_", "summary_", ii, "_raster", "<-", "output")))
#     eval(parse(text=paste0("terra::writeRaster(", index, "_", "summary_", ii, "_raster", ", ", "file=filename, overwrite=TRUE)")))
#   }
# }
# {
#   gen_summary_rasters(index="indexP")
#   #gen_summary_rasters(index="Q")
#   #gen_summary_rasters(index="V0")
# }

##################################################################
cat('## organize data for entire region... \n')

# generate data frame of monthly summary statistics (mean, median, etc.) for the entire region
get_entire_region_month_summary <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_", index, "_month_region.Rdata"))
  eval(parse(text=paste0("data", "<-", "region_", index, "_sims")))
  
  index_data <- as.matrix(dplyr::select(data, -date))
  index_mean <- rowMeans(index_data)
  index_median <- apply(index_data, MARGIN=1, median)
  index_sd <- apply(index_data, MARGIN=1, sd)
  index_l95 <- apply(index_data, MARGIN=1, FUN= function(X){ quantile(X, probs=c(0.025), na.rm=T)})
  index_u95 <- apply(index_data, MARGIN=1, FUN= function(X){ quantile(X, probs=c(0.975), na.rm=T)})
  output <- data.frame(date=data$date, mean=index_mean, median=index_median, sd=index_sd,
                       l95=index_l95, u95=index_u95)
  filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", "entire_region_monthly_summary_dataframe.Rdata")
  eval(parse(text=paste0(index, "_", "entire_region_monthly_summary_dataframe", "<-", "output")))
  eval(parse(text=paste0("save(", index, "_", "entire_region_monthly_summary_dataframe", ", ", "file=filename)")))
}

{
  get_entire_region_month_summary("indexP")
  #get_entire_region_month_summary("Q")
  #get_entire_region_month_summary("V0")
}

# generate data frame of yearly summary statistics (mean, median, etc.) for the entire region
get_entire_region_year_summary <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_", index, "_month_region.Rdata"))
  eval(parse(text=paste0("data", "<-", "region_", index, "_sims")))
  
  index_mean_list <- list()
  index_median_list <- list()
  index_sd_list <- list()
  index_l95_list <- list()
  index_u95_list <- list()
  
  years <- lubridate::year(dates) %>% unique() %>% sort()
  for (tt in years) {
    index_data <- as.matrix(data %>% mutate(year=lubridate::year(date)) %>% filter(year==tt) %>% dplyr::select(-date, -year))
    index_mean_list[[tt]] <- mean(index_data)
    index_median_list[[tt]] <- median(index_data)
    index_sd_list[[tt]] <- sd(index_data)
    index_l95_list[[tt]] <- quantile(index_data, probs=c(0.025), na.rm=T)
    index_u95_list[[tt]] <- quantile(index_data, probs=c(0.975), na.rm=T)
  }
  index_mean <- unlist(index_mean_list)
  index_median <- unlist(index_median_list)
  index_sd <- unlist(index_sd_list)
  index_l95 <- unlist(index_l95_list)
  index_u95 <- unlist(index_u95_list)
  
  output <- data.frame(year=years, mean=index_mean, median=index_median, sd=index_sd,
                       l95=index_l95, u95=index_u95)
  filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", "entire_region_yearly_summary_dataframe.Rdata")
  eval(parse(text=paste0(index, "_", "entire_region_yearly_summary_dataframe", "<-", "output")))
  eval(parse(text=paste0("save(", index, "_", "entire_region_yearly_summary_dataframe", ", ", "file=filename)")))
}

{
  get_entire_region_year_summary("indexP")
  #get_entire_region_year_summary("Q")
  #get_entire_region_year_summary("V0")
}

# generate data frame of typical year summary statistics (mean, median, etc.) for the entire region
get_entire_region_typical_year <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_", index, "_month_region.Rdata"))
  eval(parse(text=paste0("data", "<-", "region_", index, "_sims")))
  
  months_in_year <- 12
  index_data <- as.matrix(dplyr::select(data, -date))
  typical_years <- matrix(NA, nrow=months_in_year, ncol=ncol(index_data))
  for (kk in 1:months_in_year) {
    typical_years[kk, ] <- colMeans(index_data[seq(kk, nrow(index_data), months_in_year), ])
  }
  typical_year <- apply(typical_years, 1, median)
  typical_year_mean <- apply(typical_years, 1, mean)
  typical_year_sd <- apply(typical_years, 1, sd)
  typical_year_l95 <- apply(typical_years, 1, function(x) quantile(x, probs=0.025, names=FALSE))
  typical_year_u95 <- apply(typical_years, 1, function(x) quantile(x, probs=0.975, names=FALSE))
  output <- data.frame(month=1:months_in_year, mean=typical_year_mean, median=typical_year, sd=typical_year_sd, l95=typical_year_l95, u95=typical_year_u95)
  
  filename <- paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_", "entire_region_typical_year_dataframe.Rdata")
  eval(parse(text=paste0(index, "_", "entire_region_typical_year_dataframe", "<-", "output")))
  eval(parse(text=paste0("save(", index, "_", "entire_region_typical_year_dataframe", ", ", "file=filename)")))
}

{
  get_entire_region_typical_year("indexP")
  #get_entire_region_typical_year("Q")
  #get_entire_region_typical_year("V0")
}

##################################################################
cat('## deleting the raw simulation output... \n')
unlink(file.path(MVSE_OUTPUT_DATA_FOLDER, "raw_output"), recursive=TRUE)

##################################################################
step_3_time <- Sys.time() - step_3_time
cat("## step 3 time:", step_3_time, attributes(step_3_time)$units, "\n")

cat("\n")