#################################################################################################
#
# We create summary statistics for each municipality in Brazil. Summary statistics are performed 
# for each month, each year, a typical year and for the entire time frame. 
#
# The following scripts must be run and their output must be generated before this script 
# can be run:
#  1. 01_organise_population_data.R
#  2. 01_organise_adm2_sf_data.R
#  3. transform_grib_2_nc.sh
#  4. 02_generate_climate_data.R
#  5. 01_run_MVSE.R
#
# Outputs: A collection of .Rdata files that contain summary statistics for the posterior 
# distribution of Index P, as well as the monthly cases and climate data for each municipality. 
#
#################################################################################################

#################################################################################################
# Definition of required packages, directories and local data

# libraries
library(pacman)
p_load(MVSE)
p_load(tidyverse, lubridate, magrittr, pracma)

# input directories
CASE_DATA_FOLDER <- "data/bra2_cases/"
POPULATION_DATA_FOLDER <- "data/bra2_population/"
RAW_OUTPUT_DATA_FOLDER <- "output/data/raw_output/"
CLIMATE_DATA_FOLDER <- "data/bra2_climate/"

# output directories
SUMMARY_OUTPUT_DATA_FOLDER <- "output/data/summary_output/"

# regions (state, city, city_id)
regions <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")) %>% 
  select(state, city, city_id) %>% 
  filter(!grepl(pattern="Município ignorado", x=city, fixed=TRUE)) %>% 
  unique()

#################################################################################################
# compute summary statistics for the climatic variables (i.e. temperature, humidity, precipitation)
# for each municipality

summarise_climate_region <- function(regions, metric) {
  metrics <- c("temperature", "humidity", "precipitation")
  
  # summary for typical year 
  typical_year_list <- list()
  
  # summary for entire period (mean, median, sd, 95% CI)
  mean_list <- list()
  median_list <- list()
  sd_list <- list()
  l95_list <- list()
  u95_list <- list()
  
  # summary AUC of index P per year, typical year and entire period
  AUC_typical_year_list <- list()
  AUC_list <- list()
  
  # compute the basic summary statistics for each region
  cat("## Computing some basic summary statistics in each region...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    city_id <- regions[ii, "city_id"]
    
    # retrieve climate data for region 
    climate_data_file <- paste0("BR_mun_", city_id, "_climate_entire_region.csv")
    if (file.size(paste0(CLIMATE_DATA_FOLDER, climate_data_file))<8000) {
      typical_year_list[[ii]] <- NA
      
      mean_list[[ii]] <- NA
      median_list[[ii]] <- NA
      sd_list[[ii]] <- NA
      l95_list[[ii]] <- NA
      u95_list[[ii]] <- NA
      
      AUC_typical_year_list[[ii]] <- NA
      AUC_list[[ii]] <- NA
      next
    }
    climate_data <- read.csv(file=paste0(CLIMATE_DATA_FOLDER, climate_data_file)) %>% 
      arrange(desc(-year), desc(-month))
    
    metric_colnames <- c("T", "H", "R")
    names(metric_colnames) <- metrics
    data <- climate_data[[metric_colnames[metric]]]
    months_in_year <- 12
    num_years <- nrow(data)/months_in_year
    
    # summary typical year index P
    typical_year <- rep(NA, months_in_year)
    for (tt in 1:months_in_year) {
      typical_year[tt] <- mean(data[seq(tt, length(data), months_in_year)])
    }
    unadj_typical_year <- typical_year
    lowess_adj_typical_year <- lowess(x=c(1:months_in_year, months_in_year+1:months_in_year), 
                                      y=c(unadj_typical_year, unadj_typical_year), f=1/4)
    lowess_adj_typical_year <- lowess_adj_typical_year$y[c(months_in_year+1:6, 7:months_in_year)]
    typical_year <- lowess_adj_typical_year
    typical_year_list[[ii]] <- typical_year
    
    # summaries for entire time frame
    mean_list[[ii]] <- round(mean(data), 3)
    median_list[[ii]] <- round(median(data), 3)
    sd_list[[ii]] <- round(sd(data), 3)
    l95_list[[ii]] <- round(quantile(data, probs=0.025, names=FALSE), 3)
    u95_list[[ii]] <- round(quantile(data, probs=0.975, names=FALSE), 3)
    
    # summaries for AUC
    AUC_typical_year_list[[ii]] <- trapz(x=1:months_in_year, y=typical_year)
    AUC_list[[ii]] <- trapz(x=1:length(data), y=data)
    
    setTxtProgressBar(pb, ii)
  }
  months <- 1:12
  typical_year <- cbind(months, do.call(cbind, typical_year_list))
  colnames(typical_year) <- c("month", pull(regions, city_id))
  
  entire_tf_colnames <- pull(regions, city_id)
  entire_mean <- unlist(mean_list)
  entire_median <- unlist(median_list)
  entire_sd <- unlist(sd_list)
  entire_l95 <- unlist(l95_list)
  entire_u95 <- unlist(u95_list)
  names(entire_mean) <- entire_tf_colnames
  names(entire_median) <- entire_tf_colnames
  names(entire_sd) <- entire_tf_colnames
  names(entire_l95) <- entire_tf_colnames
  names(entire_u95) <- entire_tf_colnames
  
  AUC_typical_year <- unlist(AUC_typical_year_list)
  names(AUC_typical_year) <- entire_tf_colnames
  entire_AUC <- unlist(AUC_list)
  names(entire_AUC) <- entire_tf_colnames
  
  region_basic_summary_data <- list(typical_year=typical_year, 
                                    entire_mean=entire_mean, entire_median=entire_median, entire_sd=entire_sd, 
                                    entire_u95=entire_u95, entire_l95=entire_l95, 
                                    AUC_typical_year=AUC_typical_year, entire_AUC=entire_AUC)
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_basic_", metric, "_summary_data.Rdata")
  save(region_basic_summary_data, file=fileout)
}
{ summarise_climate_region(regions, metric="temperature")
  summarise_climate_region(regions, metric="humidity")
  summarise_climate_region(regions, metric="precipitation") }

#################################################################################################
# compute summary statistics for the posterior distribution of Index P estimated for each 
# municipality. 

summarise_region <- function(regions, years=2000:2014) {
  
  # summary index P for each month (mean, median, sd, 95% CI)
  mean_month_list <- list()
  median_month_list <- list()
  sd_month_list <- list()
  l95_month_list <- list()
  u95_month_list <- list()
  
  # summary cumulative index P for each year  (mean, median, sd, 95% CI)
  mean_year_list <- list()
  median_year_list <- list()
  sd_year_list <- list()
  l95_year_list <- list()
  u95_year_list <- list()
  
  # summary typical year index P (mean, median, sd, 95% CI)
  mean_typical_year_list <- list()
  median_typical_year_list <- list()
  sd_typical_year_list <- list()
  l95_typical_year_list <- list()
  u95_typical_year_list <- list()
  
  # summary normalised typical year index P (mean, median, sd, 90% CI, 75%CI)
  mean_normalised_typical_year_list <- list()
  median_normalised_typical_year_list <- list()
  sd_normalised_typical_year_list <- list()
  l90_normalised_typical_year_list <- list()
  u90_normalised_typical_year_list <- list()
  l75_normalised_typical_year_list <- list()
  u75_normalised_typical_year_list <- list()
  
  # summary index P for entire period (mean, median, sd, 95%CI)
  mean_list <- list()
  median_list <- list()
  sd_list <- list()
  l95_list <- list()
  u95_list <- list()
  
  # summary AUC of index P per year, typical year and entire period
  # (mean, median, sd, 95%CI)
  mean_AUC_per_year_list <- list()
  median_AUC_per_year_list <- list()
  sd_AUC_per_year_list <- list()
  l95_AUC_per_year_list <- list()
  u95_AUC_per_year_list <- list()
  mean_AUC_typical_year_list <- list()
  median_AUC_typical_year_list <- list()
  mean_AUC_list <- list()
  median_AUC_list <- list()
  sd_AUC_list <- list()
  
  # compute basic summary statistics for Index P
  cat("## Computing some basic summary statistics for index P in each region...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", city_id, "_SIMP.Rdata"))) {
      mean_month_list[[ii]] <- NA
      median_month_list[[ii]] <- NA
      sd_month_list[[ii]] <- NA
      l95_month_list[[ii]] <- NA
      u95_month_list[[ii]] <- NA
      
      mean_year_list[[ii]] <- NA
      median_year_list[[ii]] <- NA
      sd_year_list[[ii]] <- NA
      l95_year_list[[ii]] <- NA
      u95_year_list[[ii]] <- NA
      
      median_typical_year_list[[ii]] <- NA
      mean_typical_year_list[[ii]] <- NA
      sd_typical_year_list[[ii]] <- NA
      l95_typical_year_list[[ii]] <- NA
      u95_typical_year_list[[ii]] <- NA
      
      median_normalised_typical_year_list[[ii]] <- NA
      mean_normalised_typical_year_list[[ii]] <- NA
      sd_normalised_typical_year_list[[ii]] <- NA
      l90_normalised_typical_year_list[[ii]] <- NA
      u90_normalised_typical_year_list[[ii]] <- NA
      l75_normalised_typical_year_list[[ii]] <- NA
      u75_normalised_typical_year_list[[ii]] <- NA
      
      mean_list[[ii]] <- NA
      median_list[[ii]] <- NA
      sd_list[[ii]] <- NA
      l95_list[[ii]] <- NA
      u95_list[[ii]] <- NA
      
      mean_AUC_per_year_list[[ii]] <- NA
      median_AUC_per_year_list[[ii]] <- NA
      sd_AUC_per_year_list[[ii]] <- NA
      l95_AUC_per_year_list[[ii]] <- NA
      u95_AUC_per_year_list[[ii]] <- NA
      mean_AUC_typical_year_list[[ii]] <- NA
      median_AUC_typical_year_list[[ii]] <- NA
      mean_AUC_list[[ii]] <- NA
      median_AUC_list[[ii]] <- NA
      sd_AUC_list[[ii]] <- NA
      next
    }
    
    load(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", city_id, "_SIMP.Rdata")) #sim_p
    sim_p %<>% 
      mutate(year=lubridate::year(as.Date(date))) %>% 
      filter(year %in% years)
    data <- as.matrix(sim_p[, -1])
    months_in_year <- 12
    num_years <- nrow(data)/months_in_year
    
    # summary index P for each month
    mean_month_list[[ii]] <- round(rowMeans(data), 3)
    median_month_list[[ii]] <- round(apply(data, 1, median), 3)
    sd_month_list[[ii]] <- round(apply(data, 1, sd), 3)
    l95_month_list[[ii]] <- round(apply(data, 1, function(x) quantile(x, probs=0.025, names=FALSE)), 3)
    u95_month_list[[ii]] <- round(apply(data, 1, function(x) quantile(x, probs=0.975, names=FALSE)), 3)
    
    # summary index P for each year 
    cum_per_year <- matrix(NA, nrow=num_years, ncol=ncol(data))
    for (tt in 1:num_years) {
      cum_per_year[tt, ] <- colSums(data[seq(((tt-1)*months_in_year+1), (tt*months_in_year)), ])
    }
    mean_year_list[[ii]] <- round(rowMeans(cum_per_year), 3)
    median_year_list[[ii]] <- round(apply(cum_per_year, 1, median), 3)
    sd_year_list[[ii]] <- round(apply(cum_per_year, 1, sd), 3)
    l95_year_list[[ii]] <- round(apply(cum_per_year, 1, function(x) quantile(x, probs=0.025, names=FALSE)), 3)
    u95_year_list[[ii]] <- round(apply(cum_per_year, 1, function(x) quantile(x, probs=0.975, names=FALSE)), 3)
    
    # summary typical year index P
    typical_years <- matrix(NA, nrow=months_in_year, ncol=ncol(data))
    for (tt in 1:months_in_year) {
      typical_years[tt, ] <- colMeans(data[seq(tt, nrow(data), months_in_year), ])
    }
    for (jj in 1:ncol(typical_years)) {
      unadj_typical_year <- typical_years[, jj]
      lowess_adj_typical_year <- lowess(x=c(1:months_in_year, months_in_year+1:months_in_year), 
                                        y=c(unadj_typical_year, unadj_typical_year), f=1/4)
      lowess_adj_typical_year <- lowess_adj_typical_year$y[c(months_in_year+1:6, 7:months_in_year)]
      typical_years[, jj] <- lowess_adj_typical_year
    }
    mean_typical_year_list[[ii]] <- apply(typical_years, 1, mean)
    median_typical_year_list[[ii]] <- apply(typical_years, 1, median)
    sd_typical_year_list[[ii]] <- apply(typical_years, 1, sd)
    l95_typical_year_list[[ii]] <- apply(typical_years, 1, function(x) quantile(x, probs=0.025, names=FALSE))
    u95_typical_year_list[[ii]] <- apply(typical_years, 1, function(x) quantile(x, probs=0.975, names=FALSE))
    
    # summary of normalised typical year index P
    normalised_typical_years <- apply(typical_years, 2, function(x) (x-min(x))/(max(x)-min(x)))
    mean_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, mean)
    median_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, median)
    sd_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, sd)
    l90_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, function(x) quantile(x, probs=0.05, names=FALSE))
    u90_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, function(x) quantile(x, probs=0.95, names=FALSE))
    l75_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, function(x) quantile(x, probs=0.125, names=FALSE))
    u75_normalised_typical_year_list[[ii]] <- apply(normalised_typical_years, 1, function(x) quantile(x, probs=0.875, names=FALSE))
    
    # summaries for entire time frame
    cum_entire_period <- colSums(data)
    mean_list[[ii]] <- round(mean(cum_entire_period), 3)
    median_list[[ii]] <- round(median(cum_entire_period), 3)
    sd_list[[ii]] <- round(sd(cum_entire_period), 3)
    l95_list[[ii]] <- round(quantile(cum_entire_period, probs=0.025, names=FALSE), 3)
    u95_list[[ii]] <- round(quantile(cum_entire_period, probs=0.975, names=FALSE), 3)
    
    # AUC summaries
    AUC_per_year <- matrix(NA, nrow=num_years, ncol=ncol(data)) 
    tmp_data <- rbind(data[1, ], data)
    for (tt in 1:num_years) {
      AUC_per_year[tt, ] <- apply(tmp_data[((tt-1)*months_in_year+1):(tt*months_in_year+1), ], 2, 
                                  function(z) trapz(x=0:months_in_year, y=z)/12)
    }
    mean_AUC_per_year_list[[ii]] <- apply(AUC_per_year, 1, mean)
    median_AUC_per_year_list[[ii]] <- apply(AUC_per_year, 1, median)
    sd_AUC_per_year_list[[ii]] <- apply(AUC_per_year, 1, sd)
    l95_AUC_per_year_list[[ii]] <- apply(AUC_per_year, 1, function(x) quantile(x, probs=0.025, names=FALSE))
    u95_AUC_per_year_list[[ii]] <- apply(AUC_per_year, 1, function(x) quantile(x, probs=0.975, names=FALSE))
    
    tmp_typical_years <- rbind(typical_years[1, ], typical_years)
    AUC_typical_years <- apply(tmp_typical_years, 2, function(z) trapz(x=0:months_in_year, y=z)/12)
    mean_AUC_typical_year_list[[ii]] <- mean(AUC_typical_years)
    median_AUC_typical_year_list[[ii]] <- median(AUC_typical_years)
    
    cum_AUCs <- apply(tmp_data, 2, function(z) trapz(x=0:nrow(data), y=z)/12)
    mean_AUC_list[[ii]] <- mean(cum_AUCs)
    median_AUC_list[[ii]] <- median(cum_AUCs)
    sd_AUC_list[[ii]] <- sd(cum_AUCs)
    
    setTxtProgressBar(pb, ii)
  }
  dates <- sim_p[, 1]
  monthly_colnames <- c("date", pull(regions, city_id))
  monthly_mean <- cbind(dates, do.call(cbind, mean_month_list))
  monthly_median <- cbind(dates, do.call(cbind, median_month_list))
  monthly_sd <- cbind(dates, do.call(cbind, sd_month_list))
  monthly_l95 <- cbind(dates, do.call(cbind, l95_month_list))
  monthly_u95 <- cbind(dates, do.call(cbind, u95_month_list))
  colnames(monthly_mean) <- monthly_colnames
  colnames(monthly_median) <- monthly_colnames
  colnames(monthly_sd) <- monthly_colnames
  colnames(monthly_l95) <- monthly_colnames
  colnames(monthly_u95) <- monthly_colnames
  
  years <- sort(unique(lubridate::year(dates)))
  yearly_colnames <- c("year", pull(regions, city_id))
  yearly_mean <- cbind(years, do.call(cbind, mean_year_list))
  yearly_median <- cbind(years, do.call(cbind, median_year_list))
  yearly_sd <- cbind(years, do.call(cbind, sd_year_list))
  yearly_l95 <- cbind(years, do.call(cbind, l95_year_list))
  yearly_u95 <- cbind(years, do.call(cbind, u95_year_list))
  colnames(yearly_mean) <- yearly_colnames
  colnames(yearly_median) <- yearly_colnames
  colnames(yearly_sd) <- yearly_colnames
  colnames(yearly_l95) <- yearly_colnames
  colnames(yearly_u95) <- yearly_colnames
  
  months <- 1:12
  typical_year_colnames <- c("month", pull(regions, city_id))
  typical_year_mean <- cbind(months, do.call(cbind, mean_typical_year_list))
  typical_year_median <- cbind(months, do.call(cbind, median_typical_year_list))
  typical_year_sd <- cbind(months, do.call(cbind, sd_typical_year_list))
  typical_year_l95 <- cbind(months, do.call(cbind, l95_typical_year_list))
  typical_year_u95 <- cbind(months, do.call(cbind, u95_typical_year_list))
  colnames(typical_year_mean) <- typical_year_colnames
  colnames(typical_year_median) <- typical_year_colnames
  colnames(typical_year_sd) <- typical_year_colnames
  colnames(typical_year_l95) <- typical_year_colnames
  colnames(typical_year_u95) <- typical_year_colnames
  
  normalised_typical_year_mean <- cbind(months, do.call(cbind, mean_normalised_typical_year_list))
  normalised_typical_year_median <- cbind(months, do.call(cbind, median_normalised_typical_year_list))
  normalised_typical_year_sd <- cbind(months, do.call(cbind, sd_normalised_typical_year_list))
  normalised_typical_year_l90 <- cbind(months, do.call(cbind, l90_normalised_typical_year_list))
  normalised_typical_year_u90 <- cbind(months, do.call(cbind, u90_normalised_typical_year_list))
  normalised_typical_year_l75 <- cbind(months, do.call(cbind, l75_normalised_typical_year_list))
  normalised_typical_year_u75 <- cbind(months, do.call(cbind, u75_normalised_typical_year_list))
  colnames(normalised_typical_year_mean) <- typical_year_colnames
  colnames(normalised_typical_year_median) <- typical_year_colnames
  colnames(normalised_typical_year_sd) <- typical_year_colnames
  colnames(normalised_typical_year_l90) <- typical_year_colnames
  colnames(normalised_typical_year_u90) <- typical_year_colnames
  colnames(normalised_typical_year_l75) <- typical_year_colnames
  colnames(normalised_typical_year_u75) <- typical_year_colnames
  
  entire_tf_colnames <- pull(regions, city_id)
  entire_mean <- do.call(cbind, mean_list)
  entire_median <- do.call(cbind, median_list)
  entire_sd <- do.call(cbind, sd_list)
  entire_l95 <- do.call(cbind, l95_list)
  entire_u95 <- do.call(cbind, u95_list)
  colnames(entire_mean) <- entire_tf_colnames
  colnames(entire_median) <- entire_tf_colnames
  colnames(entire_sd) <- entire_tf_colnames
  colnames(entire_l95) <- entire_tf_colnames
  colnames(entire_u95) <- entire_tf_colnames
  
  AUC_per_year_colnames <- c("year", pull(regions, city_id))
  years <- sort(unique(lubridate::year(dates)))
  mean_AUC_per_year <- cbind(years, do.call(cbind, mean_AUC_per_year_list))
  median_AUC_per_year <- cbind(years, do.call(cbind, median_AUC_per_year_list))
  sd_AUC_per_year <- cbind(years, do.call(cbind, sd_AUC_per_year_list))
  l95_AUC_per_year <- cbind(years, do.call(cbind, l95_AUC_per_year_list))
  u95_AUC_per_year <- cbind(years, do.call(cbind, u95_AUC_per_year_list))
  colnames(mean_AUC_per_year) <- AUC_per_year_colnames
  colnames(median_AUC_per_year) <- AUC_per_year_colnames
  colnames(sd_AUC_per_year) <- AUC_per_year_colnames
  colnames(l95_AUC_per_year) <- AUC_per_year_colnames
  colnames(u95_AUC_per_year) <- AUC_per_year_colnames
  
  mean_AUC_typical_year <- do.call(cbind, mean_AUC_typical_year_list)
  median_AUC_typical_year <- do.call(cbind, median_AUC_typical_year_list)
  colnames(mean_AUC_typical_year) <- entire_tf_colnames
  colnames(median_AUC_typical_year) <- entire_tf_colnames
  
  mean_AUC <- do.call(cbind, mean_AUC_list)
  median_AUC <- do.call(cbind, median_AUC_list)
  sd_AUC <- do.call(cbind, sd_AUC_list)
  colnames(mean_AUC) <- entire_tf_colnames
  colnames(median_AUC) <- entire_tf_colnames
  colnames(sd_AUC) <- entire_tf_colnames
  
  region_basic_summary_data <- list(monthly_mean=monthly_mean, monthly_median=monthly_median,
                                    monthly_sd=monthly_sd, monthly_l95=monthly_l95,
                                    monthly_u95=monthly_u95, yearly_mean=yearly_mean, 
                                    yearly_median=yearly_median, yearly_sd=yearly_sd, 
                                    yearly_l95=yearly_l95, yearly_u95=yearly_u95,
                                    typical_year_mean=typical_year_mean, 
                                    typical_year_median=typical_year_median, typical_year_sd=typical_year_sd, 
                                    typical_year_l95=typical_year_l95, typical_year_u95=typical_year_u95, 
                                    normalised_typical_year_mean=normalised_typical_year_mean, 
                                    normalised_typical_year_median=normalised_typical_year_median, normalised_typical_year_sd=normalised_typical_year_sd, 
                                    normalised_typical_year_l90=normalised_typical_year_l90, normalised_typical_year_u90=normalised_typical_year_u90,
                                    normalised_typical_year_l75=normalised_typical_year_l75, normalised_typical_year_u75=normalised_typical_year_u75,
                                    entire_mean=entire_mean, entire_median=entire_median, entire_sd=entire_sd, 
                                    entire_u95=entire_u95, entire_l95=entire_l95, 
                                    mean_AUC_per_year=mean_AUC_per_year, median_AUC_per_year=median_AUC_per_year, 
                                    sd_AUC_per_year=sd_AUC_per_year, l95_AUC_per_year=l95_AUC_per_year, 
                                    u95_AUC_per_year=u95_AUC_per_year, mean_AUC_typical_year=mean_AUC_typical_year, 
                                    median_AUC_typical_year=median_AUC_typical_year, mean_AUC=mean_AUC, 
                                    median_AUC=median_AUC, sd_AUC=sd_AUC)
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_basic_summary_data.Rdata")
  save(region_basic_summary_data, file=fileout)
}
{ summarise_region(regions, years=2000:2014) }

#################################################################################################
# compute summary statisitcs on the timing of the peak and trough of Index P each year

compute_yearly_peak_trough <- function(regions, years=2000:2014) {
  peak_list <- list()
  trough_list <- list()
  
  MONTHS_IN_YEAR <- 12
  cat("## Computing the timing of the yearly peak and trough for index P in each region...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    state <- regions[ii, "state"]
    city <- regions[ii, "city"]
    city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", city_id, "_SIMP.Rdata"))) {
      counts <- rep(NA, MONTHS_IN_YEAR)
      names(counts) <- 1:MONTHS_IN_YEAR
      peak_list[[ii]] <- c(city_id=city_id, state=state, city=city, stat="peak", year=NA, counts)
      trough_list[[ii]] <- c(city_id=city_id, state=state, city=city, stat="trough", year=NA, counts)
      next
    }
    
    # organise Index P data
    load(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", city_id, "_SIMP.Rdata")) #sim_p
    sim_p %<>% 
      mutate(year=lubridate::year(as.Date(date))) %>% 
      filter(year %in% years)
    dates <- sim_p[, 1]
    sim_p <- mutate(sim_p, month=lubridate::month(date)) %>% 
      select(-date) %>% select(year, month, dplyr::everything())
    
    # compute the timing of peak and trough timing over the 15 years across all samples
    region_peak_list <- list()
    region_trough_list <- list()
    for (jj in seq_along(years)) {
      data <- sim_p %>% 
        filter(year==years[jj]) %>%
        gather(key="sim", value="value", 3:ncol(.)) %>%
        mutate(value=as.numeric(value)) %>%
        group_by(year, sim) %>% 
        summarise(peak_month=month[which(value==max(value))], 
                  trough_month=month[which(value==min(value))], .groups="keep") %>% 
        mutate(peak_month=factor(peak_month, levels=as.character(1:MONTHS_IN_YEAR))) %>%
        mutate(trough_month=factor(trough_month, levels=as.character(1:MONTHS_IN_YEAR)))
      
      peak_counts <- as.numeric(table(data$peak_month))
      names(peak_counts) <- 1:MONTHS_IN_YEAR
      trough_counts <- as.numeric(table(data$trough_month))
      names(trough_counts) <- 1:MONTHS_IN_YEAR
      region_peak_list[[jj]] <- c(city_id=city_id, state=state, city=city, year=years[jj], stat="peak", peak_counts)
      region_trough_list[[jj]] <- c(city_id=city_id, state=state, city=city, year=years[jj], stat="trough", trough_counts)
    }
    peak_list[[ii]] <- do.call(rbind, region_peak_list)
    trough_list[[ii]] <- do.call(rbind, region_trough_list)
    setTxtProgressBar(pb, ii)
  }
  region_yearly_peak_data <- do.call(rbind, peak_list) %>% 
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), year=as.numeric(year), 
           `1`=as.numeric(`1`),`2`=as.numeric(`2`), `3`=as.numeric(`3`), `4`=as.numeric(`4`), `5`=as.numeric(`5`),
           `6`=as.numeric(`6`), `7`=as.numeric(`7`), `8`=as.numeric(`8`), `9`=as.numeric(`9`), `10`=as.numeric(`10`), 
           `11`=as.numeric(`11`), `12`=as.numeric(`12`))
  region_yearly_trough_data <- do.call(rbind, trough_list) %>%
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), year=as.numeric(year), 
           `1`=as.numeric(`1`),`2`=as.numeric(`2`), `3`=as.numeric(`3`), `4`=as.numeric(`4`), `5`=as.numeric(`5`),
           `6`=as.numeric(`6`), `7`=as.numeric(`7`), `8`=as.numeric(`8`), `9`=as.numeric(`9`), `10`=as.numeric(`10`), 
           `11`=as.numeric(`11`), `12`=as.numeric(`12`))
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_yearly_peak_data.Rdata")
  save(region_yearly_peak_data, file=fileout)
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_yearly_trough_data.Rdata")
  save(region_yearly_trough_data, file=fileout)
}
{compute_yearly_peak_trough(regions, years=2000:2014)}

#################################################################################################
# compute summary statisitcs on the timing of the peak and trough of Index P during a typical year

compute_typical_year_peak_trough <- function(regions, years=2000:2014) {
  peak_list <- list()
  trough_list <- list()
  
  MONTHS_IN_YEAR <- 12
  cat("## Computing the timing of peak and trough for index P during a typical year in each region...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    state <- regions[ii, "state"]
    city <- regions[ii, "city"]
    city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", city_id, "_SIMP.Rdata"))) {
      counts <- rep(NA, MONTHS_IN_YEAR)
      names(counts) <- 1:MONTHS_IN_YEAR
      peak_list[[ii]] <- c(city_id=city_id, state=state, city=city, stat="peak", counts)
      trough_list[[ii]] <- c(city_id=city_id, state=state, city=city, stat="trough", counts)
      next
    }
    
    # compute typical year for each simulation
    load(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", city_id, "_SIMP.Rdata")) #sim_p
    sim_p %<>% 
      mutate(year=lubridate::year(as.Date(date))) %>% 
      filter(year %in% years)
    data <- as.matrix(sim_p[, -1])
    num_years <- nrow(data)/MONTHS_IN_YEAR
    typical_years <- matrix(NA, nrow=MONTHS_IN_YEAR, ncol=ncol(data))
    for (tt in 1:MONTHS_IN_YEAR) {
      typical_years[tt, ] <- colMeans(data[seq(tt, nrow(data), MONTHS_IN_YEAR), ])
    }
    for (jj in 1:ncol(typical_years)) {
      unadj_typical_year <- typical_years[, jj]
      lowess_adj_typical_year <- lowess(x=c(1:MONTHS_IN_YEAR, MONTHS_IN_YEAR+1:MONTHS_IN_YEAR), 
                                        y=c(unadj_typical_year, unadj_typical_year), f=1/4)
      lowess_adj_typical_year <- lowess_adj_typical_year$y[c(MONTHS_IN_YEAR+1:6, 7:MONTHS_IN_YEAR)]
      typical_years[, jj] <- lowess_adj_typical_year
    }
    
    # compute the peak and trough timing for each typical year
    data <- as.data.frame(typical_years) %>% 
      setNames(1:ncol(data)) %>% 
      mutate(month=1:MONTHS_IN_YEAR) %>% 
      select(month, dplyr::everything()) %>%
      gather(key="sim", value="value", 2:ncol(.)) %>% 
      group_by(sim) %>% 
      summarise(peak_month=month[which(value==max(value))], 
                trough_month=month[which(value==min(value))], .groups="keep") %>% 
      mutate(peak_month=factor(peak_month, levels=as.character(1:MONTHS_IN_YEAR))) %>%
      mutate(trough_month=factor(trough_month, levels=as.character(1:MONTHS_IN_YEAR)))
    peak_counts <- as.numeric(table(data$peak_month))
    names(peak_counts) <- 1:MONTHS_IN_YEAR
    trough_counts <- as.numeric(table(data$trough_month))
    names(trough_counts) <- 1:MONTHS_IN_YEAR
    peak_list[[ii]] <- c(city_id=city_id, state=state, city=city, stat="peak", peak_counts)
    trough_list[[ii]] <- c(city_id=city_id, state=state, city=city, stat="trough", trough_counts)
    setTxtProgressBar(pb, ii)
  }
  region_typical_year_peak_data <- do.call(rbind, peak_list) %>% 
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id),
           `1`=as.numeric(`1`),`2`=as.numeric(`2`), `3`=as.numeric(`3`), `4`=as.numeric(`4`), `5`=as.numeric(`5`),
           `6`=as.numeric(`6`), `7`=as.numeric(`7`), `8`=as.numeric(`8`), `9`=as.numeric(`9`), `10`=as.numeric(`10`), 
           `11`=as.numeric(`11`), `12`=as.numeric(`12`))
  region_typical_year_trough_data <- do.call(rbind, trough_list) %>%
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), 
           `1`=as.numeric(`1`),`2`=as.numeric(`2`), `3`=as.numeric(`3`), `4`=as.numeric(`4`), `5`=as.numeric(`5`),
           `6`=as.numeric(`6`), `7`=as.numeric(`7`), `8`=as.numeric(`8`), `9`=as.numeric(`9`), `10`=as.numeric(`10`), 
           `11`=as.numeric(`11`), `12`=as.numeric(`12`))
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_typical_year_peak_data.Rdata")
  save(region_typical_year_peak_data, file=fileout)
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_typical_year_trough_data.Rdata")
  save(region_typical_year_trough_data, file=fileout)
}
{ compute_typical_year_peak_trough(regions, years=2000:2014) }

#################################################################################################
# compute the number of cases during a typical year for each municipality. 

compute_typical_year_cases <- function(regions, years=2000:2014, case_metric="abs_cases", logscale=FALSE) {
  # read in the case data and population data
  case_data <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")) %>%
    filter(year %in% years)
  if (case_metric=="incidence")
    pop_data <- readRDS(paste0(POPULATION_DATA_FOLDER, "bra2_population_data_2000_2014.rds")) %>% 
      filter(year %in% years)
  
  typical_year_list <- list()
  MONTHS_IN_YEAR <- 12
  cat("## Computing the typical year for cases in each region...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    this_state <- regions[ii, "state"]
    this_city <- regions[ii, "city"]
    this_city_id <- regions[ii, "city_id"]
    
    region_case_data <- case_data %>%
      filter(state==this_state, city==this_city) %>% 
      arrange(-desc(year), -desc(month))
    if (case_metric=="incidence") {
      region_case_data %<>% 
        left_join(pop_data, by=c("state", "city", "year")) %>%
        mutate(cases=cases/n*10^5)
    }
    if (logscale) region_case_data %<>% mutate(cases=log(cases+1))
    
    # only compute for regions with complete case data
    if (any(is.na(region_case_data$cases))) {
      typical_year <- rep(NA, MONTHS_IN_YEAR)
    } else {
      region_case_data <- as.matrix(region_case_data[, "cases"])
      typical_year <- rep(NA, MONTHS_IN_YEAR)
      for (tt in 1:MONTHS_IN_YEAR) {
        typical_year[tt] <- mean(region_case_data[seq(tt, length(region_case_data), MONTHS_IN_YEAR)], na.rm=FALSE)
      }
      unadj_typical_year <- typical_year
      lowess_adj_typical_year <- lowess(x=c(1:MONTHS_IN_YEAR, MONTHS_IN_YEAR+1:MONTHS_IN_YEAR), 
                                        y=c(unadj_typical_year, unadj_typical_year), f=1/4)
      lowess_adj_typical_year <- lowess_adj_typical_year$y[c(MONTHS_IN_YEAR+1:6, 7:MONTHS_IN_YEAR)]
      typical_year <- lowess_adj_typical_year
    }
    typical_year_list[[ii]] <- data.frame(city_id=this_city_id, state=this_state, city=this_city, 
                                          month=1:MONTHS_IN_YEAR, cases=typical_year)
    setTxtProgressBar(pb, ii)
  }
  region_typical_year_cases_data <- do.call(rbind, typical_year_list) %>%
    as.data.frame() %>%
    mutate(month=as.numeric(month), cases=as.numeric(cases), city_id=as.numeric(city_id))
  case_metric_name <- case_metric
  if (logscale) case_metric_name <- paste0("log_", case_metric)
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_typical_year_", case_metric_name, "_data.Rdata")
  save(region_typical_year_cases_data, file=fileout)
  return()
}
{ compute_typical_year_cases(regions, years=2000:2014, case_metric="abs_cases")
  compute_typical_year_cases(regions, years=2000:2014, case_metric="abs_cases", logscale=TRUE)
  compute_typical_year_cases(regions, years=2000:2014, case_metric="incidence")
  compute_typical_year_cases(regions, years=2000:2014, case_metric="incidence", logscale=TRUE) }

######################################################################################################
# compute the typical year correlations between the case data and index P for varying 
# phase shifts in Index P

compute_typical_year_corr <- function(regions, years=2000:2014, case_metric="abs_cases", logscale=FALSE) {
  # read in the case data and population data
  case_data <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")) %>%
    filter(year %in% years)
  if (case_metric=="incidence")
    pop_data <- readRDS(paste0(POPULATION_DATA_FOLDER, "bra2_population_data_2000_2014.rds")) %>% 
      filter(year %in% years)
  
  # compute correlation summary statistics for the typical year data under different shifts
  corr_data_list <- list()
  MONTHS_IN_YEAR <- 12
  cat("## Computing the correlation between case data and index P for a typical year under different shifts...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    this_state <- regions[ii, "state"]
    this_city <- regions[ii, "city"]
    this_city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata"))) {
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, shift=NA,
                                mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
      next
    }
    
    # compute typical year for case data
    region_case_data <- case_data %>%
      filter(state==this_state, city==this_city) %>% 
      arrange(-desc(year), -desc(month))
    if (case_metric=="incidence") {
      region_case_data %<>% 
        left_join(pop_data, by=c("state", "city", "year")) %>%
        unique() %>%
        mutate(cases=cases/n*10^5)
    }
    if (logscale) region_case_data %<>% mutate(cases=log(cases+1))
    
    # only compute for regions with complete case data
    if (any(is.na(region_case_data$cases))) { 
      typical_year_cases <- rep(NA, MONTHS_IN_YEAR)
    } else {
      region_case_data <- as.matrix(region_case_data[, "cases"])
      typical_year_cases <- rep(NA, MONTHS_IN_YEAR)
      for (tt in 1:MONTHS_IN_YEAR) {
        typical_year_cases[tt] <- mean(region_case_data[seq(tt, length(region_case_data), MONTHS_IN_YEAR)], na.rm=FALSE)
      }
      unadj_typical_year <- typical_year_cases
      lowess_adj_typical_year <- lowess(x=c(1:MONTHS_IN_YEAR, MONTHS_IN_YEAR+1:MONTHS_IN_YEAR), 
                                        y=c(unadj_typical_year, unadj_typical_year), f=1/4)
      lowess_adj_typical_year <- lowess_adj_typical_year$y[c(MONTHS_IN_YEAR+1:6, 7:MONTHS_IN_YEAR)]
      typical_year_cases <- lowess_adj_typical_year
    }
    
    # compute typical year for each simulation
    load(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata")) #sim_p
    sim_p %<>% 
      mutate(year=lubridate::year(as.Date(date))) %>% 
      filter(year %in% years)
    data <- as.matrix(sim_p[, -1])
    num_years <- nrow(data)/MONTHS_IN_YEAR
    typical_year_indexP <- matrix(NA, nrow=MONTHS_IN_YEAR, ncol=ncol(data))
    for (tt in 1:MONTHS_IN_YEAR) {
      typical_year_indexP[tt, ] <- colMeans(data[seq(tt, nrow(data), MONTHS_IN_YEAR), ])
    }
    for (jj in 1:ncol(typical_year_indexP)) {
      unadj_typical_year <- typical_year_indexP[, jj]
      lowess_adj_typical_year <- lowess(x=c(1:MONTHS_IN_YEAR, MONTHS_IN_YEAR+1:MONTHS_IN_YEAR), 
                                        y=c(unadj_typical_year, unadj_typical_year), f=1/4)
      lowess_adj_typical_year <- lowess_adj_typical_year$y[c(MONTHS_IN_YEAR+1:6, 7:MONTHS_IN_YEAR)]
      typical_year_indexP[, jj] <- lowess_adj_typical_year
    }
    
    # correlation for each shift
    region_corr_data_list <- list()
    shifts <- seq(-5, 6)
    for (jj in seq_along(shifts)) {
      shift <- shifts[jj]
      shifted_typical_year_indexP <- typical_year_indexP
      if (shift!=0) shifted_typical_year_indexP <- rbind(tail(shifted_typical_year_indexP, shift), 
                                                         head(shifted_typical_year_indexP, -shift))
      if (length(unique(typical_year_cases))==1) {
        region_corr_data_list[[jj]] <- c(city_id=this_city_id, state=this_state, city=this_city, shift=shift,
                                         mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
      }
      else {
        corrs <- apply(shifted_typical_year_indexP, 2, function(x) cor(x, typical_year_cases, method="spearman"))
        region_corr_data_list[[jj]] <- c(city_id=this_city_id, state=this_state, city=this_city, shift=shift,
                                         mean=mean(corrs), median=median(corrs), sd=sd(corrs), 
                                         l95=quantile(corrs, probs=0.025, names=FALSE), 
                                         u95=quantile(corrs, probs=0.975, names=FALSE))
      }
    }
    region_corr_data_list <- do.call(rbind, region_corr_data_list)
    corr_data_list[[ii]] <- region_corr_data_list
    setTxtProgressBar(pb, ii)
  }
  region_typical_year_corr_data <- do.call(rbind, corr_data_list) %>%
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), shift=as.numeric(shift), mean=as.numeric(mean), 
           median=as.numeric(median), sd=as.numeric(sd), l95=as.numeric(l95), u95=as.numeric(u95))
  
  case_metric_name <- case_metric
  if (logscale) case_metric_name <- paste0("log_", case_metric)
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_typical_year_corr_", case_metric_name, "_data.Rdata")
  save(region_typical_year_corr_data, file=fileout)
}
{ compute_typical_year_corr(regions, years=2000:2014, case_metric="abs_cases", logscale=FALSE)
  compute_typical_year_corr(regions, years=2000:2014, case_metric="abs_cases", logscale=TRUE)
  compute_typical_year_corr(regions, years=2000:2014, case_metric="incidence", logscale=FALSE)
  compute_typical_year_corr(regions, years=2000:2014, case_metric="incidence", logscale=TRUE) }

#################################################################################################
# compute the typical year correlations between the climate variables and log incidence

compute_typical_year_climate_corr <- function(regions, climate_variable, years=2000:2014) {
  # read in the case data
  case_data <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")) %>% filter(year %in% years)
  pop_data <- readRDS(paste0(POPULATION_DATA_FOLDER, "bra2_population_data_2000_2014.rds")) %>% filter(year %in% years)
  
  # compute correlation summary statistics for the typical year data at zero phase shift
  corr_data_list <- list()
  MONTHS_IN_YEAR <- 12
  cat("## Computing the correlation between case data and climate variables for a typical year under different shifts...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    this_state <- regions[ii, "state"]
    this_city <- regions[ii, "city"]
    this_city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata"))) {
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, shift=0, corr=NA)
      next
    }
    
    # compute typical year for case data
    region_case_data <- case_data %>%
      filter(state==this_state, city==this_city) %>% 
      arrange(-desc(year), -desc(month)) %>% 
      left_join(pop_data, by=c("state", "city", "year")) %>%
      unique() %>%
      mutate(cases=cases/n*10^5) %>% 
      mutate(cases=log(cases+1))
    
    # only compute for regions with complete case data
    if (any(is.na(region_case_data$cases))) { 
      typical_year_cases <- rep(NA, MONTHS_IN_YEAR)
    } else {
      region_case_data <- as.matrix(region_case_data[, "cases"])
      typical_year_cases <- rep(NA, MONTHS_IN_YEAR)
      for (tt in 1:MONTHS_IN_YEAR) {
        typical_year_cases[tt] <- mean(region_case_data[seq(tt, length(region_case_data), MONTHS_IN_YEAR)], na.rm=FALSE)
      }
      unadj_typical_year <- typical_year_cases
      lowess_adj_typical_year <- lowess(x=c(1:MONTHS_IN_YEAR, MONTHS_IN_YEAR+1:MONTHS_IN_YEAR), 
                                        y=c(unadj_typical_year, unadj_typical_year), f=1/4)
      lowess_adj_typical_year <- lowess_adj_typical_year$y[c(MONTHS_IN_YEAR+1:6, 7:MONTHS_IN_YEAR)]
      typical_year_cases <- lowess_adj_typical_year
    }
    
    # compute typical year for climate variable
    climate_variable_colnames <- c("T", "H", "R")
    names(climate_variable_colnames) <- c("temperature", "humidity", "precipitation")
    climate_data <- read.csv(file.path(CLIMATE_DATA_FOLDER, paste0("BR_mun_", this_city_id, "_climate_entire_region.csv"))) %>%
      arrange(desc(-year), desc(-month), desc(-day)) %>% 
      filter(year %in% years)
    climate_data <- as.vector(climate_data[, climate_variable_colnames[climate_variable]])
    typical_year_climate <- rep(NA, MONTHS_IN_YEAR)
    for (tt in 1:MONTHS_IN_YEAR) {
      typical_year_climate[tt] <- mean(climate_data[seq(tt, length(climate_data), MONTHS_IN_YEAR)], na.rm=FALSE)
    }
    unadj_typical_year <- typical_year_climate
    lowess_adj_typical_year <- lowess(x=c(1:MONTHS_IN_YEAR, MONTHS_IN_YEAR+1:MONTHS_IN_YEAR), 
                                      y=c(unadj_typical_year, unadj_typical_year), f=1/4)
    lowess_adj_typical_year <- lowess_adj_typical_year$y[c(MONTHS_IN_YEAR+1:6, 7:MONTHS_IN_YEAR)]
    typical_year_climate <- lowess_adj_typical_year
    
    # correlation at zero phase shift
    if (length(unique(typical_year_cases))==1 || length(unique(typical_year_climate))==1) {
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, shift=0, corr=NA)
    }
    else {
      corr <- cor(typical_year_climate, typical_year_cases, method="spearman")
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, shift=0, corr=corr)
    }
    setTxtProgressBar(pb, ii)
  }
  region_typical_year_corr_data <- do.call(rbind, corr_data_list) %>%
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), shift=as.numeric(shift), corr=as.numeric(corr))
  
  fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_typical_year_corr_log_incidence_", climate_variable, "_data.Rdata")
  save(region_typical_year_corr_data, file=fileout)
}
{ compute_typical_year_climate_corr(regions, climate_variable="temperature", years=2000:2014) 
  compute_typical_year_climate_corr(regions, climate_variable="humidity", years=2000:2014) 
  compute_typical_year_climate_corr(regions, climate_variable="precipitation", years=2000:2014) } 


#################################################################################################
# compute correlation between case data and index P for each region across the entire period 

compute_entire_period_corr <- function(regions, years=2000:2014, phase_shift=FALSE, case_metric="abs_cases", logscale=FALSE) {
  # read in the case data
  case_data <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds")) %>% 
    filter(year %in% years)
  if (case_metric=="incidence")
    pop_data <- readRDS(paste0(POPULATION_DATA_FOLDER, "bra2_population_data_2000_2014.rds")) %>% 
      filter(year %in% years)
  
  # set up the output name
  case_metric_name <- case_metric
  if (logscale) case_metric_name <- paste0("log_", case_metric)
  
  # compute correlation summary statistics for the entire period
  corr_data_list <- list()
  cat("## Computing the correlation between case data and index P across entire period...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    this_state <- regions[ii, "state"]
    this_city <- regions[ii, "city"]
    this_city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata"))) { # no index P data
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
      if (phase_shift) corr_data_list[[ii]] <- c(corr_data_list[[ii]], shift=NA)
      next
    }
    
    # organise case data and index P data
    region_case_data <- case_data %>%
      filter(state==this_state, city==this_city) %>% 
      arrange(-desc(year), -desc(month))
    if (case_metric=="incidence") {
      region_case_data %<>% 
        left_join(pop_data, by=c("state", "city", "year")) %>%
        unique() %>%
        mutate(cases=cases/n*10^5)
    }
    if (logscale) region_case_data %<>% mutate(cases=log(cases+1))
    region_cases <- pull(region_case_data, cases)
    
    if (any(is.na(region_cases))) { # incomplete case data
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
      if (phase_shift) corr_data_list[[ii]] <- c(corr_data_list[[ii]], shift=NA)
      next
    }
    
    load(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata")) #sim_p
    sim_p %<>% 
      mutate(year=lubridate::year(as.Date(date))) %>% 
      filter(year %in% years)
    sim_p <- sim_p %>% arrange(-desc(date))
    dates <- sim_p[, 1]
    data <- as.matrix(sim_p[, -1])
    
    # smooth the case and index P data
    lowess_adj_region_cases <- lowess(x=1:length(region_cases), y=region_cases, f=3/length(region_cases))$y
    lowess_adj_indexP_data <- apply(data, 2, function(x) lowess(x=1:length(x), y=x, f=3/length(x))$y)

    # correct for shifts in data
    if (phase_shift) {
      shifts <- seq(-5, 6)
      shift_corr_data_list <- list()
      for (jj in seq_along(shifts)) {
        shift <- shifts[jj]
        
        if (shift>0) {
          shifted_lowess_adj_region_cases <- lowess_adj_region_cases %>% tail(-shift)
          shifted_lowess_adj_indexP_data <- head(lowess_adj_indexP_data, -shift)
          shifted_dates <- tail(dates, -shift)
        } else if (shift<0) {
          shifted_lowess_adj_region_cases <- lowess_adj_region_cases %>% head(shift)
          shifted_lowess_adj_indexP_data <- tail(lowess_adj_indexP_data, shift)
          shifted_dates <- head(dates, shift)
        } else {
          shifted_lowess_adj_region_cases <- lowess_adj_region_cases 
          shifted_lowess_adj_indexP_data <- lowess_adj_indexP_data
          shifted_dates <- dates
        }
        
        if (length(unique(shifted_lowess_adj_region_cases))==1) {
          shift_corr_data_list[[jj]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                          mean=NA, median=NA, sd=NA, l95=NA, u95=NA, shift=shift)
        } else {
          corrs <- apply(shifted_lowess_adj_indexP_data, 2, function(x) cor(x, shifted_lowess_adj_region_cases, method="spearman"))
          shift_corr_data_list[[jj]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                          mean=mean(corrs), median=median(corrs), sd=sd(corrs), 
                                          l95=quantile(corrs, probs=0.025, names=FALSE), 
                                          u95=quantile(corrs, probs=0.975, names=FALSE), 
                                          shift=shift)
        }
      }
      corr_data_list[[ii]] <- do.call(rbind, shift_corr_data_list)
    } else {
      if (length(unique(lowess_adj_region_cases))==1) {
        corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                  mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
      } else {
        corrs <- apply(lowess_adj_indexP_data, 2, function(x) cor(x, lowess_adj_region_cases, method="spearman"))
        corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                  mean=mean(corrs), median=median(corrs), sd=sd(corrs), 
                                  l95=quantile(corrs, probs=0.025, names=FALSE), 
                                  u95=quantile(corrs, probs=0.975, names=FALSE))
      }
    }
    setTxtProgressBar(pb, ii)
  }
  close(pb)
  
  region_entire_period_corr_data <- do.call(rbind, corr_data_list) %>% 
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), mean=as.numeric(mean), median=as.numeric(median), 
           sd=as.numeric(sd), l95=as.numeric(l95), u95=as.numeric(u95))
  
  if (phase_shift) {
    region_entire_period_corr_opt_shift_data <- region_entire_period_corr_data %>%
      mutate(shift=as.numeric(shift))
    fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_entire_period_corr_phase_shift_", case_metric_name, "_data.Rdata")
    save(region_entire_period_corr_opt_shift_data, file=fileout)
  }
  else {
    fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_entire_period_corr_", case_metric_name, "_data.Rdata")
    save(region_entire_period_corr_data, file=fileout)
  }
}
{ compute_entire_period_corr(regions, years=2000:2014, case_metric="abs_cases", logscale=FALSE)
  compute_entire_period_corr(regions, years=2000:2014, case_metric="abs_cases", logscale=TRUE)
  compute_entire_period_corr(regions, years=2000:2014, case_metric="incidence", logscale=FALSE)
  compute_entire_period_corr(regions, years=2000:2014, case_metric="incidence", logscale=TRUE)
  compute_entire_period_corr(regions, years=2000:2014, phase_shift=TRUE, case_metric="abs_cases", logscale=FALSE)
  compute_entire_period_corr(regions, years=2000:2014, phase_shift=TRUE, case_metric="abs_cases", logscale=TRUE)
  compute_entire_period_corr(regions, years=2000:2014, phase_shift=TRUE, case_metric="incidence", logscale=FALSE)
  compute_entire_period_corr(regions, years=2000:2014, phase_shift=TRUE, case_metric="incidence", logscale=TRUE) }

#################################################################################################
# compute correlation between Index P and cases for each year with varying phae shifts 
# (NOT USED IN ANY ANALYSIS)

compute_yearly_corr <- function(regions, opt_shift=FALSE, case_metric="abs_cases", logscale=FALSE) {
  # read in the case data
  case_data <- readRDS(paste0(CASE_DATA_FOLDER, "bra_DENV_case_data_2000_2014.rds"))
  if (case_metric=="incidence")
    pop_data <- readRDS(paste0(POPULATION_DATA_FOLDER, "bra2_population_data_2000_2014.rds"))
  
  # set up the output name
  case_metric_name <- case_metric
  if (logscale) case_metric_name <- paste0("log_", case_metric)
  
  # get the optimal shift for region
  if (opt_shift) {
    load(paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_typical_year_corr_", case_metric_name, "_data.Rdata")) #region_typical_year_corr_data
    opt_shift_data <- region_typical_year_corr_data %>% 
      group_by(city_id, state, city) %>% 
      summarise(opt_shift=ifelse(is.na(max(mean)), NA, shift[which(median==max(median))]), .groups="keep") %>% 
      as.data.frame()
  }
  
  # compute correlation summary statistics for the entire period
  corr_data_list <- list()
  cat("## Computing the correlation between case data and index P for each year...\n")
  pb <- txtProgressBar(min=1, max = nrow(regions), style = 3)
  for (ii in 1:nrow(regions)) {
    this_state <- regions[ii, "state"]
    this_city <- regions[ii, "city"]
    this_city_id <- regions[ii, "city_id"]
    
    if (!file.exists(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata"))) {
      corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, 
                                year=NA, mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
      if (opt_shift) corr_data_list[[ii]] <- c(corr_data_list[[ii]], shift=NA)
      next
    }
    
    # organise case data and index P data
    region_case_data <- case_data %>%
      filter(state==this_state, city==this_city) %>% 
      arrange(-desc(year), -desc(month))
    if (case_metric=="incidence") {
      region_case_data %<>% 
        left_join(pop_data, by=c("state", "city", "year")) %>%
        unique() %>%
        mutate(cases=cases/n*10^5)
    }
    if (logscale) region_case_data %<>% mutate(cases=log(cases+1))
    
    load(paste0(RAW_OUTPUT_DATA_FOLDER, "BR_mun_", this_city_id, "_SIMP.Rdata")) #sim_p
    sim_p <- sim_p %>% arrange(-desc(date))
    dates <- sim_p[, 1]
    years <- lubridate::year(dates) %>% unique()
    
    # correct for shifts in data
    if (opt_shift) {
      region_opt_shift <- opt_shift_data %>% filter(state==this_state, city==this_city) %>% pull(opt_shift)
      if (is.na(region_opt_shift)) {
        corr_data_list[[ii]] <- c(city_id=this_city_id, state=this_state, city=this_city, year=NA,
                                  mean=NA, median=NA, sd=NA, l95=NA, u95=NA, shift=NA)
        setTxtProgressBar(pb, ii)
        next
      }
      else if (region_opt_shift>0) {
        tmp_cases <- region_case_data %>% pull(cases) %>% head(-region_opt_shift)
        region_case_data <- tail(region_case_data, -region_opt_shift)
        region_case_data$cases <- tmp_cases
        sim_p <- tail(sim_p, -region_opt_shift)
        dates <- tail(dates, -region_opt_shift)
      }
      else if (region_opt_shift<0) {
        tmp_cases <- region_case_data %>% pull(cases) %>% tail(region_opt_shift)
        region_case_data <- head(region_case_data, region_opt_shift)
        region_case_data$cases <- tmp_cases
        sim_p <- head(sim_p, region_opt_shift)
        dates <- head(dates, region_opt_shift)
      }
    }
    
    # compute correlation for region
    region_corr_data <- list()
    for (jj in seq_along(years)) {
      climate_yearly_data <- as.matrix(sim_p[which(lubridate::year(dates)==years[jj]), -1])
      region_yearly_cases <- region_case_data %>% filter(year==years[jj]) %>% pull(cases)
      if (length(unique(region_yearly_cases))==1) {
        region_corr_data[[jj]] <- c(city_id=this_city_id, state=this_state, city=this_city, year=years[jj],
                                    mean=NA, median=NA, sd=NA, l95=NA, u95=NA)
        if (opt_shift) region_corr_data[[jj]] <- c(region_corr_data[[jj]], shift=NA)
      }
      else {
        corrs <- apply(climate_yearly_data, 2, function(x) cor(x, region_yearly_cases, method="spearman"))
        region_corr_data[[jj]] <- c(city_id=this_city_id, state=this_state, city=this_city, year=years[jj],
                                    mean=mean(corrs), median=median(corrs), sd=sd(corrs), 
                                    l95=quantile(corrs, probs=0.025, names=FALSE), 
                                    u95=quantile(corrs, probs=0.975, names=FALSE))
        if (opt_shift) region_corr_data[[jj]] <- c(region_corr_data[[jj]], shift=region_opt_shift)
      }
    }
    corr_data_list[[ii]] <- do.call(rbind, region_corr_data)
    setTxtProgressBar(pb, ii)
  }
  region_yearly_corr_data <- do.call(rbind, corr_data_list) %>% 
    as.data.frame() %>% 
    mutate(city_id=as.numeric(city_id), year=as.numeric(year), mean=as.numeric(mean), median=as.numeric(median), 
           sd=as.numeric(sd), l95=as.numeric(l95), u95=as.numeric(u95))
  
  case_metric_name <- case_metric
  if (logscale) case_metric_name <- paste0("log_", case_metric)
  if (opt_shift) {
    region_yearly_corr_opt_shift_data <- region_yearly_corr_data %>%
      mutate(shift=as.numeric(shift))
    fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_yearly_corr_opt_shift_", case_metric_name, "_data.Rdata")
    save(region_yearly_corr_opt_shift_data, file=fileout)
  }
  else {
    fileout <- paste0(SUMMARY_OUTPUT_DATA_FOLDER, "region_yearly_corr_", case_metric_name, "_data.Rdata")
    save(region_yearly_corr_data, file=fileout)
  }
}
# { compute_yearly_corr(regions, case_metric="abs_cases", logscale=FALSE)
#   compute_yearly_corr(regions, case_metric="abs_cases", logscale=TRUE)
#   compute_yearly_corr(regions, case_metric="incidence", logscale=FALSE)
#   compute_yearly_corr(regions, case_metric="incidence", logscale=TRUE)
#   compute_yearly_corr(regions, opt_shift=TRUE, case_metric="abs_cases", logscale=FALSE)
#   compute_yearly_corr(regions, opt_shift=TRUE, case_metric="abs_cases", logscale=TRUE)
#   compute_yearly_corr(regions, opt_shift=TRUE, case_metric="incidence", logscale=FALSE)
#   compute_yearly_corr(regions, opt_shift=TRUE, case_metric="incidence", logscale=TRUE) }


