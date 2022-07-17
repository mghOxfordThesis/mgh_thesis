#################################################################################################
#
#  This script is plots diagnostic figures of the estimated Index P for the region of interest. 
# 
#  This script is run after:
#  1. 00_simplify_SHP_files.R (optional)
#  2. 00_generate_MVSE_maps.R
#  3. 01_run_MVSE_per_cell.R
#  4. 01_run_MVSE_region.R
#  5. 02_run_MVSE_per_cell.R
#  6. 03_MVSE2raster.R
#
#################################################################################################

###############################################################################
options(bitmapType='cairo') # make sure it works on ADA without support for X11

# creates an appropriate color gradient given a vector of data
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(min(x),max(x),length.out=colsteps))])
}

###############################################################################
cat("## Step 4: making some maps...\n")
step_4_time <- Sys.time()

## some color pallettes
circColors <- rev(c(brewer.pal(9,"Blues")[4:9],rev(brewer.pal(9,"BuPu")[4:9]))) # circle colors
levColors <- c('black',rev(brewer.pal(11, "RdYlBu"))) # colourblind friendly

## some data
dates <- as.Date(dates, format="%Y-%m-%d")

###############################################################################
cat("## generating monthly maps...\n")

# monthly maps for each summary statistic (mean, median, sd, l95, u95)
gen_monthly_maps <- function(index, metric) {
  raster_data <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_monthly_", metric, "_rasters.tif"))
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_monthly_", metric, "_maps.pdf"), w=12, h=5)
  index_max <- -Inf
  for (mm in 1:length(dates)) {
    month_max <- max(raster_data[[mm]][], na.rm=TRUE)
    if (month_max>index_max) index_max <- month_max
  }
  for(mm in 1:length(dates)){
    rrMonth <- raster_data[[mm]]
    # g <- gplot(rrMonth) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
    #                        limits=c(0, ceiling(index_max))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste("P, time:",dates[mm]))
    # print(g)
    
    rrMonth <- rasterToPoints(rrMonth) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrMonth, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, ceiling(index_max))) + 
      labs(title=paste("P, time:",dates[mm])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_monthly_maps(index="indexP", metric="mean")
  gen_monthly_maps(index="indexP", metric="median")
  gen_monthly_maps(index="indexP", metric="sd")
  gen_monthly_maps(index="indexP", metric="l95")
  gen_monthly_maps(index="indexP", metric="u95")
  
  # gen_monthly_maps(index="Q", metric="mean")
  # gen_monthly_maps(index="Q", metric="median")
  # gen_monthly_maps(index="Q", metric="sd")
  # gen_monthly_maps(index="Q", metric="l95")
  # gen_monthly_maps(index="Q", metric="u95")
  
  # gen_monthly_maps(index="V0", metric="mean")
  # gen_monthly_maps(index="V0", metric="median")
  # gen_monthly_maps(index="V0", metric="sd")
  # gen_monthly_maps(index="V0", metric="l95")
  # gen_monthly_maps(index="V0", metric="u95")
}

# monthly summary map (l95, median and u95)
gen_monthly_summary_map1 <- function(index) {
  metrics <- c("l95", "median", "u95")
  metric_titles <- c("2.5 percentile", "median", "97.5 percentile")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_monthly_", metrics[ii], "_rasters.tif"))
  }
  
  # compute metrics maxes
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:length(dates)) {
      month_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (month_max>metric_maxes[ii]) metric_maxes[ii] <- month_max
    }
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_monthly_", "summary", "_map1.pdf"), w=12, h=5)
  for (mm in 1:length(dates)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrMonth <- raster_data[[ii]][[mm]]
      # g <- gplot(rrMonth) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", dates[mm]))
      
      rrMonth <- rasterToPoints(rrMonth) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrMonth, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\ntime: ", dates[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_monthly_summary_map1("indexP")
  # gen_monthly_summary_map1("Q")
  # gen_monthly_summary_map1("V0")
}

# monthly summary map (mean, median and sd)
gen_monthly_summary_map2 <- function(index) {
  metrics <- c("mean", "median", "sd")
  metric_titles <- c("mean", "median", "sd")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_monthly_", metrics[ii], "_rasters.tif"))
  }
  
  # compute metrics maxes
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:length(dates)) {
      month_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (month_max>metric_maxes[ii]) metric_maxes[ii] <- month_max
    }
  }
  metric_maxes[1] <- max(metric_maxes[1:2]); metric_maxes[2] <- max(metric_maxes[1:2]);
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_monthly_", "summary", "_map2.pdf"), w=12, h=5)
  for (mm in 1:length(dates)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrMonth <- raster_data[[ii]][[mm]]
      # g <- gplot(rrMonth) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", dates[mm]))
      
      rrMonth <- rasterToPoints(rrMonth) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrMonth, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\ntime: ", dates[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_monthly_summary_map2("indexP")
  # gen_monthly_summary_map2("Q")
  # gen_monthly_summary_map2("V0")
}

###############################################################################
cat("## generating yearly maps...\n")

# yearly summary map (l95, median and u95)
gen_yearly_summary_map1 <- function(index) {
  metrics <- c("l95", "median", "u95")
  metric_titles <- c("2.5 percentile", "median", "97.5 percentile")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_yearly_", metrics[ii], "_rasters.tif"))
  }
  
  # compute metrics maxes
  years <- lubridate::year(dates) %>% unique() %>% sort()
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:length(years)) {
      year_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (year_max>metric_maxes[ii]) metric_maxes[ii] <- year_max
    }
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_yearly_", "summary", "_map1.pdf"), w=12, h=5)
  for (mm in 1:length(years)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrYear <- raster_data[[ii]][[mm]]
      # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
      
      rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\ntime: ", years[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_yearly_summary_map1("indexP")
  # gen_yearly_summary_map1("Q")
  # gen_yearly_summary_map1("V0")
}

# yearly summary map (mean, median, sd)
gen_yearly_summary_map2 <- function(index) {
  metrics <- c("mean", "median", "sd")
  metric_titles <- c("mean", "median", "sd")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_yearly_", metrics[ii], "_rasters.tif"))
  }
  
  # compute metrics maxes
  years <- lubridate::year(dates) %>% unique() %>% sort()
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:length(years)) {
      year_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (year_max>metric_maxes[ii]) metric_maxes[ii] <- year_max
    }
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_yearly_", "summary", "_map2.pdf"), w=12, h=5)
  for (mm in 1:length(years)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrYear <- raster_data[[ii]][[mm]]
      # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
      
      rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\ntime: ", years[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_yearly_summary_map2("indexP")
  # gen_yearly_summary_map2("Q")
  # gen_yearly_summary_map2("V0")
}

# yearly months above 1 summary map (mean, median)
gen_yearly_months_above_one_summary_map1 <- function(index) {
  metrics <- c("mean", "median")
  metric_titles <- c("mean", "median")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_yearly_months_above_one_", metrics[ii], "_rasters.tif"))
  }
  
  
  years <- lubridate::year(dates) %>% unique() %>% sort()
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_yearly_months_above_one_", "summary", "_map1.pdf"), w=12, h=5)
  for (mm in 1:length(years)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrYear <- raster_data[[ii]][[mm]]
      # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
      
      rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, 12)) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\ntime: ", years[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_yearly_months_above_one_summary_map1("indexP")
  #gen_yearly_months_above_one_summary_map1("Q")
  #gen_yearly_months_above_one_summary_map1("V0")
}

# yearly months above 1 summary map (sd, lower bound, upper bound)
gen_yearly_months_above_one_summary_map2 <- function(index) {
  metrics <- c("sd", "l95", "u95")
  metric_titles <- c("sd", "2.5 percentile", "97.5 percentile")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_yearly_months_above_one_", metrics[ii], "_rasters.tif"))
  }
  
  
  years <- lubridate::year(dates) %>% unique() %>% sort()
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:length(years)) {
      year_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (year_max>metric_maxes[ii]) metric_maxes[ii] <- year_max
    }
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_yearly_months_above_one_", "summary", "_map2.pdf"), w=12, h=5)
  for (mm in 1:length(years)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrYear <- raster_data[[ii]][[mm]]
      # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
      
      rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\ntime: ", years[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
if (all_plots) {
  gen_yearly_months_above_one_summary_map2("indexP")
  #gen_yearly_months_above_one_summary_map1("Q")
  #gen_yearly_months_above_one_summary_map1("V0")
}

###############################################################################
cat("## generating entire period maps...\n")

# entire period summary map (l95, median and u95)
gen_entire_period_summary_map1 <- function(index) {
  metrics <- c("l95", "median", "u95")
  metric_titles <- c("2.5 percentile", "median", "97.5 percentile")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::raster(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_entire_period_", metrics[ii], "_raster.tif"))
  }
  
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    metric_maxes[ii] <- max(raster_data[[ii]][], na.rm=TRUE)
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_entire_period_", "summary", "_map1.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      #geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      #geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, ceiling(metric_maxes[ii]))) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_entire_period_summary_map1("indexP")
  # gen_entire_period_summary_map1("Q")
  # gen_entire_period_summary_map1("V0")
}

# entire period summary map (mean, median, sd)
gen_entire_period_summary_map2 <- function(index) {
  metrics <- c("mean", "median", "sd")
  metric_titles <- c("mean", "median", "sd")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::raster(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_entire_period_", metrics[ii], "_raster.tif"))
  }
  
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    metric_maxes[ii] <- max(raster_data[[ii]][], na.rm=TRUE)
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_entire_period_", "summary", "_map2.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, ceiling(metric_maxes[ii]))) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_entire_period_summary_map2("indexP")
  # gen_entire_period_summary_map2("Q")
  # gen_entire_period_summary_map2("V0")
}

###############################################################################
cat("## generating typical year maps...\n")

# typical year map for a given index
gen_typical_year_map <- function(index) {
  raster_data <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_median_rasters.tif"))
  
  # compute metrics maxes
  index_max <- -Inf
  for (ii in 1:nlayers(raster_data)) {
    month_max <- max(raster_data[[ii]][], na.rm=TRUE)
    if (month_max>index_max) index_max <- month_max
  }
  
  # generate plots
  plot_list <- list()
  month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
  for (ii in 1:nlayers(raster_data)) {
    rrMonth <- raster_data[[ii]]
    # g <- gplot(rrMonth) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
    #                        limits=c(0, ceiling(index_max))) +
    #   labs(x="longitude", y="latitude", title=paste0(index, "\nmonth: ", month_names[ii])) +
    #   theme(legend.position="bottom", legend.title=element_blank(),
    #         axis.title=element_text(size=10), plot.title=element_text(size=10),
    #         axis.text=element_text(size=10)) +
    #   coord_equal()
    
    rrMonth <- rasterToPoints(rrMonth) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrMonth, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, ceiling(index_max))) + 
      labs(title=month_names[ii]) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="bottom")
    common_legend <- get_legend(g)
    plot_list[[ii]] <- g + theme(legend.position="none")
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=2)
  g <- cowplot::plot_grid(g, common_legend, ncol=1, rel_heights=c(1, 0.1))
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_median_map.pdf"), w=12, h=5)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_map("indexP")
  # gen_typical_year_map("Q")
  # gen_typical_year_map("V0")
}

# typical year map along with uncertainty
gen_typical_year_summary_maps1 <- function(index) {
  metrics <- c("l95", "median", "u95")
  metric_titles <- c("2.5 percentile", "median", "97.5 percentile")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_", metrics[ii], "_rasters.tif"))
  }
  
  # compute metrics maxes
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:12) {
      month_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (month_max>metric_maxes[ii]) metric_maxes[ii] <- month_max
    }
  }
  
  # generate plots
  month_names <- c("January", "February", "March", "April", "May",
                   "June", "July", "August", "September", "October", "November", "December")
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_", "summary", "_maps1.pdf"), w=12, h=5)
  for (mm in seq_along(month_names)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrMonth <- raster_data[[ii]][[mm]]
      # g <- gplot(rrMonth) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
      #                        limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() +
      #   labs(x="longitude", y="latitude",
      #        title=paste0(index, "(", metric_titles[ii], ")", "\nmonth: ", month_names[mm]))
      rrMonth <- rasterToPoints(rrMonth) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrMonth, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\nmonth: ", month_names[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
{
  gen_typical_year_summary_maps1("indexP")
  # gen_typical_year_summary_maps1("Q")
  # gen_typical_year_summary_maps1("V0")
}

# typical year map along with uncertainty
gen_typical_year_summary_maps2 <- function(index) {
  metrics <- c("mean", "sd")
  metric_titles <- c("mean", "sd")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_", metrics[ii], "_rasters.tif"))
  }
  
  # compute metrics maxes
  metric_maxes <- rep(-Inf, length(metrics))
  for (ii in seq_along(metrics)) {
    for (mm in 1:12) {
      month_max <- max(raster_data[[ii]][[mm]][], na.rm=TRUE)
      if (month_max>metric_maxes[ii]) metric_maxes[ii] <- month_max
    }
  }
  
  # generate plots
  month_names <- c("January", "February", "March", "April", "May",
                   "June", "July", "August", "September", "October", "November", "December")
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_", "summary", "_maps2.pdf"), w=12, h=5)
  for (mm in seq_along(month_names)) {
    plot_list <- list()
    for (ii in seq_along(metrics)) {
      rrMonth <- raster_data[[ii]][[mm]]
      # g <- gplot(rrMonth) + geom_tile(aes(fill = value)) +
      #   scale_fill_gradientn(colours=levColors, na.value="white",
      #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
      #                        limits=c(0, ceiling(metric_maxes[ii]))) +
      #   coord_equal() +
      #   labs(x="longitude", y="latitude",
      #        title=paste0(index, "(", metric_titles[ii], ")", "\nmonth: ", month_names[mm]))
      rrMonth <- rasterToPoints(rrMonth) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
      g <- ggplot() +
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
        geom_tile(data=rrMonth, aes(x=x, y=y, fill=value)) + 
        geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
        scale_fill_gradientn(colours=levColors, na.value="white",
                             guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                             limits=c(0, ceiling(metric_maxes[ii]))) + 
        labs(title=paste0(index, ": ", metric_titles[ii], "\nmonth: ", month_names[mm])) + 
        theme_bw() + 
        theme(panel.background = element_blank(), axis.line=element_blank(), 
              axis.text=element_blank(), axis.ticks=element_blank(), 
              axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
              legend.position="right")
      plot_list[[ii]] <- g
    }
    g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
    print(g)
  }
  a <-dev.off()
}
{
  gen_typical_year_summary_maps2("indexP")
  # gen_typical_year_summary_maps2("Q")
  # gen_typical_year_summary_maps2("V0")
}

# typical year peak/trough month summary map (mean, median)
gen_typical_year_peak_month_summary_map1 <- function(index) {
  metrics <- c("mean", "median")
  metric_titles <- c("mean", "median")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_peak_month_", metrics[ii], "_rasters.tif"))
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_peak_month_", "summary", "_map1.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, 12)) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_peak_month_summary_map1("indexP")
  #gen_typical_year_peak_month_summary_map1("Q")
  #gen_typical_year_peak_month_summary_map1("V0")
}

# typical year peak/trough month summary map (sd, lb, up)
gen_typical_year_peak_month_summary_map2 <- function(index) {
  metrics <- c("sd", "lb", "ub")
  metric_titles <- c("sd", "Lower bound", "Upper bound")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_peak_month_", metrics[ii], "_rasters.tif"))
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_peak_month_", "summary", "_map2.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, 12)) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_peak_month_summary_map2("indexP")
  #gen_typical_year_peak_month_summary_map2("Q")
  #gen_typical_year_peak_month_summary_map2("V0")
}

# typical year peak/trough month summary map (mean, median)
gen_typical_year_trough_month_summary_map1 <- function(index) {
  metrics <- c("mean", "median")
  metric_titles <- c("mean", "median")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_trough_month_", metrics[ii], "_rasters.tif"))
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_trough_month_", "summary", "_map1.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, 12)) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_trough_month_summary_map1("indexP")
  #gen_typical_year_trough_month_summary_map1("Q")
  #gen_typical_year_trough_month_summary_map1("V0")
}

# typical year peak/trough month summary map (sd, lb, up)
gen_typical_year_trough_month_summary_map2 <- function(index) {
  metrics <- c("sd", "lb", "ub")
  metric_titles <- c("sd", "Lower bound", "Upper bound")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_trough_month_", metrics[ii], "_rasters.tif"))
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_trough_month_", "summary", "_map2.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, 12)) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_trough_month_summary_map2("indexP")
  #gen_typical_year_trough_month_summary_map2("Q")
  #gen_typical_year_trough_month_summary_map2("V0")
}

# typical year months above 1 summary map (mean, median)
gen_typical_year_months_above_one_summary_map1 <- function(index) {
  metrics <- c("mean", "median")
  metric_titles <- c("mean", "median")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_months_above_one_", metrics[ii], "_rasters.tif"))
  }
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_months_above_one_", "summary", "_map1.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, 12)) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_months_above_one_summary_map1("indexP")
  #gen_typical_year_months_above_one_summary_map1("Q")
  #gen_typical_year_months_above_one_summary_map1("V0")
}

# typical year months above 1 summary map (sd, lower bound, upper bound)
gen_typical_year_months_above_one_summary_map2 <- function(index) {
  metrics <- c("sd", "l95", "u95")
  metric_titles <- c("sd", "2.5 percentile", "97.5 percentile")
  raster_data <- list()
  for (ii in seq_along(metrics)) {
    raster_data[[ii]] <- raster::brick(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_typical_year_months_above_one_", metrics[ii], "_rasters.tif"))
  }
  
  
  years <- lubridate::year(dates) %>% unique() %>% sort()
  
  # generate plots
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_typical_year_months_above_one_", "summary", "_map2.pdf"), w=12, h=5)
  plot_list <- list()
  for (ii in seq_along(metrics)) {
    rrYear <- raster_data[[ii]]
    # g <- gplot(rrYear) + geom_tile(aes(fill = value)) +
    #   scale_fill_gradientn(colours=levColors, na.value="white",
    #                        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0, ceiling(metric_maxes[ii]))) +
    #   coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste0(index, "(", metric_titles[ii], ")\ntime: ", years[mm]))
    
    rrYear <- rasterToPoints(rrYear) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    g <- ggplot() +
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrYear, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
      scale_fill_gradientn(colours=levColors, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           limits=c(0, 12)) + 
      labs(title=paste0(index, ": ", metric_titles[ii])) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(), 
            axis.text=element_blank(), axis.ticks=element_blank(), 
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
            legend.position="right")
    plot_list[[ii]] <- g
  }
  g <- cowplot::plot_grid(plotlist=plot_list, nrow=1)
  print(g)
  a <-dev.off()
}
{
  gen_typical_year_months_above_one_summary_map2("indexP")
  #gen_typical_year_months_above_one_summary_map2("Q")
  #gen_typical_year_months_above_one_summary_map2("V0")
}


###############################################################################
cat("## generating entire region plots...\n")

# summary line plot with median and 95% CI for entire region
gen_entire_region_plot1 <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_entire_region_monthly_summary_dataframe.Rdata"))
  eval(parse(text=paste0("data", "<-", index, "_entire_region_monthly_summary_dataframe")))
  
  ind <- c("median"="median", "l95"="95CI", "u95"="95CI")
  index_df <- dplyr::select(data, c("date", "median", "l95", "u95")) %>%
    gather(key=percentile, value=value, 2:ncol(.)) %>%
    mutate(color=ind[percentile])
  
  # ggplot object
  lab_dates <- seq(index_df$date[1], tail(index_df$date, 1),
                   by=paste0(round(tail(index_df$date, 1)-index_df$date[1])/5, " days"))
  p <- ggplot(data=index_df) +
    geom_line(aes(x=date, y=value, group=percentile, color=color, alpha=color)) +
    theme_bw() +
    labs(y=index) +
    scale_color_manual(breaks=c("median", "95CI"), values=c("tomato3", "grey37"),
                       labels=c("median", "95% CI")) +
    scale_alpha_manual(breaks=c("median", "95CI"), values=c(1, 0.5), labels=NULL, guide="none") +
    scale_x_date(breaks=lab_dates, date_labels="%b %y") +
    theme(legend.position = c(0.9, 0.9), legend.title=element_blank(),
          legend.background = element_rect(color="black"))
  if (index=="indexP") p <- p + geom_hline(yintercept=1, color="black", alpha=0.5, linetype="dashed")
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_region_", "summary", ".pdf"), w=12, h=5)
  print(p)
  a <-dev.off()
}

{
  gen_entire_region_plot1("indexP")
  # gen_entire_region_plot1("Q")
  # gen_entire_region_plot1("V0")
}

# summary line plot with mean and +-sd for entire region
gen_entire_region_plot2 <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_entire_region_monthly_summary_dataframe.Rdata"))
  eval(parse(text=paste0("data", "<-", index, "_entire_region_monthly_summary_dataframe")))
  
  ind <- c("mean"="mean", "lower"="+/- 1 sd", "upper"="+/- 1 sd")
  index_df <- dplyr::select(data, c("date", "mean", "sd")) %>%
    mutate(lower=mean-sd, upper=mean+sd) %>%
    mutate(lower=ifelse(lower<0, 0, lower)) %>% 
    dplyr::select(date, mean, lower, upper) %>%
    gather(key=percentile, value=value, 2:ncol(.)) %>% 
    mutate(color=ind[percentile])
  
  # ggplot object
  lab_dates <- seq(index_df$date[1], tail(index_df$date, 1),
                   by=paste0(round(tail(index_df$date, 1)-index_df$date[1])/5, " days"))
  p <- ggplot(data=index_df) +
    geom_line(aes(x=date, y=value, group=percentile, color=color, alpha=color)) +
    theme_bw() +
    labs(y=index) +
    scale_color_manual(breaks=c("mean", "+/- 1 sd"), values=c("tomato3", "grey37"),
                       labels=c("mean", "+/- 1 sd")) +
    scale_alpha_manual(breaks=c("mean", "+/- 1 sd"), values=c(1, 0.5), labels=NULL, guide="none") +
    scale_x_date(breaks=lab_dates, date_labels="%b %y") +
    theme(legend.position = c(0.9, 0.9), legend.title=element_blank(),
          legend.background = element_rect(color="black"))
  if (index=="indexP") p <- p + geom_hline(yintercept=1, color="black", alpha=0.5, linetype="dashed")
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_region_", "summary2", ".pdf"), w=12, h=5)
  print(p)
  a <-dev.off()
}

{
  gen_entire_region_plot2("indexP")
  # gen_entire_region_plot2("Q")
  # gen_entire_region_plot2("V0")
}

# get typical year for region with mean and sd
gen_entire_region_typical_year1 <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_entire_region_typical_year_dataframe.Rdata")) 
  eval(parse(text=paste0("data", "<-", index, "_entire_region_typical_year_dataframe")))
  
  index_wide <- dplyr::select(data, c("month", "mean", "sd")) %>%
    mutate(lower=mean-sd, upper=mean+sd) %>%
    mutate(lower=ifelse(lower<0, 0, lower)) %>% 
    dplyr::select(month, mean, lower, upper)
  index_long <- index_wide %>%
    gather(key=percentile, value=value, 2:ncol(.))
  
  # ggplot object
  p <- ggplot() +
    geom_line(data=index_long, aes(x=month, y=value, group=percentile), color="tomato3") +
    geom_ribbon(data=index_wide, mapping=aes(x=month, ymin=lower, ymax=upper), fill="tomato3", alpha=0.7) + 
    theme_bw() +
    scale_x_continuous(breaks=seq(1, 12), limits=c(1, 12)) + 
    labs(y=index, x="Month", title="Mean and +- 1 sd")
  if (index=="indexP") p <- p + geom_hline(yintercept=1, color="black", alpha=0.5, linetype="dashed")
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_region_", "typical_year1.pdf"), w=12, h=5)
  print(p)
  a <-dev.off()
}

{
  gen_entire_region_typical_year1("indexP")
  # gen_entire_region_typical_year1("Q")
  # gen_entire_region_typical_year1("V0")
}


# typical year summary plot with median, l95 and u95
gen_entire_region_typical_year2 <- function(index) {
  load(paste0(MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, tag, "_", index, "_entire_region_typical_year_dataframe.Rdata")) 
  eval(parse(text=paste0("data", "<-", index, "_entire_region_typical_year_dataframe")))
  
  index_wide <- dplyr::select(data, c("month", "median", "l95", "u95"))
  index_long <- index_wide %>%
    gather(key=percentile, value=value, 2:ncol(.))
  
  # ggplot object
  p <- ggplot() +
    geom_line(data=index_long, aes(x=month, y=value, group=percentile), color="tomato3") +
    geom_ribbon(data=index_wide, mapping=aes(x=month, ymin=l95, ymax=u95), fill="tomato3", alpha=0.7) + 
    theme_bw() +
    scale_x_continuous(breaks=seq(1, 12), limits=c(1, 12)) + 
    labs(y=index, x="Month", title="Median and 95% CI")
  if (index=="indexP") p <- p + geom_hline(yintercept=1, color="black", alpha=0.5, linetype="dashed")
  pdf(paste0(MVSE_OUTPUT_PLOTS_FOLDER, region, tag, "_", index, "_region_", "typical_year2.pdf"), w=12, h=5)
  print(p)
  a <-dev.off()
}

{
  gen_entire_region_typical_year2("indexP")
  # gen_entire_region_typical_year2("Q")
  # gen_entire_region_typical_year2("V0")
}

###############################################################################

step_4_time <- Sys.time() - step_4_time
cat("## step 4 time:", step_4_time, attributes(step_4_time)$units, "\n")

cat("\n")