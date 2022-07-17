########################################################################################################################
# 
# We create some plots of the DENV transmission suitability (i.e. Index P) across the world. 
# The database of Index P maps for different countries is required to be able to run this code. 
# 
########################################################################################################################

########################################################################################################################
# required packages
library(raster)
library(sp)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(rgdal)
library(scales)
library(modifiedmk)

# some necessary paths
EXTERNAL_HARD_DRIVE_PATH <- ""
MVSE_OUTPUT_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "mvse_output/")

SHP_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "../climate_maps/maps_data_SHP/")
FIGURES_FOLDER <- file.path("figures")

############################################################################################################
# generate a world map of DENV transmission suitability (Index P): mean Index P for a typical year

gen_meanP_typicalYear_map <- function(regions, file_name, output=TRUE) {
  raster_data_list <- list()
  for (ii in 1:length(regions)) {
    region <- regions[ii]
    REGION_MVSE_OUTPUT_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, region, "/")
    REGION_MVSE_OUTPUT_DATA_FOLDER <- paste0(REGION_MVSE_OUTPUT_FOLDER, "data/")
    if (file.exists(paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, "_", "indexP", "_typical_year_", "mean", "_rasters.tif"))) {
      print(region)
      indexP_typical_year_mean_rasters <- raster::brick(paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/",
                                                                region, "_", "indexP", "_typical_year_", "mean", "_rasters.tif"))
    } else {
      load(paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, "_", "indexP", "_typical_year_", "mean", "_rasters.Rdata"))
    }
    raster_data <- calc(indexP_typical_year_mean_rasters, mean)
    rrRaster <- rasterToPoints(raster_data) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    raster_data_list[[ii]] <- rrRaster
  }
  
  # generate plot
  col_scale <- c(colorRampPalette(rev(brewer.pal(11, "RdYlBu"))[1:5])(5), 
                 colorRampPalette(rev(brewer.pal(11, "RdYlBu"))[6:11])(15), 
                 rep(rev(brewer.pal(11, "RdYlBu"))[11], 5))
  world_map <- map_data("world")
  g <- ggplot() + 
    geom_map(data=world_map, map=world_map, mapping=aes(x=long, y=lat, map_id=region), 
             color="black", fill="lightgrey", size=0.1)
  for (ii in seq_along(raster_data_list)) {
    g <- g + geom_tile(data=raster_data_list[[ii]], aes(x=x, y=y, fill=value))
  }
  g <- g + 
    coord_cartesian(ylim=c(-56, 50), xlim=c(-120, 170)) + 
  # cale_y_continuous(limits=c(-56, 54), expand=c(0, 0)) + 
  #   scale_x_continuous(limits=c(-120, 170)) + 
    scale_fill_gradientn(colours=rep(col_scale, each=3), 
                         na.value="white", guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), name=NULL, 
                         limits=c(0, 2.5), breaks=seq(0, 2.5, 0.5)) + 
    theme(panel.grid=element_blank(), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), 
          panel.background=element_blank(), panel.border=element_rect(colour="black", fill=NA)) + 
    theme(legend.position=c(0.05, 0.35), legend.key.size=unit(0.4, "cm"))
  
  if (output) {
    scale_factor <- 0.3
    pdf(file=file.path(FIGURES_FOLDER, file_name), w=8.27, h=11.69*scale_factor)
    print(g)
    a <- dev.off()
  } else {
    return(g)
  }
  return()
}
# examples
# regions <- list.files(MVSE_OUTPUT_FOLDER)
# gen_meanP_typicalYear_map(regions, file_name="meanP_typical_map.pdf")

############################################################################################################
# generate a regional map of DENV transmission suitability (Index P) for a particular time period: 
# mean Index P for 1981-1990 and 2010-2019

gen_meanP_map <- function(region, region_id, years, file_name, output=TRUE) {
  REGION_MVSE_OUTPUT_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, region, "/")
  REGION_MVSE_OUTPUT_DATA_FOLDER <- paste0(REGION_MVSE_OUTPUT_FOLDER, "data/")
  if (file.exists(paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", 
                         region, "_", "indexP", "_yearly_", "mean", "_rasters.tif"))) {
    indexP_yearly_mean_rasters <- raster::brick(x=paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", 
                                                         region, "_", "indexP", "_yearly_", "mean", "_rasters.tif"))
  } else {
    load(file=paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, "_", "indexP", "_yearly_", "mean", "_rasters.Rdata"))
  }
  all_years <- indexP_yearly_mean_rasters %>% names() %>% gsub("X", "", .) %>% substr(., start=1, stop=4)
  requested_years_idx <- which(sapply(all_years, function(x) x %in% as.character(years)))
  subset_indexP_yearly_mean_rasters <- subset(indexP_yearly_mean_rasters, requested_years_idx)
  raster_data <- raster::calc(subset_indexP_yearly_mean_rasters, mean)
  rrRaster <- rasterToPoints(raster_data) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
  
  
  # generate plot
  load(file=paste0(SHP_DATA_FOLDER, region, "/simplified_gadm40_", region_id, "_0.Rdata")) #crop_extent
  g <- ggplot() +  
    geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
    geom_tile(data=rrRaster, aes(x=x, y=y, fill=value)) + 
    geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=0.5) + 
    scale_fill_gradientn(colours=rev(brewer.pal(11, "Spectral")), na.value="white",
                         guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), 
                         limits=c(0, 3), breaks=seq(0, 3, 0.5), name=NULL) + 
    theme_bw() + 
    theme(panel.background = element_blank(), axis.line=element_blank(), 
          axis.text=element_blank(), axis.ticks=element_blank(), 
          axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), 
          plot.title=element_text(hjust=0.5)) + 
    theme(legend.position="right", legend.key.height=unit(2, "cm"), legend.key.width=unit(0.5, "cm"))
  
  if (output) {
    scale_factor <- 0.3
    pdf(file=file.path(FIGURES_FOLDER, file_name), w=8.27, h=11.69*scale_factor)
    print(g)
    a <- dev.off()
  } else {
    return(g)
  }
  return()
}

# examples
# gen_meanP_map(region="Kenya", region_id="KEN", years=1981:1990, file_name="meanP_1981to1990_Kenya.pdf")
# gen_meanP_map(region="Nigeria", region_id="NGA", years=1981:1990, file_name="meanP_1981to1990_Nigeria.pdf")
# gen_meanP_map(region="Kenya", region_id="KEN", years=2010:2019, file_name="meanP_2010to2019_Kenya.pdf")
# gen_meanP_map(region="Nigeria", region_id="NGA", years=2010:2019, file_name="meanP_2010to2019_Nigeria.pdf")

############################################################################################################
# generate a world map of DENV transmission suitability (Index P): mean Index P for 1981-1990, 2010-2019 and time series

gen_meanP_comparison_map <- function(region, region_id, capital_city_data, file_name, output=TRUE) {
  REGION_MVSE_OUTPUT_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, region, "/")
  REGION_MVSE_OUTPUT_DATA_FOLDER <- paste0(REGION_MVSE_OUTPUT_FOLDER, "data/")
  if (file.exists(paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", 
                         region, "_", "indexP", "_yearly_", "mean", "_rasters.tif"))) {
    indexP_yearly_mean_rasters <- raster::brick(x=paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", 
                                                         region, "_", "indexP", "_yearly_", "mean", "_rasters.tif"))
  } else {
    load(file=paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", region, "_", "indexP", "_yearly_", "mean", "_rasters.Rdata"))
  }
  all_years <- indexP_yearly_mean_rasters %>% names() %>% gsub("X", "", .) %>% substr(., start=1, stop=4)
  load(file=paste0(SHP_DATA_FOLDER, region, "/simplified_gadm40_", region_id, "_0.Rdata")) #crop_extent
  
  years_list <- list(1981:1989, 2010:2019)
  col_scale <- c(colorRampPalette(rev(brewer.pal(11, "RdYlBu"))[1:5])(5), 
                 colorRampPalette(rev(brewer.pal(11, "RdYlBu"))[6:11])(15), 
                 rep(rev(brewer.pal(11, "RdYlBu"))[11], 5))
  plot_list <- list()
  for (ii in seq_along(years_list)) {
    years <- years_list[[ii]]
    requested_years_idx <- which(sapply(all_years, function(x) x %in% as.character(years)))
    subset_indexP_yearly_mean_rasters <- subset(indexP_yearly_mean_rasters, requested_years_idx)
    raster_data <- raster::calc(subset_indexP_yearly_mean_rasters, mean)
    rrRaster <- rasterToPoints(raster_data) %>% as.data.frame() %>% setNames(c("x", "y", "value"))
    
    x_range <- c(min(c(rrRaster$x, capital_city_data$x+capital_city_data$x_adj), na.rm=TRUE)-2, 
                 max(c(rrRaster$x, capital_city_data$x+capital_city_data$x_adj), na.rm=TRUE)+1)
    
    # generate plot
    g <- ggplot() +  
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill="grey", color=NA, lwd=0.5) + 
      geom_tile(data=rrRaster, aes(x=x, y=y, fill=value)) + 
      geom_sf(data=st_as_sf(crop_extent), mapping=aes(geometry=geometry), fill=NA, color="black", lwd=1) + 
      geom_point(data=capital_city_data, mapping=aes(x=x, y=y), shape=24, fill="green", color="black") + 
      scale_fill_gradientn(colours=col_scale, na.value="white",
                           guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), 
                           breaks=seq(0, 3, 0.5), limits=c(0, 3), name=NULL) + 
      labs(title=paste(years[1], years[length(years)], sep="-")) + 
      scale_x_continuous(limits=x_range) + 
      theme_bw() + 
      theme(panel.background = element_blank(), axis.line=element_blank(),
            axis.text=element_blank(), axis.ticks=element_blank(),
            axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(),
            plot.title=element_text(hjust=0.5)) +
      theme(legend.position="none", legend.key.height=unit(2, "cm"), legend.key.width=unit(0.5, "cm")) + 
      theme(plot.margin=margin(0, 0, 0, 0, "cm"), plot.title=element_text(vjust=-2), 
            plot.background=element_rect(fill=NA, colour=NA))
    if (ii==1) g <- g + geom_text(data=capital_city_data, mapping=aes(x=x+x_adj, y=y+y_adj, label=city), size=3.5) 
    plot_list[[ii]] <- g
  }
  plot_title <- ggdraw() + draw_label(gsub("([a-z])([A-Z])","\\1 \\2", region), fontface="bold")
  p <- cowplot::plot_grid(plot_list[[1]], NULL, plot_list[[2]], rel_widths=c(1, -0.25, 1), nrow=1)
  p <- cowplot::plot_grid(plot_title, p, ncol=1, rel_heights=c(0.12, 1))
  if (output) {
    scale_factor <- 0.3
    pdf(file=file.path(FIGURES_FOLDER, file_name), w=8.27, h=11.69*scale_factor)
    print(p)
    a <- dev.off()
  } else {
    return(p)
  }
  return()
}
# example
# gen_meanP_comparison_map(region="Angola", region_id="AGO", 
#                                 capital_city_data=data.frame(y=-8.8147, x=13.2302, city="Luanda"), 
#                                 file_name="meanP_Angola_comparison_map.pdf")

############################################################################################################
# generate time series of Index P for a region of interest. 

gen_indexP_timeseries <- function(region) {
  REGION_MVSE_OUTPUT_FOLDER <- paste0(MVSE_OUTPUT_FOLDER, region, "/")
  REGION_MVSE_OUTPUT_DATA_FOLDER <- paste0(REGION_MVSE_OUTPUT_FOLDER, "data/")
  load(paste0(REGION_MVSE_OUTPUT_DATA_FOLDER, "summary_output/", 
              region, "_indexP_entire_region_monthly_summary_dataframe.Rdata")) #indexP_entire_region_monthly_summary_dataframe
  region_data <- indexP_entire_region_monthly_summary_dataframe %>% 
    mutate(lb=mean-sd, ub=mean+sd) %>% 
    mutate(lb=ifelse(lb<0, 0, lb))
  
  # perform trend analysis
  trend_data <- region_data %>%
    transmute(Month=lubridate::month(date), Year=lubridate::year(date), indexP=mean) %>% 
    arrange(-desc(Year), -desc(Month))
  trend_output <- mmkh(trend_data$indexP)
  slope <- trend_output["Sen's slope"]
  pValue <- trend_output["new P-value"]
  trend_plot_data <- data.frame(x=as.Date("2012/01/01"), y=3.4, 
                                label=paste(round(slope, 5)*12, "/year ", "(p=", round(pValue, 3), ")", sep=""))
  
  p <- ggplot(data=region_data) + 
    geom_hline(yintercept=1, color="black", linetype="dashed") + 
    geom_line(mapping=aes(x=date, y=mean), color="firebrick3") + 
    geom_line(mapping=aes(x=date, y=lowess(date, mean, f=1/2)$y), color="black") + 
    geom_text(data=trend_plot_data, mapping=aes(x=x, y=y, label=label), size=3.5) + 
    # geom_line(mapping=aes(x=date, y=lb), color="grey") + 
    # geom_line(mapping=aes(x=date, y=ub), color="grey") +
    scale_x_date(breaks=seq(as.Date("1980/1/1"), as.Date("2020/1/1"), "5 years"), 
                 date_labels = "%Y", limits=as.Date(c("1980-01-01", "2020-03-01")), expand=c(0.01, 0.7)) + 
    coord_cartesian(ylim=c(0, 3.5)) + 
    labs(y="Mean Index P") + 
    theme_bw() + 
    theme(axis.title.x=element_blank())
  return(p)
}
# example
# gen_indexP_timeseries("Angola")

############################################################################################################
# combined plot for presentation

gen_combined_plot <- function(regions, region1, region1_id, region1_capital_city_data, 
                              region2, region2_id, region2_capital_city_data, file_name) {
  meanP_typicalYear_map <- gen_meanP_typicalYear_map(regions, output=FALSE)
  region1_meanP_comparison_map <- gen_meanP_comparison_map(region=region1, region_id=region1_id, 
                                                           capital_city_data=region1_capital_city_data, 
                                                           output=FALSE)
  region2_meanP_comparison_map <- gen_meanP_comparison_map(region=region2, region_id=region2_id, 
                                                           capital_city_data=region2_capital_city_data, 
                                                           output=FALSE)
  region1_timeseries <- gen_indexP_timeseries(region=region1, region_id=region1_id)
  region2_timeseries <- gen_indexP_timeseries(region=region2, region_id=region2_id)
  
  pA <- cowplot::plot_grid(meanP_typicalYear_map, labels="A")
  pB <- cowplot::plot_grid(region1_meanP_comparison_map, region1_timeseries, ncol=1, labels=c("B", ""))
  pC <- cowplot::plot_grid(region2_meanP_comparison_map, region2_timeseries, ncol=1, labels=c("C", ""))
  pBC <- cowplot::plot_grid(pB, pC, nrow=1)
  pABC <- cowplot::plot_grid(pA, NULL, pBC, ncol=1, rel_heights=c(0.8, 0.02, 1))
  scale_factor <- 0.6
  pdf(file=file.path(FIGURES_FOLDER, file_name), w=8.27, h=11.69*scale_factor)
  print(pABC)
  a <- dev.off()
}
# example
# gen_combined_plot(regions=regions, 
#                   region1="Angola", region1_id="AGO", 
#                   region1_capital_city_data=data.frame(y=-8.8147, x=13.2302, y_adj=0.3, x_adj=-3, 
#                                                        city="Luanda"), 
#                   region2="CentralAfricanRepublic", region2_id="CAF", 
#                   region2_capital_city_data=data.frame(y=4.3947, x=18.5582, y_adj=-0.6, x_adj=2.2, 
#                                                        city="Bangui"), 
#                   file_name="indexP_combined_plots.pdf")
