#################################################################################################
#
#  This script plots some diagnostic figures for the climate variables. It can take considerable
#  time to run this script, so it should only be run when the monthly climate data want to 
#  be examined. 
# 
#  This script is run after:
#  1. 01_data2rasters.R
#  2. 02_curate_calculate.R
#  3. 03_do_relative_humidity_range.R
#
#
#################################################################################################

#################################################################################################
print("<MVSE> plotting climate maps ...")

load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag,"_", "dates.Rdata"))   # dates
dates <- as.Date(unlist(dates), format="%Y-%m-%d")

levColors <- c('black',rev(brewer.pal(11, "RdYlBu")))   # colourblind friendly

fileout<- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata")
load(fileout) # temp_rast_list
pdf(paste0(RAW_CLIMATE_PLOT_FOLDER, region, tag, '_monthly_maps_temperature.pdf'), w=10, h=10)
for(mm in 1:length(dates)){
  rrMonth<- temp_rast_list[[mm]]
  g <- gplot(rrMonth) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colours=rep(levColors,each=10), na.value="white", guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0,35)) +
    coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste("temp, time:", dates[mm]))
  print(g)
}
a<-dev.off()

fileout<- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_humrel_list.Rdata")
load(fileout) # humrel_rast_list
pdf(paste0(RAW_CLIMATE_PLOT_FOLDER, region, tag, '_monthly_maps_humidity.pdf'), w=10, h=10)
for(mm in 1:length(dates)){
  rrMonth<- humrel_rast_list[[mm]]
  g<- gplot(rrMonth) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colours=rep(levColors,each=10), na.value="white", guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(20,100)) +
    coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste("hum, time:", dates[mm]))
  print(g)
}
a<-dev.off()

fileout<- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_prec_list.Rdata")
load(fileout) # prec_rast_list
pdf(paste0(RAW_CLIMATE_PLOT_FOLDER, region, tag, '_monthly_maps_precipitation.pdf'), w=10, h=10)
for(mm in 1:length(dates)){
  rrMonth<- prec_rast_list[[mm]]
  g<- gplot(rrMonth) + geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colours=rep(levColors,each=10), na.value="white", guide=guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits=c(0,0.02)) +
    coord_equal() + xlab("longitude") + ylab("latitude") + ggtitle(paste("prec, time:",dates[mm]))
  print(g)
}
a<-dev.off()

