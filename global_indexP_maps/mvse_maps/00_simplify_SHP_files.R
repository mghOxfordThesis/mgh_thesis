#################################################################################################
#
#  This script is optional. It simplifies the geographical polygons to speed up the plotting 
#  process later. 
#
#
#################################################################################################

#################################################################################################
# load libraries
library(pacman)
pacman::p_load(rgdal, sf, shapefiles)


# folder path
EXTERNAL_HARD_DRIVE_PATH <- ""
SHP_DATA_FOLDER <- paste0(EXTERNAL_HARD_DRIVE_PATH, "../climate_maps/maps_data_SHP/")

#################################################################################################
# simplify_SHP_file takes as argument a file name for a SHP file (file) and a file path (file_path) 
# and saves an Rdata object of the simplified simple feature object. 
simplify_SHP_file <- function(file, file_path) {
  crop_extent <- readOGR(paste0(file_path, file, ".shp"))

  for(i in 1:length(crop_extent@polygons)){
    for(j in 1:length(crop_extent@polygons[[i]]@Polygons)){
      temp <- as.data.frame(crop_extent@polygons[[i]]@Polygons[[j]]@coords)
      names(temp) <- c("x", "y")
      temp2 <- dp(temp, 0.03) ##simplification threshold
      crop_extent@polygons[[i]]@Polygons[[j]]@coords <- as.matrix(cbind(temp2$x, temp2$y))
    }
  }
  crop_extent <- crop_extent %>% st_as_sf() %>% st_make_valid()
  save(crop_extent, file=paste0(file_path, "simplified_", file, ".Rdata"))
}

#################################################################################################
# here is an example of how it might be run with Haiti
{ simplify_SHP_file(file="gadm40_HTI_0", file_path=paste0(SHP_DATA_FOLDER, "Haiti/"))}

