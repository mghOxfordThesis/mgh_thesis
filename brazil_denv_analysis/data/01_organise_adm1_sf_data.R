#################################################################################################
#
# We organise the polygon data collected from GADM (https://gadm.org/) for Brazil at the 
# state level. There is an incorrectly defined administrative boundary which has to be manually 
# corrected. 
# 
# There is one output: "bra1_sf_data.rds", which is a simple feature objects of the administrative b
# boundaries of the 27 states across Brazil. 
#################################################################################################

#################################################################################################
# Definition of required packages, directories and local data

# libraries
library(pacman)
p_load("tidyverse", "sf", "sfheaders", "magrittr")

# local data
bra1 <- readRDS("bra2_sf/gadm36_BRA_1_sf.rds") %>%
  transmute(state=NAME_1)

#################################################################################################
# Itapiranga, Amazonas is part of Santa Catarina - here we fix that

amazonas_geometry <- bra1[which(bra1$state=="Amazonas"), "geometry"]
amazonas_geometry <- sfheaders::sf_remove_holes(amazonas_geometry)
bra1[which(bra1$state=="Amazonas"), "geometry"] <- amazonas_geometry

santa_catarina_geometry <- bra1[which(bra1$state=="Santa Catarina"), "geometry"]
itapiranga_geometry <- st_intersection(amazonas_geometry, santa_catarina_geometry)
santa_catarina_geometry <- st_difference(santa_catarina_geometry, amazonas_geometry)
bra1[which(bra1$state=="Santa Catarina"), "geometry"] <- santa_catarina_geometry

#################################################################################################
# save correct simple feature object
bra1 %<>% mutate(geometry=st_make_valid(geometry))

saveRDS(bra1, "bra2_sf/bra1_sf_data.rds")

