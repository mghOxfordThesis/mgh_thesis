#################################################################################################
#
#  This script organises the climate data into a format that is suitable for Index P estimation. 
#  The main output is a list of data frames which contains the climate variable time series 
#  for pixel in the region of interest. 
# 
#  This script is run after:
#  1. 01_data2rasters.R
#  2. 02_curate_calculate.R
#  3. 03_do_relative_humidity_range.R
#  4. 04_plotting_climate_maps.R
#
#################################################################################################

#################################################################################################
print("<MVSE> organizing data as MVSE expects it ...")

load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_temp_list.Rdata"))      # temp_rast_list
load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_humrel_list.Rdata"))    # humrel_rast_list
load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_prec_list.Rdata"))      # prec_rast_list

load(paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag,"_", "dates.Rdata"))  ## dates
dates$date<- as.Date(dates$date, format="%Y-%m-%d")

get_cell_THP <- function(ri, ci) {
  cell_THP <- purrr::pmap(list(temp=temp_rast_list, humrel=humrel_rast_list, prec=prec_rast_list), 
                          function(temp, humrel, prec) c(temp[ri, ci], humrel[ri, ci], prec[ri, ci])) %>%
    do.call(rbind, .)
  out <- data.frame(T=cell_THP[, 1], H=cell_THP[, 2], R=cell_THP[, 3], 
                   year=format(dates$date,"%Y"), month=format(dates$date,"%m"), day=format(dates$date,"%d"), 
                   date=dates$date)
  return(out)
}
Nrows<- nrow(temp_rast_list[[1]]) #all lists same size
Ncols<- ncol(temp_rast_list[[1]]) #all lists same size
cell_indices <- expand.grid(row=1:Nrows, col=1:Ncols) %>% as.data.frame() %>% arrange(desc(-row), desc(-col))
list_cell_THP <- purrr::map2(cell_indices[, 1], cell_indices[, 2], get_cell_THP)
outFile <- paste0(MVSE_INPUT_DATA_FOLDER, region, tag, "_list_cell_THP.Rdata")
save(list_cell_THP, file=outFile)
