#################################################################################################
#
#  This script is aggregates the posterior distribution of Index P at the pixel level into 
#  a smaller number of more manageable files. The pixels are aggregated in a specific order 
#  so that the raster files can be correctly re-constructed. 
# 
#  This script is run after:
#  1. 00_simplify_SHP_files.R (optional)
#  2. 00_generate_MVSE_maps.R
#  3. 01_run_MVSE_per_cell.R
#  4. 01_run_MVSE_region.R
#
#################################################################################################

########################################################################################
cat("## Step 2: aggregate simulated P per geo-pixel ... \n")
step_2_time <- Sys.time()

# just to get raster matrix dims
filenameIn <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_dimRaster.Rdata")
load(filenameIn) # rdim
ori_l<- rdim[1]
ori_c<- rdim[2]
N <- ori_l*ori_c

########################################################################################
# aggregate all P simulations in order
chunks <- 250
f <- function(x, index) {
  fileout<- paste0(MVSE_OUTPUT_DATA_FOLDER, "temporary/", region, tag, "_", x, "_SIM", index, ".Rdata")
  load(file=fileout) # sim_p, sim_1, or sim_v0
  if (index=="P") return(sim_p)
  else if (index=="Q") return(sim_q)
  else return(sim_v0)
}

if (N <= chunks) {
  allPs <- purrr::map(1:N, function(x) f(x, "P"))
  #allQs <- purrr::map(1:N, function(x) f(x, "Q"))
  #allV0s <- purrr::map(1:N, function(x) f(x, "V0"))
  save(allPs, file=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_indexP_month_per_cell_parallel1.Rdata"))
  #save(allQs, file=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_Q_month_per_cell_parallel1.Rdata"))
  #save(allV0s, file=paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_V0_month_per_cell_parallel1.Rdata"))
} else {
  left_endpoints <- seq(1, N, chunks)
  right_endpoints <- c(seq(chunks, N, chunks), N)
  pb <- txtProgressBar(min=1, max = length(left_endpoints), style = 3)
  for (ii in seq_along(left_endpoints)) {
    le <- left_endpoints[ii]; re <- right_endpoints[ii]
    fileout1 <- paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_indexP_month_per_cell_parallel", ii, ".Rdata")
    # fileout2 <- paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_Q_month_per_cell_parallel", ii, ".Rdata")
    # fileout3 <- paste0(MVSE_OUTPUT_DATA_FOLDER, "raw_output/", region, tag, "_V0_month_per_cell_parallel", ii, ".Rdata")
    if (!overwrite && file.exists(fileout1) && file.exists(fileout2) && file.exists(fileout3)) {
      purrr::map(le:re, function(x) file.remove(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_", x, "_SIMP.Rdata")))
      # purrr::map(le:re, function(x) file.remove(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_", x, "_SIMQ.Rdata")))
      # purrr::map(le:re, function(x) file.remove(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_", x, "_SIMV0.Rdata")))
      setTxtProgressBar(pb, ii)
      next
    }
    allPs <- purrr::map(le:re, function(x) f(x, "P"))
    # allQs <- purrr::map(le:re, function(x) f(x, "Q"))
    # allV0s <- purrr::map(le:re, function(x) f(x, "V0"))
    save(allPs, file=fileout1)
    # save(allQs, file=fileout2)
    # save(allV0s, file=fileout3)
    purrr::map(le:re, function(x) file.remove(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_", x, "_SIMP.Rdata")))
    # purrr::map(le:re, function(x) file.remove(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_", x, "_SIMQ.Rdata")))
    # purrr::map(le:re, function(x) file.remove(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_", x, "_SIMV0.Rdata")))
    setTxtProgressBar(pb, ii)
  }
}

########################################################################################
# delete the temporary files
cat("## deleting temporary data... \n")
unlink(paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary"), recursive=TRUE)

step_2_time <- Sys.time() - step_2_time
cat("## step 2 time:", step_2_time, attributes(step_2_time)$units, "\n")
cat("\n")