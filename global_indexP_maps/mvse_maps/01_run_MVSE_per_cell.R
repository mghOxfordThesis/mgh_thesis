#################################################################################################
#
#  This script is estimates the posterior distribution of the Index P time series for each pixel. 
#  Parallel executation of the computation is performed if no_clus>1. 
# 
#  This script is run after:
#  1. 00_simplify_SHP_files.R (optional)
#  2. 00_generate_MVSE_maps.R
#
#################################################################################################

#################################################################################################
#' Copy arguments into env and re-bind any function's lexical scope to bindTargetEnv .
#' See http://winvector.github.io/Parallel/PExample.html for example use.
#' Used to send data along with a function in situations such as parallel execution
#' (when the global environment would not be available). Typically called within
#' a function that constructs the worker function to pass to the parallel processes
#' (so we have a nice lexical closure to work with).
#'
#' @param bindTargetEnv environment to bind to
#' @param objNames additional names to lookup in parent environment and bind
#' @param names of functions to NOT rebind the lexical environments of
bindToEnv <- function(bindTargetEnv=parent.frame(),objNames,doNotRebind=c()) {
  # Bind the values into environment
  # and switch any list_cell_THP[[ip]] functions to this environment!
  for(var in objNames) {
    val <- get(var,envir=parent.frame())
    if(is.function(val) && (!(var %in% doNotRebind))) {
      # replace function's lexical environment with our target (DANGEROUS)
      environment(val) <- bindTargetEnv
    }
    # assign object to target environment, only after any possible alteration
    assign(var,val,envir=bindTargetEnv)
  }
}

#################################################################################################
cat("## Step 1: run MCMC per geo-pixel ...\n")
step_1_time <- Sys.time()

#################################################################################################
fileInDim <- paste0(RAW_CLIMATE_RASTER_FOLDER, region, tag, "_dimRaster.Rdata")

# just to get raster matrix dims
load(fileInDim) # rdim
ori_l<- rdim[1]
ori_c<- rdim[2]
cat(paste("## data with",ori_l,"lines","and", ori_c, "columns \n"))

#################################################################################################
cat("## loading MVSE input data \n")

# start up a parallel cluster (tips http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)
if (no_clus > ori_c) { no_clus <- ori_c; print('corrected nclus to n cols of map') }
parallelCluster <- parallel::makeCluster(no_clus, outfile="")
junk <- clusterEvalQ(parallelCluster, library(MVSE)) # load MVSE on each NODE
cat("## "); print(parallelCluster)
ArrayObjNames<- c( "thisdata",
                   "thisdatanames",
                   "mvse_model_args",
                   "fitting_args",
                   "MVSE_OUTPUT_DATA_FOLDER",
                   "tag", 
                   "overwrite", 
                   "region")


# takes the init row, the final row and the number of clusters and 
# returns how the approx. equal allocation of the rows across the clusters
buildParcelsSeq <- function(init, final, clusters){
  # corner case when (final-init+1)==clusters
  if (length(init:final)==clusters) return(lapply(init:final, function(x) c(x, x)))
  
  # when clusters < (final-init+1)
  ap <- round(seq(init, final, length.out=clusters+1))
  f <- function(x) {
    if (x==(length(ap)-1)) return(c(ap[x], ap[x+1]))
    else return(c(ap[x], ap[x+1]-1))
  }
  allPar <- lapply(1:(length(ap)-1), f)
  return(allPar)
}

# here, make sure to name all variables and functions that will be used per 'thread'
mkWorker <- function() {
  bindToEnv(objNames=ArrayObjNames) # load variables / names needed on NODE
  function(parcel) {
    pb <- txtProgressBar(min=1, max = length(parcel), style = 3)
    for (x in parcel) {
      setTxtProgressBar(pb, which(parcel==x))
      if (!overwrite) {
        fileout1 <- paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_",x,"_SIMP.Rdata")
        #fileout2 <- paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_",x,"_SIMQ.Rdata")
        #fileout3 <- paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_",x,"_SIMV0.Rdata")
        #if (file.exists(fileout1) && file.exists(fileout2) && file.exists(fileout3)) next
        if (file.exists(fileout1)) next
      }
      
      climate_data <- thisdata[[which(thisdatanames==x)]]
      if (sum(is.na(climate_data$T))==nrow(climate_data)) { #no data available
        sim_p <- NA; sim_q <- NA; sim_v0 <- NA
      } else { # construct mvse model
        if (is.null(mvse_model_args$model_category) || (mvse_model_args$model_category!="aegypti")) {
          this_mvse_model <- mvse_model(priors=mvse_model_args$priors, climate_data=climate_data, warning=FALSE)
        } else {
          this_mvse_model <- mvse_model(model_category="aegypti", climate_data=climate_data, warning=FALSE)
        }
        
        # run fitting procedure
        this_mvse_fit <- fitting(this_mvse_model, iter=fitting_args$iter, warmup=fitting_args$warmup, 
                                 seed=fitting_args$seed, init=fitting_args$init, 
                                 gauJump=fitting_args$gauJump, samples=fitting_args$samples, 
                                 verbose=FALSE)
        sim_p <- MVSE::extract(this_mvse_fit, pars="indexP")$indexP
        # sim_q <- MVSE::extract(this_mvse_fit, pars="Q")$Q
        # sim_v0 <- MVSE::extract(this_mvse_fit, pars="V0")$V0
        rm(this_mvse_model)
        rm(this_mvse_fit)
        rm(climate_data)
      }
      fileout <- paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_",x,"_SIMP.Rdata")
      save(sim_p, file=fileout)
      # fileout <- paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_",x,"_SIMQ.Rdata")
      # save(sim_q, file=fileout)
      # fileout <- paste0(MVSE_OUTPUT_DATA_FOLDER,"temporary/", region, tag,"_",x,"_SIMV0.Rdata")
      # save(sim_v0, file=fileout)
      rm(sim_p)
      # rm(sim_q)
      # rm(sim_v0)
    }
    close(pb)
  }
}

# estimate index P for each pixel over time
cat('## running mvse fitting for each pixel ... \n')
chunks <- 20000
N <- ori_l*ori_c
if (N <= chunks) {
  thisdata <- list_cell_THP
  thisdatanames <- 1:N
  clusterParcels <- split(1:N, 1:no_clus)
  parallel::parLapply(parallelCluster, X=clusterParcels, fun=mkWorker()) 
  gc(verbose=FALSE)
} else {
  left_endpoints <- seq(1, N, chunks)
  right_endpoints <- c(seq(chunks, N, chunks), N)
  for (ii in seq_along(left_endpoints)) {
    str <- left_endpoints[ii]
    end <- right_endpoints[ii]
    thisdata <- list_cell_THP[str:end]
    thisdatanames <- str:end
    clusterParcels <- split(str:end, 1:no_clus)
    parallel::parLapply(parallelCluster, X=clusterParcels, fun=mkWorker()) 
    gc(verbose=FALSE)
  }
}

stopCluster(parallelCluster)

step_1_time <- Sys.time() - step_1_time
cat("## step 1 time:", step_1_time, attributes(step_1_time)$units, "\n")
cat("\n")
