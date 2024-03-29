##########################################################################################################
#' Class \code{mvsefit}: fitted MVSE model
#'
#' An object of class \code{mvsefit} contains the output derived from fitting a
#' MVSE model as returned by the method \link{sampling}.
#'
#' @slot model_name The model name as a string.
#' @slot model_category The model category as a string.
#' @slot mvsemodel The instance of S4 class \link{mvsemodel}.
#' @slot mvse_args A list containing the arguments used for sampling (e.g. \code{iter}, \code{warmup}, \code{init}, \code{gauJumps}, etc.).
#' @slot sim A list containing simulation results including the posterior draws for the factors (rho, eta, alpha) and the generated
#' quantities (indexP, Q, V0) as well as various pieces of data used by many of the methods for \code{mvsefit} objects.
#'
#' @section Methods:
#' \subsection{Printing, plotting and summarizing:}{
#' \describe{\item{\code{show}}{Print the default summary for the model.}}}
#' \subsection{Extracting posterior draws:}{
#' \describe{\item{\code{show}}{Print the default summary for the model.}}}
#' \subsection{Simulating:}{
#' \describe{\item{\code{simulate}}{Simulate Index P, Q and V0 using the fitted model.}}}
#'
#'
#' @name mvsefit-class
#' @rdname mvsefit
#' @import ggplot2 dplyr tidyr RColorBrewer
#' @importFrom bayesplot bayesplot_grid
#' @importFrom coda spectrum0.ar
#' @include utils.R mvsemodel.R
setClass("mvsefit",
         slots = c(
           model_name = "character",
           model_category = "character",
           mvsemodel = "mvsemodel",
           mvse_args = "list",
           sim = "list"
         )
)

##########################################################################################################
#' Construct a MVSE fit
#'
#' Construct an instance of S4 class \code{mvsefit}.
#'
#' @param model_name A character string naming the model; defaults to \code{"anon_model"}.
#' @param model_category A character string naming the model category. If \code{"aegypti"} is specified,
#' the prior distributions of ecological/epidemiological human/vector parameters will automatically be set.
#' Otherwise, it defaults to \code{user-defined} and the user must specified the prior distributions.
#' @param mvsemodel An object of class \code{\linkS4class{mvsemodel}} used to generate the fit.
#' @param mvse_args A list containing the arguments used for sampling
#' (e.g. \code{iter}, \code{warmup}, \code{init}, \code{gauJumps}, etc.).
#' @param sim A list containing simulation results including the posterior draws for the factors (rho, eta, alpha) and the generated
#' quantities (indexP, Q, V0) as well as various pieces of data used by many of the methods for \code{mvsefit} objects.
#'
#' @return An instance of S4 class \code{\linkS4class{mvsefit}}.
#' @noRd
#' @export
mvse_fit <- function(model_name="anon_model", model_category="user-defined", mvsemodel, mvse_args, sim) {
  new("mvsefit", model_name=model_name, model_category=model_category, mvsemodel=mvsemodel, mvse_args=mvse_args, sim=sim)
}

##########################################################################################################
#' Print the default summary for the model
#'
#' @param object An object of class \code{\linkS4class{mvsefit}}.
#'
#' @return NULL
#' @noRd
#' @export
setMethod("show", "mvsefit",
          function(object) {
            cat("## Class:", class(object), "\n")
            cat("## Model name:", object@model_name, "\n")
            cat("## Model category:", object@model_category, "\n")
            cat(paste("## iter=", object@mvse_args$iter, "; warmup=", object@mvse_args$warmup, sep=""), "\n")
            cat(paste("## post-warmup iter=", object@mvse_args$iter*object@mvse_args$warmup, "; accepted draws=", object@sim$accepted, sep=""), "\n")
            if (object@model_category=="aegypti") vars <- c("rho", "eta", "indexP", "Q", "V0")
            else vars <- c("rho", "eta", "alpha", "indexP", "Q", "V0")
            draws <- object@sim[vars];
            all_gen_pars <- c("indexP", "Q", "V0")
            draws[all_gen_pars] <- lapply(draws[all_gen_pars], function(x) as.vector(as.matrix(x[, -1])))
            compute_summary_stats <- function(x) {
              draw <- draws[[x]]
              mean <- mean(draw)
              sd <- sd(draw)
              p <- quantile(draw, c(0.25, 0.5, 0.75))
              if (x %in% all_gen_pars) {
                n_eff <- NA; se_mean <- NA
              }
              else {
                n_eff <- length(draw)*sd^2/spectrum0.ar(draw)$spec^2
                se_mean <- sd/sqrt(n_eff)
              }
              return(round(c(mean=mean, se_mean=se_mean, sd=sd, p, n_eff=n_eff), 3))
            }
            summary_stats <- t(sapply(vars, function(x) compute_summary_stats(x)))
            print(summary_stats)
          }
)

##########################################################################################################
#' Print a customizable summary for the model.
#'
#' @param x An object of class \code{\linkS4class{mvsefit}}.
#' @param pars A character vector of parameter names. The default is all parameters
#' (i.e. \code{"rho"}, \code{"eta"}, \code{"alpha"}, \code{"indexP"}, \code{"Q"}, \code{"V0"}).
#' @param probs A numeric vector of quantiles of interest. The default is
#' \code{c(0.25, 0.5, 0.75)}.
#' @param digits_summary The number of decimal places to use when printing the summary, defaulting to 3.
#' @param ... Additional arguments passed to the \code{print} method for \code{mvsefit} objects.
#'
#' @return NULL
#' @export
setMethod("print", "mvsefit",
          function(x, pars=NULL, probs=c(0.25, 0.5, 0.75), digits_summary=3, ...) {
            cat("## Class:", class(x), "\n")
            cat("## Model name:", x@model_name, "\n")
            cat("## Model category:", x@model_category, "\n")
            cat(paste("## iter=", x@mvse_args$iter, "; warmup=", x@mvse_args$warmup, sep=""), "\n")
            cat(paste("## post-warmup iter=", x@mvse_args$iter*x@mvse_args$warmup, "; accepted draws=", x@sim$accepted, sep=""), "\n")
            if (is.null(pars)) {
              if (x@model_category=="aegypti") pars <- c("rho", "eta", "indexP", "Q", "V0")
              else pars <- c("rho", "eta", "alpha", "indexP", "Q", "V0")
            }
            draws <- x@sim[pars];
            all_gen_pars <- c("indexP", "Q", "V0")
            gen_pars <- all_gen_pars[which(all_gen_pars %in% pars)]
            if (length(gen_pars)>0) {
              draws[gen_pars] <- lapply(draws[gen_pars], function(x) as.vector(as.matrix(x[, -1])))
            }
            compute_summary_stats <- function(x) {
              draw <- draws[[x]]
              mean <- mean(draw)
              sd <- sd(draw)
              p <- quantile(draw, probs)
              if (x %in% all_gen_pars) {
                n_eff <- NA; se_mean <- NA
              }
              else {
                n_eff <- length(draw)*sd^2/spectrum0.ar(draw)$spec^2
                se_mean <- sd/sqrt(n_eff)
              }
              return(round(c(mean=mean, se_mean=se_mean, sd=sd, p, n_eff=n_eff), digits_summary))
            }
            summary_stats <- t(sapply(pars, function(x) compute_summary_stats(x)))
            print(summary_stats)
          }
)

##########################################################################################################
#' Extract samples from a fitted MVSE model
#'
#' Extract samples from a fitted model represented by an instance of class \code{\linkS4class{mvsefit}}.
#'
#' @param x An object of class \code{\linkS4class{mvsefit}}.
#' @param pars A character vector of parameter names. The default is all parameters
#' (i.e. \code{"rho"}, \code{"eta"}, \code{"alpha"}, \code{"indexP"}, \code{"Q"}, \code{"V0"}).
#'
#' @return A named list, every element of which is either a numeric vector or dataframe representing samples
#' for a parameter.
#' 
#' @examples 
#' # obtain climate data
#' data(climateFSA)
#' 
#' # define a mvse model 
#' priors <- list()
#' priors$mosq_life_exp <- list(dist="normal", pars=c(mean=12, sd=2)) # mosquito life expectancy (days)
#' priors$mosq_inc_per <- list(dist="normal", pars=c(mean=7, sd=2))   # mosquito-virus incubation period (days)
#' priors$mosq_biting_freq <- list(dist="normal", pars=c(mean=0.25, 0.01)) # mosquito biting frequency (bites/female/day)
#' priors$human_life_exp <- list(dist="normal", pars=c(mean=71.1, sd=2)) # human life expectancy (years)
#' priors$human_inc_per <- list(dist="normal", pars=c(mean=5.8, sd=1)) # human-virus incubation period (days)
#' priors$human_inf_per <- list(dist="normal", pars=c(mean=5.9, sd=1)) # human-virus infectious period (days)
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=priors)
#' 
#' # run the MCMC sampling procedure to estimate the epi-entomological parameters as well Index P, Q and V0
#' user_fit <- fitting(object=user_model, iter=10^5, warmup=0.2, seed=123, init=c(rho=1, eta=10, alpha=3), samples=1000)
#' 
#' # extract the sampled parameters
#' all_output <- extract(user_fit) # all parameters
#' indexP_output <- extract(user_fit, pars="indexP") # only Index P
#' 
#' @usage NULL
#' @name extract
#' @export
setGeneric("extract",
           def=function(object, pars=NULL)
           {standardGeneric("extract")})

#' @rdname extract
setMethod("extract", "mvsefit",
          function(object, pars=NULL) {
            if (is.null(pars)) {
              if (object@model_category=="aegypti") pars <- c("rho", "eta", "indexP", "Q", "V0")
              else pars <- c("rho", "eta", "alpha", "indexP", "Q", "V0")
            }
            return(object@sim[pars])
          }
)

##########################################################################################################
#' Plot a pairs plot of the MCMC draws of the factors rho, eta and alpha
#'
#' A square plot matrix with univariate marginal distributions along the diagonal (as histograms) and
#' bivariate distributions off the diagonal (as scatterplots).
#'
#' @param object An object of S4 class \code{\linkS4class{mvsefit}}.
#' @param filename The name of the file where the plot in PNG format is saved (optional)..
#' @param options A named list of options for the plot.
#'
#' @return Many ggplot objects organized into a grid via \code{\link{bayesplot_grid}} is returned.
#' 
#' @examples
#' # obtain climate data
#' data(climateFSA)
#' 
#' # define a mvse model 
#' priors <- list()
#' priors$mosq_life_exp <- list(dist="normal", pars=c(mean=12, sd=2)) # mosquito life expectancy (days)
#' priors$mosq_inc_per <- list(dist="normal", pars=c(mean=7, sd=2))   # mosquito-virus incubation period (days)
#' priors$mosq_biting_freq <- list(dist="normal", pars=c(mean=0.25, 0.01)) # mosquito biting frequency (bites/female/day)
#' priors$human_life_exp <- list(dist="normal", pars=c(mean=71.1, sd=2)) # human life expectancy (years)
#' priors$human_inc_per <- list(dist="normal", pars=c(mean=5.8, sd=1)) # human-virus incubation period (days)
#' priors$human_inf_per <- list(dist="normal", pars=c(mean=5.9, sd=1)) # human-virus infectious period (days)
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=priors)
#' 
#' # run the MCMC sampling procedure to estimate the epi-entomological parameters as well Index P, Q and V0
#' user_fit <- fitting(object=user_model, iter=10^5, warmup=0.2, seed=123, init=c(rho=1, eta=10, alpha=3), samples=1000)
#' 
#' # plot the pairs plot for the MCMC draws of alpha, rho, eta
#' mcmc_pairs(user_fit, options=list(max_draws=10^5))
#' 
#' @usage NULL
#' @name mcmc_pairs
#' @export
setGeneric("mcmc_pairs",
           def=function(object, filename=NULL, options=list(max_draws=10^4))
           {standardGeneric("mcmc_pairs")})

#' @rdname mcmc_pairs
setMethod("mcmc_pairs", "mvsefit",
          function(object, filename=NULL, options=list(max_draws=10^4)) {
            if (object@model_category=="aegypti") factor_names <- c("rho", "eta")
            else factor_names <- c("rho", "eta", "alpha")
            factors <- MVSE::extract(object)[factor_names]
            plot_list <- as.list(rep(NA, length(factors)^2))
            pal <- brewer.pal(4, "Blues")
            for (ii in seq_along(factor_names)) {
              max_draws <- min(options$max_draws, length(factors[[ii]]))
              this_factor <- factors[[ii]][1:max_draws]
              binwidth <- (max(this_factor)-min(this_factor))/30
              p_df <- data.frame(x=this_factor)
              p <- ggplot(p_df) +
                theme_bw() +
                geom_histogram(aes(x=x, y=..count../sum(count)),
                               binwidth=binwidth, fill=pal[2], color=pal[3]) +
                labs(title=factor_names[ii]) +
                theme(axis.title = element_blank(), axis.line.y=element_blank(), axis.text.y = element_blank(),
                      axis.ticks.y=element_blank(), panel.background=element_blank(), panel.border=element_blank(),
                      panel.grid=element_blank(), plot.title=element_text(hjust=0.5))
              plot_list[[ii+length(factors)*(ii-1)]] <- p
            }

            for (ii in seq_along(factor_names)) {
              for (jj in seq_along(factor_names)) {
                if (ii==jj) next
                max_draws <- min(options$max_draws, length(factors[[jj]]))
                x_values <- factors[[jj]][1:max_draws]
                max_draws <- min(options$max_draws, length(factors[[ii]]))
                y_values <- factors[[ii]][1:max_draws]
                p_df <- data.frame(x=x_values, y=y_values)
                p <- ggplot(data=p_df) +
                  theme_bw() +
                  geom_point(aes(x=x, y=y), color=pal[4], fill=pal[3], shape=21, size=0.5) +
                  theme(axis.title=element_blank(), axis.line=element_line(color="black"), panel.border=element_blank(),
                        panel.background=element_blank(), panel.grid=element_blank())
                plot_list[[jj+length(factors)*(ii-1)]] <- p
              }
            }
            p <- bayesplot_grid(plots=plot_list, grid_args=list(nrow=length(factors)))
            if (!is.null(filename)) {
              png(filename, width = 1024, height = 768, res=150)
              print(p)
              dev.off()
            }
            return(p)
          }
)

##########################################################################################################
#' Plot trace plots of the MCMC draws of the factors rho, eta and alpha
#'
#' Displays a separate plot of iterations vs. sampled values for each of the factors.
#'
#' @param object An object of S4 class \code{\linkS4class{mvsefit}}.
#' @param filename The name of the file where the plot in PNG format is saved (optional).
#' @param options A named list of options for the plot.
#'
#' @return Many ggplot objects organized into a grid via \code{\link{bayesplot_grid}} is returned.
#' 
#' @examples
#' # obtain climate data
#' data(climateFSA)
#' 
#' # define a mvse model 
#' aegypti_model <- mvse_model(model_name="my_aegypti_mvse_model", model_category="aegypti", climate_data=climateFSA)
#' 
#' # run the MCMC sampling procedure to estimate the epi-entomological parameters as well Index P, Q and V0
#' aegypti_fit <- fitting(object=aegypti_model, iter=10^5, warmup=0.2, seed=123, init=c(rho=1, eta=10, alpha=3), samples=1000)
#' 
#' mcmc_traceplot(aegypti_fit, options=list(max_iter=10^5))
#' @usage NULL
#' @name mcmc_traceplot
#' @export
setGeneric("mcmc_traceplot",
           def=function(object, filename=NULL, options=list(max_iter=10^4))
           {standardGeneric("mcmc_traceplot")})

#' @rdname mcmc_traceplot
setMethod("mcmc_traceplot", "mvsefit",
          function(object, filename=NULL, options=list(max_iter=10^4)) {
            if (object@model_category=="aegypti") factor_names <- c("rho", "eta")
            else factor_names <- c("rho", "eta", "alpha")
            factors <- MVSE::extract(object)[factor_names]

            plot_list <- as.list(rep(NA, length(factors)))
            pal <- brewer.pal(4, "Blues")
            for (ii in seq_along(factor_names)) {
              max_iter <- min(options$max_iter, length(factors[[ii]]))
              this_factor <- factors[[ii]]
              this_factor <- this_factor[(length(this_factor)-max_iter+1):length(this_factor)]
              p_df <- data.frame(x=(length(this_factor)-max_iter+1):length(this_factor), y=this_factor)
              p <- ggplot(p_df) +
                theme_bw() +
                geom_line(aes(x=x, y=y), color=pal[3]) +
                labs(title=factor_names[ii]) +
                theme(axis.title = element_blank(), axis.line=element_line(color="black"),
                      panel.background=element_blank(), panel.border=element_blank(),
                      panel.grid=element_blank(), plot.title=element_text(hjust=0.5))
              plot_list[[ii]] <- p
            }
            p <- bayesplot_grid(plots=plot_list, grid_args=list(nrow=1, bottom="post-warmup iteration"))
            if (!is.null(filename)) {
              png(filename, width = 1024, height = 768, res=150)
              print(p)
              dev.off()
            }
            return(p)
          }
)

##########################################################################################################
#' Plot the posterior distributions of the factors rho, eta and alpha
#'
#' Displays a separate plot (i.e. histogram) for the posterior distribution of each of the factors.
#'
#' @param x An object of S4 class \code{\linkS4class{mvsefit}}.
#' @param factors If not NULL, a character vector indicating which factors (i.e. \code{"rho"},
#' \code{"eta"} and \code{"alpha"}) to include in the plots. By default, all factors are included.
#' @param filename The name of the file where the plot in PNG format is saved.
#'
#' @return Many ggplot objects organized into a grid via \code{\link{bayesplot_grid}} is returned.
#' 
#' @examples
#' # obtain climate data
#' data(climateFSA)
#' 
#' # define a mvse model 
#' aegypti_model <- mvse_model(model_name="my_aegypti_mvse_model", model_category="aegypti", climate_data=climateFSA)
#' 
#' # run the MCMC sampling procedure to estimate the epi-entomological parameters as well Index P, Q and V0
#' aegypti_fit <- fitting(object=aegypti_model, iter=10^5, warmup=0.2, seed=123, init=c(rho=1, eta=10, alpha=3), samples=1000)
#' 
#' # plot the posterior distributions of the factors
#' mcmc_factor_dist(aegypti_fit)
#' mcmc_factor_dist(aegypti_fit, factors=c("eta"))
#' @usage NULL
#' @name mcmc_factor_dist
#' @export
setGeneric("mcmc_factor_dist",
           def=function(object, factors=NULL, filename=NULL)
           {standardGeneric("mcmc_factor_dist")})

#' @rdname mcmc_factor_dist
setMethod("mcmc_factor_dist", "mvsefit",
          function(object, factors=NULL, filename=NULL) {
            if (object@model_category=="aegypti") factor_names <- c("rho", "eta")
            else factor_names <- c("rho", "eta", "alpha")
            if (!is.null(factors)) this_factors <- MVSE::extract(object)[factors]
            else this_factors <- MVSE::extract(object)[factor_names]; factors <- factor_names

            plot_list <- as.list(rep(NA, length(this_factors)))
            pal <- brewer.pal(4, "Blues")
            ii <- 1
            for (this_factor_name in factors) {
              this_factor <- this_factors[[this_factor_name]]
              binwidth <- (max(this_factor)-min(this_factor))/30
              p_df <- data.frame(x=this_factor)
              p <- ggplot(p_df) +
                theme_bw() +
                geom_histogram(aes(x=x, y=..count../sum(count)),
                               binwidth=binwidth, fill=pal[2], color=pal[3]) +
                labs(title=this_factor_name) +
                theme(axis.title = element_blank(), axis.line=element_line(color="black"),
                      panel.background=element_blank(), panel.border=element_blank(),
                      panel.grid=element_blank(), plot.title=element_text(hjust=0.5))
              plot_list[[ii]] <- p
              ii <- ii + 1
            }
            p <- bayesplot_grid(plots=plot_list, grid_args=list(nrow=1))
            if (!is.null(filename)) {
              png(filename, width = 1024, height = 768, res=150)
              print(p)
              dev.off()
            }
            return(p)
          }
)

##########################################################################################################
#' Plot a summary time series of a transmission potential index
#'
#' Plots a summary time series of one of the three transmission potential indices
#' (i.e. index P, Q, V0). The mean time series along with the 95% credible interval are plotted.
#'
#' @param x An object of S4 class \code{\linkS4class{mvsefit}}.
#' @param index A character string specifying the index to plot. Options include \code{"indexP"},
#' \code{Q} or \code{V0}. Defaults to \code{"indexP"}.
#' @param filename The name of the file where the plot in PNG format is saved.
#' @param options A named list of options for the plot.
#'
#' @return A ggplot object is returned.
#'
#' @examples
#' # obtain climate data
#' data(climateFSA)
#' 
#' # define a mvse model 
#' aegypti_model <- mvse_model(model_name="my_aegypti_mvse_model", model_category="aegypti", climate_data=climateFSA)
#' 
#' # run the MCMC sampling procedure to estimate the epi-entomological parameters as well Index P, Q and V0
#' aegypti_fit <- fitting(object=aegypti_model, iter=10^5, warmup=0.2, seed=123, init=c(rho=1, eta=10, alpha=3), samples=1000)
#' 
#' # plot the distribution of sampled Index P
#' mcmc_index_dist(aegypti_fit, index="indexP")
#' mcmc_index_dist(aegypti_fit, index="indexP", options=list(smoothing=7))
#' @usage NULL
#' @name mcmc_index_dist
#' @export
setGeneric("mcmc_index_dist",
           def=function(object, index="indexP", filename=NULL, options=list(max_samples=1000, smoothing=NULL))
           {standardGeneric("mcmc_index_dist")})

#' @rdname mcmc_index_dist
setMethod("mcmc_index_dist", "mvsefit",
          function(object, index="indexP", filename=NULL, options=list(smoothing=NULL)) {
            data <- object@sim[[index]]

            # non-smoothed index
            index_data <- as.matrix(select(data, -date))
            index_mean <- rowMeans(index_data)
            index_l95 <- apply(index_data, MARGIN=1, FUN= function(X){ quantile(X, probs=c(0.025), na.rm=T)})
            index_u95 <- apply(index_data, MARGIN=1, FUN= function(X){ quantile(X, probs=c(0.975), na.rm=T)})
            ind <- c("mean"="mean", "l95"="95CI", "u95"="95CI")
            index_df <- data.frame(date=data$date, mean=index_mean, l95=index_l95, u95=index_u95) %>%
              gather(key=percentile, value=value, 2:ncol(.)) %>%
              mutate(color=ind[percentile])

            # smoothed index P
            smoothing <- options$smoothing
            if (!is.null(smoothing)) {
              index_smooth <- matrix(NA, nrow=length(data$date), ncol=length(smoothing))
              for (ii in 1:length(smoothing)) {
                smooth_data <- apply(index_data, 2, function(x) .smoothUDSeries(x, smoothing[ii]))
                index_smooth_mean <- rowMeans(smooth_data)
                index_smooth[, ii] <- index_smooth_mean
              }
              index_smooth <- as.data.frame(index_smooth) %>% mutate(date=data$date) %>% select(date, dplyr::everything())
              index_smooth <- gather(index_smooth, key=variable, value=value, 2:ncol(indexP_smooth))
            }

            # ggplot object
            lab_dates <- seq(index_df$date[1], tail(index_df$date, 1),
                             by=paste0(round(tail(index_df$date, 1)-index_df$date[1])/5, " days"))
            p <- ggplot(data=index_df) +
              geom_line(aes(x=date, y=value, group=percentile, color=color, alpha=color)) +
              theme_bw() +
              labs(y=index) +
              scale_color_manual(breaks=c("mean", "95CI"), values=c("tomato3", "grey77"),
                                 labels=c("mean", "95% CI")) +
              scale_alpha_manual(breaks=c("mean", "95CI"), values=c(1, 0.5), labels=NULL, guide="none") +
              scale_x_date(breaks=lab_dates, date_labels="%b %y") +
              coord_cartesian(ylim=c(0, 10)) + 
              theme(legend.position = c(0.9, 0.9), legend.title=element_blank(),
                    legend.background = element_rect(color="black"))
            if (index=="indexP") p <- p + geom_hline(yintercept=1, color="black", alpha=0.5, linetype="dashed")
            if (!is.null(smoothing))
              p <- p + geom_line(data=index_smooth, aes(x=date, y=value, group=variable), color="black", size=1.2)
            if (!is.null(filename)) {
              png(filename, width = 1024, height = 768, res=150)
              print(p)
              dev.off()
            }
            return(p)
          }
)
