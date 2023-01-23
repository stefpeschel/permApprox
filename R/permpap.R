#' @title Approximate permutation p-values
#'
#' @description
#'   Approximate the p-value for a one-sided permutation hypothesis test,
#'   where the p-value is defined as the probability of observing a test
#'   statistic (in absolute terms) greater than or equal to the observed value
#'   under the null. Approximation methods include the empirical p-value,
#'   Gamma approximation, and GPD tail approximation.
#'
#' @param tObs numeric value giving the observed test statistic
#' @param tPerm numeric vector with permutation test statistics
#' @param includeObs logical. Indicates whether the observed test statistic
#'   should be included in the permutation distribution.
#' @param fitThresh threshold for initial p-value, above which a GPD or Gamma
#'   distribution is fitted for p-value approximation (e.g. 2 times alpha).
#'   Defaults is 0.1.
#' @param nonZero logical. If \code{TRUE} (default), a p-value of zero is
#'   strictly avoided.
#' @param method character indicating the method used for p-value estimation.
#'   Possible values are: "gpd" (default), "gamma", "empirical".
#'
#'
#' @param threshMethod method for threshold detection.
#'   Possible values are:
#'   "fix", "ftr", "minPR", "PRbelowAlpha", "fwdStop"
#' @param threshMethodPar fix threshold/number of exceedances or start value for
#'   threshold / number of exceedances, depending on threshMethod
#' @param nExceedMin minimum number of exceedances

#' @param stepSize numeric giving the value by which either the threshold or the
#'   number of exceedances is increased in each step, depending on threshMethod
#' @param ffwd if TRUE, large steps are made at the beginning until H0 is accepted
#'   for the first time

#' @param fitMethod ...
#' @param fitInformation ...
#' @param optimMethod ...
#' @param shape0,scale0 start value of the shape and scale parameter
#' @param shapeMin,scaleMin minimum value of the shape and scale parameter
#' @param shapeMax,scaleMax maximum value of the shape and scale parameter
#' @param gofTest goodness-of-fit test method
#' @param gofAlpha significance level of the GOF test
#' @param seed ...
#' @param cores
#' @param plotHist
#' @param histBreaks
#' @param histXlim
#' @param gpdEstimate
#' ...
#'
#' @export
#' @import foreach

# tPerm  <- gpddata
# tObs  <- max(tPerm)-0.01
# method = "gpd"
# threshVec = NULL
# threshMethod = "PRbelowAlpha"
# thresh0 = NULL
# nExceed0 = NULL
# nExceedMin = 1
# fitThresh = 0.1
# stepSize = 1
# includeObs = FALSE
# fitMethod = "ml"
# fitInformation = "observed"
# optimMethod = NULL
# shape0 = NULL
# scale0 = NULL
# shapeMin = -Inf
# scaleMin = -Inf
# shapeMax = Inf
# scaleMax = Inf
# gofTest = "ad"
# gofAlpha = 0.05
# gofTailRMMeth = "allrej"
# gofTailRMPar = NULL
# seed = NULL
# cores = 4L
# plotHist = TRUE
# histBreaks = 100
# histXlim = NULL
# histYlim = NULL
# plotPvals = TRUE
# plotPvalsxvar = "thresh"
# plotTitle = ""
# maintext = ""
# verbose = TRUE
# gpdEstimate = NULL


permpap <- function(tPerm,
                    tObs,
                    useAllPerm = FALSE,
                    method = "gpd",
                    fitThresh = 0.2,
                    includeObs = TRUE,
                    fitMethod = "MLE1D",
                    constraint = "tObs",
                    tol = 1e-8,
                    eps = "quantile",
                    epsVal = 0.9,
                    threshMethod = "PRbelowAlpha",
                    thresh0 = NULL,
                    nExceed0 = NULL,
                    nExceedMin = 10,
                    stepSize = 1,
                    gofTest = "ad",
                    gofAlpha = 0.05,
                    gofTailRMMeth = "allrej",
                    gofTailRMPar = NULL,
                    seed = NULL,
                    cores = 1L,
                    verbose = FALSE,
                    gpdEstimate = NULL,
                    plotHist = F,
                    histBreaks = 100,
                    histXlim = NULL,
                    histYlim = NULL,
                    plotPvals = TRUE,
                    plotPvalsxvar = "thresh",
                    plotTitle = "",
                    maintext = "",
                    jpoint = 0.2,
                    usegrid = "none",
                    alpha = 0.05,
                    multAdj = "adaptBH",
                    trueNullMethod = "convest",
                    pTrueNull = NULL,
                    ...) {

  # possible threshold detection methods are:
  # - "minimum proportion of rejected GOF-tests"
  # - "failure-to-reject"
  # threshMethod <- match.arg(threshMethod, choices = c("fix",
  #                                                     "ftr",
  #                                                     "minPR",
  #                                                     "PRbelowAlpha",
  #                                                     "fwdStop"))

  method <- match.arg(method, choices = c("gpd", "gamma", "empirical"))

  fitMethod <- match.arg(fitMethod, choices = c("LME", "MLE1D", "MLE2D", "MOM",
                                                "NLS2", "WNLLSM", "ZSE"))

  constraint <- match.arg(constraint, choices = c("none", "shapePos", "tObs",
                                                  "tObsMax"))
  #
  # gofTailRMMeth <- match.arg(gofTailRMMeth, choices = c("none", "num", "prop",
  #                                                       "allrej", "changepoint"))

  if (constraint == "shapePos" && !fitMethod %in% c("MLE1D", "MLE2D", "NLS2")) {
    stop("Constraint \"shapePos\" only available for methods ",
         "MLE1D, MLE2D, and NLS2.")
  }

  if (length(tObs) == 1) {
    tPerm <- matrix(tPerm, ncol = length(tPerm), nrow = 1)
  }

  # Number of permutations
  nPerm <- ncol(tPerm)

  #  Overall number of permutation test statistics
  ntPerm <- length(tPerm)

  # Number of tests
  nTest <- length(tObs)

  # Take absolute values (we test H1: tObs > t)
  tPerm <- abs(tPerm)
  tObs <- abs(tObs)

  pemp <- get_pvals_emp(tObs = tObs, tPerm = tPerm, nTest = nTest, nPerm = nPerm,
                        ntPerm = ntPerm, useAllPerm = useAllPerm)
  pvals <- pemp$pvals

  if (includeObs) tPerm <- rbind(tObs, tPerm)

  #-----------------------------------------------------------------------------
  # Empirical p-value(s) (observed test statistic is always included)


  if (method == "empirical") {

    return(list(pvals = pvals))

  } else if (all(pvals > fitThresh)) { # no distribution is fitted

    # out <- list(pvals = pvals)

    # if (method == "gpd") {
    #   out$shape <- out$scale <- out$gofPval <- out$thresh <- out$excessPerm <-
    #     out$excessObs <- out$nExceed <- out$shapeVec <- out$scaleVec <-
    #     out$gofPvalVec <- out$nExceedVec <- NA
    # } else if (method == "gamma") {
    #   out$shape <- out$rate <- out$gofPval <- NA
    # }

    return(list(pvals = pvals))
  }

  pvalsEmp <- pvals

  #-----------------------------------------------------------------------------
  # Fitting a Gamma distribution

  if (method == "gamma") {

    output <- get_pvals_gamma(pvals = pvals,
                              tPerm = tPerm,
                              tObs = tObs,
                              nTest = nTest,
                              fitThresh = fitThresh)

    output$pvalsEmp <- pvalsEmp

    return(output)
  }

  #-----------------------------------------------------------------------------
  # Tail approximation using the GPD

  args <- as.list(environment())

  res <- do.call(get_pvals_gpd, args)

  pvals <- res$pvals

  #-----------------------------------------------------------------------------
  # Multiple testing adjustment

  if (multAdj == "YB") {

    pAdj <- NULL

  } else if (multAdj == "Lee") {

    args <- as.list(environment())
    args$res <- NULL
    args$cores <- 1

    #---------------------------------------------------------------------------
    # Initialize parallel stuff

    if (verbose) {
      # Create progress bar:
      pb <- utils::txtProgressBar(0, nPerm, style=3)

      # Function for progress bar
      progress <- function(n) {
        utils::setTxtProgressBar(pb, n)
      }
    }

    if (cores > 1) {
      if (parallel::detectCores() < cores) cores <- parallel::detectCores()

      cl <- parallel::makeCluster(cores, outfile = "")
      doSNOW::registerDoSNOW(cl)
      '%do_or_dopar%' <- get('%dopar%')

    } else {
      '%do_or_dopar%' <- get('%do%')
    }

    if (verbose) {
      opts <- list(progress = progress)
    } else {
      opts <- list()
    }

    #---------------------------------------------------------------------------
    # Compute p-values

    loopres <-
      foreach(i = 1:nPerm,
              .export = c("get_gpd_thresh", "fit_gpd", "get_thresh_idx",
                          "get_pvals_gpd", "get_pvals_emp", ".est_gpd_params",
                          "gpdAd_adapt", "gpdCvm_adapt",
                          "gpd_LME", "gpd_MLE1D", ".MLE1D_fk", ".MLE1D_fp",
                          "gpd_MLE2D", ".MLE2D_negloglik", "gpd_MOM",
                          "gpd_NLS2", ".NLS2_gpdf", ".NLS2_ecdf", ".NLS2_gpdf2",
                          ".NLS2_ecdf2", ".NLS2_rss1", ".NLS2_rss2",
                          "gpd_WNLLSM", ".WNLLSM_sum_i", ".WNLLSM_WLLS1",
                          ".WNLLSM_WLLS", ".WNLSM_WLS1", ".WNLSM_WLS",
                          "gpd_ZSE", ".ZSE_lx"),
              .packages = c("permpap", "foreach"),
              .options.snow = opts) %do_or_dopar% {


                if (verbose) progress(i)

                args$tObs <- tPerm[, i]
                args$tPerm <- cbind(tObs, tPerm[, -i])

                res <- do.call(get_pvals_gpd, args)
                res$pvals
              }

    if (verbose) {
      # Close progress bar
      close(pb)
    }

    # Stop cluster
    if (cores > 1) parallel::stopCluster(cl)

    #---------------------------------------------------------------------------

    pPerm <- matrix(unlist(loopres), nrow = nTest, ncol = nPerm)

    minP <- apply(pPerm, 2, min)

    # pAdjust would be 0 or 1 with the original method
    pAdj <- sapply(1:nTest, function(i) sum(minP <= pvals[i]) / nPerm)

    #descdist(minP, discrete = FALSE)

    minPfit <- fitdistrplus::fitdist(data = minP, distr = "beta")
    minPfit$estimate
    #plot(minPfit)

    pvalsAdj <- sapply(1:nTest, function(i)  pbeta(pvals[i],
                                               minPfit$estimate[1],
                                               minPfit$estimate[2]))

  } else if (multAdj == "none") {
    pvalsAdj <- NULL

  } else {
    pvalsAdj <- multAdjust(pvals,
                           adjust = multAdj,
                           trueNullMethod = trueNullMethod,
                           pTrueNull = pTrueNull,
                           verbose = verbose)
  }


  #-----------------------------------------------------------------------------
  callArgs <- mget(names(formals()),sys.frame(sys.nframe()))
  callArgs$gpdEstimate <- NULL

  output <- list(pvals = pvals,
                 pvalsAdj = pvalsAdj,
                 pvalsEmp = pvalsEmp
                 #fdr.output = fdr.output,
                 #sp = sp,
                 #fdr.pa = fdr.pa
                 )

  return(output)
}


