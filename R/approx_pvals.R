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


approx_pvals <- function(tPerm,
                         tObs,
                         method = "gpd",
                         fitThresh = 0.1,
                         includeObs = TRUE,
                         fitMethod = "MLE1D",
                         constraint = "tObs",
                         tol = 1e-8,
                         eps = NULL,
                         factor = 1,
                         finalOnly = FALSE,
                         nonZero = TRUE,
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
    tPerm <- as.matrix(tPerm)
  }

  if (includeObs) tPerm <- rbind(tObs, tPerm)

  # Number of permutations
  nPerm <- nrow(tPerm)

  # Number of tests
  nTest <- ncol(tPerm)

  # Take absolute values (we test H1: tObs > t)
  tPerm <- abs(tPerm)
  tObs <- abs(tObs)


  # Number of permutation test statistics being greater than or equal to tObs
  nlarger <- sapply(1:nTest, function(i) {
    sum(tPerm[, i] >= tObs[i])
  })

  #-----------------------------------------------------------------------------
  # Empirical p-value(s) (observed test statistic is always included)
  pval <- (nlarger + 1) / (nPerm + 1)

  if (method == "empirical") {

    return(list(pval = pval))

  } else if (all(pval > fitThresh)) {
    # no distribution is fitted

    # out <- list(pvals = pvals)

    # if (method == "gpd") {
    #   out$shape <- out$scale <- out$gofPval <- out$thresh <- out$excessPerm <-
    #     out$excessObs <- out$nExceed <- out$shapeVec <- out$scaleVec <-
    #     out$gofPvalVec <- out$nExceedVec <- NA
    # } else if (method == "gamma") {
    #   out$shape <- out$rate <- out$gofPval <- NA
    # }

    return(list(pval = pval))
  }

  # Indices of p-values above threshold
  idxFit <- which(pval <= fitThresh)

  #-----------------------------------------------------------------------------
  # Fitting a Gamma distribution

  if (method == "gamma") {

    shape <- rate <- gofPval <- rep(NA, nTest)

    # capture.output(gammafit <- fitdistrplus::fitdist(data = tPerm[, 1],
    #                                                  distr = "gamma",
    #                                                  method = "mle"),
    #                file='NULL')

    for (i in idxFit) {
      gammafit <- fitdistrplus::fitdist(data = tPerm[, i],
                                        distr = "gamma",
                                        method = "mle")

      shape[i] <- as.numeric(gammafit$estimate["shape"])
      rate[i] <- as.numeric(gammafit$estimate["rate"])

      pval[i] <- pgamma(q = tObs[i], shape = shape[i], rate = rate[i],
                        lower.tail = FALSE)

      cvmtest <- goftest::cvm.test(x = tPerm[, i], null = "gamma",
                                   shape = shape[i], rate = rate[i])
      gofPval[i] <- cvmtest$p.value
    }

    output <- list(pval = pval, shape = shape, rate = rate, gofPval = gofPval)

    return(output)
  }


  #=============================================================================
  # Tail approximation using the GPD

  # Maximum value at which the GPD density must be positive
  if (constraint == "tObs") {
    tMax <- tObs

  } else if (constraint == "tObsMax") {
    tMax <- rep(max(tObs), length(tObs))

  } else {
    tMax <- NULL
  }

  #-----------------------------------------------------------------------------
  # Initialize parallel stuff

  if (verbose) {
    # Create progress bar:
    pb <- utils::txtProgressBar(0, length(idxFit), style=3)

    # Function for progress bar
    progress <- function(n) {
      utils::setTxtProgressBar(pb, n)
    }
  }

  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
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

  loopres <-
    foreach(i = seq_along(idxFit),
            .export = c("get_gpd_thresh", "fit_gpd", "get_thresh_idx",
                        ".est_gpd_params",
                        "gpdAd_adapt", "gpdCvm_adapt",
                        "gpd_LME", "gpd_MLE1D", "gpd_MLE2D",
                        "gpd_MOM", "gpd_NLS2", "gpd_WNLLSM", "gpd_ZSE"),
            .packages = "permpap",
            .combine='comb', .multicombine=TRUE,
            .init=list(list(), list(), list(), list(), list(), list(),
                       list(), list(), list()),
            .options.snow = opts) %do_or_dopar% {

              if (verbose) progress(i)

              out <- list()

              threshList <- get_gpd_thresh(tPerm = tPerm[, idxFit[i]],
                                           tObs = tObs[idxFit[i]],
                                           tMax = tMax,
                                           tol = tol,
                                           eps = eps,
                                           threshVec = threshVec,
                                           threshMethod = threshMethod,
                                           thresh0 = thresh0,
                                           nExceed0 = nExceed0,
                                           nExceedMin = nExceedMin,
                                           stepSize = stepSize,
                                           includeObs = includeObs,
                                           fitMethod = fitMethod,
                                           constraint = constraint,
                                           fitInformation = fitInformation,
                                           optimMethod = optimMethod,
                                           shape0 = shape0,
                                           scale0 = scale0,
                                           shapeMin = shapeMin,
                                           scaleMin = scaleMin,
                                           shapeMax = shapeMax,
                                           scaleMax = scaleMax,
                                           gofTest = gofTest,
                                           gofAlpha = gofAlpha,
                                           gofTailRMMeth = gofTailRMMeth,
                                           gofTailRMPar = gofTailRMPar,
                                           seed = seed,
                                           cores = 1,
                                           verbose = F)

              thresh <- threshList$thresh
              out$thresh <- thresh
              out$nExceed <- threshList$nExceed
              out$fitFailed <- FALSE
              out$zeroRepl <- FALSE

              if (is.na(thresh)) {

                if (verbose) {
                  message("In test ", idxFit[i],
                          ": No iteration led to a good fit ",
                          "(GOF test rejected for ",
                          "all thresholds). Empirical p-value used.")
                }

                gofPval <- shape <- scale <- thresh <- excessPerm <- nExceed <-
                  excessObs <- NA

                out$shape <- out$scale <- out$gofPval <- NA

                out$pval <- (nlarger[idxFit[i]] + 1) / (nPerm + 1)

                out$fitFailed <- TRUE

              } else {
                # Fit and test the GPD distribution
                fittestres <- fit_gpd(data = tPerm[, idxFit[i]],
                                      thresh = thresh,
                                      fitMethod = fitMethod,
                                      tol = tol,
                                      eps = eps,
                                      factor = factor,
                                      constraint = constraint,
                                      maxVal = tMax[idxFit[i]],
                                      gofTest = gofTest,
                                      ...)

                out$shape <- fittestres$shape
                out$scale <- fittestres$scale
                out$gofPval <- fittestres$pval

                out$pval <- (out$nExceed / nPerm) *
                  pgpd_upper_tail(q = tObs[idxFit[i]] - thresh,
                                  loc = 0,
                                  scale = out$scale,
                                  shape = out$shape)

                if (out$pval == 0) {
                  if (verbose) {
                    message("In test ", idxFit[i],
                            ": GPD approximation led to a p-value of zero. ",
                            "Empirical p-value used.")
                  }

                  out$pval <- (nlarger[idxFit[i]] + 1) / (nPerm + 1)
                  out$zeroRepl <- TRUE
                }
              }

              out$names <- names(out)

              out
            }

  if (verbose) {
    # Close progress bar
    close(pb)
  }

  # Stop cluster
  if (cores > 1) parallel::stopCluster(cl)

  names(loopres) <- loopres[[length(loopres)]][[1]]
  loopres[[length(loopres)]] <- NULL

  threshVec <- nExceedVec <- shapeVec <- scaleVec <- gofPvalVec <-
    fitFailed <- zeroRepl <- rep(NA, nTest)

  threshVec[idxFit] <- unlist(loopres$thresh)
  nExceedVec[idxFit] <- unlist(loopres$nExceed)
  shapeVec[idxFit] <- unlist(loopres$shape)
  scaleVec[idxFit] <- unlist(loopres$scale)
  gofPvalVec[idxFit] <- unlist(loopres$gofPval)
  pval[idxFit] <- unlist(loopres$pval)
  fitFailed[idxFit] <- unlist(loopres$fitFailed)
  zeroRepl[idxFit] <- unlist(loopres$zeroRepl)

  callArgs <- mget(names(formals()),sys.frame(sys.nframe()))
  callArgs$gpdEstimate <- NULL

  output <- list(pval = pval,
                 idxFit = idxFit,
                 shape = shapeVec,
                 scale = scaleVec,
                 thresh = threshVec,
                 nExceed = nExceedVec,
                 gofPval = gofPvalVec,
                 tMax = tMax,
                 eps = eps,
                 fitFailed = fitFailed,
                 zeroRepl = zeroRepl,
                 callArgs = callArgs)

  return(output)
}


