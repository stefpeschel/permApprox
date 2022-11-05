
get_gpd_estimates <- function(tPerm,
                              tObs,
                              tMax,
                              tol,
                              threshVec = NULL,
                              threshMethod = "PRbelowAlpha",
                              thresh0 = NULL,
                              nExceed0 = NULL,
                              nExceedMin = 1,
                              stepSize = 1,
                              includeObs = FALSE,
                              fitMethod = "MLE",
                              constraint = "none",
                              fitInformation = "observed",
                              optimMethod = NULL,
                              shape0 = NULL, scale0 = NULL,
                              shapeMin = -Inf, scaleMin = -Inf,
                              shapeMax = Inf, scaleMax = Inf,
                              gofTest = "ad",
                              seed = NULL,
                              cores = 1L,
                              verbose = FALSE,
                              ...) {

  if (!is.null(seed)) set.seed(seed)

  nPerm <- length(tPerm)

  # Sort permutation test statistics in increasing order
  tSort <- sort(tPerm, decreasing = FALSE)
  nPerm_gpd <- nPerm


  #-----------------------------------------------------------------------------
  # Define vector with possible thresholds
  if (is.null(threshVec)) {

    if (threshMethod == "fix") {
      if (!is.null(thresh0)) {
        threshVec <- thresh0
      } else if (!is.null(nExceed0)) {
        threshVec <- tSort[tSort < max(tObs, tPerm)]
        threshVec <- c(0, threshVec)
        tmp <- sort(tSort, decreasing = TRUE)[nExceed0 + 1]
        threshVec <- threshVec[threshVec == tmp]

      } else {
        stop("Value for thresh0 or nExceed0 must be set")
      }

    } else {
      # Thresholds must be smaller than the observed test statistic and the largest
      # tPerm to ensure at least 1 exceedance
      threshVec <- tSort[tSort < min(tObs, max(tPerm))]
      threshVec <- c(0, threshVec)

      #-----------------
      # Adapt threshold vector to start threshold or start number of exceedances
      if (!is.null(thresh0) && !is.null(nExceed0)) {
        stop("Either thresh0 or nExceed0 must be set to NULL")
      }

      if (!is.null(thresh0)) {
        threshVec <- threshVec[threshVec > thresh0]
        threshVec <- c(thresh0, threshVec)
      }

      if (!is.null(nExceed0)) {
        tmp <- sort(tSort, decreasing = TRUE)[nExceed0 + 1]
        threshVec <- threshVec[threshVec >= tmp]
      }

      #-----------------
      # Adapt threshold vector according to stepSize
      threshVec <- threshVec[c(TRUE, rep(FALSE, stepSize - 1))]

      # Adapt threshold vector to minimum number of exceedances
      if (nExceedMin > 1) {
        tmp <- sort(tSort, decreasing = TRUE)[nExceedMin]
        threshVec <- threshVec[threshVec < tmp]
      }

      # Make thresholds unique
      threshVec <- unique(threshVec)
    }
  }

  # Number of iterations
  niter <- length(threshVec)

  #-----------------------------------------------------------------------------
  # Load functions

  # path <- "F:/Nextcloud/project2/permpap/R/"
  # files.sources <- list.files(path = path)
  # files.sources <- files.sources[!files.sources %in%
  #                                  c("approx_pval.R",
  #                                    "est_gpd_params.R",
  #                                    "fit_gpd.R",
  #                                    "get_gpd_estimates.R",
  #                                    "plothist.R")]
  # sapply(paste0(path, files.sources), source)

  #-----------------------------------------------------------------------------
  # Initialize parallel stuff

  if (verbose) {
    # Create progress bar:
    pb <- utils::txtProgressBar(0, niter, style=3)

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

  loopres <- foreach(i = seq_along(threshVec),
                     .export = c("gpd_fit_test", "gpd_optim", "gpdAd_adapt",
                                 "gpdCvm_adapt", "pgpd_upper_tail",
                                 "gpd_LME", "gpd_MLE1D", "gpd_MLE2D",
                                 "gpd_MOM", "gpd_NLS2", "gpd_WNLLSM", "gpd_ZSE"),
                     .combine='comb', .multicombine=TRUE,
                     .init=list(list(), list(), list(), list(), list(),
                                list()),
                     .options.snow = opts) %do_or_dopar% {
                       if (verbose) progress(i)

                       thresh <- threshVec[i]

                       out <- list(idx = i)

                       # exceedPerm are the exceedances (test statistics above the threshold)
                       exceedPerm.tmp <- tSort[tSort > thresh]

                       # number of exceedances
                       nExceed.tmp <- out$nExceed <- length(exceedPerm.tmp)

                       # Fit and test the GPD distribution
                       fittestres <-
                         fit_gpd(data = tSort,
                                 thresh = thresh,
                                 tol = tol,
                                 nExceed = nExceed.tmp,
                                 exceedPerm = exceedPerm.tmp,
                                 fitMethod = fitMethod,
                                 maxVal = tMax,
                                 constraint = constraint,
                                 gofTest = gofTest,
                                 ...)

                       shape.tmp <- out$shape <- fittestres$shape
                       scale.tmp <- out$scale <- fittestres$scale
                       #out$negLogLik <- fittestres$negLogLik
                       out$gofPval <- fittestres$pval

                       excessObs <- tObs - thresh

                       if (is.na(shape.tmp)) {
                         out$pval <- NA
                       } else {
                         # out$pval <- (nExceed.tmp / nPerm_gpd) * (VGAM::pgpd(q = excessObs,
                         #                                                         loc = 0,
                         #                                                         scale = scale.tmp,
                         #                                                         shape = shape.tmp,
                         #                                                         lower.tail = FALSE))

                         out$pval <- (nExceed.tmp / nPerm_gpd) * pgpd_upper_tail(q = excessObs,
                                                                                 loc = 0,
                                                                                 scale = scale.tmp,
                                                                                 shape = shape.tmp)
                       }

                       out
                     }

  if (verbose) {
    # Close progress bar
    close(pb)
  }

  # Stop cluster
  if (cores > 1) parallel::stopCluster(cl)

  out <- list()

  out$idxVec <- unlist(loopres[[1]])
  out$nExceedVec <- unlist(loopres[[2]])
  out$shapeVec <- unlist(loopres[[3]])
  out$scaleVec <- unlist(loopres[[4]])
  out$negLogLikVec <- NA
  out$gofPvalVec <- unlist(loopres[[5]])
  out$pvalVec <- unlist(loopres[[6]])
  out$threshVec <- threshVec

  return(out)
}
