#' @title Compute p-values via GPD approximation
#' @keywords internal

get_pvals_gpd <- function(tObs,
                          tPerm,
                          nTest,
                          nPerm,
                          ntPerm,
                          useAllPerm,
                          constraint,
                          fitThresh,
                          includeObs,
                          fitMethod,
                          tol,
                          eps,
                          epsVal,
                          threshMethod,
                          thresh0,
                          nExceed0,
                          nExceedMin,
                          stepSize,
                          gofTest,
                          gofAlpha,
                          gofTailRMMeth,
                          gofTailRMPar,
                          cores,
                          verbose,
                          ...) {

  # Maximum value at which the GPD density must be positive
  if (constraint == "tObs") {
    tMax <- tObs

  } else if (constraint == "tObsMax") {
    tMax <- rep(max(tObs), length(tObs))

  } else {
    tMax <- NULL
  }

  pemp <- get_pvals_emp(tObs = tObs, tPerm = tPerm, nTest = nTest,
                        nPerm = nPerm, ntPerm = ntPerm, useAllPerm = useAllPerm)
  pvals <- pemp$pvals
  nlarger <- pemp$nlarger

  # Indices of p-values below threshold (only these are fitted)
  idxFit <- which(pvals <= fitThresh)

  fitted <- rep(FALSE, nTest)
  fitted[idxFit] <- TRUE

  #---------------------------------------------------------------------------
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

  #---------------------------------------------------------------------------
  # Compute p-values

  loopres <-
    foreach(i = seq_along(idxFit),
            .export = c("get_gpd_thresh", "fit_gpd", "get_thresh_idx",
                        ".est_gpd_params", "get_pvals_emp",
                        "gpdAd_adapt", "gpdCvm_adapt",
                        "gpd_LME", "gpd_MLE1D", ".MLE1D_fk", ".MLE1D_fp",
                        "gpd_MLE2D", ".MLE2D_negloglik", "gpd_MOM",
                        "gpd_NLS2", ".NLS2_gpdf", ".NLS2_ecdf", ".NLS2_gpdf2",
                        ".NLS2_ecdf2", ".NLS2_rss1", ".NLS2_rss2",
                        "gpd_WNLLSM", ".WNLLSM_sum_i", ".WNLLSM_WLLS1",
                        ".WNLLSM_WLLS", ".WNLSM_WLS1", ".WNLSM_WLS",
                        "gpd_ZSE", ".ZSE_lx"),
            .packages = "permpap",
            .combine='comb', .multicombine=TRUE,
            .init=list(list(), list(), list(), list(), list(), list(),
                       list(), list(), list()),
            .options.snow = opts) %do_or_dopar% {

              if (verbose) progress(i)

              out <- list()

              tPermTmp <- if (useAllPerm) {
                as.vector(tPerm)
              } else {
                tPerm[idxFit[i], ]
              }

              threshList <- get_gpd_thresh(tPerm = tPermTmp,
                                           tObs = tObs[idxFit[i]],
                                           useAllPerm = useAllPerm,
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
                fittestres <- fit_gpd(data = tPermTmp,
                                      thresh = thresh,
                                      fitMethod = fitMethod,
                                      tol = tol,
                                      eps = eps,
                                      epsVal = epsVal,
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

  #---------------------------------------------------------------------------

  names(loopres) <- loopres[[length(loopres)]][[1]]
  loopres[[length(loopres)]] <- NULL

  threshVec <- nExceedVec <- shapeVec <- scaleVec <- gofPvalVec <-
    fitFailed <- zeroRepl <- rep(NA, nTest)

  threshVec[idxFit] <- unlist(loopres$thresh)
  nExceedVec[idxFit] <- unlist(loopres$nExceed)
  shapeVec[idxFit] <- unlist(loopres$shape)
  scaleVec[idxFit] <- unlist(loopres$scale)
  gofPvalVec[idxFit] <- unlist(loopres$gofPval)
  pvals[idxFit] <- unlist(loopres$pval)
  fitFailed[idxFit] <- unlist(loopres$fitFailed)
  zeroRepl[idxFit] <- unlist(loopres$zeroRepl)

  output <- list(pvals = pvals,
                 fitted = fitted,
                 shape = shapeVec,
                 scale = scaleVec,
                 thresh = threshVec,
                 nExceed = nExceedVec,
                 gofPval = gofPvalVec,
                 tMax = tMax,
                 eps = eps,
                 fitFailed = fitFailed,
                 zeroRepl = zeroRepl)

  return(output)
}
