#' @title Compute p-values via GPD approximation
#' @keywords internal

get_pvals_gpd <- function(tObs,
                          tPerm,
                          alternative,
                          pEmp,
                          nExtreme,
                          nTest,
                          nPerm,
                          constraint,
                          fitThresh,
                          gammaOnFail,
                          gofTestGamma,
                          includeObs,
                          fitMethod,
                          tol,
                          eps,
                          epsType,
                          threshMethod,
                          thresh0,
                          threshPoss,
                          exceed0,
                          exceedMin,
                          stepSize,
                          gofTest,
                          gofAlpha,
                          cores,
                          verbose,
                          ...) {
  
  if (is.null(thresh0) & is.null(exceed0)) {
    message("exceed0 set to used number of permutations.")
    # Will be set by get_gpd_thresh
  }

  # Maximum value at which the GPD density must be positive
  if (constraint == "tObsMax") {
    tMax <- rep(max(tObs), length(tObs))

  } else {
    tMax <- tObs
  }

  pvals <- pEmp

  # Indices of p-values below threshold (only these are fitted)
  idxFit <- which(pvals <= fitThresh)
  fitted <- rep(FALSE, nTest)
  fitted[idxFit] <- TRUE

  # Type of p-value estimation
  approxType <- rep(NA, nTest)

  approxType[!seq_along(tObs) %in% idxFit] <- "empirical"

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
                        "gpd_NLS2", "pgpd_upper_tail",
                        ".NLS2_gpdf", ".NLS2_ecdf", ".NLS2_gpdf2",
                        ".NLS2_ecdf2", ".NLS2_rss1", ".NLS2_rss2",
                        "gpd_WNLLSM", ".WNLLSM_sum_i", ".WNLLSM_WLLS1",
                        ".WNLLSM_WLLS", ".WNLSM_WLS1", ".WNLSM_WLS",
                        "gpd_ZSE", ".ZSE_lx"),
            #.packages = "permAprox",
            .combine='comb', .multicombine=TRUE,
            .init=list(list(), list(), list(), list(), list(), list(),
                       list(), list(), list(), list()),
            .options.snow = opts) %do_or_dopar% {
              
              if (verbose) progress(i)

              out <- list()

              if (alternative == "less") {
                tPermUsed <- tPerm[idxFit[i], ]
                tPermUsed <- abs(tPermUsed[tPermUsed < 0])

              } else if (alternative == "greater") {
                tPermUsed <- tPerm[idxFit[i], ]
                tPermUsed <- tPermUsed[tPermUsed > 0]

              } else {
                tPermUsed <- abs(tPerm[idxFit[i], ])
                tPermUsed <- tPermUsed[tPermUsed > 0]
              }

              out$tPermUsed <- tPermUsed

              threshList <- get_gpd_thresh(tPerm = tPermUsed,
                                           tObs = abs(tObs[idxFit[i]]),
                                           tMax = abs(tMax[idxFit[i]]),
                                           tol = tol,
                                           threshPoss = threshPoss,
                                           threshMethod = threshMethod,
                                           thresh0 = thresh0,
                                           exceed0 = exceed0,
                                           exceedMin = exceedMin,
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
              out$zeroRepl <- FALSE

              if (is.na(thresh)) {

                if (verbose) {
                  message("In test ", idxFit[i],
                          ": No iteration led to a good fit ",
                          "(GOF test rejected for ",
                          "all thresholds).")
                }

                gofPval <- shape <- scale <- thresh <- excessPerm <- nExceed <-
                  excessObs <- NA

                out$shape <- out$scale <- out$gofPval <- NA

                if (gammaOnFail) {
                  # Fit Gamma distribution (warning about NANs is suppressed)
                  suppressWarnings(gammafit <- fitdistrplus::fitdist(data = tPermUsed,
                                                                     distr = "gamma",
                                                                     method = "mle"))

                  shape <- as.numeric(gammafit$estimate["shape"])
                  rate <- as.numeric(gammafit$estimate["rate"])

                  nUsed <- length(tPermUsed)
                  pval_gamma <- (nUsed / nPerm) * pgamma(q = abs(tObs[idxFit[i]]),
                                                         shape = shape,
                                                         rate = rate,
                                                         lower.tail = FALSE)

                  # Goodness-of-fit test
                  cvmtest <- goftest::cvm.test(x = tPermUsed,
                                               null = "gamma",
                                               shape = shape,
                                               rate = rate)

                  out$gofPval <- cvmtest$p.value

                  if (gofTestGamma && (out$gofPval <= gofAlpha)) {
                    if (verbose) {
                      message(" Empirical p-value used.")
                    }
                    out$approxType <- "empirical"
                    out$pval <- pvals[i]

                  } else {
                    if (verbose) {
                      message(" Gamma approximation used.")
                    }
                    out$approxType <- "gamma"
                    out$pval <- pval_gamma
                  }

                } else {
                  out$pval <- pvals[i]
                }

              } else {
                # Fit and test the GPD distribution
                fittestres <- fit_gpd(data = tPermUsed,
                                      thresh = thresh,
                                      fitMethod = fitMethod,
                                      tol = tol,
                                      eps = eps,
                                      epsType = epsType,
                                      constraint = constraint,
                                      maxVal = abs(tMax[idxFit[i]]),
                                      gofTest = gofTest,
                                      ...)

                out$shape <- fittestres$shape
                out$scale <- fittestres$scale
                out$gofPval <- fittestres$pval
                out$approxType <- "gpd"

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

                  #out$pval <- (nlarger[idxFit[i]] + 1) / (nPerm + 1)
                  out$pval <- pvals[i]
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
    zeroRepl <- rep(NA, nTest)

  approxType <- rep("empirical", nTest)

  threshVec[idxFit] <- unlist(loopres$thresh)
  nExceedVec[idxFit] <- unlist(loopres$nExceed)
  shapeVec[idxFit] <- unlist(loopres$shape)
  scaleVec[idxFit] <- unlist(loopres$scale)
  gofPvalVec[idxFit] <- unlist(loopres$gofPval)
  pvals[idxFit] <- unlist(loopres$pval)
  approxType[idxFit] <- unlist(loopres$approxType)
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
                 approxType = approxType,
                 zeroRepl = zeroRepl,
                 tPermUsed = loopres$tPermUsed)

  return(output)
}
