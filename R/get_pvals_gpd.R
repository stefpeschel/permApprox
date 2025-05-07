#' @title Compute p-values via GPD approximation
#' @keywords internal

get_pvals_gpd <- function(p_empirical,
                          obs_stats,
                          perm_stats,
                          n_test,
                          fit_thresh,
                          alternative,
                          control,
                          ...) {

  pvals <- p_empirical

  # Assign control arguments
  fit_method <- control$fit_method
  include_obs <- control$include_obs
  constraint <- control$constraint
  eps <- control$eps
  eps_type <- control$eps_type
  tol <- control$tol
  thresh_method <- control$thresh_method
  thresh0 <- control$thresh0
  thresh_step <- control$thresh_step
  exceed0 <- control$exceed0
  exceed_min <- control$exceed_min
  gof_test <- control$gof_test
  gof_alpha <- control$gof_alpha
  cores <- control$cores
  verbose <- control$verbose
browser()
  # Include observed test statistic
  if (include_obs) perm_stats <- cbind(obs_stats, perm_stats)

  # Maximum value at which the GPD density must be positive
  if (constraint == "support_at_max") {
    max_stats <- rep(max(obs_stats), length(obs_stats))

  } else {
    max_stats <- obs_stats
  }

  # Indices of p-values below threshold (only these are fitted)
  idx_fit <- which(pvals <= fit_thresh)

  # Tests for which a GPD fit is performed
  fitted <- rep(FALSE, n_test)
  fitted[idx_fit] <- TRUE

  # Method that is finally used
  method_used <- rep("empirical", n_test)
  method_used[idx_fit] <- "gpd"

  #---------------------------------------------------------------------------
  # Initialize parallel stuff

  if (verbose) {
    # Create progress bar:
    pb <- utils::txtProgressBar(0, length(idx_fit), style=3)

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
    foreach(i = seq_along(idx_fit),
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

              perm_stats_row <- perm_stats[i, ]
              perm_stats_used <- get_perm_stats_used(perm_stats_row, alternative)

              n_used <- length(perm_stats_used)

              out$perm_stats_used <- perm_stats_used

              threshList <- get_gpd_thresh(perm_stats = perm_stats_used,
                                           obs_stats = abs(obs_stats[idx_fit[i]]),
                                           tol = tol,
                                           threshPoss = threshPoss,
                                           thresh_method = thresh_method,
                                           thresh0 = thresh0,
                                           exceed0 = exceed0,
                                           exceed_min = exceed_min,
                                           thresh_step = thresh_step,
                                           includeObs = includeObs,
                                           fit_method = fit_method,
                                           constraint = constraint,
                                           fitInformation = fitInformation,
                                           optimMethod = optimMethod,
                                           shape0 = shape0,
                                           scale0 = scale0,
                                           shapeMin = shapeMin,
                                           scaleMin = scaleMin,
                                           shapeMax = shapeMax,
                                           scaleMax = scaleMax,
                                           gof_test = gof_test,
                                           gof_alpha = gof_alpha,
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
                  message("In test ", idx_fit[i],
                          ": No iteration led to a good fit ",
                          "(GOF test rejected for ",
                          "all thresholds).")
                }

                gof_p_value <- shape <- scale <- thresh <- excessPerm <- nExceed <-
                  excessObs <- NA

                out$shape <- out$scale <- out$gof_p_value <- NA

                if (gammaOnFail) {
                  # Fit Gamma distribution (warning about NANs is suppressed)
                  suppressWarnings(gammafit <- fitdistrplus::fitdist(data = perm_statsUsed,
                                                                     distr = "gamma",
                                                                     method = "mle"))

                  shape <- as.numeric(gammafit$estimate["shape"])
                  rate <- as.numeric(gammafit$estimate["rate"])

                  nUsed <- length(perm_statsUsed)
                  pval_gamma <- (nUsed / n_perm) * pgamma(q = abs(obs_stats[idx_fit[i]]),
                                                         shape = shape,
                                                         rate = rate,
                                                         lower.tail = FALSE)

                  # Goodness-of-fit test
                  cvmtest <- gof_test::cvm.test(x = perm_statsUsed,
                                               null = "gamma",
                                               shape = shape,
                                               rate = rate)

                  out$gof_p_value <- cvmtest$p.value

                  if (gof_testGamma && (out$gof_p_value <= gof_alpha)) {
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
                fittestres <- fit_gpd(data = perm_statsUsed,
                                      thresh = thresh,
                                      fit_method = fit_method,
                                      tol = tol,
                                      eps = eps,
                                      eps_type = eps_type,
                                      constraint = constraint,
                                      maxVal = abs(max_stats[idx_fit[i]]),
                                      gof_test = gof_test,
                                      ...)

                out$shape <- fittestres$shape
                out$scale <- fittestres$scale
                out$gof_p_value <- fittestres$pval
                out$approxType <- "gpd"

                out$pval <- (out$nExceed / n_perm) *
                  pgpd_upper_tail(q = obs_stats[idx_fit[i]] - thresh,
                                  loc = 0,
                                  scale = out$scale,
                                  shape = out$shape)

                if (out$pval == 0) {
                  if (verbose) {
                    message("In test ", idx_fit[i],
                            ": GPD approximation led to a p-value of zero. ",
                            "Empirical p-value used.")
                  }

                  #out$pval <- (nlarger[idx_fit[i]] + 1) / (n_perm + 1)
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

  threshVec <- nExceedVec <- shapeVec <- scaleVec <- gof_p_value_vec <-
    zeroRepl <- rep(NA, n_test)

  approxType <- rep("empirical", n_test)

  threshVec[idx_fit] <- unlist(loopres$thresh)
  nExceedVec[idx_fit] <- unlist(loopres$nExceed)
  shapeVec[idx_fit] <- unlist(loopres$shape)
  scaleVec[idx_fit] <- unlist(loopres$scale)
  gof_p_value_vec[idx_fit] <- unlist(loopres$gof_p_value)
  pvals[idx_fit] <- unlist(loopres$pval)
  approxType[idx_fit] <- unlist(loopres$approxType)
  zeroRepl[idx_fit] <- unlist(loopres$zeroRepl)

  output <- list(pvals = pvals,
                 fitted = fitted,
                 shape = shapeVec,
                 scale = scaleVec,
                 thresh = threshVec,
                 nExceed = nExceedVec,
                 gof_p_value = gof_p_value_vec,
                 max_stats = max_stats,
                 eps = eps,
                 approxType = approxType,
                 zeroRepl = zeroRepl,
                 perm_statsUsed = loopres$perm_statsUsed)

  return(output)
}
