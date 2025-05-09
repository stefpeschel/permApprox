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
  n_perm <- ncol(perm_stats)

  # Unpack control
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

  # Transform test statistics for tail modeling
  transformed <- lapply(seq_len(n_test), function(i) {
    transform_stats(perm_stats = perm_stats[i, ],
                    obs_stats = obs_stats[i],
                    alternative = alternative)
  })

  # Extract transformed obs_stats and perm_stats into vectors/matrices
  trans_obs <- sapply(transformed, function(x) x$obs_stats)
  trans_perm <- lapply(transformed, function(x) x$perm_stats)

  # Determine which tests to fit
  idx_fit <- which(pvals < fit_thresh & trans_obs > 0)

  # Tests for which a Gamma fit is performed
  fitted <- logical(n_test)
  fitted[idx_fit] <- TRUE

  # Method that is finally used
  method_used <- rep("empirical", n_test)
  method_used[idx_fit] <- "gamma"

  # Maximum value at which the GPD density must be positive
  if (constraint == "support_at_max") {
    max_stats <- rep(max(obs_stats), length(obs_stats))

  } else if (constraint == "support_at_obs") {
    max_stats <- obs_stats

  } else {
    max_stats <- NULL
  }

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

              obs <- trans_obs[idx_fit[i]]
              perm <- trans_perm[[idx_fit[i]]]

              out$perm_stat_used <- perm

              threshList <- get_gpd_thresh(perm_stats = perm,
                                           obs_stats = obs,
                                           tol = tol,
                                           thresh_method = thresh_method,
                                           thresh0 = thresh0,
                                           exceed0 = exceed0,
                                           exceed_min = exceed_min,
                                           thresh_step = thresh_step,
                                           gof_test = gof_test,
                                           gof_alpha = gof_alpha)

              thresh <- threshList$thresh
              out$thresh <- thresh
              out$n_exceed <- threshList$n_exceed
              out$zero_replaced <- FALSE

              if (is.na(thresh)) {

                if (verbose) {
                  message("In test ", idx_fit[i],
                          ": No iteration led to a good fit ",
                          "(GOF test rejected for ",
                          "all thresholds).")
                }

                gof_p_value <- shape <- scale <- thresh <- n_exceed <- NA

                out$shape <- out$scale <- out$gof_p_value <- NA

                out$pval <- p_empirical[idx_fit[i]]

              } else {

                if (is.null(max_stats)) {
                  maxVal <- NULL
                } else {
                  maxVal <- max_stats[idx_fit[i]]
                }

                # Fit and test the GPD distribution
                fittestres <- fit_gpd(data = perm,
                                      thresh = thresh,
                                      fit_method = fit_method,
                                      tol = tol,
                                      eps = eps,
                                      eps_type = eps_type,
                                      constraint = constraint,
                                      maxVal = maxVal,
                                      gof_test = gof_test,
                                      ...)

                out$shape <- fittestres$shape
                out$scale <- fittestres$scale
                out$gof_p_value <- fittestres$pval
                out$method_used <- "gpd"

                out$p_value <- (out$n_exceed / n_perm) *
                  pgpd_upper_tail(q = obs_stats[idx_fit[i]] - thresh,
                                  loc = 0,
                                  scale = out$scale,
                                  shape = out$shape)

                if (out$p_value == 0) {
                  if (verbose) {
                    message("In test ", idx_fit[i],
                            ": GPD approximation led to a p-value of zero. ",
                            "Empirical p-value used.")
                  }

                  out$p_value <- p_empirical[idx_fit[i]]
                  out$zero_replaced <- TRUE
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
  # Evaluate loop results

  names(loopres) <- loopres[[length(loopres)]][[1]]
  loopres[[length(loopres)]] <- NULL

  threshVec <- n_exceedVec <- shapeVec <- scaleVec <- gof_p_value_vec <-
    zero_replaced <- rep(NA, n_test)

  # List to store the test statistics used for the fit
  perm_stats_fitted <- vector("list", length = n_test)

  method_used <- rep("empirical", n_test)

  threshVec[idx_fit] <- unlist(loopres$thresh)
  n_exceedVec[idx_fit] <- unlist(loopres$n_exceed)
  shapeVec[idx_fit] <- unlist(loopres$shape)
  scaleVec[idx_fit] <- unlist(loopres$scale)
  gof_p_value_vec[idx_fit] <- unlist(loopres$gof_p_value)
  pvals[idx_fit] <- unlist(loopres$p_value)
  method_used[idx_fit] <- unlist(loopres$method_used)
  zero_replaced[idx_fit] <- unlist(loopres$zero_replaced)
  perm_stats_fitted[idx_fit] <- loopres$perm_stat_used

  output <- list(p_values = pvals,
                 fitted = fitted,
                 shape = shapeVec,
                 scale = scaleVec,
                 thresh = threshVec,
                 n_exceed = n_exceedVec,
                 gof_p_value = gof_p_value_vec,
                 max_stats = max_stats,
                 eps = eps,
                 method_used = method_used,
                 zero_replaced = zero_replaced,
                 perm_stats_used = perm_stats_fitted)

  return(output)
}


