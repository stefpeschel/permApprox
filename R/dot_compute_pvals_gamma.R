#' @title Compute p-values via Gamma approximation
#' 
#' @inheritParams compute_p_values
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pgamma
#' @keywords internal
.compute_pvals_gamma <- function(
    obs_stats,
    perm_stats,
    idx_fit,
    control,
    cores,
    parallel_min,
    verbose
) {
  # Numbers of tests and permutations (before including obs)
  n_test <- length(obs_stats)
  n_perm <- nrow(perm_stats)
  
  # Early exit if nothing to fit
  if (length(idx_fit) == 0L) {
    return(list(
      p_values          = rep(NA_real_, n_test),
      shape             = rep(NA_real_, n_test),
      rate              = rep(NA_real_, n_test),
      gof_p_value       = rep(NA_real_, n_test),
      gof_rejected      = rep(FALSE, n_test),
      method_used       = rep("empirical", n_test),
      perm_stats_fitted = vector("list", length = n_test)
    ))
  }
  
  # Control arguments
  include_obs <- isTRUE(control$include_obs)
  gof_test    <- control$gof_test
  gof_alpha   <- control$gof_alpha
  fit_args    <- control$fit_args %||% list()  # optional fitting arguments
  
  # Optionally include observed statistics in the fit
  if (include_obs) {
    # obs_stats is a row vector, perm_stats has rows = permutations
    perm_stats <- rbind(obs_stats, perm_stats)
  }
  
  # Initialize output containers (full length = n_test)
  shape        <- rep(NA_real_, n_test)
  rate         <- rep(NA_real_, n_test)
  gof_p_value  <- rep(NA_real_, n_test)
  gof_rejected <- rep(FALSE, n_test)
  method_used  <- rep("empirical", n_test)
  perm_stats_fitted <- vector("list", length = n_test)
  pvals        <- rep(NA_real_, n_test)
  
  # Helper: fit Gamma for a single test j
  fit_one <- function(j) {
    obs  <- obs_stats[j]
    perm <- perm_stats[, j]
    
    # Store used permutation stats
    perm_fit <- perm
    
    # Build argument list for fitdist:
    # data and distr are fixed, user can override e.g. method/start/...
    args_list <- c(
      list(
        data  = perm,
        distr = "gamma"
      ),
      fit_args
    )
    
    # Fit Gamma distribution (suppress NaN warnings)
    gamma_fit <- tryCatch(
      suppressWarnings(
        do.call(fitdistrplus::fitdist, args_list)
      ),
      error = function(e) NULL
    )
    
    if (is.null(gamma_fit)) {
      return(list(
        p_value      = NA_real_,
        shape        = NA_real_,
        rate         = NA_real_,
        gof_p_value  = NA_real_,
        gof_rejected = FALSE,
        method_used  = "empirical",
        perm_fit     = perm_fit
      ))
    }
    
    shape_j <- as.numeric(gamma_fit$estimate["shape"])
    rate_j  <- as.numeric(gamma_fit$estimate["rate"])
    
    # Tail probability under fitted Gamma
    factor_len <- if (include_obs) length(perm) / n_perm else 1
    p_val <- factor_len * stats::pgamma(
      q          = abs(obs),
      shape      = shape_j,
      rate       = rate_j,
      lower.tail = FALSE
    )
    
    gof_p  <- NA_real_
    rej    <- FALSE
    method <- "gamma"
    
    if (identical(gof_test, "cvm")) {
      cvmtest <- tryCatch(
        goftest::cvm.test(
          x    = perm,
          null = "gamma",
          shape = shape_j,
          rate  = rate_j
        ),
        error = function(e) NULL
      )
      
      if (!is.null(cvmtest)) {
        gof_p <- cvmtest$p.value
        if (!is.null(gof_alpha) && gof_p <= gof_alpha) {
          rej    <- TRUE
          method <- "empirical"  # reject Gamma fit, fall back to empirical
        }
      }
    }
    
    list(
      p_value      = p_val,
      shape        = shape_j,
      rate         = rate_j,
      gof_p_value  = gof_p,
      gof_rejected = rej,
      method_used  = method,
      perm_fit     = perm_fit
    )
  }
  
  # Decide if we run in parallel
  n_workers <- .choose_workers(cores, length(idx_fit), parallel_min)
  run_parallel <- (n_workers > 1L)
  
  n_tasks <- length(idx_fit)
  
  if (isTRUE(verbose)) {
    mode_txt <- if (run_parallel) "parallel" else "sequential"
    writeLines(sprintf(
      "Fitting Gamma distribution for %d tests (%s) ...",
      n_tasks, mode_txt
    ))
  }
  
  res_list <- progressr::with_progress(
    handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
    {
      if (verbose) {
        p <- progressr::progressor(along = idx_fit)
      } else {
        # dummy progressor that does nothing
        p <- function(...) NULL
      }
      
      if (run_parallel) {
        strategy <- if (.Platform$OS.type != "windows" && 
                        future::supportsMulticore())
          future::multicore else future::multisession
        
        old_plan <- future::plan(strategy, workers = n_workers)
        on.exit(future::plan(old_plan), add = TRUE)
        
        future.apply::future_lapply(
          idx_fit,
          function(j) {
            res <- fit_one(j)
            p()   # progress update
            res
          },
          future.packages = c("permApprox", "fitdistrplus", "goftest", "progressr"),
          future.seed     = TRUE
        )
        
      } else {
        # Sequential loop with the SAME progress bar
        lapply(idx_fit, function(j) {
          res <- fit_one(j)
          p()   # progress update
          res
        })
      }
    }
  )
  
  if (isTRUE(verbose)) writeLines("Done.")
  
  # Fill results back into full-length vectors
  for (k in seq_along(idx_fit)) {
    j   <- idx_fit[k]
    res <- res_list[[k]]
    
    pvals[j]              <- res$p_value
    shape[j]              <- res$shape
    rate[j]               <- res$rate
    gof_p_value[j]        <- res$gof_p_value
    gof_rejected[j]       <- res$gof_rejected
    method_used[j]        <- res$method_used
    perm_stats_fitted[[j]] <- res$perm_fit
  }
  
  output <- list(
    p_values          = pvals,
    shape             = shape,
    rate              = rate,
    gof_p_value       = gof_p_value,
    gof_rejected      = gof_rejected,
    method_used       = method_used,
    perm_stats_fitted = perm_stats_fitted
  )
  
  return(output)
}
