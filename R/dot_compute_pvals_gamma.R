#' Compute p-values via Gamma approximation
#' 
#' @description
#' Internal workhorse for computing permutation p-values using a
#' Gamma distribution approximation. For each test in \code{idx_fit}, the
#' function:
#' \enumerate{
#'   \item performs a discreteness screening of the permutation distribution,
#'   \item fits a Gamma distribution to the (optionally augmented) permutation
#'         statistics using \pkg{fitdistrplus},
#'   \item optionally evaluates a goodness-of-fit test (currently Cramér–von
#'         Mises via \pkg{goftest}),
#'   \item computes tail probabilities under the fitted Gamma, scaled to the
#'         effective permutation sample size, and
#'   \item returns Gamma-based p-values and diagnostics.
#' }
#'
#' @param idx_fit Integer vector of test indices (subset of
#'   \code{1:length(obs_stats)}) for which Gamma approximation should be
#'   attempted. Typically defined in \code{\link{perm_approx}} as
#'   those tests with empirical p-values below \code{fit_thresh}.
#'
#' @param control A control list created by \code{\link{make_gamma_ctrl}}.
#'   
#' @inheritParams perm_approx
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pgamma
#' @keywords internal
#'
#' @return
#' A named \code{list} with one element per output quantity. Except for
#' \code{perm_stats_fit}, each element is a vector of length \code{n_test}:
#' \describe{
#'   \item{\code{p_value}}{Gamma-based tail p-values (or \code{NA_real_} if no
#'         valid Gamma fit was obtained or if the goodness-of-fit test rejects
#'         the Gamma approximation). These are combined with empirical
#'         p-values in \code{\link{perm_approx}}.}
#'   \item{\code{shape}}{Estimated Gamma shape parameter
#'         (\code{NA_real_} if no fit was obtained).}
#'   \item{\code{rate}}{Estimated Gamma rate parameter
#'         (\code{NA_real_} if no fit was obtained).}
#'   \item{\code{gof_p_value}}{Goodness-of-fit p-value (Cramér–von Mises test
#'         if \code{control$gof_test == "cvm"}, otherwise \code{NA_real_}).}
#'   \item{\code{status}}{Factor describing the per-test fitting outcome with
#'         levels \code{"not_selected"}, \code{"success"}, \code{"discrete"},
#'         \code{"no_threshold"}, \code{"gof_reject"}, \code{"fit_failed"}.
#'         For Gamma approximation, the \code{"no_threshold"} level is not
#'         used.}
#'   \item{\code{discrete}}{Logical flag indicating whether the permutation
#'         distribution was classified as too discrete for reliable Gamma
#'         approximation.}
#'   \item{\code{gof_rejected}}{Logical flag indicating whether the
#'         goodness-of-fit test rejected the Gamma approximation at
#'         \code{control$gof_alpha}.}
#'   \item{\code{method_used}}{Character vector indicating which p-value
#'         source should be used downstream for each test (typically
#'         \code{"gamma"} for successful Gamma fits and \code{"empirical"}
#'         otherwise).}
#'   \item{\code{perm_stats_fit}}{List of length \code{n_test} containing the
#'         permutation statistics used in the Gamma fit for each test (empty
#'         numeric vectors for tests not fitted or not successful).}
#' }

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
  
  ## -------------------------------------------------------------------------
  ## Prepare result containers
  ## -------------------------------------------------------------------------
  status_levels <- c("not_selected", "discrete", "fit_failed", "gof_reject", 
                     "success")
  
  p_value      <- rep(NA_real_, n_test)
  shape        <- rep(NA_real_, n_test)
  rate         <- rep(NA_real_, n_test)
  gof_p_value  <- rep(NA_real_, n_test)
  gof_rejected <- rep(FALSE,     n_test)
  method_used  <- rep("empirical", n_test)
  discrete     <- rep(FALSE,     n_test)
  status       <- factor(rep("not_selected", n_test), levels = status_levels)
  perm_stats_fit <- vector("list", length = n_test)
  
  .pack_result <- function() {
    list(
      p_value       = p_value,
      shape         = shape,
      rate          = rate,
      gof_p_value   = gof_p_value,
      status        = status,
      discrete      = discrete,
      gof_rejected  = gof_rejected,
      method_used   = method_used,
      perm_stats_fit = perm_stats_fit
    )
  }
  
  # Early exit if nothing to fit
  if (length(idx_fit) == 0L) {
    return(.pack_result())
  }
  
  ## -------------------------------------------------------------------------
  ## Discreteness screening
  ## -------------------------------------------------------------------------
  uniq_counts <- vapply(
    idx_fit,
    function(i) length(unique(perm_stats[, i])),
    integer(1)
  )
  is_discrete_sub <- uniq_counts < max(1L, floor(n_perm * 0.02))
  discrete[idx_fit] <- is_discrete_sub
  status[idx_fit[is_discrete_sub]] <- "discrete"
  
  idx_non_dis <- idx_fit[!is_discrete_sub]
  
  # Early exit if all tests are discrete
  if (length(idx_non_dis) == 0L) {
    return(.pack_result())
  }
  
  ## -------------------------------------------------------------------------
  ## Helper: fit Gamma for a single test j
  ## -------------------------------------------------------------------------
  fit_one <- function(j) {
    obs  <- obs_stats[j]
    perm <- perm_stats[, j]
    
    # Store used permutation stats
    perm_fit <- perm
    
    # Build argument list for fitdist:
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
        ok           = FALSE,
        gof_fail     = FALSE,
        perm_fit     = perm_fit
      ))
    }
    
    shape_j <- as.numeric(gamma_fit$estimate["shape"])
    rate_j  <- as.numeric(gamma_fit$estimate["rate"])
    
    # Tail probability under fitted Gamma
    factor_len <- if (include_obs) length(perm) / n_perm else 1
    p_val_raw <- factor_len * stats::pgamma(
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
          x     = perm,
          null  = "gamma",
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
    
    # Only accept p-values that are finite and in [0, 1]
    p_val_valid <- is.finite(p_val_raw) && p_val_raw >= 0 && p_val_raw <= 1
    ok <- (!rej && p_val_valid)
    
    list(
      p_value      = if (ok) p_val_raw else NA_real_,
      shape        = shape_j,
      rate         = rate_j,
      gof_p_value  = gof_p,
      gof_rejected = rej,
      method_used  = method,
      ok           = ok,
      gof_fail     = rej,
      perm_fit     = perm_fit
    )
  }
  
  ## -------------------------------------------------------------------------
  ## Decide if we run in parallel
  ## -------------------------------------------------------------------------
  n_workers   <- .choose_workers(cores, length(idx_non_dis), parallel_min)
  run_parallel <- (n_workers > 1L)
  n_tasks      <- length(idx_non_dis)
  
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
        p <- progressr::progressor(along = idx_non_dis)
      } else {
        p <- function(...) NULL
      }
      
      if (run_parallel) {
        strategy <- if (.Platform$OS.type != "windows" && 
                        future::supportsMulticore())
          future::multicore else future::multisession
        
        old_plan <- future::plan(strategy, workers = n_workers)
        on.exit(future::plan(old_plan), add = TRUE)
        
        future.apply::future_lapply(
          idx_non_dis,
          function(j) {
            res <- fit_one(j)
            p()
            res
          },
          future.packages = c("permApprox", "fitdistrplus", "goftest", "progressr"),
          future.seed     = TRUE
        )
        
      } else {
        lapply(idx_non_dis, function(j) {
          res <- fit_one(j)
          p()
          res
        })
      }
    }
  )
  
  if (isTRUE(verbose)) writeLines("Done.")
  
  ## -------------------------------------------------------------------------
  ## Fill results back into full-length vectors
  ## -------------------------------------------------------------------------
  for (k in seq_along(idx_non_dis)) {
    j   <- idx_non_dis[k]
    res <- res_list[[k]]
    
    p_value[j]        <- res$p_value
    shape[j]          <- res$shape
    rate[j]           <- res$rate
    gof_p_value[j]    <- res$gof_p_value
    gof_rejected[j]   <- res$gof_rejected
    method_used[j]    <- res$method_used
    perm_stats_fit[[j]] <- res$perm_fit
  }
  
  # Derive status from ok / gof_fail (only for non-discrete, attempted fits)
  ok_vec   <- vapply(res_list, `[[`, logical(1), "ok")
  gof_fail <- vapply(res_list, `[[`, logical(1), "gof_fail")
  
  status[idx_non_dis[ ok_vec]]               <- "success"
  status[idx_non_dis[!ok_vec &  gof_fail]]   <- "gof_reject"
  status[idx_non_dis[!ok_vec & !gof_fail]]   <- "fit_failed"
  
  .pack_result()
}
