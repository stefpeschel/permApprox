#' Compute p-values via GPD tail approximation (constrained or unconstrained)
#'
#' @description
#' Internal workhorse for computing permutation p-values using the
#' Generalized Pareto Distribution (GPD) tail approximation.  
#' For each test, the function:
#' \enumerate{
#'   \item receives a set of tests selected for GPD fitting (indices in
#'         \code{idx_fit}, typically based on empirical p-values below a
#'         threshold in the calling function),
#'   \item performs discreteness screening,
#'   \item detects a GPD threshold using \code{\link{.find_gpd_thresh}},
#'   \item evaluates the appropriate support boundary
#'         (depending on the constraint in \code{control}),
#'   \item computes \eqn{\varepsilon}-values via the epsilon function
#'         defined in \code{control} (\code{eps_fun}, \code{eps_tune},
#'         \code{eps_args}),
#'   \item fits the GPD tail using \code{\link{fit_gpd}},
#'   \item optionally performs adaptive epsilon refinement
#'         (zero-guard) to avoid numerical underflow, and
#'   \item returns GPD-based p-values and diagnostics.
#' }
#'
#' The function supports parallel execution for threshold detection and
#' GPD fitting through the \pkg{future} framework.
#'
#' @param obs_stats Numeric vector of observed test statistics
#'   (\code{length = n_test}).
#'
#' @param perm_stats Numeric matrix of permutation statistics with dimension
#'   \code{B x n_test}.  
#'   If \code{control$include_obs = TRUE}, the vector \code{obs_stats} is
#'   appended as an additional row prior to threshold detection and fitting.
#'   
#' @param idx_fit Integer vector of test indices (subset of
#'   \code{1:length(obs_stats)}) for which GPD approximation should be
#'   attempted. Typically defined in \code{\link{perm_approx}} as
#'   those tests with empirical p-values below \code{approx_thresh}.
#'
#' @param control A control list created by \code{\link{make_gpd_ctrl}}.
#'   It contains GPD fitting options, threshold detection settings,
#'   constraint mode, epsilon function (\code{eps_fun}),
#'   epsilon tuning parameter (\code{eps_tune}),
#'   optional epsilon arguments (\code{eps_args}),
#'   and zero-guard configuration.
#'
#' @inheritParams perm_approx
#'
#' @param ... Additional arguments passed downstream to internal fitting
#'   functions such as \code{fit_gpd}.
#'
#' @return
#' A named \code{list} with one element per output quantity. Except for
#' \code{perm_stats_fit}, each element is a vector of length \code{n_test}:
#' \describe{
#'   \item{\code{p_value}}{GPD-based tail p-values (or \code{NA_real_} if no
#'         valid GPD fit was obtained for a test). These are combined with
#'         empirical p-values in \code{\link{perm_approx}}.}
#'   \item{\code{thresh}}{Selected GPD threshold (or \code{NA_real_} if none
#'         found).}
#'   \item{\code{n_exceed}}{Number of exceedances used for the GPD fit.}
#'   \item{\code{shape}}{Estimated GPD shape parameter (\code{NA_real_} if no
#'         fit was obtained).}
#'   \item{\code{scale}}{Estimated GPD scale parameter (\code{NA_real_} if no
#'         fit was obtained).}
#'   \item{\code{epsilon}}{Final epsilon used in the GPD fit (after refinement
#'         if zero-guard was active).}
#'   \item{\code{gof_p_value}}{Goodness-of-fit p-value (Anderson–Darling or
#'         Cramér–von Mises depending on \code{control$gof_test}).}
#'   \item{\code{status}}{Factor describing the per-test fitting outcome:
#'     \code{"not_selected"}, \code{"success"}, \code{"discrete"},
#'     \code{"no_threshold"}, \code{"gof_reject"}, \code{"fit_failed"}.}
#'   \item{\code{discrete}}{Logical flag indicating whether the permutation
#'         distribution was too discrete for reliable tail fitting.}
#'   \item{\code{perm_stats_fit}}{List of length \code{n_test} containing the
#'         exceedances used for each successful GPD fit (empty numeric vectors
#'         for tests not fitted or not successful).}
#' }
#'
#' @keywords internal


.compute_pvals_gpd <- function(
    obs_stats,
    perm_stats,
    idx_fit,
    control,
    cores,
    parallel_min,
    verbose,
    ...
) {
  
  # Numbers of tests and permutations
  n_test <- length(obs_stats)
  n_perm <- nrow(perm_stats)
  
  if (isTRUE(control$include_obs)) {
    stopifnot(length(obs_stats) == ncol(perm_stats))
    perm_stats <- rbind(perm_stats, obs_stats)
  }
  
  ## -------------------------------------------------------------------------
  ## Prepare result vectors
  ## -------------------------------------------------------------------------
  status_levels <- c("not_selected", "discrete", "no_threshold",
                     "fit_failed", "gof_reject", "success")
  
  p_value     <- rep(NA_real_,    n_test)
  thresh      <- rep(NA_real_,    n_test)
  n_exceed    <- rep(NA_integer_, n_test)
  shape       <- rep(NA_real_,    n_test)
  scale       <- rep(NA_real_,    n_test)
  epsilon     <- rep(NA_real_,    n_test)
  gof_p_value <- rep(NA_real_,    n_test)
  status      <- factor(rep("not_selected", n_test), levels = status_levels)
  discrete    <- rep(FALSE, n_test)
  
  ## If nothing to fit, attach empty perm_stats_fit and return
  perm_stats_fit <- vector("list", n_test)
  
  .pack_result <- function() {
    list(
      p_value       = p_value,
      thresh        = thresh,
      n_exceed      = n_exceed,
      shape         = shape,
      scale         = scale,
      epsilon       = epsilon,
      gof_p_value   = gof_p_value,
      status        = status,
      discrete      = discrete,
      perm_stats_fit = perm_stats_fit
    )
  }
  
  ## If nothing to fit, return immediately
  if (length(idx_fit) == 0L) {
    return(.pack_result())
  }
  
  ## -------------------------------------------------------------------------
  ## Discreteness screening only on idx_fit
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
  
  ## Early return if all tests are discrete
  if (length(idx_non_dis) == 0L) {
    return(.pack_result())
  }
  
  ## -------------------------------------------------------------------------
  ## Threshold detection per test (subset)
  ## -------------------------------------------------------------------------
  find_thresh_one <- function(j) {
    .find_gpd_thresh(
      perm_stats    = perm_stats[, j],
      obs_stat      = obs_stats[j],
      tol           = control$tol,
      thresh_method = control$thresh_method,
      thresh0       = control$thresh0,
      exceed0       = control$exceed0,
      exceed0_min   = control$exceed0_min,
      exceed_min    = control$exceed_min,
      thresh_step   = control$thresh_step,
      gof_test      = control$gof_test,
      gof_alpha     = control$gof_alpha,
      ...
    )
  }
  
  if (isTRUE(verbose)) writeLines("Run threshold detection ...")
  
  n_workers <- .choose_workers(cores, length(idx_non_dis), parallel_min)
  run_parallel_thresh <- (n_workers > 1L)
  
  if (run_parallel_thresh) {
    strategy <- if (.Platform$OS.type != "windows" && future::supportsMulticore())
      future::multicore else future::multisession
    old_plan <- future::plan(strategy, workers = n_workers)
    on.exit(future::plan(old_plan), add = TRUE)
    
    res_list_thr <- progressr::with_progress(
      handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
      {
        if (verbose) p <- progressr::progressor(along = idx_non_dis)
        
        future.apply::future_lapply(
          idx_non_dis,
          function(j) {
            res <- find_thresh_one(j)
            if (verbose) p()
            res
          },
          future.packages = c("permApprox","progressr"),
          future.seed = TRUE
        )
      }
    )
    
  } else {
    res_list_thr <- progressr::with_progress(
      handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
      {
        if (verbose) p <- progressr::progressor(along = idx_non_dis)
        
        lapply(idx_non_dis, function(j) {
          res <- find_thresh_one(j)
          if (verbose) p()
          res
        })
      }
    )
  }
  
  thresh[idx_non_dis]   <- vapply(res_list_thr, `[[`, numeric(1),  "thresh")
  n_exceed[idx_non_dis] <- vapply(res_list_thr, `[[`, integer(1), "n_exceed")
  
  idx_no_thr <- idx_non_dis[is.na(thresh[idx_non_dis])]
  status[idx_no_thr] <- "no_threshold"
  
  if (isTRUE(verbose)) writeLines("Done.")
  
  idx_valid <- idx_non_dis[!is.na(thresh[idx_non_dis])]
  
  ## Early return if there is no test with valid threshold
  if (length(idx_valid) == 0L) {
    return(.pack_result())
  }
  
  ## -------------------------------------------------------------------------
  ## Support boundary vector according to constraint (subset)
  ## -------------------------------------------------------------------------
  
  support_boundaries <- rep(NA_real_, n_test)
  sb <- switch(
    control$constraint,
    "support_at_max" = rep(max(obs_stats, na.rm = TRUE), length(idx_valid)),
    "support_at_obs" = obs_stats[idx_valid],
    "unconstrained"  = rep(NA_real_, length(idx_valid))
  )
  support_boundaries[idx_valid] <- sb
  
  ## -------------------------------------------------------------------------
  ## Helper to run GPD fit with given epsilon vector on an index subset
  ## -------------------------------------------------------------------------
  
  .fit_one_gpd <- function(i, eps) {
    obs_i    <- obs_stats[i]
    perm_i   <- perm_stats[, i]
    eps_i    <- eps[i]
    thr_i    <- thresh[i]
    exceed_i <- perm_i[perm_i > thr_i]
    bound_i  <- support_boundaries[i]
    
    fit_res <- fit_gpd(
      data             = exceed_i,
      thresh           = thr_i,
      fit_method       = control$fit_method,
      tol              = control$tol,
      epsilon          = eps_i,
      constraint       = control$constraint,
      support_boundary = bound_i,
      gof_test         = control$gof_test,
      gof_alpha        = control$gof_alpha,
      ...
    )
    
    gof_fail <- is.finite(fit_res$p_value) && (fit_res$p_value < control$gof_alpha)
    
    n_perm_eff <- if (isTRUE(control$include_obs)) n_perm + 1L else n_perm
    p_tail     <- (n_exceed[i] / n_perm_eff) * .pgpd_upper_tail(
      q        = obs_i - thr_i,
      location = 0,
      scale    = fit_res$scale,
      shape    = fit_res$shape
    )
    p_tail_valid <- is.finite(p_tail) && p_tail >= 0 && p_tail <= 1
    
    list(
      p_value     = if (!gof_fail && p_tail_valid) p_tail else NA_real_,
      shape       = fit_res$shape,
      scale       = fit_res$scale,
      gof_p_value = fit_res$p_value,
      epsilon     = fit_res$epsilon,
      ok          = (!gof_fail && p_tail_valid),
      gof_fail    = gof_fail
    )
  }
  
  ## -------------------------------------------------------------------------
  ## Initial epsilon definition and fit (subset)
  ## -------------------------------------------------------------------------
  epsilons <- .define_eps(
    perm_stats   = perm_stats,
    obs_stats    = obs_stats,
    sample_size  = control$sample_size,
    constraint   = control$constraint,
    eps_fun      = control$eps_fun,
    eps_tune     = control$eps_tune,
    eps_args     = control$eps_args
  )
  
  if (isTRUE(verbose)) writeLines("Run GPD fit ...")
  
  n_workers_fit <- .choose_workers(cores, length(idx_valid), parallel_min)
  run_parallel_fit_all <- (n_workers_fit > 1L)
  
  if (run_parallel_fit_all) {
    strategy <- if (.Platform$OS.type != "windows" && future::supportsMulticore())
      future::multicore else future::multisession
    old_plan <- future::plan(strategy, workers = n_workers_fit)
    on.exit(future::plan(old_plan), add = TRUE)
    
    res_list <- progressr::with_progress(
      handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
      {
        if (verbose) p <- progressr::progressor(along = idx_valid)
        future.apply::future_lapply(
          idx_valid,
          function(j) {
            res <- .fit_one_gpd(j, eps = epsilons)
            if (verbose) p()
            res
          },
          future.packages = c("permApprox", "progressr"),
          future.seed = TRUE
        )
      }
    )
    
  } else {
    res_list <- progressr::with_progress(
      handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
      {
        if (verbose) p <- progressr::progressor(along = idx_valid)
        lapply(idx_valid, function(j) {
          res <- .fit_one_gpd(j, eps = epsilons)
          if (verbose) p()
          res
        })
      }
    )
  }
  
  if (isTRUE(verbose)) writeLines("Done.")
  
  ## -------------------------------------------------------------------------
  ## Optional: Adaptive epsilon refinement (zero-guard)
  ## -------------------------------------------------------------------------
  
  supp_constr <- control$constraint %in% c("support_at_obs", "support_at_max")
  
  if (isTRUE(supp_constr) && isTRUE(control$zero_guard)) {
    
    # Index of zero (or underflow) p-values
    .get_zero_idx <- function(pv, idx) {
      min_pos <- .Machine$double.xmin
      z <- idx[which(!is.na(pv[idx]) & pv[idx] <= min_pos)]
      as.integer(z)
    }
    
    # Extract p-values from results data frame
    .get_pvals_from_res <- function(res) {
      vapply(res, `[[`, numeric(1), "p_value")
    } 
    
    # Identify current number of zero p-values
    pvals_tmp <- p_value
    pvals_tmp[idx_valid] <- .get_pvals_from_res(res_list)
    
    # Indices of zero p-values
    zero_idx <- .get_zero_idx(pvals_tmp, idx_valid)
    
    if (length(zero_idx) > 0L) {
      
      ## --- Helpers ---------------------------------------------------------
      
      # Run GPD fit for a subset without progress bar
      .run_gpd_fit <- function(epsilons, idx_subset = idx_valid) {
        
        n_workers_fit <- .choose_workers(cores, length(idx_subset), parallel_min)
        run_parallel_fit <- (n_workers_fit > 1L)
        
        if (run_parallel_fit) {
          strategy <- if (.Platform$OS.type != "windows" && future::supportsMulticore()) 
            future::multicore else future::multisession
          old_plan <- future::plan(strategy, workers = n_workers_fit)
          on.exit(future::plan(old_plan), add = TRUE)
          
          res_list <- future.apply::future_lapply(
            idx_subset, 
            function(i) .fit_one_gpd(i, eps = epsilons), 
            future.packages = c("permApprox","progressr"),
            future.seed = TRUE
          )
        } else {
          res_list <- lapply(idx_subset,
                             function(i) .fit_one_gpd(i, eps = epsilons))
        }
        res_list
      }
      
      # Count the number of zeros 
      # (or values below the machine's smalles possible floating-point number)
      .count_zeros <- function(pv) {
        min_pos <- .Machine$double.xmin
        sum(is.finite(pv) & pv <= min_pos)
      }
      
      ## Digits to show for the tuning parameter based on bisect tolerance
      .digits_from_tol <- function(tol) {
        if (!is.finite(tol) || tol <= 0) return(6L)
        d <- ceiling(-log10(tol)) + 1L  # +1 to make changes at the tolerance visible
        d <- max(0L, min(16L, d))       # clamp to sane bounds
        as.integer(d)
      }
      
      ## Evaluate a candidate tuning parameter on a subset only
      .eval_tp_subset <- function(tp, idx_subset) {
        eps_try <- .define_eps(
          perm_stats   = perm_stats,
          obs_stats    = obs_stats,
          sample_size  = control$sample_size,
          constraint   = control$constraint,
          eps_fun      = control$eps_fun,
          eps_tune     = tp,
          eps_args     = control$eps_args
        )
        res_sub <- .run_gpd_fit(eps_try, idx_subset = idx_subset)
        pv_sub  <- vapply(res_sub, `[[`, numeric(1), "p_value")
        list(res = res_sub, pvals = pv_sub)
      }
      
      ## --- Configuration parameters ----------------------------------------
      tp0     <- control$eps_tune
      step    <- control$eps_retry$step_init
      grow    <- control$eps_retry$grow
      max_exp <- control$eps_retry$max_expand_iter
      max_bis <- control$eps_retry$bisect_iter_max
      bis_tol <- control$eps_retry$bisect_tol
      
      tp_digits <- .digits_from_tol(bis_tol)
      tp_fmt    <- function(x) sprintf(paste0("%.", tp_digits, "f"), x)
      
      # Work only on the subset that currently underflows
      need_idx <- zero_idx
      
      if (isTRUE(verbose)) {
        message(sprintf(
          "Zero-guard: start with tp=%s; zeros: %d (of %d selected tests)",
          tp_fmt(tp0), length(zero_idx), length(idx_valid)
        ))
      }
      
      # --- EXPANSION PHASE (large increasing steps) -------------------------
      if (isTRUE(verbose)) message("Zero-guard: expanding tuning parameter ...")
      tp_curr <- tp0
      tp_lo   <- tp0            # last tp that still produces zeros
      tp_hi   <- NA_real_       # first tp with no zeros in need_idx
      exp_it  <- 0L
      
      repeat {
        exp_it <- exp_it + 1L
        tp_try <- tp_curr + step
        ev <- .eval_tp_subset(tp_try, need_idx)
        zeros_now <- .count_zeros(ev$pvals)
        if (isTRUE(verbose)) {
          message(sprintf("  expand %02d: tp=%s; zeros in subset: %d",
                          exp_it, tp_fmt(tp_try), zeros_now))
        }
        
        if (zeros_now == 0L) {
          tp_hi <- tp_try   # bracket found
          if (isTRUE(verbose)) {
            message(sprintf("  bracket found at tp=%s (subset zeros=0)", tp_fmt(tp_hi)))
          }
          break
        } else {
          # still zeros → move forward and inflate step
          tp_lo   <- tp_try
          tp_curr <- tp_try
          step    <- step * grow
          
          # shrink subset to those still zero to save runtime
          need_idx_old <- need_idx
          need_idx <- need_idx[which(is.finite(ev$pvals) & 
                                       ev$pvals <= .Machine$double.xmin)]
          
          if (isTRUE(verbose) && length(need_idx) != length(need_idx_old)) {
            message(sprintf("    refining subset to the %d zero(s)", length(need_idx)))
          }
          if (exp_it >= max_exp || length(need_idx) == 0L) {
            if (isTRUE(verbose)) {
              message("  expansion stop: reached limit or no subset left")
            }
            break
          }
        }
      }
      
      # If expansion never cleared zeros, we stop here with best effort tp_curr
      final_tp <- tp_curr
      
      # --- BISECTION PHASE (minimal safe tp) --------------------------------
      if (!is.na(tp_hi)) {
        if (isTRUE(verbose)) 
          message("Zero-guard: bisecting tuning parameter ...")
        tp_left  <- tp_lo   # zeros occur
        tp_right <- tp_hi   # no zeros
        bis_it   <- 0L
        ref_idx  <- zero_idx  # start from original zero set; subset fits only
        
        while ((tp_right - tp_left > bis_tol) &&
               bis_it < max_bis &&
               length(ref_idx) > 0L) {
          bis_it <- bis_it + 1L
          tp_mid <- 0.5 * (tp_left + tp_right)
          ev <- .eval_tp_subset(tp_mid, ref_idx)
          zeros_mid <- .count_zeros(ev$pvals)
          if (isTRUE(verbose)) {
            message(sprintf("  bisect %02d: tp=%s; zeros in subset: %d",
                            bis_it, tp_fmt(tp_mid), zeros_mid))
          }
          
          if (zeros_mid > 0L) {
            tp_left <- tp_mid
            # keep only those still zero at tp_mid to reduce work further
            ref_idx_old <- ref_idx
            ref_idx <- ref_idx[which(is.finite(ev$pvals) & 
                                       ev$pvals <= .Machine$double.xmin)]
            if (isTRUE(verbose) && length(ref_idx) != length(ref_idx_old)) {
              message(sprintf("    refining subset to the %d zero(s)", 
                              length(ref_idx)))
            }
          } else {
            tp_right <- tp_mid
            # no need to persist interim fits here; we will refit all at the end
          }
        }
        final_tp <- tp_right
      }
      
      # --- Final refit for ALL selected tests with final_tp -----------------
      eps_final <- .define_eps(
        perm_stats   = perm_stats,
        obs_stats    = obs_stats,
        sample_size  = control$sample_size,
        constraint   = control$constraint,
        eps_fun      = control$eps_fun,
        eps_tune     = final_tp,
        eps_args     = control$eps_args
      )
      if (isTRUE(verbose)) {
        message("Zero-guard: refitting all selected tests with tuning parameter = ", 
                tp_fmt(final_tp))
      }
      
      if (run_parallel_fit_all) {
        strategy <- if (.Platform$OS.type != "windows" && 
                        future::supportsMulticore())
          future::multicore else future::multisession
        old_plan <- future::plan(strategy, workers = n_workers_fit)
        on.exit(future::plan(old_plan), add = TRUE)
        
        res_list <- progressr::with_progress(
          handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
          {
            if (verbose) p <- progressr::progressor(along = idx_valid)
            future.apply::future_lapply(
              idx_valid,
              function(j) {
                res <- .fit_one_gpd(j, eps = eps_final)
                if (verbose) p()
                res
              },
              future.packages = c("permApprox", "progressr"),
              future.seed = TRUE
            )
          }
        )
        
      } else {
        res_list <- progressr::with_progress(
          handlers = list(progressr::handler_txtprogressbar(clear = TRUE)),
          {
            if (verbose) p <- progressr::progressor(along = idx_valid)
            lapply(idx_valid, function(j) {
              res <- .fit_one_gpd(j, eps = eps_final)
              if (verbose) p()
              res
            })
          }
        )
      }
      
      pvals_tmp[idx_valid] <- vapply(res_list, `[[`, numeric(1), "p_value")
      
      if (isTRUE(verbose)) {
        zeros_final <- .count_zeros(pvals_tmp[idx_valid])
        message(sprintf("Zero-guard: %d zeros (of %d tests) after refit", 
                        zeros_final, length(idx_valid)))
      }
    }
  }
  
  ## -------------------------------------------------------------------------
  ## Collect final results
  ## -------------------------------------------------------------------------
  p_value[idx_valid]     <- vapply(res_list, `[[`, numeric(1), "p_value")
  shape[idx_valid]       <- vapply(res_list, `[[`, numeric(1), "shape")
  scale[idx_valid]       <- vapply(res_list, `[[`, numeric(1), "scale")
  gof_p_value[idx_valid] <- vapply(res_list, `[[`, numeric(1), "gof_p_value")
  epsilon[idx_valid]     <- vapply(res_list, `[[`, numeric(1), "epsilon")
  
  ok_vec   <- vapply(res_list, `[[`, logical(1), "ok")
  gof_fail <- vapply(res_list, `[[`, logical(1), "gof_fail")
  status[idx_valid[ ok_vec]]             <- "success"
  status[idx_valid[!ok_vec &  gof_fail]] <- "gof_reject"
  status[idx_valid[!ok_vec & !gof_fail]] <- "fit_failed"
  
  for (i in idx_valid) {
    if (identical(as.character(status[i]), "success")) {
      perm_stats_fit[[i]] <- perm_stats[perm_stats[, i] > thresh[i], i]
    }
  }
  
  .pack_result()
}
