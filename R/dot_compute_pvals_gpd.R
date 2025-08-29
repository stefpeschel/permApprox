#' Compute p-values via GPD tail approximation (constrained or unconstrained)
#'
#' @keywords internal
#' @param obs_stats Numeric vector of observed statistics (length = n_test).
#' @param perm_stats Numeric matrix of permutation statistics with shape B x n_test.
#'   If `control$include_obs` is TRUE, `obs_stats` will be row-bound to `perm_stats`
#'   before thresholding / fitting.
#' @param p_empirical Numeric vector of empirical p-values (length = n_test).
#' @param fit_thresh Numeric scalar; threshold under which GPD fit is attempted.
#' @param control List created by `make_gpd_ctrl()`.
#' @return A data.frame with columns p_empirical, p_value (final), thresh, n_exceed,
#'   shape, scale, epsilon, gof_p_value, status (factor), discrete (logical).
#'   The list of exceedances used for fitting is attached as attribute `perm_stats_fit`.
#'
.compute_pvals_gpd <- function(obs_stats,
                               perm_stats,
                               p_empirical,
                               fit_thresh,
                               control,
                               cores,
                               verbose,
                               ...) {
  ## Sizes
  n_test <- ncol(perm_stats)
  
  if (isTRUE(control$include_obs)) {
    stopifnot(length(obs_stats) == ncol(perm_stats))
    perm_stats <- rbind(perm_stats, obs_stats)
  }
  
  n_perm <- nrow(perm_stats)
  
  ## -------------------------------------------------------------------------
  ## Select tests to fit
  ## -------------------------------------------------------------------------
  idx_fit <- which(!is.na(p_empirical) & 
                     (p_empirical < fit_thresh) & 
                     (obs_stats > 0))
  
  ## -------------------------------------------------------------------------
  ## Prepare result data.frame (start with empirical p-values)
  ## -------------------------------------------------------------------------
  status_levels <- c("not_selected", "success", "discrete", "no_threshold",
                     "gof_reject", "fit_failed")
  
  res_df <- data.frame(
    p_value     = as.numeric(p_empirical),
    p_empirical = as.numeric(p_empirical),
    thresh      = rep(NA_real_,    n_test),
    n_exceed    = rep(NA_integer_, n_test),
    shape       = rep(NA_real_,    n_test),
    scale       = rep(NA_real_,    n_test),
    epsilon     = rep(NA_real_,    n_test),
    gof_p_value = rep(NA_real_,    n_test),
    status      = factor(rep("not_selected", n_test), 
                         levels = status_levels, ordered = TRUE),
    discrete    = rep(FALSE, n_test),
    stringsAsFactors = FALSE
  )
  
  ## If nothing to fit, attach empty perm_stats_fit and return
  perm_stats_fit <- vector("list", n_test)
  if (length(idx_fit) == 0L) {
    attr(res_df, "perm_stats_fit") <- perm_stats_fit
    return(res_df)
  }
  
  ## -------------------------------------------------------------------------
  ## Discreteness screening only on idx_fit
  ## -------------------------------------------------------------------------
  uniq_counts <- vapply(idx_fit, function(i) length(unique(perm_stats[, i])), integer(1))
  is_discrete_sub <- uniq_counts < max(1L, floor(n_perm * 0.02))
  res_df$discrete[idx_fit] <- is_discrete_sub
  res_df$status[idx_fit[is_discrete_sub]] <- "discrete"
  idx_non_dis <- idx_fit[!is_discrete_sub]
  
  ## Early return if all tests are discrete
  if (length(idx_non_dis) == 0L) {
    attr(res_df, "perm_stats_fit") <- perm_stats_fit
    return(res_df)
  }
  
  ## -------------------------------------------------------------------------
  ## Threshold detection per test (subset)
  ## -------------------------------------------------------------------------
  find_thresh_one <- function(j) {
    .find_gpd_thresh(
      perm_stats    = perm_stats[, j],
      obs_stats     = obs_stats[j],
      tol           = control$tol,
      thresh_method = control$thresh_method,
      thresh0       = control$thresh0,
      exceed0       = control$exceed0,
      exceed_min    = control$exceed_min,
      thresh_step   = control$thresh_step,
      gof_test      = control$gof_test,
      gof_alpha     = control$gof_alpha,
      ...
    )
  }
  
  if (isTRUE(verbose)) message("Run threshold detection ...")
  
  n_workers <- min(max(1L, cores), max(1L, length(idx_non_dis)))
  run_parallel_thresh <- (n_workers > 1L)
  
  if (run_parallel_thresh) {
    strategy <- if (.Platform$OS.type != "windows" && 
                    future::supportsMulticore())
      future::multicore else future::multisession
    
    old_plan <- future::plan(strategy, 
                             workers = min(cores, length(idx_non_dis)))
    on.exit(future::plan(old_plan), add = TRUE)
    
    res_list_thr <- progressr::with_progress({
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
    })
  } else {
    res_list_thr <- with_progress({
      if (verbose) p <- progressr::progressor(along = idx_non_dis)
      lapply(idx_non_dis, function(j) {
        res <- find_thresh_one(j)
        if (verbose) p()
        res
      })
    })
  }
  
  thresh_vec   <- res_df$thresh
  n_exceed_vec <- res_df$n_exceed
  thresh_vec[idx_non_dis]   <- vapply(res_list_thr, `[[`, numeric(1),  "thresh")
  n_exceed_vec[idx_non_dis] <- vapply(res_list_thr, `[[`, integer(1), "n_exceed")
  
  res_df$thresh   <- thresh_vec
  res_df$n_exceed <- n_exceed_vec
  
  idx_no_thr <- idx_non_dis[is.na(thresh_vec[idx_non_dis])]
  res_df$status[idx_no_thr] <- "no_threshold"
  
  if (isTRUE(verbose)) message("Done.")
  
  idx_valid <- idx_non_dis[!is.na(thresh_vec[idx_non_dis])]
  
  ## Early return if there is no test with valid threshold
  if (length(idx_valid) == 0L) {
    attr(res_df, "perm_stats_fit") <- perm_stats_fit
    return(res_df)
  }
  
  ## -------------------------------------------------------------------------
  ## Support boundary vector according to constraint (subset)
  ## -------------------------------------------------------------------------
  
  support_boundaries <- rep(NA_real_, n_test)
  sb <- switch(control$constraint,
               "support_at_max" = rep(max(obs_stats, na.rm = TRUE), length(idx_valid)),
               "support_at_obs" = obs_stats[idx_valid],
               "unconstrained"  = rep(NA_real_, length(idx_valid)))
  support_boundaries[idx_valid] <- sb
  
  ## -------------------------------------------------------------------------
  ## Helper to run GPD fit with given epsilon vector on an index subset
  ## -------------------------------------------------------------------------
  
  .fit_one_gpd <- function(i, eps) {
    obs_i    <- obs_stats[i]
    perm_i   <- perm_stats[, i]
    eps_i    <- eps[i]
    thr_i    <- thresh_vec[i]
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
    
    p_tail <- (n_exceed_vec[i] / n_perm) * .pgpd_upper_tail(
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
    eps_rule     = control$eps_rule,
    eps_par      = control$eps_par
  )
  
  if (isTRUE(verbose)) message("Run GPD fit ...")
  
  n_workers_fit <- if (cores > 1L && length(idx_valid) > cores) {
    min(cores, length(idx_valid))
  } else {
    1L
  }
  run_parallel_fit_all <- (n_workers_fit > 1L)
  
  if (run_parallel_fit_all) {
    strategy <- if (.Platform$OS.type != "windows" && future::supportsMulticore()) 
      future::multicore else future::multisession
    old_plan <- future::plan(strategy, workers = n_workers_fit)
    on.exit(future::plan(old_plan), add = TRUE)
    
    res_list <- progressr::with_progress({
      if (verbose) p <- progressr::progressor(along = idx_valid)
      
      future.apply::future_lapply(
        idx_valid, 
        function(j) {
          res <- .fit_one_gpd(j, eps = epsilons)
          if (verbose) p()
          res
        },
        future.packages = c("permApprox","progressr"),
        future.seed = TRUE
      )
    })
    
  } else {
    res_list <- progressr::with_progress({
      if (verbose) p <- progressr::progressor(along = idx_valid)
      lapply(idx_valid, function(j) {
        res <- .fit_one_gpd(j, eps = epsilons)
        if (verbose) p()
        res
      })
    })
  }
  
  if (isTRUE(verbose)) message("Done.")
  
  ## -------------------------------------------------------------------------
  ## Optional: Adaptive epsilon refinement
  ## -------------------------------------------------------------------------
  
  constr_fit <- control$constraint != "unconstrained"
  
  if (isTRUE(constr_fit) && isTRUE(control$zero_guard)) {
    
    # Index of zero (or underflow) p-values
    .get_zero_idx <- function(pv, idx) {
      min_pos <- .Machine$double.xmin
      z <- idx[which(!is.na(pv[idx]) & pv[idx] <= min_pos)]
      as.integer(z)
    }
    
    # Extract p-values from results data frame
    .get_pvals_from_res <- function(res, idxs) {
      vapply(res, `[[`, numeric(1), "p_value")
    } 
    
    # Identify current number of zero p-values
    pvals_tmp <- res_df$p_value
    pvals_tmp[idx_valid] <- .get_pvals_from_res(res_list, idx_valid)
    
    # Indices of zero p-values
    zero_idx <- .get_zero_idx(pvals_tmp, idx_valid)
    
    if (length(zero_idx) > 0L) {
      
      ## --- Helpers ---------------------------------------------------------
      
      # Run GPD fit for a subset without progress bar
      .run_gpd_fit <- function(epsilons, idx_subset = idx_valid) {
        
        n_workers_fit <- if (cores > 1L && length(idx_subset) > cores) {
          min(cores, length(idx_subset))
        } else {
          1L
        }
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
      
      ## Digits to show for target_factor based on bisect tolerance
      .digits_from_tol <- function(tol) {
        if (!is.finite(tol) || tol <= 0) return(6L)
        d <- ceiling(-log10(tol)) + 1L  # +1 to make changes at the tolerance visible
        d <- max(0L, min(16L, d))       # clamp to sane bounds
        as.integer(d)
      }
      
      ## Evaluate a candidate target_factor on a subset only
      .eval_tf_subset <- function(tf, idx_subset) {
        eps_try <- .define_eps(
          perm_stats   = perm_stats,
          obs_stats    = obs_stats,
          sample_size  = control$sample_size,
          constraint   = control$constraint,
          eps_rule     = control$eps_rule,
          eps_par      = tf
        )
        res_sub <- .run_gpd_fit(eps_try, idx_subset = idx_subset)
        pv_sub  <- vapply(res_sub, `[[`, numeric(1), "p_value")
        list(res = res_sub, pvals = pv_sub)
      }
      
      ## --- Configuration parameters ----------------------------------------
      tf0     <- control$eps_par
      step    <- control$eps_retry$step_init
      grow    <- control$eps_retry$grow
      max_exp <- control$eps_retry$max_expand_iter
      max_bis <- control$eps_retry$bisect_iter_max
      bis_tol <- control$eps_retry$bisect_tol
      
      tf_digits <- .digits_from_tol(bis_tol)
      tf_fmt    <- function(x) sprintf(paste0("%.", tf_digits, "f"), x)
      
      # Work only on the subset that currently underflows
      need_idx <- zero_idx
      
      if (isTRUE(verbose)) {
        message(sprintf(
          "Zero-guard: start with tf=%s; zeros: %d (of %d selected tests)",
          tf_fmt(tf0), length(zero_idx), length(idx_valid)
        ))
      }
      
      # --- EXPANSION PHASE (large increasing steps) -------------------------
      if (isTRUE(verbose)) message("Zero-guard: expanding target_factor ...")
      tf_curr <- tf0
      tf_lo   <- tf0            # last tf that still produces zeros
      tf_hi   <- NA_real_       # first tf with no zeros in need_idx
      exp_it  <- 0L
      
      repeat {
        exp_it <- exp_it + 1L
        tf_try <- tf_curr + step
        ev <- .eval_tf_subset(tf_try, need_idx)
        zeros_now <- .count_zeros(ev$pvals)
        if (isTRUE(verbose)) {
          message(sprintf("  expand %02d: tf=%s; zeros in subset: %d",
                          exp_it, tf_fmt(tf_try), zeros_now))
        }
        
        if (zeros_now == 0L) {
          tf_hi <- tf_try   # bracket found
          if (isTRUE(verbose)) {
            message(sprintf("  bracket found at tf=%s (subset zeros=0)", tf_fmt(tf_hi)))
          }
          break
        } else {
          # still zeros → move forward and inflate step
          tf_lo   <- tf_try
          tf_curr <- tf_try
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
      
      # If expansion never cleared zeros, we stop here with best effort tf_curr
      final_tf <- tf_curr
      
      # --- BISECTION PHASE (minimal safe tf) --------------------------------
      if (!is.na(tf_hi)) {
        if (isTRUE(verbose)) 
          message("Zero-guard: bisecting target_factor ...")
        tf_left  <- tf_lo   # zeros occur
        tf_right <- tf_hi   # no zeros
        bis_it   <- 0L
        ref_idx  <- zero_idx  # start from original zero set; subset fits only
        
        while ((tf_right - tf_left > bis_tol) && bis_it < max_bis && length(ref_idx) > 0L) {
          bis_it <- bis_it + 1L
          tf_mid <- 0.5 * (tf_left + tf_right)
          ev <- .eval_tf_subset(tf_mid, ref_idx)
          zeros_mid <- .count_zeros(ev$pvals)
          if (isTRUE(verbose)) {
            message(sprintf("  bisect %02d: tf=%s; zeros in subset: %d",
                            bis_it, tf_fmt(tf_mid), zeros_mid))
          }
          
          if (zeros_mid > 0L) {
            tf_left <- tf_mid
            # keep only those still zero at tf_mid to reduce work further
            ref_idx_old <- ref_idx
            ref_idx <- ref_idx[which(is.finite(ev$pvals) & 
                                       ev$pvals <= .Machine$double.xmin)]
            if (isTRUE(verbose) && length(ref_idx) != length(ref_idx_old)) {
              message(sprintf("    refining subset to the %d zero(s)", 
                              length(ref_idx)))
            }
          } else {
            tf_right <- tf_mid
            # no need to persist interim fits here; we will refit all at the end
          }
        }
        final_tf <- tf_right
      }
      
      # --- Final refit for ALL selected tests with final_tf -----------------
      eps_final <- .define_eps(
        perm_stats   = perm_stats,
        obs_stats    = obs_stats,
        sample_size  = control$sample_size,
        constraint   = control$constraint,
        eps_rule     = control$eps_rule,
        eps_par      = final_tf
      )
      if (isTRUE(verbose)) {
        message("Zero-guard: refitting all selected tests with target_factor = ", 
                tf_fmt(final_tf))
      }
      
      if (run_parallel_fit_all) {
        strategy <- if (.Platform$OS.type != "windows" && future::supportsMulticore()) 
          future::multicore else future::multisession
        old_plan <- future::plan(strategy, workers = n_workers_fit)
        on.exit(future::plan(old_plan), add = TRUE)
        
        res_list <- progressr::with_progress({
          if (verbose) p <- progressr::progressor(along = idx_valid)
          
          future.apply::future_lapply(
            idx_valid, 
            function(j) {
              res <- .fit_one_gpd(j, eps = eps_final)
              if (verbose) p()
              res
            },
            future.packages = c("permApprox","progressr"),
            future.seed = TRUE
          )
        })
        
      } else {
        res_list <- progressr::with_progress({
          if (verbose) p <- progressr::progressor(along = idx_valid)
          lapply(idx_valid, function(j) {
            res <- .fit_one_gpd(j, eps = eps_final)
            if (verbose) p()
            res
          })
        })
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
  ## Collect results into data.frame
  ## -------------------------------------------------------------------------
  res_df$p_value[idx_valid]      <- vapply(res_list, `[[`, numeric(1), "p_value")
  res_df$shape[idx_valid]        <- vapply(res_list, `[[`, numeric(1),  "shape")
  res_df$scale[idx_valid]        <- vapply(res_list, `[[`, numeric(1),  "scale")
  res_df$gof_p_value[idx_valid]  <- vapply(res_list, `[[`, numeric(1),  "gof_p_value")
  res_df$epsilon[idx_valid]      <- vapply(res_list, `[[`, numeric(1),  "epsilon")
  
  # Replace NAs
  idx_na <- which(is.na(res_df$p_value))
  res_df$p_value[idx_na] <- res_df$p_empirical[idx_na]
  
  ok_vec   <- vapply(res_list, `[[`, logical(1), "ok")
  gof_fail <- vapply(res_list, `[[`, logical(1), "gof_fail")
  res_df$status[idx_valid[ ok_vec]]             <- "success"
  res_df$status[idx_valid[!ok_vec &  gof_fail]] <- "gof_reject"
  res_df$status[idx_valid[!ok_vec & !gof_fail]] <- "fit_failed"
  
  for (i in idx_valid) {
    if (identical(as.character(res_df$status[i]), "success")) {
      perm_stats_fit[[i]] <- perm_stats[perm_stats[, i] > res_df$thresh[i], i]
    }
  }
  
  attr(res_df, "perm_stats_fit") <- perm_stats_fit
  res_df
}