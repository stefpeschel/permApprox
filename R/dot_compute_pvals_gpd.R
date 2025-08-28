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
  
  if (isTRUE(control$verbose)) message("Run threshold detection ...")
  
  n_workers <- min(max(1L, control$cores), max(1L, length(idx_non_dis)))
  run_parallel_thresh <- (n_workers > 1L)
  
  if (run_parallel_thresh) {
    strategy <- if (.Platform$OS.type != "windows" && 
                    future::supportsMulticore())
      future::multicore else future::multisession
    
    old_plan <- future::plan(strategy, 
                             workers = min(control$cores, length(idx_non_dis)))
    on.exit(future::plan(old_plan), add = TRUE)
    
    res_list_thr <- with_progress({
      if (control$verbose) p <- progressor(along = idx_non_dis)
      future.apply::future_lapply(
        idx_non_dis,
        function(j) {
          res <- find_thresh_one(j)
          if (control$verbose) p()
          res
        },
        future.packages = "permApprox",
        future.seed = TRUE
      )
    })
  } else {
    res_list_thr <- with_progress({
      if (control$verbose) p <- progressor(along = idx_non_dis)
      lapply(idx_non_dis, function(j) {
        res <- find_thresh_one(j)
        if (control$verbose) p()
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
  
  if (isTRUE(control$verbose)) message("Done.")
  
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
  run_gpd_fit <- function(epsilons, idx_subset = idx_valid) {
    fit_one_gpd <- function(i) {
      obs_i    <- obs_stats[i]
      perm_i   <- perm_stats[, i]
      eps_i    <- epsilons[i]
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
    
    n_workers_fit <- if (control$cores > 1L && length(idx_subset) > control$cores) {
      min(control$cores, length(idx_subset))
    } else {
      1L
    }
    run_parallel_fit <- (n_workers_fit > 1L)
    
    if (run_parallel_fit) {
      strategy <- if (.Platform$OS.type != "windows" && future::supportsMulticore()) 
        future::multicore else future::multisession
      old_plan <- future::plan(strategy, workers = n_workers_fit)
      on.exit(future::plan(old_plan), add = TRUE)
      res_list <- future.apply::future_lapply(idx_subset, function(i) 
        fit_one_gpd(i), future.seed = TRUE)
    } else {
      pb <- utils::txtProgressBar(min = 0, max = length(idx_subset), style = 3)
      on.exit(close(pb), add = TRUE)
      k <- 0L
      res_list <- lapply(idx_subset, 
                         function(i) { 
                           r <- fit_one_gpd(i)
                           k <<- k + 1L
                           utils::setTxtProgressBar(pb, k); r 
                         })
    }
    res_list
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
  
  if (isTRUE(control$verbose)) message("Run GPD fit ...")
  res_list <- run_gpd_fit(epsilons)
  if (isTRUE(control$verbose)) message("Done.")
  browser()
  ## -------------------------------------------------------------------------
  ## Adaptive epsilon refinement (Fibonacci-ish bracketing: golden growth + bisection)
  ## -------------------------------------------------------------------------
  get_pvals_from_res <- function(res, idxs) vapply(res, `[[`, numeric(1), "p_value")
  
  pvals_tmp <- res_df$p_value
  pvals_tmp[idx_valid] <- get_pvals_from_res(res_list, idx_valid)
  
  .get_zero_idx <- function(pv, idx) {
    min_pos <- .Machine$double.xmin
    z <- idx[which(!is.na(pv[idx]) & pv[idx] <= min_pos)]
    as.integer(z)
  }
  
  zero_idx <- .get_zero_idx(pvals_tmp, idx_valid)
  
  if (isTRUE(control$zero_guard) && length(zero_idx) > 0L) {
    ## Retry controls with sensible defaults
    er <- control$eps_retry
    step_init        <- er$step_init
    grow_phi         <- er$grow
    max_expand_iter  <- er$max_expand_iter
    bisect_iter_max  <- er$bisect_iter_max
    bisect_tol       <- er$bisect_tol
    
    eps_base_scalar  <- as.numeric(control$eps_par)
    if (!is.finite(eps_base_scalar)) eps_base_scalar <- 0
    
    ## ------------------ Expansion (golden growth) ------------------
    step_cur <- step_init
    eps_low  <- eps_base_scalar
    eps_high <- NA_real_
    
    expand_iter <- 0L
    still_zeros <- length(zero_idx) > 0L
    while (still_zeros && expand_iter < max_expand_iter) {
      expand_iter <- expand_iter + 1L
      prev_scalar <- eps_low
      cand_scalar <- eps_low + step_cur
      eps_par_try <- cand_scalar
      
      epsilons <- .define_eps(
        perm_stats   = perm_stats,
        obs_stats    = obs_stats,
        sample_size  = control$sample_size,
        constraint   = control$constraint,
        eps_rule     = control$eps_rule,
        eps_par      = eps_par_try
      )
      
      res_list_zero <- run_gpd_fit(epsilons, idx_subset = zero_idx)
      
      ## Update p-values
      pvals_tmp[zero_idx] <- vapply(res_list_zero, `[[`, numeric(1), "p_value")
      
      zero_idx_new <- .get_zero_idx(pvals_tmp, idx_valid)
      if (length(zero_idx_new) == 0L) {
        ## bracket found after this probe
        eps_high <- cand_scalar
        eps_low  <- prev_scalar
        
        ## Optional immediate backoff: try half the step to catch overshoot
        mid_scalar   <- 0.5 * (eps_low + eps_high)
        eps_par_try2 <- mid_scalar
        
        epsilons2 <- .define_eps(
          perm_stats   = perm_stats,
          obs_stats    = obs_stats,
          sample_size  = control$sample_size,
          constraint   = control$constraint,
          eps_rule     = control$eps_rule,
          eps_par      = eps_par_try2
        )
        
        ## Refit only the *current zero set* for speed.
        res_list_zero2 <- run_gpd_fit(epsilons2, idx_subset = zero_idx)
        
        ## Update p-values
        pvals_tmp[zero_idx] <- vapply(res_list_zero2, `[[`, numeric(1), "p_value")
        
        zero_idx_new2 <- .get_zero_idx(pvals_tmp, idx_valid)
        if (length(zero_idx_new2) == 0L) {
          ## Half-step still clears zeros → tighten upper bracket (go further back).
          eps_high <- mid_scalar
          ## Note: zero_idx remains the last failing set for subsequent bisection.
        } else {
          ## Half-step brought zeros back → move low up and keep these as the active zero set.
          eps_low  <- mid_scalar
          zero_idx <- zero_idx_new2
        }
        break
      } else {
        ## continue expanding from the new baseline
        eps_low   <- cand_scalar
        zero_idx  <- zero_idx_new
        step_cur  <- step_cur * grow_phi
      }
    }
    
    ## ------------------ Bisection within [eps_low, eps_high] ------------------
    if (is.finite(eps_high) && (eps_high - eps_low) > bisect_tol) {
      b_it <- 0L
      while (b_it < bisect_iter_max && length(zero_idx) > 0L && (eps_high - eps_low) > bisect_tol) {
        b_it <- b_it + 1L
        mid <- 0.5 * (eps_low + eps_high)
        eps_par_try <- .set_eps_scalar(control$eps_par, mid)
        
        epsilons <- .define_eps(
          perm_stats   = perm_stats,
          obs_stats    = obs_stats,
          sample_size  = control$sample_size,
          constraint   = control$constraint,
          eps_rule     = control$eps_rule,
          eps_par      = eps_par_try
        )
        
        res_list_zero <- run_gpd_fit(epsilons, idx_subset = zero_idx)
        pvals_tmp[zero_idx] <- vapply(res_list_zero, `[[`, numeric(1), "p_value")
        
        zero_idx_new <- .get_zero_idx(pvals_tmp, idx_valid)
        if (length(zero_idx_new) == 0L) {
          eps_high <- mid
        } else {
          eps_low  <- mid
        }
        zero_idx <- zero_idx_new
      }
      ## Adopt the minimal epsilon that clears zeros (eps_high if finite, else eps_low)
      final_scalar <- if (is.finite(eps_high)) eps_high else eps_low
      control$eps_par <- .set_eps_scalar(control$eps_par, final_scalar)
    } else if (length(zero_idx) == 0L) {
      ## Zeros already cleared during expansion; use eps_low->eps_high bracket high
      final_scalar <- if (is.finite(eps_high)) eps_high else eps_low
      control$eps_par <- .set_eps_scalar(control$eps_par, final_scalar)
    } else {
      ## No bracket found (hit max expansion); keep last scalar (eps_low)
      control$eps_par <- .set_eps_scalar(control$eps_par, eps_low)
    }
    
    ## Refit of all valid tests at the final epsilon
    if (isTRUE(control$verbose)) 
      message("Refitting all tests with final epsilon after zero-guard ...")
    
    epsilons <- .define_eps(
      perm_stats   = perm_stats,
      obs_stats    = obs_stats,
      sample_size  = control$sample_size,
      constraint   = control$constraint,
      eps_rule     = control$eps_rule,
      eps_par      = control$eps_par
    )
    res_list <- run_gpd_fit(epsilons)
    pvals_tmp[idx_valid] <- vapply(res_list, `[[`, numeric(1),  "p_value")
  }
  
  ## -------------------------------------------------------------------------
  ## Collect results into data.frame
  ## -------------------------------------------------------------------------
  res_df$p_value[idx_valid]      <- pvals_tmp[idx_valid]
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
