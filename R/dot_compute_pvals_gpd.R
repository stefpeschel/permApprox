#' @title Compute p-values via GPD tail approximation
#'
#' @import progressr future future.apply goftest
#' @keywords internal
.compute_pvals_gpd <- function(p_empirical,
                               obs_stats,
                               perm_stats,
                               n_test,
                               fit_thresh,
                               alternative,
                               control,
                               ...) {
  
  # Unpack control arguments
  cores <- control$cores

  # Empirical p-values and number of permutations
  p_values  <- p_empirical
  n_perm <- ncol(perm_stats)

  # Transform test statistics for tail modeling
  transformed <- lapply(seq_len(n_test), function(i) {
    .transform_stats(perm_stats = perm_stats[i, ],
                     obs_stats  = obs_stats[i],
                     alternative = alternative)
  })
  trans_obs  <- vapply(transformed, `[[`, numeric(1), "obs_stats")
  trans_perm <- lapply(transformed, `[[`, "perm_stats")

  # Determine which tests to fit
  idx_fit <- which(p_values < fit_thresh & trans_obs > 0)

  # Setup future plan
  if (cores > 1) {
    # multisession works across platforms; adjust if you prefer multicore
    future::plan(future::multisession, workers = cores)
  } else {
    future::plan(future::sequential)
  }
  
  # Define the fit function for one test
  fit_one <- function(i, control, trans_obs, trans_perm) {
    out <- list()

    obs_i  <- trans_obs[i]
    perm_i <- trans_perm[[i]]

    # Include observed stat if requested
    if (control$include_obs) perm_i <- c(obs_i, perm_i)
    
    if (length(unique(perm_i)) < n_perm * 0.02) { # Check for discreteness
      if (control$verbose) {
        message(paste0("In test ", i, ": Data too discrete to fit a GPD. ", 
                       "Empirical p-value used."))
      }
      
      out$thresh       <- NA
      out$n_exceed     <- NA
      out$zero_replaced <- FALSE
      out$method_used   <- "empirical"
      out$p_value       <- p_empirical[i]
      out$shape         <- NA
      out$scale         <- NA
      out$gof_p_value   <- NA
      out$epsilon       <- NA
      
    } else {
      # Compute threshold
      thresh_res <- .find_gpd_thresh(
        perm_stats    = perm_i,
        obs_stats     = obs_i,
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
      
      thresh    <- thresh_res$thresh
      n_exceed  <- thresh_res$n_exceed
      out$thresh       <- thresh
      out$n_exceed     <- n_exceed
      out$zero_replaced <- FALSE
      
      # No valid threshold → fallback to empirical
      if (is.na(thresh)) {
        out$method_used   <- "empirical"
        out$p_value       <- p_empirical[i]
        out$shape         <- NA
        out$scale         <- NA
        out$gof_p_value   <- NA
        
      } else {
        # Possibly constrain support
        support_boundary <- switch(control$constraint,
                                   support_at_max = max(trans_obs, na.rm = TRUE),
                                   support_at_obs = obs_i,
                                   NULL)
        
        # Fit GPD
        fit_res <- fit_gpd(
          data             = perm_i[perm_i > thresh],
          thresh           = thresh,
          fit_method       = control$fit_method,
          tol              = control$tol,
          eps_fun          = control$eps_fun,
          eps_par          = control$eps_par,
          constraint       = control$constraint,
          support_boundary = support_boundary,
          gof_test         = control$gof_test,
          gof_alpha        = control$gof_alpha,
          ...
        )
        
        out$shape       <- fit_res$shape
        out$scale       <- fit_res$scale
        out$gof_p_value <- fit_res$p_value
        out$epsilon     <- fit_res$epsilon
        out$method_used <- "gpd"
        
        # Compute upper‐tail probability (p-value)
        p_gpd <- (n_exceed / n_perm) *
          .pgpd_upper_tail(q        = obs_i - thresh,
                           location = 0,
                           scale    = out$scale,
                           shape    = out$shape)
        
        if (p_gpd == 0 | is.na(p_gpd)) {
          out$p_value       <- p_empirical[i]
          out$zero_replaced <- TRUE
          out$method_used   <- "empirical"
        } else {
          out$p_value <- p_gpd
        }
      }
    }

    out
  }
  
  if (control$verbose) {
    # Set progress bar handler
    handlers("progress")
  }
  
  # Run fits in parallel
  results_list <- with_progress({
    if (control$verbose) {
      p <- progressor(along = idx_fit)
    }
    
    future.apply::future_lapply(
      idx_fit,
      FUN = function(i) {
        res <- fit_one(i, control = control, 
                       trans_obs = trans_obs, trans_perm = trans_perm)
        if (control$verbose) {
          p()  # update progress bar
        }
        res
      },
      future.packages = "permApprox",
      future.seed = TRUE
    )
  })
  
  # Shut down the future plan (back to sequential)
  future::plan(future::sequential)
  
  # Assemble outputs
  # prepare vectors
  fitted         <- logical(n_test)
  method_used    <- rep("empirical", n_test)
  thresh_vec     <- numeric(n_test)
  n_exceed_vec   <- integer(n_test)
  shape_vec      <- numeric(n_test)
  scale_vec      <- numeric(n_test)
  gof_pval_vec   <- numeric(n_test)
  zero_replaced  <- logical(n_test)
  epsilon_vec    <- numeric(n_test)
  perm_stats_fit <- vector("list", n_test)

  for (j in seq_along(idx_fit)) {
    i <- idx_fit[j]
    res <- results_list[[j]]

    p_values[i]      <- res$p_value
    fitted[i]        <- TRUE
    method_used[i]   <- res$method_used
    thresh_vec[i]    <- res$thresh
    n_exceed_vec[i]  <- res$n_exceed
    shape_vec[i]     <- res$shape
    scale_vec[i]     <- res$scale
    gof_pval_vec[i]  <- res$gof_p_value
    zero_replaced[i] <- res$zero_replaced
    epsilon_vec[i]   <- res$epsilon

    # store the permuted stats used for fitting
    # (those > threshold, or NULL if no fit)
    perm_stats_fit[[i]] <- if (!is.na(res$thresh)) {
      perm_stats[i, perm_stats[i, ] > res$thresh]
    } else {
      NULL
    }
  }

  list(
    p_values       = p_values,
    fitted         = fitted,
    method_used    = method_used,
    thresh         = thresh_vec,
    n_exceed       = n_exceed_vec,
    shape          = shape_vec,
    scale          = scale_vec,
    epsilon        = epsilon_vec,
    gof_p_value    = gof_pval_vec,
    zero_replaced  = zero_replaced,
    perm_stats_fit = perm_stats_fit
  )
}
