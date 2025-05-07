#' @title Compute p-values via Gamma approximation
#' @keywords internal

get_pvals_gamma <- function(p_empirical, perm_stats, obs_stats, n_test,
                            fit_thresh, alternative, control) {
browser()
  pvals <- p_empirical

  # Assign control arguments
  include_obs <- control$include_obs
  gof_test <- control$gof_test
  gof_alpha <- control$gof_alpha

  # Include observed test statistic
  if (include_obs) perm_stats <- cbind(obs_stats, perm_stats)

  # Number of permutations
  n_perm <- ncol(perm_stats)

  # Initialize vectors for parameter output
  shape <- rate <- gof_p_value <- rep(NA, n_test)

  # Indices of p-values below threshold (only these are fitted)
  idx_fit <- which(pvals <= fit_thresh)

  # Tests for which a Gamma fit is performed
  fitted <- rep(FALSE, n_test)
  fitted[idx_fit] <- TRUE

  # Method that is finally used
  method_used <- rep("empirical", n_test)
  method_used[idx_fit] <- "gamma"

  # List to store the test statistics used for the fit
  perm_stats_used_list <- vector("list", length = n_test)

  # Transform test statistics for tail modeling
  transformed <- lapply(seq_len(n_test), function(i) {
    transform_stats(perm_stats = perm_stats[i, ],
                    obs_stats = obs_stats[i],
                    alternative = alternative)
  })

  # Extract transformed obs_stats and perm_stats into vectors/matrices
  trans_obs <- sapply(transformed, function(x) x$obs_stats)
  trans_perm <- lapply(transformed, function(x) x$perm_stats)


  for (i in idx_fit) {

    obs <- trans_obs[i]
    perm <- trans_perm[[i]]

    # Skip if obs is not in the support (Gamma only supports > 0)
    if (obs <= 0 || length(perm) < control$exceed_min) {
      method_used <- "empirical"
      next
    }

    # Store used permutation stats
    perm_stats_used_list[[i]] <- perm

    # Fit Gamma distribution (warning about NANs is suppressed)
    suppressWarnings(gamma_fit <- fitdistrplus::fitdist(data = perm_stats_used,
                                                        distr = "gamma",
                                                        method = "mle"))

    shape[i] <- as.numeric(gamma_fit$estimate["shape"])
    rate[i] <- as.numeric(gamma_fit$estimate["rate"])

    pvals[i] <- (n_used / n_perm) * pgamma(q = abs(obs_stats[i]),
                                           shape = shape[i],
                                           rate = rate[i],
                                           lower.tail = FALSE)

    if (gof_test == "cvm") {
      # Goodness-of-fit test
      cvmtest <- goftest::cvm.test(x = perm_stats_used,
                                   null = "gamma",
                                   shape = shape[i],
                                   rate = rate[i])

      gof_p_value[i] <- cvmtest$p.value
    }
  }

  if (gof_test == "cvm") {
    gof_rejected <- which(gof_p_value <= gof_alpha)
    pvals[gof_rejected] <- p_empirical[gof_rejected]
    method_used[gof_rejected] <- "empirical"
  }

  output <- list(p_values = pvals,
                 fitted = fitted,
                 shape = shape,
                 rate = rate,
                 gof_p_value = gof_p_value,
                 gof_rejected = gof_rejected,
                 method_used = method_used,
                 perm_stats_used = perm_stats_used_list
  )

  return(output)
}
