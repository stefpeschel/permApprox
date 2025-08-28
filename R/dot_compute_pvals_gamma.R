#' @title Compute p-values via Gamma approximation
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pgamma
#' @import goftest
#'
#' @keywords internal

.compute_pvals_gamma <- function(perm_stats, obs_stats, n_test,
                                 alternative, control) {
  # Assign control arguments
  include_obs <- control$include_obs
  gof_test <- control$gof_test
  gof_alpha <- control$gof_alpha
  
  # Include observed test statistic
  if (include_obs) perm_stats <- rbind(obs_stats, perm_stats)
  
  # Number of permutations
  n_perm <- nrow(perm_stats)
  
  # Initialize vectors for parameter output
  shape <- rate <- gof_p_value <- rep(NA, n_test)
  gof_rejected <- rep(FALSE, n_test)
  
  # Method that is finally used
  method_used <- rep("gamma", n_test)
  
  # List to store the test statistics used for the fit
  perm_stats_fitted <- vector("list", length = n_test)
  pvals <- numeric(n_test)
  
  for (i in seq_along(obs_stats)) {
    
    obs <- obs_stats[i]
    perm <- perm_stats[, i]
    
    # Store used permutation stats
    perm_stats_fitted[[i]] <- perm
    
    # Fit Gamma distribution (warning about NANs is suppressed)
    suppressWarnings(gamma_fit <- fitdistrplus::fitdist(data = perm,
                                                        distr = "gamma",
                                                        method = "mle"))
    
    shape[i] <- as.numeric(gamma_fit$estimate["shape"])
    rate[i] <- as.numeric(gamma_fit$estimate["rate"])
    
    pvals[i] <- (length(perm) / n_perm) * stats::pgamma(q = abs(obs_stats[i]),
                                                        shape = shape[i],
                                                        rate = rate[i],
                                                        lower.tail = FALSE)
    
    if (gof_test == "cvm") {
      # Goodness-of-fit test
      cvmtest <- goftest::cvm.test(x = perm,
                                   null = "gamma",
                                   shape = shape[i],
                                   rate = rate[i])
      
      gof_p_value[i] <- cvmtest$p.value
      
      if (cvmtest$p.value <= gof_alpha) {
        gof_rejected[i] <- TRUE
        method_used[i] <- "empirical"
      }
    }
  }
  
  output <- list(p_values = pvals,
                 shape = shape,
                 rate = rate,
                 gof_p_value = gof_p_value,
                 gof_rejected = gof_rejected,
                 method_used = method_used,
                 perm_stats_fitted = perm_stats_fitted
  )
  
  return(output)
}
