#' @title Transform test statistics for tail modeling
#'
#' @description
#' Transforms observed and permutation test statistics according to the specified
#' alternative hypothesis. For one-sided tests, this involves sign-flipping to
#' ensure that values to be fitted (e.g., to a Gamma or GPD) are positive and lie
#' in the appropriate tail of the null distribution.
#'
#' @keywords internal
.transform_stats <- function(perm_stats, obs_stats, alternative) {
  if (alternative == "less") {
    perm_stats <- perm_stats[perm_stats < 0]
    perm_stats <- - perm_stats
    obs_stats <- -obs_stats
    
  } else if (alternative == "greater") {
    perm_stats <- perm_stats[perm_stats > 0]
    
  } else {
    perm_stats <- abs(perm_stats)
    obs_stats <- abs(obs_stats)
  }
  return(list(perm_stats = perm_stats, obs_stats = obs_stats))
}



#' @title Quantile function for the upper tail of the GPD
#'
#' @description
#' Enables the computation of very small p-values.
#'
#' @param q quantile
#' @param location location parameter of the GPD
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @keywords internal

.pgpd_upper_tail <- function(q, location = 0, shape, scale) {
  if (!is.numeric(scale) || any(scale <= 0)) stop("scale must be > 0")
  z <- (q - location) / scale
  z <- pmax(z, 0)  # support lower bound: q < location => survival = 1
  
  tiny <- 1e-12
  res <- numeric(length(z))
  
  # Valid region mask (avoid log1p() on invalid inputs)
  if (shape < 0) {
    mask <- (1 + shape * z) > 0  # within support
  } else {
    mask <- rep(TRUE, length(z)) # no finite upper endpoint
  }
  
  if (abs(shape) < tiny) {
    # Exponential limit
    res[mask] <- exp(-z[mask])
  } else {
    # Stable computation via log1p, only on valid indices
    res[mask] <- exp((-1/shape) * log1p(shape * z[mask]))
  }
  
  # Outside support (shape < 0 and beyond endpoint): survival = 0
  res[!mask] <- 0
  
  res
}
