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

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
#' @title Define number of workers for parallel processing
#' 
#' @keywords internal
.choose_workers <- function(cores, n, parallel_min) {
  if (cores > 1L && n >= parallel_min) {
    min(cores, n)
  } else {
    1L
  }
}

#-------------------------------------------------------------------------------
#' @title Helper for 'default or user'
#' @name percent-or-or-percent
#' 
#' @keywords internal
#' 
`%||%` <- function(x, y) if (is.null(x)) y else x


#-------------------------------------------------------------------------------
#' Tabulate status counts for perm_approx fits
#'
#' @param status A factor or character vector of status labels.
#'
#' @return A named integer vector with counts per status (possibly empty
#'   if \code{status} is \code{NULL} or all \code{NA}).
#'
#' @keywords internal
.perm_approx_status_counts <- function(status) {
  if (is.null(status)) return(integer(0L))
  status <- status[!is.na(status)]
  if (!length(status)) return(integer(0L))
  tab <- table(status)
  as.integer(tab) |> stats::setNames(names(tab))
}

#-------------------------------------------------------------------------------
#' Numeric summary helper for permApprox
#'
#' @param x Numeric vector.
#' @return Named numeric vector with min, median, mean, max (or NULL if no finite values).
#' @keywords internal
.perm_approx_num_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  c(
    min    = min(x),
    median = stats::median(x),
    mean   = mean(x),
    max    = max(x)
  )
}
