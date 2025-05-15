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

.pgpd_upper_tail <- function(q, location = 0, shape, scale){
  zedd <- (q - location)/scale
  use.zedd <- pmax(zedd, 0)
  pmax(1 + shape * use.zedd, 0)^(-1/shape)
}
