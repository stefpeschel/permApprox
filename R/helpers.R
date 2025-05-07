#' @title Transform test statistics for tail modeling
#'
#' @description
#' Transforms observed and permutation test statistics according to the specified
#' alternative hypothesis. For one-sided tests, this involves sign-flipping to
#' ensure that values to be fitted (e.g., to a Gamma or GPD) are positive and lie
#' in the appropriate tail of the null distribution.
#'
#' @keywords internal
transform_stats <- function(perm_stats, obs_stats, alternative) {
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
