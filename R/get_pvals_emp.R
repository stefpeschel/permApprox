#' @title Compute empirical permutation p-values
#' @keywords internal

get_pvals_emp <- function(obs_stats, perm_stats, nTest, nPerm, alternative) {

  # Number of permutation test statistics being more extreme than obs_stats

  # Compute the number of permutation statistics more extreme than obs_stats
  if (alternative == "less") {
    nExtreme <- sapply(1:nTest, function(i) {
      sum(perm_stats[i, ] <= obs_stats[i])
    })

  } else if (alternative == "greater") {
    nExtreme <- sapply(1:nTest, function(i) {
      sum(perm_stats[i, ] >= obs_stats[i])
    })

  } else {
    nExtreme <- sapply(1:nTest, function(i) {
      sum(abs(perm_stats[i, ]) >= abs(obs_stats[i]))
    })
  }

  pvals <- (nExtreme + 1) / (nPerm + 1)

  return(list(pvals = pvals, nExtreme = nExtreme))
}
