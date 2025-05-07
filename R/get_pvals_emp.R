#' @title Compute empirical permutation p-values
#' @keywords internal

get_pvals_emp <- function(obs_stats, perm_stats, n_test, n_perm, alternative) {

  perm_stats <- cbind(obs_stats, perm_stats)

  # Compute the number of permutation statistics at least as extreme as obs_stats
  if (alternative == "less") {
    n_perm_exceeding <- sapply(1:n_test, function(i) {
      sum(perm_stats[i, ] <= obs_stats[i])
    })

  } else if (alternative == "greater") {
    n_perm_exceeding <- sapply(1:n_test, function(i) {
      sum(perm_stats[i, ] >= obs_stats[i])
    })

  } else {
    n_perm_exceeding <- sapply(1:n_test, function(i) {
      sum(abs(perm_stats[i, ]) >= abs(obs_stats[i]))
    })
  }

  pvals <- (n_perm_exceeding + 1) / (n_perm + 1)

  return(list(pvals = pvals, n_perm_exceeding = n_perm_exceeding))
}
