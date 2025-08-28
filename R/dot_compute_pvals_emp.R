#' @title Compute empirical permutation p-values
#' @keywords internal

.compute_pvals_emp <- function(obs_stats, perm_stats, n_test, n_perm, alternative) {
  
  if (alternative == "less") {
    n_perm_exceeding <- vapply(seq_len(n_test), function(i) {
      sum(perm_stats[, i] <= obs_stats[i])
    }, integer(1))
    
  } else if (alternative == "greater") {
    n_perm_exceeding <- vapply(seq_len(n_test), function(i) {
      sum(perm_stats[, i] >= obs_stats[i])
    }, integer(1))
    
  } else { # two_sided
    n_perm_exceeding <- vapply(seq_len(n_test), function(i) {
      sum(abs(perm_stats[, i]) >= abs(obs_stats[i]))
    }, integer(1))
  }
  
  pvals <- (n_perm_exceeding + 1) / (n_perm + 1)
  
  return(list(pvals = pvals, n_perm_exceeding = n_perm_exceeding))
}
