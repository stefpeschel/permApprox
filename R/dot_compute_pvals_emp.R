#' @title Compute empirical permutation p-values
#' @keywords internal
.compute_pvals_emp <- function(obs_stats, perm_stats) {
  
  n_test <- length(obs_stats)
  n_perm <- nrow(perm_stats)
  
  n_perm_exceed <- vapply(seq_len(n_test), function(i) {
    sum(perm_stats[, i] >= obs_stats[i])
  }, integer(1))
  
  pvals <- (n_perm_exceed + 1) / (n_perm + 1)
  
  list(pvals = pvals, n_perm_exceed = n_perm_exceed)
}