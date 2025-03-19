#' @title Compute empirical permutation p-values
#' @keywords internal

get_pvals_emp <- function(tObs, tPerm, nTest, nPerm, alternative) {

  # Number of permutation test statistics being more extreme than tObs

  # Compute the number of permutation statistics more extreme than tObs
  if (alternative == "less") {
    nExtreme <- sapply(1:nTest, function(i) {
      sum(tPerm[i, ] <= tObs[i])
    })

  } else if (alternative == "greater") {
    nExtreme <- sapply(1:nTest, function(i) {
      sum(tPerm[i, ] >= tObs[i])
    })

  } else {
    nExtreme <- sapply(1:nTest, function(i) {
      sum(abs(tPerm[i, ]) >= abs(tObs[i]))
    })
  }

  pvals <- (nExtreme + 1) / (nPerm + 1)

  return(list(pvals = pvals, nExtreme = nExtreme))
}
