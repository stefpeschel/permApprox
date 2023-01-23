#' @title Compute empirical permutation p-values
#' @keywords internal

get_pvals_emp <- function(tObs, tPerm, nTest, nPerm, ntPerm, useAllPerm) {

  # Number of permutation test statistics being greater than or equal to tObs
  if (useAllPerm) {
    nlarger <- sapply(1:nTest, function(i) {
      sum(tPerm >= tObs[i])
    })

    nPerm <- ntPerm

  } else {
    nlarger <- sapply(1:nTest, function(i) {
      sum(tPerm[i, ] >= tObs[i])
    })
  }

  pvals <- (nlarger + 1) / (nPerm + 1)

  return(list(pvals = pvals, nlarger = nlarger))
}
