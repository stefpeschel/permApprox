#' @title Compute p-values via Gamma approximation
#' @keywords internal

get_pvals_gamma <- function(pEmp, perm_stats, obs_stats, nTest, fit_thresh, alternative,
                            gof_testGamma, gof_alpha) {

  # Number of permutations
  nPerm <- ncol(perm_stats)

  pvals <- pEmp

  # Initialize vectors for parameter output
  shape <- rate <- gofPval <- rep(NA, nTest)

  # Indices of p-values below threshold (only these are fitted)
  idxFit <- which(pvals <= fit_thresh)

  # Indicates, for which tests a gamma fit is done
  fitted <- rep(FALSE, nTest)
  fitted[idxFit] <- TRUE

  approxType <- rep("empirical", nTest)
  approxType[idxFit] <- "gamma"

  # List to store only the test statistics used for the fit
  perm_statsUsedList <- list()

  for (i in idxFit) {

    if (alternative == "less") {
      perm_statsUsed <- perm_stats[i, ]
      perm_statsUsed <- abs(perm_statsUsed[perm_statsUsed < 0])

    } else if (alternative == "greater") {
      perm_statsUsed <- perm_stats[i, ]
      perm_statsUsed <- perm_statsUsed[perm_statsUsed > 0]

    } else {
      perm_statsUsed <- abs(perm_stats[i, ])
      perm_statsUsed <- perm_statsUsed[perm_statsUsed > 0]
    }

    perm_statsUsedList[[i]] <- perm_statsUsed

    nUsed <- length(perm_statsUsed)

    # Fit Gamma distribution (warning about NANs is suppressed)
    suppressWarnings(gammafit <- fitdistrplus::fitdist(data = perm_statsUsed,
                                                       distr = "gamma",
                                                       method = "mle"))

    shape[i] <- as.numeric(gammafit$estimate["shape"])
    rate[i] <- as.numeric(gammafit$estimate["rate"])

    pvals[i] <- (nUsed / nPerm) * pgamma(q = abs(obs_stats[i]),
                       shape = shape[i],
                       rate = rate[i],
                       lower.tail = FALSE)

    # Goodness-of-fit test
    cvmtest <- gof_test::cvm.test(x = perm_statsUsed,
                                 null = "gamma",
                                 shape = shape[i],
                                 rate = rate[i])

    gofPval[i] <- cvmtest$p.value
  }

  if (gof_testGamma) {
    idxFailed <- which(gofPval <= gof_alpha)
    pvals[idxFailed] <- pEmp[idxFailed]
    approxType[idxFailed] <- "empirical"
  }

  # Remove NULL entries from perm_statsUsedList
  perm_statsUsedList <- perm_statsUsedList[!sapply(perm_statsUsedList, is.null)]

  output <- list(pvals = pvals,
                 fitted = fitted,
                 shape = shape,
                 rate = rate,
                 gofPval = gofPval,
                 fitted = fitted,
                 approxType = approxType,
                 perm_statsUsedList = perm_statsUsedList)

  return(output)
}
