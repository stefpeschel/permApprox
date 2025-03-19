#' @title Compute p-values via Gamma approximation
#' @keywords internal

get_pvals_gamma <- function(pEmp, tPerm, tObs, nTest, fitThresh, alternative,
                            gofTestGamma, gofAlpha) {

  # Number of permutations
  nPerm <- ncol(tPerm)

  pvals <- pEmp

  # Initialize vectors for parameter output
  shape <- rate <- gofPval <- rep(NA, nTest)

  # Indices of p-values below threshold (only these are fitted)
  idxFit <- which(pvals <= fitThresh)

  # Indicates, for which tests a gamma fit is done
  fitted <- rep(FALSE, nTest)
  fitted[idxFit] <- TRUE

  approxType <- rep("empirical", nTest)
  approxType[idxFit] <- "gamma"

  # List to store only the test statistics used for the fit
  tPermUsedList <- list()

  for (i in idxFit) {

    if (alternative == "less") {
      tPermUsed <- tPerm[i, ]
      tPermUsed <- abs(tPermUsed[tPermUsed < 0])

    } else if (alternative == "greater") {
      tPermUsed <- tPerm[i, ]
      tPermUsed <- tPermUsed[tPermUsed > 0]

    } else {
      tPermUsed <- abs(tPerm[i, ])
      tPermUsed <- tPermUsed[tPermUsed > 0]
    }

    tPermUsedList[[i]] <- tPermUsed

    nUsed <- length(tPermUsed)

    # Fit Gamma distribution (warning about NANs is suppressed)
    suppressWarnings(gammafit <- fitdistrplus::fitdist(data = tPermUsed,
                                                       distr = "gamma",
                                                       method = "mle"))

    shape[i] <- as.numeric(gammafit$estimate["shape"])
    rate[i] <- as.numeric(gammafit$estimate["rate"])

    pvals[i] <- (nUsed / nPerm) * pgamma(q = abs(tObs[i]),
                       shape = shape[i],
                       rate = rate[i],
                       lower.tail = FALSE)

    # Goodness-of-fit test
    cvmtest <- goftest::cvm.test(x = tPermUsed,
                                 null = "gamma",
                                 shape = shape[i],
                                 rate = rate[i])

    gofPval[i] <- cvmtest$p.value
  }

  if (gofTestGamma) {
    idxFailed <- which(gofPval <= gofAlpha)
    pvals[idxFailed] <- pEmp[idxFailed]
    approxType[idxFailed] <- "empirical"
  }

  # Remove NULL entries from tPermUsedList
  tPermUsedList <- tPermUsedList[!sapply(tPermUsedList, is.null)]

  output <- list(pvals = pvals,
                 fitted = fitted,
                 shape = shape,
                 rate = rate,
                 gofPval = gofPval,
                 fitted = fitted,
                 approxType = approxType,
                 tPermUsedList = tPermUsedList)

  return(output)
}
