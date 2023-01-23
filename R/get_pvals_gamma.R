#' @title Compute p-values via Gamma approximation
#' @keywords internal

get_pvals_gamma <- function(pvals, tPerm, tObs, nTest, fitThresh) {

  shape <- rate <- gofPval <- rep(NA, nTest)

  # Indices of p-values below threshold (only these are fitted)
  idxFit <- which(pvals <= fitThresh)

  fitted <- rep(FALSE, nTest)
  fitted[idxFit] <- TRUE

  for (i in idxFit) {

    # Fit Gamma distribution (warning about NANs is suppressed)
    suppressWarnings(gammafit <- fitdistrplus::fitdist(data = tPerm[i, ],
                                                       distr = "gamma",
                                                       method = "mle"))

    shape[i] <- as.numeric(gammafit$estimate["shape"])
    rate[i] <- as.numeric(gammafit$estimate["rate"])

    pvals[i] <- pgamma(q = tObs[i],
                       shape = shape[i],
                       rate = rate[i],
                       lower.tail = FALSE)

    # Goodness-of-fit test
    cvmtest <- goftest::cvm.test(x = tPerm[i, ],
                                 null = "gamma",
                                 shape = shape[i],
                                 rate = rate[i])

    gofPval[i] <- cvmtest$p.value
  }

  output <- list(pvals = pvals,
                 fitted = fitted,
                 shape = shape,
                 rate = rate,
                 gofPval = gofPval,
                 fitted = fitted)

  return(output)
}
