#' @title Approximate permutation p-values
#'
#' @description
#'  Approximate p-values of one-tailed permutation hypothesis tests, where the
#'  p-value is defined as the probability that a test statistic (in absolute
#'  terms) is greater than or equal to the observed value under the null.
#'  Approximation methods include the empirical p-value, the gamma
#'  approximation, and the GPD tail approximation.
#'
#' @param obs_stats Numeric vector of observed test statistic(s).
#'
#' @param perm_stats Numeric vector or matrix of permutation test statistics.
#'   In the multiple testing setting, \code{perm_stats} must be a matrix where
#'   each row contains the permutation test statistics of one test and thus
#'   belongs to an entry in \code{obs_stats}.
#'
#' @param method character indicating the method used for p-value approximation.
#'   Possible values are: "gpd" (default; fitting the GPD distribution),
#'   "gamma" (fitting a Gamma distribution), "empirical" (use the empirical
#'   p-value).
#'
#' @param fit_thresh p-value threshold above which a GPD or gamma distribution is
#'   fitted to approximate the p-value (e.g., 2 times alpha). Defaults to 0.1.
#'
#' @param control_gpd description
#'
#' @param control_gamma description
#'
#' @param control_mult_adjust description
#'
#' @export
#'
#' @import foreach

permaprox <- function(obs_stats,
                      perm_stats,
                      alternative = "twoSided",
                      nullValue = 0,
                      method = "gpd",
                      fit_thresh = 0.2,
                      control_gpd = NULL,
                      control_gamma = NULL,
                      control_mult_adjust = NULL,
                      ...) {

  alternative <- match.arg(alternative, choices = c("greater", "less", "twoSided"))

  method <- match.arg(method, choices = c("gpd", "gamma", "empirical"))


  if (length(obs_stats) == 1) {
    perm_stats <- matrix(perm_stats, ncol = length(perm_stats), nrow = 1)
  }

  # Number of permutations and tests
  nPerm <- ncol(perm_stats)
  nTest <- length(obs_stats)

  # Determine centering vector
  if (is.numeric(nullValue)) {
    center_vec <- nullValue
  } else {
    center_vec <- switch(
      nullValue,
      mean = if (nTest == 1) mean(perm_stats) else rowMeans(perm_stats),
      median = if (nTest == 1) median(perm_stats) else apply(perm_stats, 1, median),
      stop("Invalid 'nullValue': must be numeric, 'mean', or 'median'")
    )
  }

  # Center permutation statistics
  if (nTest == 1) {
    perm_stats <- perm_stats - center_vec
  } else {
    perm_stats <- sweep(perm_stats, 1, center_vec, FUN = "-")
  }

  # Center observed statistics
  obs_stats <- obs_stats - center_vec

  pEmpList <- get_pvals_emp(obs_stats = obs_stats, perm_stats = perm_stats, nTest = nTest,
                            nPerm = nPerm, alternative = alternative)

  pEmp <- pEmpList$pvals
  nExtreme <- pEmpList$nExtreme

  if (includeObs) perm_stats <- cbind(obs_stats, perm_stats)
  #-----------------------------------------------------------------------------
  # Empirical p-value(s) (observed test statistic is always included)

  if (method == "empirical" | all(pEmp > fit_thresh)) {

    pvals <- pEmp

    gammaFit <- gpdFit <- NULL

    approxType <- "empirical"

  } else if (method == "gamma") {
    gammaFit <- get_pvals_gamma(pEmp = pEmp,
                                perm_stats = perm_stats,
                                obs_stats = obs_stats,
                                nTest = nTest,
                                fit_thresh = fit_thresh,
                                alternative = alternative,
                                gof_testGamma,
                                gof_alpha)

    pvals <- gammaFit$pvals
    approxType <- gammaFit$approxType
    gammaFit$pvals <- NULL
    gpdFit <- NULL

  } else if (method == "gpd") {

    # Tail approximation using the GPD

    args <- as.list(environment())

    gpdFit <- do.call(get_pvals_gpd, args)

    pvals <- gpdFit$pvals
    approxType <- gpdFit$approxType
    gammaFit <- NULL
  }

  #-----------------------------------------------------------------------------
  # Multiple testing adjustment

  if (multAdj == "none") {
    pAdjust <- NULL
    adjustRes <- NULL

  } else {
    adjustRes <- mult_adjust(pvals,
                             perm_stats = perm_stats,
                             pPerm = pPerm,
                             adjust = multAdj,
                             trueNullMethod = trueNullMethod,
                             pTrueNull = pTrueNull,
                             nseq = nseq,
                             verbose = verbose)
    pAdjust <- adjustRes$pAdjust
  }

  #-----------------------------------------------------------------------------
  callArgs <- mget(names(formals()),sys.frame(sys.nframe()))
  callArgs$gpdEstimate <- NULL

  output <- list(p = pvals,
                 pAdjust = pAdjust,
                 pEmp = pEmp,
                 gammaFit = gammaFit,
                 gpdFit = gpdFit,
                 approxType = approxType,
                 adjustRes = adjustRes,
                 perm_stats = perm_stats,
                 obs_stats = obs_stats
                 #fdr.output = fdr.output,
                 #sp = sp,
                 #fdr.pa = fdr.pa
  )

  return(output)
}


