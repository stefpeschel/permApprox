#' @title Approximate permutation p-values
#'
#' @description
#'  Approximate p-values of one-tailed permutation hypothesis tests, where the
#'  p-value is defined as the probability that a test statistic (in absolute
#'  terms) is greater than or equal to the observed value under the null.
#'  Approximation methods include the empirical p-value, the gamma
#'  approximation, and the GPD tail approximation.
#'
#' @param tObs numeric vector giving the observed test statistic(s).
#' @param tPerm numeric vector or matrix with permutation test statistics. In
#'   the multiple testing setting, \code{tPerm} must be a matrix where each row
#'   contains the permutation test statistics of a test and thus belongs to an
#'   entry in \code{tObs}.
#' @param method character indicating the method used for p-value approximation.
#'   Possible values are: "gpd" (default; fitting the GPD distribution),
#'   "gamma" (fitting a Gamma distribution), "empirical" (use the empirical
#'   p-value).
#' @param fitThresh p-value threshold above which a GPD or gamma distribution is
#'   fitted to approximate the p-value (e.g., 2 times alpha). Defaults to 0.1.
#' @param includeObs logical. Indicates whether the observed test statistic
#'   **tObs** should be included in the fitting process. Note that **tObs** is
#'   always included when computing the empirical p-value.
#' @param fitMethod character giving the method used for fitting the GPD
#'   parameters. Possible methods are: "LME", "MLE1D", "MLE2D", "MOM", "NLS2",
#'   "WNLLSM", and "ZSE".
#' @param constraint character defining the constraint, under which the GPD is
#'   fitted. Possible values are:
#'   \describe{
#'   \item{\code{"none"}}{No constraint.}
#'   \item{\code{"shapePos"}}{GPD shape parameter must be positive.}
#'   \item{\code{"tObs"}}{The upper bound of f(x) (=-scale/shape) must be
#'   greater than **tObs**.}
#'   \item{\code{"tObsMax"}}{The upper bound of f(x) (=-scale/shape) must be
#'   greater than the maximum **tObs** of all tests (in the multiple testing
#'   setting).}
#'   }
#' @param tol convergence tolerance (used for fitting GPD parameters)
#' @param eps small numeric value or proportion (used for constraint)
#' @param epsType character defining the type of epsilon. Default is "quantile".
#' @param threshMethod method for threshold detection.
#'   Possible values are: "fix", "ftr", "minPR", "PRbelowAlpha", "fwdStop"
#' @param thresh0 numeric value giving the initial threshold.
#' @param exceed0 numeric value giving the initial number of exceedances. Either
#'   \code{thresh0} or \code{exceed0} must be given.
#' @param exceedMin numeric giving the minimum number of exceedances.
#' @param stepSize numeric giving the value by which either the threshold or the
#'   number of exceedances is increased in each step, depending on threshMethod.
#' @param gofTest goodness-of-fit test method. Possible values are "ad"
#'   (Anderson-Darling-Test; default) or "cvm" (Cramér–von Mises criterion).
#' @param gofAlpha significance level of the GOF test.
#' @param multAdj character specifying the method used to adjust for multiple
#'   testing. Possible methods are "lfdr" (local FDR), "adaptBH" (default;
#'   adaptive Benjamini Hochberg), "rbFDR" (resampling-based FDR), as well as
#'   the methods contained in \code{p.adjust.methods}.
#' @param alpha significance level.
#' @param trueNullMethod method used for identifying the proportion of
#'   true null hypotheses. Only used if \code{multAdj} is set to "adaptBH".
#' @param pTrueNull proportion of true null hypotheses.
#' @param seed random seed.
#' @param cores number of cores used for parallelization.
#' @param verbose logical. If TRUE, messages returned by internal functions are
#'   printed.
#' @param gpdEstimate not used.
#' @param nseq length of the sequence of possible cut points. Used for
#'   resampling-based multiple testing adjustment.
#' @param pPerm description
#'
#' @export
#' @import foreach

permaprox <- function(tObs,
                      tPerm,
                      alternative = "twoSided",
                      nullValue = 0,
                      method = "gpd",
                      fitThresh = 0.2,
                      gammaOnFail = TRUE,
                      includeObs = TRUE,
                      fitMethod = "MLE1D",
                      constraint = "tObsMax",
                      tol = 1e-8,
                      eps = 0.8,
                      epsType = "quantile",
                      threshMethod = "PRbelowAlpha",
                      thresh0 = NULL,
                      threshPoss = NULL,
                      exceed0 = NULL,
                      exceedMin = 0.1,
                      stepSize = 1,
                      gofTest = "ad",
                      gofAlpha = 0.05,
                      gofTestGamma = TRUE,
                      multAdj = "adaptBH",
                      alpha = 0.05,
                      trueNullMethod = "convest",
                      pTrueNull = NULL,
                      seed = NULL,
                      cores = 1L,
                      verbose = FALSE,
                      gpdEstimate = NULL,
                      nseq = 100,
                      pPerm = NULL,
                      ...) {
  
  # possible threshold detection methods are:
  # - Failure to reject
  # - Failure to reject (min 5 subsequent acceptances)
  # - Proportion of rejections below alpha
  # - Forward Stop
  # - GOF change point
  threshMethod <- match.arg(threshMethod, choices = c("fix",
                                                      "ftr",
                                                      "ftrMin5",
                                                      "PRbelowAlpha",
                                                      "fwdStop",
                                                      "gofCP"))

  alternative <- match.arg(alternative, choices = c("greater", "less", "twoSided"))

  method <- match.arg(method, choices = c("gpd", "gamma", "empirical"))

  fitMethod <- match.arg(fitMethod, choices = c("LME", "MLE1D", "MLE2D", "MOM",
                                                "NLS2", "WNLLSM", "ZSE"))

  constraint <- match.arg(constraint, choices = c("none", "shapePos", "tObs",
                                                  "tObsMax"))

  gofTest <- match.arg(gofTest, choices = c("ad", "cvm"))

  if (constraint == "shapePos" && !fitMethod %in% c("MLE1D", "MLE2D", "NLS2")) {
    stop("Constraint \"shapePos\" only available for methods ",
         "MLE1D, MLE2D, and NLS2.")
  }
  
  if (length(tObs) == 1) {
    tPerm <- matrix(tPerm, ncol = length(tPerm), nrow = 1)
  }
  
  # Number of permutations and tests
  nPerm <- ncol(tPerm)
  nTest <- length(tObs)
  
  # Determine centering vector
  if (is.numeric(nullValue)) {
    center_vec <- nullValue
  } else {
    center_vec <- switch(
      nullValue,
      mean = if (nTest == 1) mean(tPerm) else rowMeans(tPerm),
      median = if (nTest == 1) median(tPerm) else apply(tPerm, 1, median),
      stop("Invalid 'nullValue': must be numeric, 'mean', or 'median'")
    )
  }
  
  # Center permutation statistics
  if (nTest == 1) {
    tPerm <- tPerm - center_vec
  } else {
    tPerm <- sweep(tPerm, 1, center_vec, FUN = "-")
  }
  
  # Center observed statistics
  tObs <- tObs - center_vec
  
  pEmpList <- get_pvals_emp(tObs = tObs, tPerm = tPerm, nTest = nTest,
                            nPerm = nPerm, alternative = alternative)
  
  pEmp <- pEmpList$pvals
  nExtreme <- pEmpList$nExtreme
  
  if (includeObs) tPerm <- cbind(tObs, tPerm)
  #-----------------------------------------------------------------------------
  # Empirical p-value(s) (observed test statistic is always included)

  if (method == "empirical" | all(pEmp > fitThresh)) {

    pvals <- pEmp

    gammaFit <- gpdFit <- NULL

    approxType <- "empirical"

  } else if (method == "gamma") {
    gammaFit <- get_pvals_gamma(pEmp = pEmp,
                                tPerm = tPerm,
                                tObs = tObs,
                                nTest = nTest,
                                fitThresh = fitThresh,
                                alternative = alternative,
                                gofTestGamma,
                                gofAlpha)

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
                             tPerm = tPerm,
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
                 tPerm = tPerm,
                 tObs = tObs
                 #fdr.output = fdr.output,
                 #sp = sp,
                 #fdr.pa = fdr.pa
  )

  return(output)
}


