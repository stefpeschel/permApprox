#' @title Approximate permutation p-values
#'
#' @description
#'   Approximate the p-value for a one-sided permutation hypothesis test,
#'   where the p-value is defined as the probability of observing a test
#'   statistic (in absolute terms) greater than or equal to the observed value
#'   under the null. Approximation methods include the empirical p-value,
#'   Gamma approximation, and GPD tail approximation.
#'
#' @param tObs numeric value giving the observed test statistic
#' @param tPerm numeric vector with permutation test statistics
#' @param includeObs logical. Indicates whether the observed test statistic
#'   should be included in the permutation distribution.
#' @param fitThresh threshold for initial p-value, above which a GPD or Gamma
#'   distribution is fitted for p-value approximation (e.g. 2 times alpha).
#'   Defaults is 0.1.
#' @param nonZero logical. If \code{TRUE} (default), a p-value of zero is
#'   strictly avoided.
#' @param method character indicating the method used for p-value estimation.
#'   Possible values are: "gpd" (default), "gamma", "empirical".
#'
#'
#' @param threshMethod method for threshold detection.
#'   Possible values are:
#'   "fix", "ftr", "minPR", "PRbelowAlpha", "fwdStop"
#' @param threshMethodPar fix threshold/number of exceedances or start value for
#'   threshold / number of exceedances, depending on threshMethod
#' @param exceedMin minimum number of exceedances

#' @param stepSize numeric giving the value by which either the threshold or the
#'   number of exceedances is increased in each step, depending on threshMethod
#' @param ffwd if TRUE, large steps are made at the beginning until H0 is accepted
#'   for the first time

#' @param fitMethod ...
#' @param fitInformation ...
#' @param optimMethod ...
#' @param shape0,scale0 start value of the shape and scale parameter
#' @param shapeMin,scaleMin minimum value of the shape and scale parameter
#' @param shapeMax,scaleMax maximum value of the shape and scale parameter
#' @param gofTest goodness-of-fit test method
#' @param gofAlpha significance level of the GOF test
#' @param seed ...
#' @param cores
#' @param plotHist
#' @param histBreaks
#' @param histXlim
#' @param gpdEstimate
#' ...
#'
#' @export
#' @import foreach

# tPerm  <- gpddata
# tObs  <- max(tPerm)-0.01
# method = "gpd"
# threshVec = NULL
# threshMethod = "PRbelowAlpha"
# thresh0 = NULL
# exceed0 = NULL
# exceedMin = 1
# fitThresh = 0.1
# stepSize = 1
# includeObs = FALSE
# fitMethod = "ml"
# fitInformation = "observed"
# optimMethod = NULL
# shape0 = NULL
# scale0 = NULL
# shapeMin = -Inf
# scaleMin = -Inf
# shapeMax = Inf
# scaleMax = Inf
# gofTest = "ad"
# gofAlpha = 0.05
# gofTailRMMeth = "allrej"
# gofTailRMPar = NULL
# seed = NULL
# cores = 4L
# plotHist = TRUE
# histBreaks = 100
# histXlim = NULL
# histYlim = NULL
# plotPvals = TRUE
# plotPvalsxvar = "thresh"
# plotTitle = ""
# maintext = ""
# verbose = TRUE
# gpdEstimate = NULL

permaprox <- function(tPerm,
                    tObs,
                    useAllPerm = FALSE,
                    method = "gpd",
                    fitThresh = 0.2,
                    includeObs = TRUE,
                    fitMethod = "MLE1D",
                    constraint = "tObs",
                    tol = 1e-8,
                    eps = "quantile",
                    epsVal = 0.9,
                    threshMethod = "PRbelowAlpha",
                    thresh0 = NULL,
                    exceed0 = NULL,
                    exceedMin = 0.1,
                    stepSize = 1,
                    gofTest = "ad",
                    gofAlpha = 0.05,
                    gofTailRMMeth = "allrej",
                    gofTailRMPar = NULL,
                    seed = NULL,
                    cores = 1L,
                    verbose = FALSE,
                    gpdEstimate = NULL,
                    plotHist = F,
                    histBreaks = 100,
                    histXlim = NULL,
                    histYlim = NULL,
                    plotPvals = TRUE,
                    plotPvalsxvar = "thresh",
                    plotTitle = "",
                    maintext = "",
                    jpoint = 0.2,
                    alpha = 0.05,
                    multAdj = "adaptBH",
                    trueNullMethod = "convest",
                    pTrueNull = NULL,
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

  method <- match.arg(method, choices = c("gpd", "gamma", "empirical"))

  fitMethod <- match.arg(fitMethod, choices = c("LME", "MLE1D", "MLE2D", "MOM",
                                                "NLS2", "WNLLSM", "ZSE"))

  constraint <- match.arg(constraint, choices = c("none", "shapePos", "tObs",
                                                  "tObsMax"))

  if (constraint == "shapePos" && !fitMethod %in% c("MLE1D", "MLE2D", "NLS2")) {
    stop("Constraint \"shapePos\" only available for methods ",
         "MLE1D, MLE2D, and NLS2.")
  }

  if (length(tObs) == 1) {
    tPerm <- matrix(tPerm, ncol = length(tPerm), nrow = 1)
  }

  # Number of permutations
  nPerm <- ncol(tPerm)

  #  Overall number of permutation test statistics
  ntPerm <- length(tPerm)

  # Number of tests
  nTest <- length(tObs)

  # Take absolute values (we test H1: tObs > t)
  tPerm <- abs(tPerm)
  tObs <- abs(tObs)

  pEmpList <- get_pvals_emp(tObs = tObs, tPerm = tPerm, nTest = nTest, nPerm = nPerm,
                            ntPerm = ntPerm, useAllPerm = useAllPerm)
  pEmp <- pEmpList$pvals

  if (includeObs) tPerm <- cbind(tObs, tPerm)
  #-----------------------------------------------------------------------------
  # Empirical p-value(s) (observed test statistic is always included)

  if (method == "empirical" | all(pEmp > fitThresh)) {

    pvals <- pEmp

    gammaFit <- gpdFit <- NULL

  } else if (method == "gamma") {

    gammaFit <- get_pvals_gamma(pvals = pEmp,
                                tPerm = tPerm,
                                tObs = tObs,
                                nTest = nTest,
                                fitThresh = fitThresh)

    pvals <- gammaFit$pvals

    gammaFit$pvals <- NULL

    gpdFit <- NULL

  } else if (method == "gpd") {

    # Tail approximation using the GPD

    args <- as.list(environment())

    gpdFit <- do.call(get_pvals_gpd, args)

    pvals <- gpdFit$pvals

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
                 adjustRes = adjustRes
                 #fdr.output = fdr.output,
                 #sp = sp,
                 #fdr.pa = fdr.pa
  )

  return(output)
}


