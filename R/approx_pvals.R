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
#' @param appMethod character indicating the method used for p-value estimation.
#'   Possible values are: "gpd" (default), "gamma", "empirical".
#'
#'
#' @param threshMethod method for threshold detection.
#'   Possible values are:
#'   "fix", "ftr", "minPR", "PRbelowAlpha", "fwdStop"
#' @param threshMethodPar fix threshold/number of exceedances or start value for
#'   threshold / number of exceedances, depending on threshMethod
#' @param nExceedMin minimum number of exceedances

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
#' @param gpd_estimate
#' ...
#'
#' @export
#' @import foreach

# tPerm  <- gpddata
# tObs  <- max(tPerm)-0.01
# appMethod = "gpd"
# threshVec = NULL
# threshMethod = "PRbelowAlpha"
# thresh0 = NULL
# nExceed0 = NULL
# nExceedMin = 1
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
# gpd_estimate = NULL


approx_pvals <- function(tPerm,
                        tObs,
                        tol = 1e-8,
                        eps = NULL,
                        includeObs = TRUE,
                        fitMethod = "MLE",
                        finalOnly = FALSE,
                        constraint = "tObs",
                        fitThresh = 0.1,
                        nonZero = TRUE,
                        appMethod = "gpd",
                        threshVec = NULL,
                        threshMethod = "PRbelowAlpha",
                        thresh0 = NULL,
                        nExceed0 = NULL,
                        nExceedMin = 1,
                        stepSize = 1,
                        gofTest = "ad",
                        gofAlpha = 0.05,
                        gofTailRMMeth = "allrej",
                        gofTailRMPar = NULL,
                        seed = NULL,
                        cores = 1L,
                        verbose = FALSE,
                        gpd_estimate = NULL,
                        plotHist = F,
                        histBreaks = 100,
                        histXlim = NULL,
                        histYlim = NULL,
                        plotPvals = TRUE,
                        plotPvalsxvar = "thresh",
                        plotTitle = "",
                        maintext = "") {

  # possible threshold detection methods are:
  # - "minimum proportion of rejected GOF-tests"
  # - "failure-to-reject"
  # threshMethod <- match.arg(threshMethod, choices = c("fix",
  #                                                     "ftr",
  #                                                     "minPR",
  #                                                     "PRbelowAlpha",
  #                                                     "fwdStop"))

  appMethod <- match.arg(appMethod, choices = c("gpd", "gamma", "empirical"))

  fitMethod <- match.arg(fitMethod, choices = c("LME", "MLE1D", "MLE2D", "MOM",
                                                "NLS2", "WNLLSM", "ZSE"))

  constraint <- match.arg(constraint, choices = c("none", "shapePos", "tObs"))
  #
  # gofTailRMMeth <- match.arg(gofTailRMMeth, choices = c("none", "num", "prop",
  #                                                       "allrej", "changepoint"))

  if (length(tObs) == 1) {
    tPerm <- as.matrix(tPerm)
  }

  if (includeObs) tPerm <- rbind(tObs, tPerm)

  # Number of permutations
  nPerm <- nrow(tPerm)

  # Number of tests
  nTest <- ncol(tPerm)

  # Take absolute values (we test H1: tObs > t)
  tPerm <- abs(tPerm)
  tObs <- abs(tObs)


  # Number of permutation test statistics being greater than or equal to tObs
  nlarger <- sapply(1:nTest, function(i) {
    sum(tPerm[, i] >= tObs[i])
  })

  #-----------------------------------------------------------------------------
  # Empirical p-value

  pvals <- nlarger / nPerm

  if (appMethod == "empirical") {

    if (any(pvals == 0)) {
      warning("One or more estimated p-values are zero. Consider including ",
              "the observed test statistic.")
    }

    return(list(pvals = pvals))

  } else if (all(pvals > fitThresh)) {
    # no distribution is fitted

    # out <- list(pvals = pvals)

    # if (appMethod == "gpd") {
    #   out$shape <- out$scale <- out$gofPval <- out$thresh <- out$excessPerm <-
    #     out$excessObs <- out$nExceed <- out$shapeVec <- out$scaleVec <-
    #     out$gofPvalVec <- out$nExceedVec <- NA
    # } else if (appMethod == "gamma") {
    #   out$shape <- out$rate <- out$gofPval <- NA
    # }

    return(list(pvals = pvals))
  }

  #-----------------------------------------------------------------------------
  # Fitting a Gamma distribution

  if (appMethod == "gamma") {

    capture.output(gammafit <- fitdistrplus::fitdist(tPerm,
                                                     distr = "gamma",
                                                     method = "mle"),
                   file='NUL')

    shape <- as.numeric(gammafit$estimate["shape"])
    rate <- as.numeric(gammafit$estimate["rate"])

    params <- list(shape = shape, rate = rate)

    pvals <- pgamma(q = tObs, shape = shape, rate = rate, lower.tail = FALSE)

    cvmtest <- goftest::cvm.test(tPerm, "gamma", shape = shape, rate = rate)
    gofPval <- cvmtest$p.value

    output <- list(pvals = pvals, params = params, gofPval = gofPval)

    return(output)
  }

  #-----------------------------------------------------------------------------
  # Tail approximation using the GPD

  # Indices of p-values above threshold
  indFit <- which(pvals <= fitThresh)

  # Maximum test statistic (needed for constraint)
  tMax <- max(tObs)

  # Add epsilon
  if (is.null(eps)) {
    eps <- tMax / 10
  }

  tMax <- tMax + eps

  if (is.null(gpd_estimate)) {
    #fitMethod.tmp <- ifelse(finalOnly, "ZSE", fitMethod)

    gpdEstList <- list()

    for (i in seq_along(indFit)) {
      gpdEstList[[i]] <- get_gpd_estimates(tPerm = tPerm[, indFit[i]],
                                           tObs = tObs[indFit[i]],
                                           tMax = tMax,
                                           tol = tol,
                                           threshVec = threshVec,
                                           threshMethod = threshMethod,
                                           thresh0 = thresh0,
                                           nExceed0 = nExceed0,
                                           nExceedMin = nExceedMin,
                                           stepSize = stepSize,
                                           includeObs = includeObs,
                                           fitMethod = fitMethod,
                                           constraint = constraint,
                                           fitInformation = fitInformation,
                                           optimMethod = optimMethod,
                                           shape0 = shape0, scale0 = scale0,
                                           shapeMin = shapeMin, scaleMin = scaleMin,
                                           shapeMax = shapeMax, scaleMax = scaleMax,
                                           gofTest = gofTest,
                                           seed = seed,
                                           cores = cores,
                                           verbose = verbose)
      gpdEstList[[i]]$tObs <- tObs[indFit[i]]
    }
  }

  getFinalPval <- function(gpdEst) {

    toassign <- c("idxVec", "nExceedVec", "shapeVec", "scaleVec", "negLogLikVec",
                  "gofPvalVec", "pvalVec", "threshVec", "tObs")

    for (i in seq_along(toassign)) {
      assign(toassign[i], gpdEst[[toassign[i]]])
    }

    if (threshMethod == "fix") {
      idxUse <- 1

      tailRem <- 0
      propTailRem <- 0
      rmTailIdx <- NA

    } else {
      threshIdxList <- getThreshIdx(threshMethod = threshMethod,
                                    idxVec = idxVec,
                                    shapeVec = shapeVec,
                                    gofPvalVec = gofPvalVec,
                                    gofAlpha = gofAlpha,
                                    nExceedMin = nExceedMin,
                                    gofTailRMMeth = gofTailRMMeth,
                                    gofTailRMPar = gofTailRMPar)

      idxUse <- threshIdxList$idxUse

      tailRem <- ifelse(!is.null(threshIdxList$tailRem), threshIdxList$tailRem, 0)
      propTailRem <- ifelse(!is.null(threshIdxList$propTailRem),
                            threshIdxList$propTailRem, 0)
      rmTailIdx <- ifelse(!is.null(threshIdxList$rmTailIdx), threshIdxList$rmTailIdx, NA)
    }

    if (is.null(idxUse)) {

      if (verbose) {
        message(paste0("No iteration led to a good fit (GOF test rejected for ",
                       "all thresholds). Empirical p-value used."))
      }

      gofPval <- shape <- scale <- thresh <- excessPerm <- nExceed <-
        excessObs <- NA

      fitFailed <- TRUE

      pvals <- (nlarger + 1) / (nPerm + 1)

    } else {
      if (finalOnly) {
        # Final fit with selected threshold
        finalfit <- get_gpd_estimates(tPerm = tPerm,
                                      tObs = tObs,
                                      threshVec = threshVec[idxUse],
                                      threshMethod = "fix",
                                      thresh0 = thresh0,
                                      nExceed0 = nExceed0,
                                      nExceedMin = nExceedMin,
                                      stepSize = stepSize,
                                      includeObs = includeObs,
                                      fitMethod = fitMethod,
                                      constraint = constraint,
                                      fitInformation = fitInformation,
                                      optimMethod = optimMethod,
                                      shape0 = shape0, scale0 = scale0,
                                      shapeMin = shapeMin, scaleMin = scaleMin,
                                      shapeMax = shapeMax, scaleMax = scaleMax,
                                      gofTest = gofTest,
                                      seed = seed,
                                      cores = 1,
                                      verbose = verbose)

        gofPval <- finalfit$gofPvalVec
        shape <- finalfit$shapeVec
        scale <- finalfit$scaleVec
        pvals <- finalfit$pvalVec

      } else {
        gofPval <- gofPvalVec[idxUse]
        shape <- shapeVec[idxUse]
        scale <- scaleVec[idxUse]
        pval <- pvalVec[idxUse]
      }

      tObs <- tObs
      thresh <- threshVec[idxUse]
      exceedPerm <- tPerm[tPerm > thresh]
      nExceed <- nExceedVec[idxUse]

      excessPerm <- exceedPerm - thresh
      excessObs <- tObs - thresh

      fitFailed <- FALSE
    }

    if (plotPvals && (!threshMethod == "fix")) {
      plot_gof_pvals(gofPvalVec = gofPvalVec, gofAlpha = gofAlpha,
                     thresh = thresh,
                     threshMethod = threshMethod,
                     propReject = threshIdxList$propReject,
                     xvar = plotPvalsxvar,
                     idxVec = idxVec, threshVec = threshVec,
                     title = "GOF p-values",
                     mar = NULL, cexaxis = NULL,
                     cexlab = NULL, cexmain = NULL, cexlegend = NULL,
                     plotPropReject = TRUE)
    }

    if (plotHist && !fitFailed) {
      plot_excess_hist(excessPerm = excessPerm, excessObs = excessObs,
                       nExceed = nExceed, shape = shape, scale = scale,
                       gofPval = gofPval, pval = pval, breaks = histBreaks,
                       maintext = maintext, main = NULL,
                       mar = NULL, xlim = NULL, ylim = NULL,
                       col = NULL, cexaxis = NULL,
                       cexlab = NULL, cexmain = NULL, cexlegend = NULL)
    }

    if (pval == 0) {
      if (verbose) {
        message("GPD approximation led to a p-value of zero. Empirical p-value used.")
      }

      pval <- (nlarger + 1) / (nPerm + 1)
    }

    out <- list(pval = pval,
                tObs = tObs,
                gofPval = gofPval,
                shape = shape,
                scale = scale,
                thresh = thresh,
                exceedPerm = exceedPerm,
                nExceed = nExceed,
                excessPerm = excessPerm,
                excessObs = excessObs,
                fitFailed = fitFailed)

    return(out)
  }

  finallist <- lapply(gpdEstList, getFinalPval)

  excessObs <- gofPval <- shape <- scale <- thresh <- nExceed <-
    fitFailed <- rep(NA, nTest)

  #exceedPerm <- excessPerm <- list()

  for (i in seq_along(indFit)) {
    pvals[indFit[i]] <- finallist[[i]]$pval
    excessObs[indFit[i]] <- finallist[[i]]$excessObs
    gofPval[indFit[i]] <- finallist[[i]]$gofPval
    shape[indFit[i]] <- finallist[[i]]$shape
    scale[indFit[i]] <- finallist[[i]]$scale
    thresh[indFit[i]] <- finallist[[i]]$thresh
    nExceed[indFit[i]] <- finallist[[i]]$nExceed
    fitFailed[[indFit[i]]] <- finallist[[i]]$fitFailed

    #exceedPerm[[i]] <- finallist[[i]]$exceedPerm
    #excessPerm[[i]] <- finallist[[i]]$excessPerm
  }

  callArgs <- mget(names(formals()),sys.frame(sys.nframe()))
  callArgs$gpd_estimate <- NULL

  output <- list(pvals = pvals,
                 shape = shape,
                 scale = scale,
                 nExceed = nExceed,
                 gofPval = gofPval,
                 #propReject = threshIdxList$propReject,
                 thresh = thresh,
                 #excessPerm = excessPerm,
                 excessObs = excessObs,
                 #tailRem = tailRem,
                 #pvalVec = pvalVec,
                 #shapeVec = shapeVec,
                 #scaleVec = scaleVec,
                 #gofPvalVec = gofPvalVec,
                 #threshVec = threshVec,
                 #nExceedVec = nExceedVec,
                 #idxVec = idxVec,
                 tMax = tMax,
                 eps = eps,
                 callArgs = callArgs)

  return(output)
}


