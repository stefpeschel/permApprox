
get_gpd_thresh <- function(tPerm,
                           tObs,
                           useAllPerm,
                           tMax,
                           tol,
                           eps = NULL,
                           threshPoss = NULL,
                           threshMethod = "PRbelowAlpha",
                           thresh0 = NULL,
                           nExceed0 = NULL,
                           nExceedMin = 1,
                           stepSize = 1,
                           includeObs = FALSE,
                           fitMethod = "MLE",
                           constraint = "none",
                           fitInformation = "observed",
                           optimMethod = NULL,
                           shape0 = NULL, scale0 = NULL,
                           shapeMin = -Inf, scaleMin = -Inf,
                           shapeMax = Inf, scaleMax = Inf,
                           gofTest = "ad",
                           gofAlpha = 0.05,
                           gofTailRMMeth = "allrej",
                           gofTailRMPar = NULL,
                           seed = NULL,
                           cores = 1L,
                           verbose = FALSE,
                           ...) {

  if (!is.null(seed)) set.seed(seed)

  nPerm <- length(tPerm)

  # Sort permutation test statistics in increasing order
  tSort <- sort(tPerm, decreasing = FALSE)
  nPerm_gpd <- nPerm

  #-----------------------------------------------------------------------------

  if (threshMethod == "fix") {
    if (!is.null(thresh0)) {
      thresh <- thresh0

    } else if (!is.null(nExceed0)) {
      threshTmp <- c(0, tSort)
      thresh <- sort(threshTmp, decreasing = TRUE)[nExceed0 + 1]

    } else {
      stop("Either thresh0 or nExceed0 must be set.")
    }

    # if ((tObs - thresh) < 0) {
    #   stop("Threshold must be smaller than tObs.")
    # }

    nExceed <- sum(tPerm > thresh)

    return(list(thresh = thresh, nExceed = nExceed))
  }

  if (is.null(threshPoss)) {
    # Define vector with possible thresholds
    #  (thresholds must be smaller than the observed test statistic and the
    #  largest tPerm to ensure at least 1 exceedance)
    threshPoss <- tSort[tSort < min(tObs, max(tPerm))]
    threshPoss <- c(0, threshPoss)

    #-----------------
    # Adapt threshold vector to start threshold or start number of exceedances
    if (!is.null(thresh0) && !is.null(nExceed0)) {
      stop("Either thresh0 or nExceed0 must be set to NULL")
    }

    if (!is.null(thresh0)) {
      threshPoss <- threshPoss[threshPoss > thresh0]
      threshPoss <- c(thresh0, threshPoss)
    }

    if (!is.null(nExceed0)) {
      tmp <- sort(tSort, decreasing = TRUE)[nExceed0 + 1]
      threshPoss <- threshPoss[threshPoss >= tmp]
    }

    # Adapt threshold vector according to stepSize
    threshPoss <- threshPoss[c(TRUE, rep(FALSE, stepSize - 1))]

    # Adapt threshold vector to minimum number of exceedances
    if (nExceedMin > 1) {
      tmp <- sort(tSort, decreasing = TRUE)[nExceedMin]
      threshPoss <- threshPoss[threshPoss < tmp]
    }

    # Make thresholds unique
    threshPoss <- unique(threshPoss)
  }

  # Number of iterations
  niter <- length(threshPoss)

  #-----------------------------------------------------------------------------
  idxVec <- nExceedVec <- shapeVec <- scaleVec <- gofPvalVec <-
    numeric(length(threshPoss))

  for (i in seq_along(threshPoss)) {
    idxVec[i] <- i

    thresh <- threshPoss[i]

    # exceedPerm are the exceedances (test statistics above the threshold)
    exceedPerm.tmp <- tSort[tSort > thresh]

    # number of exceedances
    nExceedVec[i] <- length(exceedPerm.tmp)

    # Fit and test the GPD distribution
    fittestres <- fit_gpd(data = tSort,
                          thresh = thresh,
                          fitMethod = "ZSE",
                          tol = 1e-8,
                          eps = "fix",
                          epsVal = 0,
                          factor = 1,
                          constraint = "none",
                          maxVal = NULL,
                          gofTest = gofTest,
                          ...)

    shapeVec[i] <- fittestres$shape
    scaleVec[i] <- fittestres$scale
    gofPvalVec[i] <- fittestres$pval
  }

  #-----------------------------------------------------------------------------
  threshIdxList <- get_thresh_idx(threshMethod = threshMethod,
                                idxVec = idxVec,
                                shapeVec = shapeVec,
                                gofPvalVec = gofPvalVec,
                                gofAlpha = gofAlpha,
                                nExceedMin = nExceedMin,
                                gofTailRMMeth = gofTailRMMeth,
                                gofTailRMPar = gofTailRMPar)

  idxUse <- threshIdxList$idxUse

  #tailRem <- ifelse(!is.null(threshIdxList$tailRem), threshIdxList$tailRem, 0)
  #propTailRem <- ifelse(!is.null(threshIdxList$propTailRem),
  #                      threshIdxList$propTailRem, 0)
  #rmTailIdx <- ifelse(!is.null(threshIdxList$rmTailIdx), threshIdxList$rmTailIdx, NA)

  if (is.null(idxUse)) {
    thresh <- nExceed <- NA

  } else {
    thresh <- threshPoss[idxUse]
    nExceed <- nExceedVec[idxUse]
  }

  return(list(thresh = thresh, nExceed = nExceed))
}
