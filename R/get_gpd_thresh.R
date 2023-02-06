
get_gpd_thresh <- function(tPerm,
                           tObs,
                           useAllPerm,
                           tMax,
                           tol,
                           eps = NULL,
                           threshPoss = NULL,
                           threshMethod = "PRbelowAlpha",
                           thresh0 = NULL,
                           exceed0 = NULL,
                           exceedMin = 1,
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
                           doPlot = TRUE,
                           ...) {

  if (!is.null(seed)) set.seed(seed)

  nPerm <- length(tPerm)

  # Sort permutation test statistics in increasing order
  tSort <- sort(tPerm, decreasing = FALSE)
  nPerm_gpd <- nPerm

  if (!is.null(thresh0) && !is.null(exceed0)) {
    stop("Either thresh0 or exceed0 must be set to NULL")
  }

  if (exceed0 <= 1) {
    exceed0 <- floor(nPerm * exceed0)
  }

  if (exceedMin < 1) {
    exceedMin <- floor(nPerm * exceedMin)
  }

  #-----------------------------------------------------------------------------
  if (threshMethod == "fix") {
    if (!is.null(thresh0)) {
      thresh <- thresh0

    } else if (!is.null(exceed0)) {

      threshTmp <- c(0, tSort)

      thresh <- sort(threshTmp, decreasing = TRUE)[exceed0 + 1]

    } else {
      stop("Either thresh0 or exceed0 must be set.")
    }

    # if ((tObs - thresh) < 0) {
    #   stop("Threshold must be smaller than tObs.")
    # }

    nExceed <- sum(tPerm > thresh)

    return(list(thresh = thresh, nExceed = nExceed))
  }

  if (is.null(threshPoss)) {

    # Maximum threshold to ensure the minimum number of exceedances
    threshMax <- sort(tPerm, decreasing = TRUE)[exceedMin]

    # Define vector with possible thresholds
    #  (threshold must be smaller than the observed test statistic and the
    #  maximum threshold defined before)
    threshPoss <- tSort[tSort < min(tObs, threshMax)]
    threshPoss <- c(0, threshPoss)

    # Adapt threshold vector to thresh0 or exceed0

    if (!is.null(thresh0)) {
      threshPoss <- threshPoss[threshPoss > thresh0]
      threshPoss <- c(thresh0, threshPoss)
    }

    if (!is.null(exceed0)) {
      tmp <- sort(c(tPerm, 0), decreasing = TRUE)[exceed0 + 1]
      threshPoss <- threshPoss[threshPoss >= tmp]
    }

    # Adapt threshold vector according to stepSize
    threshPoss <- threshPoss[c(TRUE, rep(FALSE, stepSize - 1))]

    # Make thresholds unique
    threshPoss <- unique(threshPoss)
  }

  # Number of iterations
  niter <- length(threshPoss)

  #-----------------------------------------------------------------------------
  idxVec <- nExceedVec <- shapeVec <- scaleVec <- gofPvalVec <-
    numeric(length(threshPoss))

  idxUse <- NA

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

    if (threshMethod == "ftr" && is.na(idxUse) && fittestres$pval > gofAlpha) {
      idxUse <- i
      #break
    } else if (threshMethod == "ftrMin5" && i >= 5 && is.na(idxUse) &&
               all(gofPvalVec[(i-5):i] > gofAlpha)) {
      #break
      idxUse <- i-5
    }
  }

  if (is.na(idxUse)) {
    idxUse <- i
  }

  #-----------------------------------------------------------------------------
  threshIdxList <- get_thresh_idx(threshMethod = threshMethod,
                                idxVec = idxVec,
                                shapeVec = shapeVec,
                                gofPvalVec = gofPvalVec,
                                gofAlpha = gofAlpha,
                                exceedMin = exceedMin,
                                gofTailRMMeth = gofTailRMMeth,
                                gofTailRMPar = gofTailRMPar)

  idxUse <- threshIdxList$idxUse

  if (is.null(idxUse) || is.na(idxUse)) {
    idxUse <- length(threshPoss)
  }

  thresh <- threshPoss[idxUse]
  nExceed <- nExceedVec[idxUse]

  if (doPlot) {
    #tmp <- gofPvalVec[(idxUse-50):(idxUse+100)]
    #thtmp <- threshPoss[(idxUse-50):(idxUse+100)]
    #thtmp <- seq(threshPoss[1], rev(threshPoss)[1], length = 10)

    plot(gofPvalVec ~ threshPoss, pch = 20,
         ylab = "AD pvalue", xlab = "threshold")
    #abline(v = thtmp, col = "lightgray")
    grid(50, NA, lwd = 1, lty = 1)
    abline(h = gofAlpha)
    abline(v = thresh, col = "red")
    points(gofPvalVec ~ threshPoss, pch = 20)
    legend("topleft",
           legend = c("AD p-values", "AD alpha", "selected threshold"),
           col = c(1, 1, 2), pch = c(20, NA, NA), lty = c(NA, 1, 1))
  }

  return(list(thresh = thresh, nExceed = nExceed))
}
