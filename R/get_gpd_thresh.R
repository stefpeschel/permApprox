
get_gpd_thresh <- function(perm_stats,
                           obs_stats,
                           tol,
                           threshPoss = NULL,
                           thresh_method = "PRbelowAlpha",
                           thresh0 = NULL,
                           exceed0 = NULL,
                           exceed_min = 1,
                           thresh_step = 1,
                           includeObs = FALSE,
                           fit_method = "MLE",
                           constraint = "none",
                           fitInformation = "observed",
                           optimMethod = NULL,
                           shape0 = NULL, scale0 = NULL,
                           shapeMin = -Inf, scaleMin = -Inf,
                           shapeMax = Inf, scaleMax = Inf,
                           gof_test = "ad",
                           gof_alpha = 0.05,
                           seed = NULL,
                           cores = 1L,
                           verbose = FALSE,
                           doPlot = FALSE,
                           ...) {

  if (!is.null(seed)) set.seed(seed)

  n_perm <- length(perm_stats)

  # Sort permutation test statistics in increasing order
  tSort <- sort(perm_stats, decreasing = FALSE)
  n_perm_gpd <- n_perm

  if (is.null(thresh0) & is.null(exceed0)) {
    exceed0 <- n_perm
  }

  if (!is.null(thresh0) && !is.null(exceed0)) {
    stop("Either thresh0 or exceed0 must be set to NULL")
  }

  if (exceed0 <= 1) {
    exceed0 <- floor(n_perm * exceed0)
  }

  if (exceed_min < 1) {
    exceed_min <- floor(n_perm * exceed_min)
  }

  #-----------------------------------------------------------------------------
  if (thresh_method == "fix") {
    if (!is.null(thresh0)) {
      thresh <- thresh0

    } else if (!is.null(exceed0)) {

      if (exceed0 > n_perm) {
        stop("'exceed0' larger than number of permutations in use (", n_perm, ").")
      }

      threshTmp <- c(0, tSort)

      thresh <- sort(threshTmp, decreasing = TRUE)[exceed0 + 1]

    } else {
      stop("Either thresh0 or exceed0 must be set.")
    }

    # if ((obs_stats - thresh) < 0) {
    #   stop("Threshold must be smaller than obs_stats.")
    # }

    nExceed <- sum(perm_stats > thresh)

    return(list(thresh = thresh, nExceed = nExceed))
  }

  if (is.null(threshPoss)) {

    # Maximum threshold to ensure the minimum number of exceedances
    threshMax <- sort(perm_stats, decreasing = TRUE)[exceed_min]

    # Define vector with possible thresholds
    #  (threshold must be smaller than the observed test statistic and the
    #  maximum threshold defined before)
    threshPoss <- tSort[tSort < min(obs_stats, threshMax)]
    threshPoss <- c(0, threshPoss)

    # Adapt threshold vector to thresh0 or exceed0

    if (!is.null(thresh0)) {
      if (max(threshPoss) <= thresh0) {
        stop("Argument 'thresh0' is larger then the maximum possible threshold: ",
                    round(max(threshPoss), 3))
      }
      threshPoss <- threshPoss[threshPoss > thresh0]
      threshPoss <- c(thresh0, threshPoss)
    }

    if (!is.null(exceed0)) {
      thresh_tmp <- sort(c(perm_stats, 0), decreasing = TRUE)[exceed0 + 1]

      if (max(threshPoss) <= thresh_tmp) {
        stop("Argument 'exceed0' is too low (minimum number of exceedances is ",
        sum(perm_stats > max(threshPoss)))
      }

      threshPoss <- threshPoss[threshPoss >= thresh_tmp]
    }

    # Adapt threshold vector according to thresh_step
    threshPoss <- threshPoss[c(TRUE, rep(FALSE, thresh_step - 1))]

    # Make thresholds unique
    threshPoss <- unique(threshPoss)
  }

  # Number of iterations
  niter <- length(threshPoss)

  #-----------------------------------------------------------------------------
  idxVec <- nExceedVec <- shapeVec <- scaleVec <- gof_p_value_vec <-
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
                          fit_method = "ZSE",
                          tol = 1e-8,
                          eps = 0,
                          eps_type = "fix",
                          factor = 1,
                          constraint = "none",
                          maxVal = NULL,
                          gof_test = gof_test)#,
                          #...)

    shapeVec[i] <- fittestres$shape
    scaleVec[i] <- fittestres$scale
    gof_p_value_vec[i] <- fittestres$pval

    if (thresh_method == "ftr" && is.na(idxUse) && fittestres$pval > gof_alpha) {
      idxUse <- i
      #break
    } else if (thresh_method == "ftrMin5" && i > 5 && is.na(idxUse) &&
               all(gof_p_value_vec[(i-5):i] > gof_alpha)) {
      #break
      idxUse <- i-5
    }
  }

  # if (thresh_method %in% c("ftr", "ftrMin5") && is.na(idxUse)) {
  #   idxUse <- i
  # }

  if (is.na(idxUse)) {
    threshIdxList <- get_thresh_idx(thresh_method = thresh_method,
                                    gof_p_value_vec = gof_p_value_vec,
                                    gof_alpha = gof_alpha)

    idxUse <- threshIdxList$idxUse
  }

  if (is.null(idxUse) || is.na(idxUse)) {
    thresh <- nExceed <- NA

  } else {
    thresh <- threshPoss[idxUse]
    nExceed <- nExceedVec[idxUse]
  }

  if (doPlot) {
    #tmp <- gof_p_value_vec[(idxUse-50):(idxUse+100)]
    #thtmp <- threshPoss[(idxUse-50):(idxUse+100)]
    #thtmp <- seq(threshPoss[1], rev(threshPoss)[1], length = 10)

    plot(gof_p_value_vec ~ threshPoss, pch = 20,
         ylab = "AD pvalue", xlab = "threshold")
    #abline(v = thtmp, col = "lightgray")
    grid(50, NA, lwd = 1, lty = 1)
    abline(h = gof_alpha)
    abline(v = thresh, col = "red")
    points(gof_p_value_vec ~ threshPoss, pch = 20)
    legend("topleft",
           legend = c("AD p-values", "AD alpha", "selected threshold"),
           col = c(1, 1, 2), pch = c(20, NA, NA), lty = c(NA, 1, 1))
  }

  return(list(thresh = thresh, nExceed = nExceed))
}
