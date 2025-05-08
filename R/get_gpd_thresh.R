
get_gpd_thresh <- function(perm_stats,
                           obs_stats,
                           tol,
                           thresh_method = "PRbelowAlpha",
                           thresh0 = NULL,
                           exceed0 = NULL,
                           exceed_min = 1,
                           thresh_step = 1,
                           gof_test = "ad",
                           gof_alpha = 0.05,
                           seed = NULL,
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

    n_exceed <- sum(perm_stats > thresh)

    return(list(thresh = thresh, n_exceed = n_exceed))
  }


  # Maximum threshold to ensure the minimum number of exceedances
  threshMax <- sort(perm_stats, decreasing = TRUE)[exceed_min]

  # Define vector with possible thresholds
  #  (threshold must be smaller than the observed test statistic and the
  #  maximum threshold defined before)
  thresh_poss <- tSort[tSort < min(obs_stats, threshMax)]
  thresh_poss <- c(0, thresh_poss)

  # Adapt threshold vector to thresh0 or exceed0

  if (!is.null(thresh0)) {
    if (max(thresh_poss) <= thresh0) {
      stop("Argument 'thresh0' is larger then the maximum possible threshold: ",
           round(max(thresh_poss), 3))
    }
    thresh_poss <- thresh_poss[thresh_poss > thresh0]
    thresh_poss <- c(thresh0, thresh_poss)
  }

  if (!is.null(exceed0)) {
    thresh_tmp <- sort(c(perm_stats, 0), decreasing = TRUE)[exceed0 + 1]

    if (max(thresh_poss) <= thresh_tmp) {
      stop("Argument 'exceed0' is too low (minimum number of exceedances is ",
           sum(perm_stats > max(thresh_poss)))
    }

    thresh_poss <- thresh_poss[thresh_poss >= thresh_tmp]
  }

  # Adapt threshold vector according to thresh_step
  thresh_poss <- thresh_poss[c(TRUE, rep(FALSE, thresh_step - 1))]

  # Make thresholds unique
  thresh_poss <- unique(thresh_poss)

  # Number of iterations
  niter <- length(thresh_poss)

  #-----------------------------------------------------------------------------
  idxVec <- n_exceed_vec <- shapeVec <- scaleVec <- gof_p_value_vec <-
    numeric(length(thresh_poss))

  idxUse <- NA

  for (i in seq_along(thresh_poss)) {
    idxVec[i] <- i

    thresh <- thresh_poss[i]

    # exceedPerm are the exceedances (test statistics above the threshold)
    exceedPerm.tmp <- tSort[tSort > thresh]

    # number of exceedances
    n_exceed_vec[i] <- length(exceedPerm.tmp)

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
    thresh <- n_exceed <- NA

  } else {
    thresh <- thresh_poss[idxUse]
    n_exceed <- n_exceed_vec[idxUse]
  }

  if (doPlot) {
    #tmp <- gof_p_value_vec[(idxUse-50):(idxUse+100)]
    #thtmp <- thresh_poss[(idxUse-50):(idxUse+100)]
    #thtmp <- seq(thresh_poss[1], rev(thresh_poss)[1], length = 10)

    plot(gof_p_value_vec ~ thresh_poss, pch = 20,
         ylab = "AD pvalue", xlab = "threshold")
    #abline(v = thtmp, col = "lightgray")
    grid(50, NA, lwd = 1, lty = 1)
    abline(h = gof_alpha)
    abline(v = thresh, col = "red")
    points(gof_p_value_vec ~ thresh_poss, pch = 20)
    legend("topleft",
           legend = c("AD p-values", "AD alpha", "selected threshold"),
           col = c(1, 1, 2), pch = c(20, NA, NA), lty = c(NA, 1, 1))
  }

  return(list(thresh = thresh, n_exceed = n_exceed))
}
