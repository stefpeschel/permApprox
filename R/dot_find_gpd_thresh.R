#' @title Adapted version of gpdAd from eva package
#'
#' @import stats graphics
#'
#' @keywords internal
#'
.find_gpd_thresh <- function(perm_stats,
                             obs_stats,
                             tol,
                             thresh_method = "pr_below_alpha",
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
  
  idx_use <- NA
  
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
                          fit_method = "LME",
                          constraint = "unconstrained",
                          gof_test = gof_test)#,
    #...)
    
    shapeVec[i] <- fittestres$shape
    scaleVec[i] <- fittestres$scale
    gof_p_value_vec[i] <- fittestres$p_value
    
    if (thresh_method == "ftr" && is.na(idx_use) && fittestres$p_value > gof_alpha) {
      idx_use <- i
      if (!doPlot) {
        break
      }
    } else if (thresh_method == "ftr_min5" && i > 5 && is.na(idx_use) &&
               all(gof_p_value_vec[(i-5):i] > gof_alpha)) {
      idx_use <- i-5
      if (!doPlot) {
        break
      }
    }
  }
  
  # if (thresh_method %in% c("ftr", "ftr_min5") && is.na(idx_use)) {
  #   idx_use <- i
  # }
  
  if (is.na(idx_use)) {
    
    # Find threshold index
    if(all(gof_p_value_vec <= gof_alpha)){
      
      idx_use <- NA
      
    } else {
      idx_H0_accept <- which(gof_p_value_vec > gof_alpha)
      idx_H0_reject <- which(gof_p_value_vec <= gof_alpha)
      
      if(thresh_method == "pr_below_alpha"){
        
        # Actual number of iterations
        n <- length(gof_p_value_vec)
        
        # Proportion of rejected GOF tests for all thresholds
        prop_reject <- sapply(1:n, function(i){
          sum(gof_p_value_vec[i:n] <= gof_alpha) / (n - i + 1)
        })
        
        # Ensure that H0 is accepted at the chosen threshold
        prop_reject2 <- prop_reject
        prop_reject2[idx_H0_reject] <- 1
        
        # Select the first threshold with a PR below alpha
        idx_use <- which(prop_reject2 <= gof_alpha)[1]
        
        if (is.na(idx_use)) {
          # Use the minimum if no threshold leads to a PR below alpha
          idx_use <- which.min(prop_reject2)
        }
        
      } else if(thresh_method == "fwd_stop"){
        
        # Log-transformed p-values
        y <- -log(1 - gof_p_value_vec)
        
        # Transform log-p-values as described in Barder et. al 2018
        ysum <- sapply(1:length(gof_p_value_vec), function(k) {
          1 / k * sum(y[1:k])
        })
        
        if (all(ysum <= gof_alpha)) {
          idx_use <- length(gof_p_value_vec)
          
        } else if (all(ysum > gof_alpha)) {
          idx_use <- 1
          
        } else{
          # Select the first value above alpha
          idx_use <- which(ysum > gof_alpha)[1]
        }
        
      } else if (thresh_method == "gof_cp") {
        
        # Add 100 fake p-values (sampled from U(0, 0.01)) to ensure a correct
        # estimate if (nearly) all hypotheses are true
        gof_p_value_tmp <- c(stats::runif(100, min = 0, max = 0.01),
                             gof_p_value_vec)
        
        # Changepoint detection
        cp <- changepoint::cpt.meanvar(gof_p_value_tmp)
        idxSel <- cp@cpts[1] - 100
        
        # Find the next larger index for which the AD test is accepted
        idx_use <- idx_H0_accept[idx_H0_accept >= idxSel][1]
        
      } else {
        stop("Threshold method not supported.")
      }
    }
  }
  
  if (is.null(idx_use) || is.na(idx_use)) {
    thresh <- n_exceed <- NA
    
  } else {
    thresh <- thresh_poss[idx_use]
    n_exceed <- n_exceed_vec[idx_use]
  }
  
  if (doPlot) {
    #tmp <- gof_p_value_vec[(idx_use-50):(idx_use+100)]
    #thtmp <- thresh_poss[(idx_use-50):(idx_use+100)]
    #thtmp <- seq(thresh_poss[1], rev(thresh_poss)[1], length = 10)
    
    plot(gof_p_value_vec ~ thresh_poss, pch = 20,
         ylab = "AD pvalue", xlab = "threshold")
    #abline(v = thtmp, col = "lightgray")
    graphics::grid(50, NA, lwd = 1, lty = 1)
    graphics::abline(h = gof_alpha)
    graphics::abline(v = thresh, col = "red")
    graphics::points(gof_p_value_vec ~ thresh_poss, pch = 20)
    graphics::legend("topleft",
                     legend = c("AD p-values", "AD alpha", "selected threshold"),
                     col = c(1, 1, 2), pch = c(20, NA, NA), lty = c(NA, 1, 1))
  }
  
  n_exceed <- as.integer(n_exceed)
  
  return(list(thresh = thresh, n_exceed = n_exceed))
}
