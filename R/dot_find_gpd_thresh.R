#' GPD threshold detection
#'
#' Select a threshold for Generalized Pareto Distribution (GPD) tail fitting
#' based on permutation test statistics and a goodness-of-fit driven rule.
#' The function evaluates a grid of candidate thresholds and returns the
#' selected threshold together with the corresponding number of exceedances.
#'
#' @param perm_stats Numeric vector of permutation test statistics.
#' @param obs_stat Numeric scalar. Observed test statistic. The threshold is
#'   restricted to values strictly smaller than \code{obs_stat} so that the
#'   observed statistic is included in the exceedance region.
#' @param tol Numeric scalar. Tolerance parameter reserved for internal use
#'   (currently not used).
#' @param thresh_method Character string specifying the threshold selection
#'   rule. Supported methods are
#'   \itemize{
#'     \item \code{"fix"}: Use a fixed threshold \code{thresh0} or a fixed
#'       number of exceedances \code{exceed0}.
#'     \item \code{"ftr"}: Failure-to-reject rule. Select the first
#'       threshold at which the GOF null hypothesis is not rejected
#'       (\code{p\_value > gof_alpha}).
#'     \item \code{"rob_ftr"} (default): 
#'       Robust variant of \code{"ftr"} requiring that the 
#'       GOF test is not rejected for at least 3 consecutive thresholds.
#'     \item \code{"fwd_stop"}: ForwardStop rule based on transformed GOF
#'       p-values as proposed by Bader et al. (2018).
#'     \item \code{"gof_cp"}: Changepoint detection on GOF p-values (with 100
#'       small pseudo p-values at the start) to identify a shift from non-GPD to 
#'       GPD behavior.
#'     \item \code{"shape_var"}: Among thresholds where the GOF test is
#'       accepted, choose the one with minimal rolling variance of the shape
#'       estimate over a window of length \code{shape_var_window}.
#'     \item \code{"rpc"}: "Rejection proportion control" rule. 
#'       Among candidate thresholds for
#'       which the GOF test is not rejected, choose the first one such that
#'       the proportion of subsequent rejections is below \code{gof_alpha}.
#'   }
#' @param thresh0 Numeric scalar or \code{NULL}. Initial (upper) threshold
#'   candidate. If supplied, candidate thresholds are restricted to be at least
#'   as large as \code{thresh0}. For \code{thresh_method == "fix"}, this value
#'   is used as the threshold if \code{exceed0} is \code{NULL}.
#' @param exceed0 Integer or numeric scalar, or \code{NULL}. Target number of
#'   exceedances at the initial threshold. If \code{0 < exceed0 <= 1}, it is
#'   interpreted as a proportion of \code{length(perm_stats)}; otherwise it is
#'   taken as an absolute count. For \code{thresh_method == "fix"}, this
#'   determines the threshold via the order statistic if \code{thresh0} is
#'   \code{NULL}. Exactly one of \code{thresh0} and \code{exceed0} must be
#'   non-\code{NULL}.
#' @param exceed_min Integer or numeric scalar. Minimal number of exceedances
#'   allowed at the selected threshold. If \code{0 < exceed_min < 1}, it is
#'   interpreted as a proportion of \code{length(perm_stats)}; otherwise it is
#'   taken as an absolute count.
#' @param thresh_step Integer scalar. Step size for thinning the grid of
#'   candidate thresholds. Only every \code{thresh_step}-th candidate is
#'   evaluated, starting from the first.
#' @param gof_test Character string specifying the goodness-of-fit test to be
#'   used in \code{\link{fit_gpd}}. Typically \code{"ad"} for the
#'   Anderson–Darling test (default).
#' @param gof_alpha Numeric scalar in \eqn{(0, 1)}. Significance level used to
#'   decide whether the GPD GOF test is accepted or rejected.
#' @param shape_var_window Integer scalar. Window length used by the
#'   \code{"shape_var"} threshold selection method to compute the rolling
#'   variance of the estimated GPD shape parameter over accepted thresholds.
#' @param seed Integer scalar or \code{NULL}. Optional random seed to ensure
#'   reproducibility for methods that involve randomness (e.g., changepoint
#'   detection with added pseudo p-values).
#' @param doPlot Logical. If \code{TRUE}, produce a diagnostic plot of GOF
#'   p-values against candidate thresholds, marking the selected threshold.
#' @param ... Additional arguments reserved for internal use (currently
#'   ignored).
#'
#' @details
#' Candidate thresholds are constructed from the sorted permutation statistics
#' and are constrained to be smaller than \code{obs_stat}. For each candidate,
#' a GPD is fitted via \code{\link{fit_gpd}} using the LME method without
#' boundary constraints, and a GOF p-value is obtained. Depending on
#' \code{thresh_method}, these p-values and the corresponding shape estimates
#' are combined to select a single threshold.
#'
#' If no threshold satisfying the method-specific criterion can be found, the
#' function returns \code{NA} for both the threshold and the number of
#' exceedances.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item \code{thresh}: Numeric scalar giving the selected threshold, or
#'     \code{NA_real_} if no suitable threshold was found.
#'   \item \code{n_exceed}: Integer giving the number of permutation statistics
#'     exceeding \code{thresh}, or \code{NA_integer_} if no suitable threshold
#'     was found.
#' }
#'
#' @import stats graphics
#'
#' @keywords internal

.find_gpd_thresh <- function(perm_stats,
                             obs_stat,
                             tol,
                             thresh_method = "rpc",
                             thresh0 = NULL,
                             exceed0 = NULL,
                             exceed0_min = NULL,
                             exceed_min = 1,
                             thresh_step = 1,
                             gof_test = "ad",
                             gof_alpha = 0.05,
                             shape_var_window = 7L,
                             seed = NULL,
                             doPlot = FALSE,
                             ...) {
  if (!is.null(seed)) set.seed(seed)
  
  n_perm <- length(perm_stats)
  
  # Sort permutation test statistics in increasing order													
  t_sort <- sort(perm_stats, decreasing = FALSE)
  
  if (is.null(thresh0) & is.null(exceed0)) {
    stop("Either thresh0 or exceed0 must be set.")
  }
  
  if (!is.null(thresh0) && !is.null(exceed0)) {
    stop("Either thresh0 or exceed0 must be set to NULL")
  }
  
  if (!is.null(exceed0)) {
    exceed0 <- if (is.finite(exceed0) && exceed0 <= 1) {
      max(1L, floor(n_perm * exceed0))
    } else {
      as.integer(exceed0)
    }
    
    exceed0 <- max(exceed0, exceed0_min)
    
    if (exceed0 > n_perm) {
      stop("'exceed0' larger than number of permutations in use (", n_perm, ").")
    }
  }
  
  exceed_min <- if (is.finite(exceed_min) && exceed_min < 1) {
    max(1L, floor(n_perm * exceed_min))
  } else {
    as.integer(exceed_min)
  }
  
  # Initial threshold vector
  thresh_init <- c(0, t_sort)
  
  #-----------------------------------------------------------------------------
  # Fix threshold or number of exceedances
  
  if (thresh_method == "fix") {
    if (!is.null(thresh0)) {
      thresh <- thresh0
    } else {
      thresh <- sort(thresh_init, decreasing = TRUE)[exceed0 + 1]
    }
    
    n_exceed <- sum(perm_stats > thresh)
    
    if (n_exceed < exceed_min | thresh >= obs_stat) {
      return(list(thresh = NA_real_, n_exceed = NA_integer_))
    } else {
      return(list(thresh = thresh, n_exceed = n_exceed))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Translate exceed0 and exceed_min into test statistics
  
  thresh_ex_min <- sort(thresh_init, decreasing = TRUE)[exceed_min + 1]
  
  if (!is.null(exceed0)) {
    thresh0 <- sort(thresh_init, decreasing = TRUE)[exceed0 + 1]
  } else {
    thresh0 <- thresh0
  }
  
  # Largest threshold so that the observed statistic is just included
  thresh_incl_obs <- max(thresh_init[thresh_init < obs_stat])
  
  #-----------------------------------------------------------------------------
  # Define possible thresholds
  
  # Lower and upper bound for thresholds
  thresh_low <- min(thresh0, thresh_incl_obs)
  thresh_upp <- min(thresh_ex_min, thresh_incl_obs)
  
  # Possible thresholds
  thresh_poss <- thresh_init[thresh_init >= thresh_low & thresh_init <= thresh_upp]
  
  # Adapt threshold vector according to thresh_step
  thresh_poss <- thresh_poss[c(TRUE, rep(FALSE, thresh_step - 1))]
  
  # Make thresholds unique
  thresh_poss <- unique(thresh_poss)
  
  # Number of iterations
  niter <- length(thresh_poss)
  
  #-----------------------------------------------------------------------------
  idx_vec <- n_exceed_vec <- shape_vec <- scale_vec <- gof_p_value_vec <-
    numeric(length(thresh_poss))
  
  idx_use <- NA
  
  for (i in seq_along(thresh_poss)) {
    idx_vec[i] <- i
    
    thresh <- thresh_poss[i]
    
    # exceedPerm are the exceedances (test statistics above the threshold)
    exceed_perm_tmp <- t_sort[t_sort > thresh]
    
    # number of exceedances					   
    n_exceed_vec[i] <- length(exceed_perm_tmp)
    
    # Fit and test the GPD distribution
    res <- fit_gpd(
      data = t_sort,
      thresh = thresh,
      fit_method = "LME",
      constraint = "unconstrained",
      gof_test = gof_test
    )
    
    shape_vec[i] <- res$shape
    scale_vec[i] <- res$scale
    gof_p_value_vec[i] <- res$p_value
    
    # Early stopping for FTR variants
    if (thresh_method == "ftr" && is.na(idx_use) && res$p_value > gof_alpha) {
      idx_use <- i
      if (!doPlot) break
      
    } else if (thresh_method == "rob_ftr" && i > min(niter, 2) && is.na(idx_use) &&
               all(gof_p_value_vec[(i-min(niter, 2)):i] > gof_alpha)) {
      idx_use <- i - min(niter, 2)
      if (!doPlot) break
    }
  }
  
  if (is.na(idx_use)) {
    
    # Find threshold index					  
    if (all(gof_p_value_vec <= gof_alpha)) {
      
      idx_use <- NA
      
    } else {
      idx_H0_accept <- which(gof_p_value_vec > gof_alpha)
      idx_H0_reject <- which(gof_p_value_vec <= gof_alpha)
      
      if (thresh_method == "rpc") {
        
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
        
      } else if (thresh_method == "fwd_stop") {
        
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
          
        } else {
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
        
      } else if (thresh_method == "shape_var") {
        
        acc <- idx_H0_accept
        if (length(acc) == 0) {
          idx_use <- NA
        } else if (length(acc) == 1) {
          idx_use <- acc[1]
        } else {
          w <- as.integer(shape_var_window)
          if (!is.finite(w) || w < 2) w <- 7L
          half <- floor((w - 1L) / 2L)
          
          # rolling variance over accepted indices only
          roll_var <- rep(Inf, length(acc))
          for (j in seq_along(acc)) {
            lo <- max(1L, j - half)
            hi <- min(length(acc), j + half)
            idx_win <- acc[lo:hi]
            sh <- shape_vec[idx_win]
            roll_var[j] <- if (length(sh) >= 2) stats::var(sh, na.rm = TRUE) else Inf
          }
          
          # choose minimal variance; tie-breaker: largest threshold (more conservative)
          if (all(!is.finite(roll_var))) {
            idx_use <- acc[which.max(thresh_poss[acc])]
          } else {
            j_best <- which(roll_var == min(roll_var, na.rm = TRUE))
            if (length(j_best) > 1) {
              j_best <- j_best[which.max(thresh_poss[acc[j_best]])]
            }
            idx_use <- acc[j_best[1]]
          }
        }
        
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