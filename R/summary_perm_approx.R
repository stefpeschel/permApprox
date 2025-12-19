#' Summary method for permApprox results
#'
#' @description
#' Provides a more detailed summary of a \code{"perm_approx"} object than
#' \code{print()}, including fit status counts, summaries of fitted
#' parameters, goodness-of-fit p-values, and empirical vs. approximated
#' p-value summaries.
#'
#' @param object An object of class \code{"perm_approx"}, as returned by
#'   \code{\link{perm_approx}}.
#' @param digits Integer. Number of significant digits to display for
#'   numeric summaries. Default is \code{3}.
#' @param alpha Numeric scalar in \code{(0, 1)}. Significance level used
#'   to count rejections based on the final (possibly adjusted) p-values
#'   in the summary output. Default is \code{0.05}.
#' @param ... Further arguments passed to or from other methods (currently
#'   ignored).
#'
#' @return The input object \code{object}, invisibly.
#'
#' @export
summary.perm_approx <- function(object, digits = 3L, alpha = 0.05, ...) {
  
  if (!inherits(object, "perm_approx")) {
    stop("Object 'object' is not of class 'perm_approx'.")
  }
  
  # Helpers ---------------------------------------------------------------
  fmt <- function(z) formatC(z, digits = digits, format = "fg")
  
  summarize_and_print <- function(label, x) {
    s <- .perm_approx_num_summary(x)
    if (is.null(s)) {
      cat("  ", label, ": no finite values\n", sep = "")
    } else {
      cat("  ", label, ":\n", sep = "")
      cat("    min = ", fmt(s["min"]),
          ", median = ", fmt(s["median"]),
          ", mean = ", fmt(s["mean"]),
          ", max = ", fmt(s["max"]), "\n", sep = "")
    }
  }
  
  # Basic quantities ------------------------------------------------------
  p_final      <- object$p_values
  p_unadj      <- object$p_unadjusted
  p_emp        <- object$p_empirical
  n_test       <- length(p_final)
  fit_method   <- object$fit_method %||% "empirical"
  fit_res      <- object$fit_result
  adjust_res   <- object$adjust_result
  ctrl         <- object$control %||% list()
  approx_thresh   <- ctrl$approx_thresh
  adjust_method <- ctrl$adjust_method %||% if (!is.null(adjust_res)) "unspecified" else "none"
  
  # Header ----------------------------------------------------------------
  cat("Summary of permApprox result\n")
  cat("----------------------------\n")
  cat("Number of tests             : ", n_test, "\n", sep = "")
  
  if (fit_method == "gpd") {
    cat("Approximation method        : GPD tail approximation\n")
  } else if (fit_method == "gamma") {
    cat("Approximation method        : Gamma approximation\n")
  } else {
    cat("Approximation method        : empirical (no parametric tail approximation)\n")
  }
  
  if (!is.null(approx_thresh)) {
    cat("Approximation threshold     : p-values <", fmt(approx_thresh), "\n", sep = "")
  }
  
  if (adjust_method == "none" || is.null(adjust_res)) {
    cat("Multiple testing adjustment : none\n")
  } else {
    cat("Multiple testing adjustment : ", adjust_method, "\n", sep = "")
  }
  
  # Fit status counts -----------------------------------------------------
  cat("\nFit status counts:\n")
  if (!is.null(fit_res) && !is.null(fit_res$status)) {
    status_counts <- .perm_approx_status_counts(fit_res$status)
    if (length(status_counts)) {
      # Order statuses in a logical order if present
      status_order <- c("Successful fits"    = "success", 
                        "GOF rejections"     = "gof_reject", 
                        "Fit failed"         = "fit_failed",
                        "No threshold found" = "no_threshold",
                        "Discrete distributions"   = "discrete",
                        "Not selected for fitting" = "not_selected")
      
      status_order <- status_order[status_order %in% names(status_counts)]
      status_counts <- status_counts[status_order]
      status_names <- names(status_order)
      names(status_counts) <- status_names
      w <- max(nchar(status_names))
      for (nm in status_names) {
        cat("  ", sprintf(paste0("%-", w, "s : %d"), nm, status_counts[[nm]]), "\n", sep = "")
      }
    } else {
      cat("  (no status information available)\n")
    }
  } else {
    cat("  (no parametric fit information available)\n")
  }
  
  # Parameter summaries ---------------------------------------------------
  if (!is.null(fit_res) && !is.null(fit_res$status)) {
    status_vec <- fit_res$status
    idx_success <- which(status_vec == "success")
    
    if (length(idx_success) > 0L) {
      cat("\n")
      if (fit_method == "gpd") {
        cat("GPD parameter summary (successful fits)\n")
        cat("--------------------------------------\n")
        summarize_and_print("shape",    fit_res$shape[idx_success])
        summarize_and_print("scale",    fit_res$scale[idx_success])
        summarize_and_print("n_exceed", fit_res$n_exceed[idx_success])
        
      } else if (fit_method == "gamma") {
        cat("Gamma parameter summary (successful fits)\n")
        cat("----------------------------------------\n")
        summarize_and_print("shape", fit_res$shape[idx_success])
        summarize_and_print("rate",  fit_res$rate[idx_success])
      }
    }
    
    # GOF p-values (all fitted tests with GOF info)
    if (!is.null(fit_res$gof_p_value)) {
      gof_vals <- fit_res$gof_p_value
      gof_vals <- gof_vals[is.finite(gof_vals)]
      if (length(gof_vals)) {
        cat("\nGoodness-of-fit p-values (all fitted tests)\n")
        cat("------------------------------------------\n")
        summarize_and_print("GOF p-values", gof_vals)
      }
    }
  }
  
  # P-value summaries -----------------------------------------------------
  cat("\nP-value summary\n")
  cat("---------------\n")
  
  cat("Empirical p-values:\n")
  summarize_and_print("empirical", p_emp)
  
  cat("\nFinal p-values (unadjusted):\n")
  summarize_and_print("unadjusted", p_unadj)
  
  if (!is.null(adjust_res)) {
    cat("\nFinal p-values (adjusted, ", adjust_method, "):\n", sep = "")
    summarize_and_print("adjusted", p_final)
    
    # Rejections at alpha
    if (is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
        alpha > 0 && alpha < 1 && any(is.finite(p_final))) {
      n_rej <- sum(p_final <= alpha, na.rm = TRUE)
      cat("  Rejections at alpha = ", fmt(alpha), ": ", n_rej, "\n", sep = "")
    }
  }
  
  invisible(object)
}
