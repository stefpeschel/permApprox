#' Print method for permApprox results
#'
#' @param x An object of class \code{"perm_approx"}.
#' @param digits Number of digits for numeric summaries.
#' @param ... Ignored.
#'
#' @export
print.perm_approx <- function(x, digits = 3L, ...) {
  
  if (!inherits(x, "perm_approx"))
    stop("Object 'x' is not of class 'perm_approx'.")
  
  # Helpers
  fmt <- function(z) formatC(z, digits = digits, format = "fg")
  
  pvals <- x$p_values
  n_test <- length(pvals)
  fit_method <- x$fit_method %||% "empirical"
  fit_res    <- x$fit_result
  
  ctrl          <- x$control %||% list()
  fit_thresh    <- ctrl$fit_thresh
  adjust_method <- ctrl$adjust_method %||% 
    if (!is.null(x$adjust_result)) "unspecified" else "none"
  
  cat("permApprox result\n")
  cat("-----------------\n")
  cat("Number of tests              : ", n_test, "\n", sep = "")
  
  # Method
  if (fit_method == "gpd") {
    cat("Approximation method         : GPD tail approximation\n")
  } else if (fit_method == "gamma") {
    cat("Approximation method         : Gamma approximation\n")
  } else {
    cat("Approximation method         : empirical (no approximation)\n")
  }
  
  # Threshold + adjustment
  if (!is.null(fit_thresh))
    cat("Approximation threshold      : p-values <", fmt(fit_thresh), "\n", sep = "")
  
  if (adjust_method == "none")
    cat("Multiple testing adjustment  : none\n")
  else
    cat("Multiple testing adjustment  : ", adjust_method, "\n", sep = "")
  
  # Fit summary (compact)
  if (!is.null(fit_res) && !is.null(fit_res$status)) {
    status <- fit_res$status
    
    success_count  <- sum(status == "success", na.rm = TRUE)
    discrete_count <- sum(status == "discrete", na.rm = TRUE)
    gof_count      <- sum(status == "gof_reject", na.rm = TRUE)
    
    cat("\nSuccessful fits              : ", success_count,  "\n", sep = "")
    cat("Discrete distributions       : ", discrete_count, "\n", sep = "")
    cat("GOF rejections               : ", gof_count,      "\n", sep = "")
  }
  
  # Final p-values
  if (any(is.finite(pvals))) {
    cat("\nFinal p-values:\n")
    cat("  min = ", fmt(min(pvals, na.rm = TRUE)), ", ", 
        "median = ", fmt(median(pvals, na.rm = TRUE)), ", ",
        "max = ", fmt(max(pvals, na.rm = TRUE)), "\n", sep = "")
  }
  
  cat("\nUse summary() for detailed fit diagnostics.\n")
  invisible(x)
}
