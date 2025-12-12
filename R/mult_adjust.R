#' Multiple testing adjustment
#'
#' The functions adjusts a vector of p-values for multiple testing.
#'
#' @param p_values Numeric. Vector with p-values.
#'
#' @param method Character. Method for p-value adjustment. Options include:
#'   \itemize{
#'     \item \code{"lfdr"}: local false discovery rates via \code{fdrtool}.
#'     \item \code{"adapt_BH"}: adaptive Benjamini-Hochberg (requires estimation
#'     of the proportion of true nulls).
#'     \item Any method supported by \code{stats::p.adjust} (e.g., "holm",
#'     "BH", "BY").
#'   }
#'
#' @param true_null_method Character. Method to estimate the proportion of true
#'   null hypotheses when \code{method = "adapt_BH"}. Options:
#'   \itemize{
#'     \item \code{"convest"} (default), \code{"lfdr"}, \code{"mean"}, \code{"hist"}:
#'       estimated via \code{limma::propTrueNull} (requires \pkg{limma}).
#'     \item \code{"farco"}: Farcomeni (2007) iterative plug-in method.
#'   }
#'
#' @param p_true_null Numeric or NULL. Pre-specified proportion of true nulls
#'   for \code{adapt_BH}. If \code{NULL}, it is estimated using
#'   \code{true_null_method}. Default: \code{NULL}.
#'
#' @param verbose Logical. If \code{TRUE}, progress messages are printed.
#'   Default: \code{FALSE}.
#'
#' @importFrom fdrtool fdrtool
#' @importFrom stats p.adjust
#' @export
mult_adjust <- function(p_values,
                        method = "BH",
                        true_null_method = "convest",
                        p_true_null = NULL,
                        verbose = FALSE) {
  
  # Check input arguments ------------------------------------------------------
  
  if (!is.numeric(p_values)) {
    stop('Argument "p_values" must be a numeric vector.')
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop('Argument "verbose" must be a single logical value.')
  }
  
  method <- match.arg(method, c(p.adjust.methods, "lfdr", "adapt_BH"))
  true_null_method <- match.arg(true_null_method,
                                c("farco", "lfdr", "mean", "hist", "convest"))
  
  #-----------------------------------------------------------------------------
  
  if (method == "lfdr") {
    
    if (verbose) {
      message("")
      message("Execute fdrtool() ...")
    }
    
    fdr_out <- fdrtool::fdrtool(
      p_values, statistic = "pvalue", plot = FALSE, verbose = verbose
    )
    
    p_adjusted <- fdr_out$lfdr
    q_values <- fdr_out$qval
    names(p_adjusted) <- names(q_values) <- names(p_values)
    
    out <- list(p_adjusted = p_adjusted, q_values = q_values)
    
  } else if (method == "adapt_BH") {
    
    m <- length(p_values)
    ind <- m:1
    o <- order(p_values, decreasing = TRUE)
    ro <- order(o)
    
    if (is.null(p_true_null)) {
      
      if (true_null_method == "farco") {
        
        R <- 0
        iter <- TRUE
        
        while (iter) {
          p_true_null <- 1 - (R / m) # proportion of true null hypotheses
          p_adjusted <- pmin(1, cummin(m * p_true_null / ind * p_values[o]))[ro]
          R_new <- sum(p_adjusted < 0.05, na.rm = TRUE)
          iter <- (R_new != R)
          R <- R_new
        }
        
      } else {
        
        # limma is optional (Suggests). Only needed for these pi0 estimators.
        if (!requireNamespace("limma", quietly = TRUE)) {
          stop(
            'For method = "adapt_BH" and true_null_method = "', true_null_method,
            '", the Bioconductor package "limma" is required.\n',
            'Install it with:\n',
            '  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n',
            '  BiocManager::install("limma")\n',
            'Alternatively, set true_null_method = "farco" or use method = "BH".'
          )
        }
        
        p_true_null <- limma::propTrueNull(p_values, method = true_null_method)
      }
      
      if (verbose) {
        message("\n Proportion of true null hypotheses: ", round(p_true_null, 2))
      }
    }
    
    p_adjusted <- pmin(1, cummin(m * p_true_null / ind * p_values[o]))[ro]
    names(p_adjusted) <- names(p_values)
    
    out <- list(p_adjusted = p_adjusted, p_true_null = p_true_null)
    
  } else {
    
    p_adjusted <- stats::p.adjust(p_values, method)
    out <- list(p_adjusted = p_adjusted)
  }
  
  out
}
