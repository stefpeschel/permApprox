#' Multiple testing adjustment
#'
#' The functions adjusts a vector of p-values for multiple testing.
#'
#' @param p_values Numeric. Vector with p-values.
#'
#' @param method Character. Method for p-value adjustment. Options include:
#'   \itemize{
#'     \item \code{"lfdr"}: local false discovery rates via \code{fdrtool}.
#'     \item \code{"adaptBH"}: adaptive Benjamini-Hochberg (requires estimation of the proportion of true nulls).
#'     \item Any method supported by \code{stats::p.adjust} (e.g., "holm", "BH", "BY").
#'   }
#'
#' @param true_null_method Character. Method to estimate the proportion of true null hypotheses
#'   when \code{method = "adaptBH"}. Options (for \code{limma::propTrueNull}):
#'   \itemize{
#'     \item \code{"convest"} (default)
#'     \item \code{"lfdr"}
#'     \item \code{"mean"}
#'     \item \code{"hist"}
#'     \item \code{"farco"}: Farcomeni (2007) iterative plug-in method.
#'   }
#'
#' @param p_true_null Numeric or NULL. Pre-specified proportion of true nulls for \code{adaptBH}.
#'   If \code{NULL}, it is estimated using \code{true_null_method}. Default: \code{NULL}.
#'
#' @param verbose if \code{TRUE}, progress messages are returned.
#'
#' @param seq_length Integer. Number of grid points when computing FDR curves in
#'   \code{"rbFDR"} method. Default: \code{100}.
#'
#' @param perm_stats Permutation test statistics.
#'
#' @param p_perm Permutation p-values.
#'
#' @param cores Integer. Number of CPU cores to use for parallel computations
#'   (e.g., in \code{"rbFDR"}). Must be >=1. Default: \code{1}.
#'
#' @param verbose Logical. If \code{TRUE}, progress messages are printed.
#'   Default: \code{FALSE}.
#'
#' @importFrom fdrtool fdrtool
#' @importFrom stats p.adjust
#' @export

mult_adjust <- function(p_values,
                        method = "adaptBH",
                        true_null_method = "convest",
                        p_true_null = NULL,
                        seq_length = 100,
                        perm_stats = NULL,
                        p_perm = NULL,
                        cores = 1,
                        verbose = FALSE) {

  # Check input arguments

  if (!is.numeric(p_values)) {
    stop('Argument "p_values" must be a numeric vector.')
  }

  method <- match.arg(method, c(p.adjust.methods, "lfdr", "adaptBH", "rbFDR"))

  true_null_method <- match.arg(true_null_method, c("farco", "lfdr", "mean",
                                                "hist", "convest"))

  if (!is.logical(verbose)) {
    stop('Argument "verbose" must be logical.')
  }

  #-----------------------------------------------------------------------------

  if (method == "lfdr") {

    if (verbose) {
      message("")
      message("Execute fdrtool() ...")
    }

    fdrout <- fdrtool::fdrtool(p_values, statistic = "pvalue", plot = FALSE,
                                verbose = verbose)

    p_adjusted <- fdrout$lfdr
    qValues <- fdrout$qval
    names(p_adjusted) <- names(qValues) <- names(p_values)

    out <- list(p_adjusted = p_adjusted, qValues = qValues)

  } else if (method == "adaptBH") {

    m <- length(p_values)
    ind <- m:1
    o <- order(p_values, decreasing = TRUE)
    ro <- order(o)

    if (is.null(p_true_null)) {
      if (true_null_method == "farco") {
        R <- 0
        iter <- TRUE

        while(iter) {
          p_true_null <- 1- (R/m)  # proportion of true null hypotheses
          p_adjusted <- pmin(1, cummin(m * p_true_null / ind * p_values[o]))[ro]
          R_new <- length(which(p_adjusted < 0.05))
          iter <- R_new != R  # stop iteration if R_new==R
          R <- R_new
        }

        if (verbose) {
          message("\n Proportion of true null hypotheses: ", round(p_true_null, 2))
        }

      } else {
        # true_null_method must be one of "lfdr", "mean", "hist", or "convest"
        p_true_null <- limma::propTrueNull(p_values, method = true_null_method)
        if (verbose) {
          message("\n Proportion of true null hypotheses: ", round(p_true_null, 2))
        }
      }
    }

    p_adjusted <- pmin(1, cummin(m * p_true_null / ind * p_values[o]))[ro]

    names(p_adjusted) <- names(p_values)

    out <- list(p_adjusted = p_adjusted, p_true_null = p_true_null)

  } else if (method == "rbFDR") {

    n_perm <- ncol(perm_stats)
    n_test <- nrow(perm_stats)

    #---------------------------------------------------------------------------
    # Compute p-values of the permutation test statistics
    if (is.null(p_perm)) {

      # Initialize parallel stuff

      if (verbose) {
        # Create progress bar:
        pb <- utils::txtProgressBar(0, n_perm, style=3)

        # Function for progress bar
        progress <- function(n) {
          utils::setTxtProgressBar(pb, n)
        }
      }

      if (cores > 1) {
        if (parallel::detectCores() < cores) cores <- parallel::detectCores()

        cl <- parallel::makeCluster(cores, outfile = "")
        doSNOW::registerDoSNOW(cl)
        '%do_or_dopar%' <- get('%dopar%')

      } else {
        '%do_or_dopar%' <- get('%do%')
      }

      if (verbose) {
        opts <- list(progress = progress)
      } else {
        opts <- list()
      }

      loopres <-
        foreach(i = 1:n_perm,
                .export = c(".find_gpd_thresh", "fit_gpd", ".find_gpd_thresh_idx",
                            ".est_gpd_params", "get_p_values_emp",
                            ".gof_gpd_ad", ".gof_gpd_cvm",
                            ".fit_gpd_lme", ".fit_gpd_mle1d", ".MLE1D_fk", ".MLE1D_fp",
                            ".fit_gpd_mle2d", ".MLE2D_negloglik", ".fit_gpd_mom",
                            ".fit_gpd_nls2", ".NLS2_gpdf", ".NLS2_ecdf", ".NLS2_gpdf2",
                            ".NLS2_ecdf2", ".NLS2_rss1", ".NLS2_rss2",
                            ".fit_gpd_wnllsm", ".WNLLSM_sum_i", ".WNLLSM_WLLS1",
                            ".WNLLSM_WLLS", ".WNLSM_WLS1", ".WNLSM_WLS",
                            ".fit_gpd_zse", ".ZSE_lx"),
                .packages = "permAprox",
                .options.snow = opts) %do_or_dopar% {

                  if (verbose) progress(i)

                  p_valuesPermRes <- permaprox(perm_stats = perm_stats[,-i],
                                          obs_stats = perm_stats[, i],
                                          tol = 1e-8,
                                          method = "gpd",
                                          eps = "quantile",
                                          epsVal = 0.95,
                                          constraint = "obs_statsMax",
                                          thresh_method = "fix",
                                          fit_method = "ZSE",
                                          exceed0 = 250,
                                          thresh_step = 1,
                                          includeObs = FALSE,
                                          cores = 1,
                                          multAdj = "none")

                  p_valuesPermRes$p
                }

      if (verbose) {
        # Close progress bar
        close(pb)
      }

      # Stop cluster
      if (cores > 1) parallel::stopCluster(cl)

      p_perm <- matrix(unlist(loopres), nrow = n_test, ncol = n_perm)

    }
    #---------------------------------------------------------------------------
    # True possible cut-points
    ord <- order(p_values, decreasing = TRUE)
    c_poss_true <- p_values[ord]

    # Sequence of possible cut points
    c_poss <- seq(max(p_values) * 1.1, 0, length.out = seq_length)

    # Get proportion of rejected hypotheses given H_0 is true
    get_v <- function(c, p_perm, pEst) {
      B <- ncol(p_perm)
      (sum(p_perm <= c) + sum(pEst <= c)) / (B + 1)
    }

    # Get proportion of rejected hypotheses
    get_r <- function(c, pEst) {
      sum(pEst <= c)
    }

    # V vector of the sequence
    v_hat_seq <- sapply(c_poss, get_v, p_perm = p_perm, pEst = p_values)

    # Linear interpolation to get V for the original p-values
    v_hat <- approx(c_poss, v_hat_seq, xout = c_poss_true)$y

    # R vector for the original p-values
    r_hat <- sapply(c_poss_true, get_r, pEst = p_values)

    # FDR for the original p-values
    fdr <- v_hat / r_hat

    # Reorder
    p_adjusted <- fdr[order(ord)]

    out <- list(p_adjusted = p_adjusted, p_perm = p_perm)

  } else {
    p_adjusted <- stats::p.adjust(p_values, method)

    out <- list(p_adjusted = p_adjusted)
  }

  return(out)
}
