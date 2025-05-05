#' Multiple testing adjustment
#'
#' The functions adjusts a vector of p-values for multiple testing
#'
#' @param pvals Numeric. Vector with p-values.
#'
#' @param method Character. Method for p-value adjustment. Options include:
#'   \itemize{
#'     \item \code{"lfdr"}: local false discovery rates via \code{fdrtool}.
#'     \item \code{"adaptBH"}: adaptive Benjamini-Hochberg (requires estimation of the proportion of true nulls).
#'     \item Any method supported by \code{stats::p.adjust} (e.g., "holm", "BH", "BY").
#'   }
#'
#' @param trueNullMethod Character. Method to estimate the proportion of true null hypotheses
#'   when \code{method = "adaptBH"}. Options (for \code{limma::propTrueNull}):
#'   \itemize{
#'     \item \code{"convest"} (default)
#'     \item \code{"lfdr"}
#'     \item \code{"mean"}
#'     \item \code{"hist"}
#'     \item \code{"farco"}: Farcomeni (2007) iterative plug-in method.
#'   }
#'
#' @param pTrueNull Numeric or NULL. Pre-specified proportion of true nulls for \code{adaptBH}.
#'   If \code{NULL}, it is estimated using \code{trueNullMethod}. Default: \code{NULL}.
#'
#' @param verbose if \code{TRUE}, progress messages are returned.
#'
#' @param nseq Integer. Number of grid points when computing FDR curves in \code{"rbFDR"} method.
#'   Default: \code{100}.
#'
#' @param perm_stats ToDo! Permutation test statistics.
#' @param pPerm ToDo! Permutation p-values.
#'
#' @param cores Integer. Number of CPU cores to use for parallel computations (e.g., in \code{"rbFDR"}).
#'   Must be >=1. Default: \code{1}.
#'
#' @param verbose Logical. If \code{TRUE}, progress messages are printed. Default: \code{FALSE}.
#'
#'
#' @importFrom fdrtool fdrtool
#' @importFrom stats p.adjust
#' @export

mult_adjust <- function(pvals,
                        method = "adaptBH",
                        trueNullMethod = "convest",
                        pTrueNull = NULL,
                        nseq = 100,
                        perm_stats = NULL,
                        pPerm = NULL,
                        cores = 1,
                        verbose = FALSE) {

  # Check input arguments

  if (!is.numeric(pvals)) {
    stop('Argument "pvals" must be a numeric vector.')
  }

  method <- match.arg(method, c(p.adjust.methods, "lfdr", "adaptBH", "rbFDR"))

  trueNullMethod <- match.arg(trueNullMethod, c("farco", "lfdr", "mean",
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

    fdrout <- fdrtool::fdrtool(pvals, statistic = "pvalue", plot = FALSE,
                                verbose = verbose)

    pAdjust <- fdrout$lfdr
    qValues <- fdrout$qval
    names(pAdjust) <- names(qValues) <- names(pvals)

    out <- list(pAdjust = pAdjust, qValues = qValues)

  } else if (method == "adaptBH") {

    m <- length(pvals)
    ind <- m:1
    o <- order(pvals, decreasing = TRUE)
    ro <- order(o)

    if (is.null(pTrueNull)) {
      if (trueNullMethod == "farco") {
        R <- 0
        iter <- TRUE

        while(iter) {
          pTrueNull <- 1- (R/m)  # proportion of true null hypotheses
          pAdjust <- pmin(1, cummin(m * pTrueNull / ind * pvals[o]))[ro]
          R_new <- length(which(pAdjust < 0.05))
          iter <- R_new != R  # stop iteration if R_new==R
          R <- R_new
        }

        if (verbose) {
          message("\n Proportion of true null hypotheses: ", round(pTrueNull, 2))
        }

      } else {
        # trueNullMethod must be one of "lfdr", "mean", "hist", or "convest"
        pTrueNull <- limma::propTrueNull(pvals, method = trueNullMethod)
        if (verbose) {
          message("\n Proportion of true null hypotheses: ", round(pTrueNull, 2))
        }
      }
    }

    pAdjust <- pmin(1, cummin(m * pTrueNull / ind * pvals[o]))[ro]

    names(pAdjust) <- names(pvals)

    out <- list(pAdjust = pAdjust, pTrueNull = pTrueNull)

  } else if (method == "rbFDR") {

    nPerm <- ncol(perm_stats)
    nTest <- nrow(perm_stats)

    #---------------------------------------------------------------------------
    # Compute p-values of the permutation test statistics
    if (is.null(pPerm)) {

      # Initialize parallel stuff

      if (verbose) {
        # Create progress bar:
        pb <- utils::txtProgressBar(0, nPerm, style=3)

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
        foreach(i = 1:nPerm,
                .export = c("get_gpd_thresh", "fit_gpd", "get_thresh_idx",
                            ".est_gpd_params", "get_pvals_emp",
                            "gpdAd_adapt", "gpdCvm_adapt",
                            "gpd_LME", "gpd_MLE1D", ".MLE1D_fk", ".MLE1D_fp",
                            "gpd_MLE2D", ".MLE2D_negloglik", "gpd_MOM",
                            "gpd_NLS2", ".NLS2_gpdf", ".NLS2_ecdf", ".NLS2_gpdf2",
                            ".NLS2_ecdf2", ".NLS2_rss1", ".NLS2_rss2",
                            "gpd_WNLLSM", ".WNLLSM_sum_i", ".WNLLSM_WLLS1",
                            ".WNLLSM_WLLS", ".WNLSM_WLS1", ".WNLSM_WLS",
                            "gpd_ZSE", ".ZSE_lx"),
                .packages = "permAprox",
                .options.snow = opts) %do_or_dopar% {

                  if (verbose) progress(i)

                  pvalsPermRes <- permaprox(perm_stats = perm_stats[,-i],
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

                  pvalsPermRes$p
                }

      if (verbose) {
        # Close progress bar
        close(pb)
      }

      # Stop cluster
      if (cores > 1) parallel::stopCluster(cl)

      pPerm <- matrix(unlist(loopres), nrow = nTest, ncol = nPerm)

    }
    #---------------------------------------------------------------------------
    # True possible cut-points
    ord <- order(pvals, decreasing = TRUE)
    c_poss_true <- pvals[ord]

    # Sequence of possible cut points
    c_poss <- seq(max(pvals) * 1.1, 0, length.out = nseq)

    # Get proportion of rejected hypotheses given H_0 is true
    get_v <- function(c, pPerm, pEst) {
      B <- ncol(pPerm)
      (sum(pPerm <= c) + sum(pEst <= c)) / (B + 1)
    }

    # Get proportion of rejected hypotheses
    get_r <- function(c, pEst) {
      sum(pEst <= c)
    }

    # V vector of the sequence
    v_hat_seq <- sapply(c_poss, get_v, pPerm = pPerm, pEst = pvals)

    # Linear interpolation to get V for the original p-values
    v_hat <- approx(c_poss, v_hat_seq, xout = c_poss_true)$y

    # R vector for the original p-values
    r_hat <- sapply(c_poss_true, get_r, pEst = pvals)

    # FDR for the original p-values
    fdr <- v_hat / r_hat

    # Reorder
    pAdjust <- fdr[order(ord)]

    out <- list(pAdjust = pAdjust, pPerm = pPerm)

  } else {
    pAdjust <- stats::p.adjust(pvals, method)

    out <- list(pAdjust = pAdjust)
  }

  return(out)
}
