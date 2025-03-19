#' Multiple testing adjustment
#'
#' The functions adjusts a vector of p-values for multiple testing
#'
#' @param pvals numeric vector with p-values
#' @param adjust character specifying the method used for adjustment.
#'   Can be \code{"lfdr"}, \code{"adaptBH"}, or one of the methods provided by
#'   \code{\link[stats]{p.adjust}}.
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}). Accepts the provided options of the
#'   \code{method} argument of \code{\link[limma]{propTrueNull}}:
#'   \code{"convest"}(default), \code{"lfdr"}, \code{"mean"}, and \code{"hist"}.
#'   Can alternatively be \code{"farco"} for
#'   the "iterative plug-in method" proposed by \cite{Farcomeni (2007)}.
#' @param pTrueNull proportion of true null hypothesis used for the adaptBH
#'   method. If \code{NULL}, the proportion is computed using the method
#'   defined via \code{trueNullMethod}.
#' @param verbose if \code{TRUE}, progress messages are returned.
#'
#' @importFrom fdrtool fdrtool
#' @importFrom stats p.adjust
#' @export

mult_adjust <- function(pvals,
                        tPerm = NULL,
                        pPerm = NULL,
                        adjust = "adaptBH",
                        trueNullMethod = "convest",
                        pTrueNull = NULL,
                        nseq = 100,
                        cores = 1,
                        verbose = FALSE) {

  # Check input arguments

  if (!is.numeric(pvals)) {
    stop('Argument "pvals" must be a numeric vector.')
  }

  adjust <- match.arg(adjust, c(p.adjust.methods, "lfdr", "adaptBH", "rbFDR"))

  trueNullMethod <- match.arg(trueNullMethod, c("farco", "lfdr", "mean",
                                                "hist", "convest"))

  if (!is.logical(verbose)) {
    stop('Argument "verbose" must be logical.')
  }

  #-----------------------------------------------------------------------------

  if (adjust == "lfdr") {

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

  } else if (adjust == "adaptBH") {

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

  } else if (adjust == "rbFDR") {

    nPerm <- ncol(tPerm)
    nTest <- nrow(tPerm)

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

                  pvalsPermRes <- permaprox(tPerm = tPerm[,-i],
                                          tObs = tPerm[, i],
                                          tol = 1e-8,
                                          method = "gpd",
                                          eps = "quantile",
                                          epsVal = 0.95,
                                          constraint = "tObsMax",
                                          threshMethod = "fix",
                                          fitMethod = "ZSE",
                                          exceed0 = 250,
                                          stepSize = 1,
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
    pAdjust <- stats::p.adjust(pvals, adjust)

    out <- list(pAdjust = pAdjust)
  }

  return(out)
}
