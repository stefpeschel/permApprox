#' @title Fit and test the GPD
#'
#' @param data data vector
#' @param thresh numeric defining the GPD threshold. Data above this threshold
#'   are the exceedances.
#' @param nExceed integer giving the number of exceedances (data above the
#'   threshold).
#' @param fit_method !!! character giving the method used for fitting the GPD
#'   distribution.
#' @param maxVal maximum value for which the probability density function must
#'   be positive. If \code{NULL}, the maximum of the data vector is used.
#' @param ... arguments passed to the function for GPD parameter estimation.
#'
#' @details
#' The maximum value (\code{maxVal}) is only relevant for light-tailed
#' distributions, i.e. if the shape parameter is negative.
#'
#' @export

fit_gpd <- function(data,
                    thresh = NULL,
                    fit_method = "MLE1D",
                    tol = 1e-8,
                    eps = 0.8,
                    eps_type = "quantile",
                    constraint = "none",
                    maxVal = NULL,
                    gof_test = "ad",
                    ...) {

  stopifnot(is.vector(data) & is.numeric(data))
  stopifnot(is.numeric(thresh) & length(thresh) == 1)

  fit_method <- match.arg(fit_method, choices = c("LME", "MLE1D", "MLE2D", "MOM",
                                                "NLS2", "WNLLSM", "ZSE"))

  constraint <- match.arg(constraint, choices = c("none", "shapePos", "obs_stats",
                                                  "obs_statsMax"))

  gof_test <- match.arg(gof_test, choices = c("ad", "cvm"))

  #-----------------------------------------------------------------------------

  # Estimate GPD parameters
  fit <- try(.est_gpd_params(data = data,
                             thresh = thresh,
                             maxVal = maxVal,
                             constraint = constraint,
                             tol = tol,
                             eps = eps,
                             eps_type = eps_type,
                             fit_method = fit_method,
                             ...),
             silent = TRUE)

  if ("try-error" %in% class(fit)) {
    shape <- scale <- negLogLik <- NA
    pval <- 0

  } else {
    shape <- fit$shape
    scale <- fit$scale
    names(shape) <- names(scale) <- NULL

    #negLogLik <- fit$negLogLik

    tSort <- sort(data, decreasing = FALSE)

    # exceedances (test statistics above the threshold)
    exceedPerm <- tSort[tSort > thresh]

    if (gof_test == "ad") {
      #testres <- try(eva::gpdAd(exceedPerm-thresh), silent = TRUE)
      testres <- try(gpdAd_adapt(exceedPerm-thresh,
                                 scale = scale, shape = shape), silent = TRUE)
    } else {
      #testres <- try(eva::gpdCvm(exceedPerm-thresh), silent = TRUE)
      testres <- try(gpdCvm_adapt(exceedPerm-thresh,
                                  scale = scale, shape = shape), silent = TRUE)
    }

    pval <- ifelse("try-error" %in% class(testres), 0, testres$p.value)

    # if ("try-error" %in% class(testres)) {
    #   pval <- 0
    #
    # } else {
    #   pval <- testres$p.value
    # }

  }

  return(list(shape = shape, scale = scale, pval = pval))#, negLogLik = negLogLik))
}

