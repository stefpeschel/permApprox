#' @title Fit and test the GPD
#'
#' @param data data vector
#' @param thresh numeric defining the GPD threshold. Data above this threshold
#'   are the exceedances.
#' @param nExceed integer giving the number of exceedances (data above the
#'   threshold).
#' @param fitMethod !!! character giving the method used for fitting the GPD
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
                    fitMethod = "MLE1D",
                    tol = 1e-8,
                    eps = NULL,
                    factor = 1,
                    constraint = "none",
                    maxVal = NULL,
                    gofTest = "ad",
                    ...) {

  stopifnot(is.vector(data) & is.numeric(data))
  stopifnot(is.numeric(thresh) & length(thresh) == 1)

  fitMethod <- match.arg(fitMethod, choices = c("LME", "MLE1D", "MLE2D", "MOM",
                                                "NLS2", "WNLLSM", "ZSE"))

  constraint <- match.arg(constraint, choices = c("none", "shapePos", "tObs",
                                                  "tObsMax"))

  gofTest <- match.arg(gofTest, choices = c("ad", "cvm"))

  #-----------------------------------------------------------------------------

  # Estimate GPD parameters
  fit <- try(.est_gpd_params(data = data,
                             thresh = thresh,
                             maxVal = maxVal,
                             constraint = constraint,
                             tol = tol,
                             eps = eps,
                             factor = factor,
                             fitMethod = fitMethod,
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

    if (gofTest == "ad") {
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

