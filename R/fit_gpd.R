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

  fit_method <- match.arg(fit_method,
                          choices = c("LME",
                                      "MLE1D",
                                      "MLE2D",
                                      "MOM",
                                      "NLS2",
                                      "WNLLSM",
                                      "ZSE"))

  constraint <- match.arg(constraint,
                          choices = c("unconstrained",
                                      "shape_nonneg",
                                      "support_at_obs",
                                      "support_at_max"))

  stopifnot(is.numeric(eps))
  stopifnot(eps_type %in% c("quantile", "fix"))

  gof_test <- match.arg(gof_test, choices = c("ad", "cvm", "none"))

  #-----------------------------------------------------------------------------
  # Fit GPD to the data
  #-----------------------------------------------------------------------------

  # Exceedances (test statistics above threshold)
  exceedances <- data[data > thresh]
  excess <- exceedances - thresh

  # Consider maximum value at which the GPD density must be positive
  if (!is.null(maxVal)) {
    maxExcess <- maxVal - thresh
    maxOrig <- maxExcess

    if (eps_type == "quantile") {
      eps <- quantile(data, eps)
    }

    maxExcess <- maxExcess + eps

  } else {
    maxExcess <- maxOrig <- NULL
  }

  # Contraint of a positive shape parameter
  if (constraint == "shape_nonneg") {
    if (!fit_method %in% c("MLE1D", "MLE2D", "NLS2")) {
      stop("Constraint \"shape_nonneg\" only available for methods ",
           "MLE1D, MLE2D, and NLS2.")
    }

    shapeMin <- 0
    shapePos <- TRUE
  } else {
    shapeMin <- -Inf
    shapePos <- FALSE
  }

  # GPD fit
  if (fit_method == "MLE2D") {
    fit <- try(.fit_gpd_mle2d(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                             tol = tol, shapeMin = shapeMin),
                  silent = TRUE)

  } else if (fit_method == "MLE1D") {
    fit <- try(.fit_gpd_mle1d(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                             tol = tol, shapePos = shapePos),
                  silent = TRUE)

  } else if (fit_method == "LME") {
    fit <- try(.fit_gpd_lme(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                           tol = tol),
                  silent = TRUE)

  } else if (fit_method == "NLS2") {
    fit <- try(.fit_gpd_nls2(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                            tol = tol, shapeMin = shapeMin),
                  silent = TRUE)

  } else if (fit_method == "WNLLSM") {
    fit <- try(.fit_gpd_wnllsm(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                              tol = tol),
                  silent = TRUE)

  } else if (fit_method == "ZSE") {
    fit <- try(.fit_gpd_zse(x = excess, maxX = maxExcess, maxXOrig = maxOrig),
                  silent = TRUE)

  } else if (fit_method == "MOM") {
    fit <- try(.fit_gpd_mom(x = excess),
                  silent = TRUE)
  }

  if ("try-error" %in% class(fit)) {
    shape <- scale <- negLogLik <- NA
    p_value <- 0

  } else {
    shape <- fit$shape
    scale <- fit$scale
    names(shape) <- names(scale) <- NULL

    #negLogLik <- fit$negLogLik

    if (gof_test == "none") {
      p_value <- NULL

    } else {
      if (gof_test == "ad") {
        testres <- try(.gof_gpd_ad(excess,
                                   scale = scale, shape = shape), silent = TRUE)
      } else if (gof_test == "cvm") {
        testres <- try(.gof_gpd_cvm(excess,
                                    scale = scale, shape = shape), silent = TRUE)
      }

      p_value <- ifelse("try-error" %in% class(testres), 0, testres$p.value)
    }

  }

  return(list(shape = shape, scale = scale, p_value = p_value))#, negLogLik = negLogLik))
}

