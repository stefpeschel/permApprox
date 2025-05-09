# Function for estimating GPD parameters from a data vector

.est_gpd_params <- function(data,
                            thresh,
                            maxVal,
                            constraint,
                            tol,
                            eps,
                            eps_type,
                            fit_method,
                            ...) {

  data <- as.numeric(data)

  exceedances <- data[data > thresh]

  excess <- exceedances - thresh

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
    gpdfit <- .fit_gpd_mle2d(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                        tol = tol, shapeMin = shapeMin)

  } else if (fit_method == "MLE1D") {
    gpdfit <- .fit_gpd_mle1d(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                        tol = tol, shapePos = shapePos)

  } else if (fit_method == "LME") {
    gpdfit <- .fit_gpd_lme(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                      tol = tol)

  } else if (fit_method == "NLS2") {
    gpdfit <- .fit_gpd_nls2(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                       tol = tol, shapeMin = shapeMin)

  } else if (fit_method == "WNLLSM") {
    gpdfit <- .fit_gpd_wnllsm(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                         tol = tol)

  } else if (fit_method == "ZSE") {
    gpdfit <- .fit_gpd_zse(x = excess, maxX = maxExcess, maxXOrig = maxOrig)

  } else if (fit_method == "MOM") {
    gpdfit <- .fit_gpd_mom(x = excess)
  }

  return(gpdfit)
}
