# Function for estimating GPD parameters from a data vector

.est_gpd_params <- function(data,
                            thresh,
                            maxVal,
                            constraint,
                            tol,
                            eps,
                            epsType,
                            fitMethod,
                            ...) {

  data <- as.numeric(data)

  exceedances <- data[data > thresh]

  excess <- exceedances - thresh

  if (!is.null(maxVal)) {
    maxExcess <- maxVal - thresh
    maxOrig <- maxExcess

    if (epsType == "quantile") {
      eps <- quantile(data, eps)
    }

    maxExcess <- maxExcess + eps

  } else {
    maxExcess <- maxOrig <- NULL
  }

  if (constraint == "shapePos") {
    if (!fitMethod %in% c("MLE1D", "MLE2D", "NLS2")) {
      stop("Constraint \"shapePos\" only available for methods ",
           "MLE1D, MLE2D, and NLS2.")
    }

    shapeMin <- 0
    shapePos <- TRUE
  } else {
    shapeMin <- -Inf
    shapePos <- FALSE
  }

  # GPD fit
  if (fitMethod == "MLE2D") {
    gpdfit <- gpd_MLE2D(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                        tol = tol, shapeMin = shapeMin)

  } else if (fitMethod == "MLE1D") {
    gpdfit <- gpd_MLE1D(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                        tol = tol, shapePos = shapePos)

  } else if (fitMethod == "LME") {
    gpdfit <- gpd_LME(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                      tol = tol)

  } else if (fitMethod == "NLS2") {
    gpdfit <- gpd_NLS2(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                       tol = tol, shapeMin = shapeMin)

  } else if (fitMethod == "WNLLSM") {
    gpdfit <- gpd_WNLLSM(x = excess, maxX = maxExcess, maxXOrig = maxOrig,
                         tol = tol)

  } else if (fitMethod == "ZSE") {
    gpdfit <- gpd_ZSE(x = excess, maxX = maxExcess, maxXOrig = maxOrig)

  } else if (fitMethod == "MOM") {
    gpdfit <- gpd_MOM(x = excess)
  }

  return(gpdfit)
}
