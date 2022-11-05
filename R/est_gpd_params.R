# Function for estimating GPD parameters from a data vector

.est_gpd_params <- function(data,
                            thresh,
                            maxVal = NULL,
                            tol = 1e-8,
                            eps = 0,
                            fitMethod = "MLE1D",
                            ...) {

  data <- as.numeric(data)

  exceedances <- data[data > thresh]

  excess <- exceedances - thresh

  if (is.null(maxVal)) {
    maxValEx <- NULL
  } else {
    maxValEx <- maxVal - thresh
  }

  # GPD fit
  if (fitMethod == "MLE2D") {
    gpdfit <- gpd_MLE2D(x = excess, maxX = maxValEx, eps = eps,
                        tol = tol, ...)

  } else if (fitMethod == "MLE1D") {
    gpdfit <- gpd_MLE1D(x = excess, maxX = maxValEx, tol = tol, eps = eps, ...)

  } else if (fitMethod == "LME") {
    gpdfit <- gpd_LME(x = excess, maxX = maxValEx, tol = tol, eps = eps, ...)

  } else if (fitMethod == "NLS2") {
    gpdfit <- gpd_NLS2(x = excess, maxX = maxValEx, tol = tol, eps = eps, ...)

  } else if (fitMethod == "WNLLSM") {
    gpdfit <- gpd_WNLLSM(x = excess, maxX = maxValEx, tol = tol, eps = eps, ...)

  } else if (fitMethod == "ZSE") {
    gpdfit <- gpd_ZSE(x = excess, maxX = maxValEx, eps = eps, ...)

  } else if (fitMethod == "MOM") {
    gpdfit <- gpd_MOM(x = excess)
  }

  return(gpdfit)
}
