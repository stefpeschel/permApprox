#' @title WNLSM and WNLLSM parameter estimation for the GPD
#'
#' @description
#'   Weighted nonlinear least squares moments (WNLSM) estimation and
#'   likelihood WNLSM (WNLLSM) estimation for the two-parameter Generalized
#'   Pareto Distribution (GPD) proposed by \cite{Zhao et al. (2019)}.
#'
#' @param x data vector
#' @param q threshold value. E.g., set to q=0.9 to use the 90% percentile.
#'   Default is 0, i.e., all data are used.
#' @param maxX numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param eps numeric. Optional small epsilon, which is added to \code{maxX}
#'   so that the GPD density is non-zero at maxX + eps.
#'   \code{eps} must be larger than zero.
#' @param tol numeric giving the desired accuracy of the optimization
#'   process.
#' @param shapeIni,scaleIni initial values for the shape and scale parameters,
#'   i.e., values where the optimization starts.
#' @param shapeMin,scaleMin lower bound on the shape and scale parameters.
#'   If != -Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @param shapeMax,scaleMax upper bound on the shape and scale parameters.
#'   If != Inf, \code{optimMethod} is set to "L-BFGS-B".
#'
#' @details A simplified version of the R code has been provided by
#'   \cite{Zhao et al. (2019)}. The function has been extended by the arguments:
#'   \code{maxX}, \code{eps}, and \code{tol} to ensure that the GPD
#'   density is non-zero at maxX + eps or max(x) + eps.
#'
#' @references
#'   \insertRef{Zhao2019new}{permAprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal

.fit_gpd_wnllsm <- function(x,
                       maxX = NULL,
                       maxXOrig = NULL,
                       method = "WNLLSM",
                       tol =  1e-8,
                       shapeIni = 0.01,
                       scaleIni = 0.01,
                       F0 = 0,
                       plotopt = FALSE,
                       ...) {

  # Parameters:
  # xi: shape
  # sigma: scale
  # b: xi/sigma
  # beta: (xi, b)

  #stopifnot(eps > 0)

  x <- sort(x)
  n <- length(x)

  bIni <- shapeIni / scaleIni
  betaIni <- c(shapeIni, bIni)

  nu0 <- floor(n * F0)

  # Get threshold u
  if (nu0==0){
    u = 0
  } else {
    u <- x[nu0]
  }

  # Empirical Cumulative Distribution Function
  y <- ecdf(x)
  cf<- y(x)

  # Exceedances
  z <- x[x>u]

  if (method == "WNLLSM"){

    #---------------------------------------------------------------------------
    # Original code proposed by Zhao et al.:
    # b_WLLS1<- optim(0.1, WLLS1, method="Brent", lower=0, upper=10)$par
    # b_WLLS<- optim(b_WLLS1, WLLS, method="Brent", lower=0, upper=10)$par
    # xi_WLLS<- mean(log(1+b_WLLS*z))
    # bets_WLLS <- c(xi_WLLS,b_WLLS)
    #---------------------------------------------------------------------------

    # Lower limit depending on the maximum value, at which the GPD density must
    # be positive (edited by SP)

    lower <- -1 / max(c(x, maxX))

    if (is.null(maxX)) {
      maxX <- maxXOrig <- max(x)
    }

    # Increase lower limit because b must be larger than 'lower'
    lower <- lower


    b_WLLS1 <- optim(par = bIni, fn = .WNLLSM_WLLS1, method = "Brent",
                     lower = lower, upper = 10,
                     control = list(reltol = tol),
                     #...,
                     z = z, x = x, nu0 = nu0, n = n, F0 = F0, u = u)$par

    b_WLLS <- optim(par = b_WLLS1, fn = .WNLLSM_WLLS, method = "Brent",
                    lower = lower, upper = 10,
                    control = list(reltol = tol),
                    #...,
                    z = z, x = x, nu0 = nu0, n = n, F0 = F0, u = u)$par

    b <- b_WLLS
    shape <- mean(log(1 + b * z))
    scale <- shape / b

    if (plotopt){
      vals <- data.frame(par = seq(from = lower, to = 1, by = 0.01), f = NA)
      vals$f <-  sapply(1:nrow(vals), function(i)
        .WNLLSM_WLLS(b = vals[i, "par"], z = z, x = x, nu0 = nu0,
                     n = n, F0 = F0, u = u))

      for (i in 1:nrow(vals)){
        vals[i, "f"] <- .WNLLSM_WLLS(b = vals[i, "par"], z = z, x = x, nu0 = nu0,
                                     n = n, F0 = F0, u = u)
      }

      plot(vals$par, vals$f, type = "l")
      abline(v = b_WLLS, col = "blue")
    }


  } else {
    # Weighted Nonlinear Least Squares Moments (WNLSM) Estimation

    #---------------------------------------------------------------------------
    # Original code proposed by Zhao et al.:
    # beta_WLS1<-optim(beta,WLS1)$par
    # beta_WLS<-optim(root_WLS1,WLS)$par
    #---------------------------------------------------------------------------

    # Edited by SP:

    # Set lower limit for b and shape
    bMin <- -1 / max(c(x, maxX))
    shapeMin <- -Inf

    if (is.null(maxX)) {
      maxX <- maxXOrig <- max(x)
    }

    # increase lower limit because b must be larger than lower
    #bMin <- bMin + sqrt(.Machine$double.eps)
    bMin <- bMin + eps

    # Set upper limit for b and shape
    bMax <- 10
    shapeMax <- Inf


    # First step
    beta_WLS1 <- optim(par = betaIni, fn = .WNLSM_WLS1,
                       control = list(reltol = tol),
                       #...,
                       x = x, nu0 = nu0, n = n, F0 = F0, u = u)$par

    # Fitting the first step with constraint not needed:
    # beta_WLS1 <- optim(par = betaIni, fn = .WNLSM_WLS1,
    #                    method = "L-BFGS-B",
    #                    lower = c(shapeMin, bMin),
    #                    upper = c(shapeMax, bMax),
    #                    control = list(factr = tol),
    #                    #...,
    #                    x = x, nu0 = nu0, n = n, F0 = F0, u = u)$par

    #------------------------------------
    # Second step
    # Fit always with constraint to ensure a positive density at max(x) and maxX

    # Estimate beta
    beta_WLS <- optim(par = beta_WLS1, fn = .WNLSM_WLS,
                      method = "L-BFGS-B",
                      lower = c(shapeMin, bMin),
                      upper = c(shapeMax, bMax),
                      control = list(factr = tol),
                      #...,
                      x = x, nu0 = nu0, n = n, F0 = F0, u = u)$par

    shape <- beta_WLS[1]

    # scale = xi / b
    scale <- shape / beta_WLS[2]
    b <- beta_WLS[2]
  }

  bound <- - scale / shape

  densMax <- VGAM::dgpd(maxXOrig, scale = scale, shape = shape)

  out <- list(shape = shape,
              scale = scale,
              bound = bound,
              densMax = densMax,
              b = b)

  class(out) <- "GPDest"

  return(out)
}


#-------------------------------------------------------------------------------
#WNLLSM helper functions

.WNLLSM_sum_i <- function(i, n){
  j<- 1:i
  sum(1/(n-j+1))
}

.WNLLSM_WLLS1 <- function(b, z, x, nu0, n, F0, u){
  if (b == 0){
    sum2 <- 1e+6

  } else {
    sum1 <- mean(log(1+b*z))  #xi
    sum2 <- 0
    for (i in (nu0+1):n){
      sum2 <- sum2 + (.WNLLSM_sum_i(i, n)+log(1-F0)-log(1+b*(x[i]-u))/sum1)^2
    }
  }

  sum2
}

.WNLLSM_WLLS <- function(b, z, x, nu0, n, F0, u){
  if (b == 0){
    sum2 <- 1e+6

  } else {
    sum1 <- mean(log(1+b*z))  #xi
    sum2 <- 0
    for (i in (nu0+1):n)
    {
      sum2 <- sum2 + ((n+1)^2*(n+2)/(i*(n-i+1)))*((i-n-1)/((n+1)*(1-F0))+(1+b*(x[i]-u))^(-1/sum1))^2
    }
  }
  sum2
}

#-------------------------------------------------------------------------------
#WNLSM (Least square estimator) helper functions

.WNLSM_WLS1 <- function(beta1, x, nu0, n, F0, u){
  if (beta1[1] == 0 & beta1[2] == 0){
    sum2 <- 1e+6

  } else {
    sum2 <- 0
    for (i in (nu0 + 1):n)
    {
      sum2 <- sum2 + (.WNLLSM_sum_i(i, n)+log(1 - F0) -
                        log(1 + beta1[2] * (x[i] - u)) / beta1[1])^2
    }
    sum2
  }
}

.WNLSM_WLS <- function(beta1, x, nu0, n, F0, u){
  if (beta1[1] == 0 & beta1[2] == 0){
    sum2 <- 1e+6

  } else {
    sum2 <- 0
    for (i in (nu0 + 1):n)
    {
      sum2 <- sum2 + ((n + 1)^2 * (n + 2) / (i * (n - i + 1))) *
        ((i - n - 1) / ((n + 1) * (1 - F0)) +
           (1 + beta1[2] * (x[i] - u))^(-1 / beta1[1]))^2
    }
    sum2
  }
}


# #-------------------------------------------------------------------------------
# # positive shape
# scale <- 1
# shape <- 0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = NULL, method = "WNLLSM")
# gpdfit$shape
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = 8,  method = "WNLLSM")
# gpdfit$shape
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = NULL, method = "WNLSM")
# gpdfit$shape
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = 8,  method = "WNLSM")
# gpdfit$shape
#
# # -> fit is not effected for positive shapes
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = NULL, method = "WNLLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = 8,  method = "WNLLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = 8, eps = 1e-5,  method = "WNLLSM")
# gpdfit$shape
# gpdfit$densMax
# # density depends on eps
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = NULL, method = "WNLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = 8, eps = 1e-8, method = "WNLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, maxX = 8, eps = 1e-5, method = "WNLSM")
# gpdfit$shape
# gpdfit$densMax
#
# -gpdfit$scale / gpdfit$shape

