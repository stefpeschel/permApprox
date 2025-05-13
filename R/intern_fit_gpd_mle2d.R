#' @title Two-dimensional MLE for the GPD
#'
#' @description
#'   Two-dimensional Maximum Likelihood Estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD).
#'
#' @param x data vector
#' @param maxX numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param optimMethod character defining the optimization method.
#' @param tol numeric giving the desired accuracy of the optimization
#'   process.
#' @param shapeIni,scaleIni initial values for the shape and scale parameters,
#'   i.e., values where the optimization starts.
#' @param shapeMin,scaleMin lower bound on the shape and scale parameters.
#'   If != -Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @param shapeMax,scaleMax upper bound on the shape and scale parameters.
#'   If != Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @param ... Further arguments passed to \code{optim()}.
#'
#' @keywords internal

.fit_gpd_mle2d <- function(x,
                      maxX = NULL,
                      maxXOrig = NULL,
                      optimMethod = "Nelder-Mead",
                      tol = 1e-8,
                      shapeIni = NULL,
                      scaleIni = NULL,
                      shapeMin = -Inf,
                      scaleMin = -Inf,
                      shapeMax = Inf,
                      scaleMax = Inf,
                      ...) {

  # x must be numeric
  x <- as.numeric(x)

  # Use box constraint estimation if parameters are constraint
  if (!(shapeMin == -Inf && scaleMin == -Inf &&
       shapeMax == Inf && scaleMax == Inf)) {
    optimMethod <- "L-BFGS-B"
  }

  if (optimMethod == "L-BFGS-B") {
    optcontr <- list(factr = tol)
  } else {
    optcontr <- list(reltol = tol)
  }

  # Method of moments estimator
  momEst <- .fit_gpd_mom(x)

  # Set start values for shape and scale (use method of moments)
  if (is.null(shapeIni)) {
    if (momEst$shape < 0) {
      shapeIni <- -0.1
    } else {
      shapeIni <- 0.1
    }
  }

  if (is.null(scaleIni)) {
    scaleIni <- 1
  }

  params <- c(shapeIni, scaleIni)

  # Actual maximum value (GPD density must be non-zero at this value)
  maxXact <- max(c(x, maxX))

  if (is.null(maxX)) {
    maxX <- maxXOrig <- maxXact
  }

  fit <- optim(par = params, fn = .MLE2D_negloglik, method = optimMethod,
               lower = c(shapeMin, scaleMin), upper = c(shapeMax, scaleMax),
               control = optcontr,
               ...,
               x = x, maxX = maxXact)

  if (fit$convergence != 0) {
    warning("GPD fit may not have succeeded.")
  }

  shape <- fit$par[1]
  scale <- fit$par[2]
  bound <- - scale / shape

  densMax <- VGAM::dgpd(maxXOrig, scale = scale, shape = shape)

  out <- list(shape = shape,
              scale = scale,
              bound = bound,
              densMax = densMax,
              negLogLik = fit$value,
              optimRes = fit)

  class(out) <- "GPDest"

  return(out)
}


.MLE2D_negloglik <- function(params, x, maxX) {
  shape <- params[1]
  scale <- params[2]

  cond1 <- scale <= 0
  cond2 <- (shape <= 0) && (maxX > (-scale/shape))

  if (cond1 || cond2) {
    nll <- 1e+6

  } else {
    y <- 1 / shape * log(1 + (shape * x) / scale)
    ll <- -length(x) * log(scale) - (1 + shape) * sum(y)
    nll <- -ll
  }
  return(nll)
}


#-------------------------------------------------------------------------------
# # positive shape
# scale <- 1
# shape <- 0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_mle2d(x, maxX = NULL)
# gpdfit$shape
#
# gpdfit <- .fit_gpd_mle2d(x, maxX = 5)
# gpdfit$shape
#
# gpdfit <- .fit_gpd_mle2d(x, maxX = 5, shapeMin = 0)
# gpdfit$shape
#
# # -> fit is nearly independent of maxX and shapePos constraint for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_mle2d(x, maxX = NULL)
# gpdfit
#
# gpdfit <- .fit_gpd_mle2d(x, maxX = 5)
# gpdfit
# -gpdfit$scale / gpdfit$shape
#
# gpdfit <- .fit_gpd_mle2d(x, maxX = 5, shapeMin = 0)
# gpdfit
#
# # -> works perfectly
