#' @title Likelihood Moment Estimation for the GPD
#'
#' @description Likelihood Moment Estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD) proposed by \cite{Zhang (2007)}.
#'
#' @param x data vector
#' @param maxX numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param r numeric defining the r constant as part of the algorithm that must
#'   be smaller than 1/2 and different from zero. Default is -1/2.
#'   See details.
#' @param eps numeric. Optional small epsilon, which is added to \code{maxX}
#'   so that the GPD density is non-zero at maxX + eps.
#' @param tol tolerance used as convergence criterion. The iteration stops if
#'   |f(b)/f'(b)/b| < tol. Defaults to 1e-8.
#'
#' @details
#'   According to \cite{Zhang (2007)}, the constant r must be smaller than 1/2
#'   and different from zero (r < 1/2 and r != 0).
#'   r close to -shape leads to estimates close to the MLEs.
#'   Without knowledge about the true shape, r = -1/2 (default)
#'   "produces LMEs with high asymptotic efficiency" \cite{(Zhang 2007)}.
#'   \cr\cr
#'
#' @references
#'   \insertRef{Zhang2007lme}{permpap}
#'
#' @importFrom Rdpack reprompt
#' @export

gpd_LME <- function(x, maxX = NULL, maxXOrig = NULL, r = -1/2, tol = 1e-8) {

  # Starting value for b
  b <- -1

  # Maximum value (GPD density must be non-zero at this value)
  xn <- max(c(x, maxX)) # Edited by SP

  if (is.null(maxX)) {
    maxX <- maxXOrig <- xn
  }

  # Starting value for error estimate
  err <- 0

  for (i in 1:100) {
    B <- 1 - b * x
    K <- log(1 - b * x)
    k <- -mean(K)

    # Function to minimize
    gb <- mean(B^(-r/k)) - 1/(1 - r)

    # Derivative of gb
    gd <- mean(B^(-r/k) * (x/B * k + K * mean(x/B))) * r/k/k

    # Estimate of b
    b <- min(b - gb/gd, (1 - 1/2^i)/xn)

    # Current error estimate
    errTmp <- abs(gb/gd/b)  # Added by SP

    # The loop breaks if the error estimate is smaller than "tol"
    # or if the error estimates do not change anymore
    if (errTmp < tol || abs((errTmp - err)) < tol) {  # Edited by SP
      break
    }

    err <- errTmp
  }

  sigma <- k/b

  #-----------------------------------------------------------------------------
  # The following has been added by Stefanie Peschel:

  shape <- -k
  scale <- sigma
  bound <- -scale / shape

  # GPD density at maxX
  densMax <- VGAM::dgpd(maxXOrig, scale = scale, shape = shape)

  out <- list(shape = shape,
              scale = scale,
              bound = bound,
              densMax = densMax,
              nLoops = i)

  class(out) <- "GPDest"

  return(out)
}


# #-------------------------------------------------------------------------------
# # positive shape
# scale <- 1
# shape <- 0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- gpd_LME(x, maxX = NULL)
# gpdfit
#
# gpdfit <- gpd_LME(x, maxX = 8)
# gpdfit
#
# # -> fit is independent of maxX for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- gpd_LME(x, maxX = NULL, r = 0.25)
# gpdfit
#
# gpdfit <- gpd_LME(x, maxX = 8)
# gpdfit
#
# -gpdfit$scale / gpdfit$shape
#
# # -> works perfectly



