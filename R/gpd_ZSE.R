#' @title Parameter estimation for the GPD proposed by Zhang & Stephens
#'
#' @description
#'   Parameter estimation for the two-parameter Generalized Pareto Distribution
#'   (GPD) proposed by \cite{Zhang & Stephens (2009)}.
#'
#' @param x data vector
#' @param q threshold value. E.g., set to q=0.9 to use the 90% percentile.
#'   Default is 0, i.e., all data are used.
#' @param maxX numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param eps numeric. Optional small epsilon, which is added to \code{maxX}
#'   so that the GPD density is non-zero at maxX + eps.
#' @param optimTol numeric giving the desired accuracy of the optimization
#'   process.
#' @param shapeIni,scaleIni initial values for the shape and scale parameters,
#'   i.e., values where the optimization starts.
#' @param shapeMin,scaleMin lower bound on the shape and scale parameters.
#'   If != -Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @param shapeMax,scaleMax upper bound on the shape and scale parameters.
#'   If != Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @references
#'   \insertRef{Zhang2009new}{permAprox}
#'
#' @importFrom Rdpack reprompt
#' @export


gpd_ZSE <- function(x, maxX = NULL, maxXOrig = NULL, m = NULL) {

  # constr: "none", "shapePos", "maxX"
  # shapePos doesn't work

  n <- length(x)
  x <- sort(x)

  if (is.null(m)) {
    m <- 20 + floor(sqrt(n))
  } else {
    stopifnot(m > 20)
  }

  # Actual maximum value (GPD density must be non-zero at this value)
  maxXact <- max(c(x, maxX))

  if (is.null(maxX)) {
    maxX <- maxXOrig <- maxXact
  }

  b <- w <- L <-
    1 / maxXact + (1 - sqrt(m / (1:m - 0.5))) / 3 / x[floor(n / 4 + 0.5)]

  for (i in 1:m){
    L[i] <- n * .ZSE_lx(b[i], x)
  }

  for (i in 1:m){
    w[i]<- 1/sum(exp(L-L[i]))
  }

  theta <- sum(b * w)

  k <- -mean(log(1 - theta * x))

  sigma <- k / theta

  shape <- -k
  scale <- sigma
  bound <- - scale / shape

  densMax <- VGAM::dgpd(maxXOrig, scale = scale, shape = shape)

  out <- list(shape = shape,
              scale = scale,
              bound = bound,
              densMax = densMax)

  class(out) <- "GPDest"

  return(out)
}


.ZSE_lx <- function(b, x) {
  k <- -mean(log(1 - b * x))
  log(b / k) + k - 1
}


# library(POT)
#
# # positive shape
# scale <- 1
# shape <- 0.25
#
# set.seed(123456)
# x <- rgpd(n = 100, scale = scale, shape = shape)
#
# gpdfit <- ZSE(x, constr = "none", maxX = NULL)
# gpdfit
#
# gpdfit <- ZSE(x, constr = "maxX", maxX = 8)
# gpdfit
#
# # -> fit is independent of maxX for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- rgpd(n = 100, scale = scale, shape = shape)
#
# gpdfit <- ZSE(x, constr = "none", maxX = NULL)
# gpdfit
#
# gpdfit <- ZSE(x, constr = "shapePos", maxX = NULL)
# gpdfit
# # shapePos constraint doesn't work
#
# gpdfit <- ZSE(x, constr = "maxX", maxX = 8)
# gpdfit
#
# -gpdfit$scale / gpdfit$shape
# # -> much larger than maxX



