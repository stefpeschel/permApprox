#' @title One-dimensional MLE for the GPD
#'
#' @description
#'   One-dimensional Maximum Likelihood Estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD) proposed by
#'   \cite{Castillo & Serra (2015)}.
#'
#' @param x data vector
#' @param maxX numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param eps numeric. Optional small epsilon, which is added to \code{maxX}
#'   so that the GPD density is non-zero at maxX + eps.
#' @param shapePos logical. If \code{TRUE}, the optimization is done under the
#'   constraint of a positive shape parameter.
#' @param optimTol numeric giving the desired accuracy of the optimization
#'   process.
#' @references
#'   \insertRef{Castillo2015likelihood}{permpap}
#'
#' @importFrom Rdpack reprompt
#' @export

gpd_MLE1D <- function(x,
                      maxX = NULL,
                      shapePos = FALSE,
                      optimTol = 1e-8,
                      eps = 0) {

  # Positive-shape constraint
  if (shapePos) {
    int <- c(-100*max(x), 0)
  } else {
    int <- c(-100*max(x), 100*max(x))
  }

  # Actual maximum value (GPD density must be non-zero at this value)
  maxXact <- max(c(x, maxX + eps))

  if (is.null(maxX)) {
    maxX <- maxXact
  }

  # Optimization
  sigma <- optimize(.MLE1D_fp, interval=int, maxX = maxXact, maximum = FALSE,
                    tol = optimTol)$minimum

  shape <- -.MLE1D_fk(sigma)
  scale <- -shape * sigma

  # GPD density at maxX
  densMax <- VGAM::dgpd(maxX, scale = scale, shape = shape)

  out <- list(shape = shape,
              scale = scale,
              bound = sigma,
              densMax = densMax)

  class(out) <- "GPDest"

  return(out)
}


# Function for calculating shape from sigma
.MLE1D_fk <- function(sigma) {
  -mean(log(1-x/sigma))
}


# Function to minimize
.MLE1D_fp <- function(sigma, maxX) {
  if ((sigma > 0) && (maxX > sigma)) {
    out <- 1e+6
  } else {
    out <- -length(x)*(-log(.MLE1D_fk(sigma)*sigma)+.MLE1D_fk(sigma)-1)
  }

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
# gpdfit <- gpd_MLE1D(x, maxX = NULL, shapePos = FALSE)
# gpdfit$shape
#
# gpdfit <- gpd_MLE1D(x, maxX = 5, shapePos = FALSE)
# gpdfit$shape
#
# gpdfit <- gpd_MLE1D(x, maxX = 5, shapePos = TRUE)
# gpdfit$shape
#
# # -> fit is independent of xobs and shapePos constraint for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- gpd_MLE1D(x, maxX = NULL, shapePos = FALSE)
# gpdfit
#
# gpdfit <- gpd_MLE1D(x, maxX = 5, shapePos = FALSE)
# gpdfit
# -gpdfit$scale / gpdfit$shape
#
# gpdfit <- gpd_MLE1D(x, maxX = 5, shapePos = TRUE)
# gpdfit
#
# # -> works perfectly



