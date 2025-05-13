#' @title Parameter estimation via two-step NLS for the GPD
#'
#' @description
#'   Parameter estimation for the two-parameter Generalized Pareto Distribution
#'   (GPD) using a two-step Nonlinear Least Squares (NLS) approach proposed by
#'   \cite{Song & Song (2012)}.
#'
#' @param x data vector
#' @param q threshold value. E.g., set to q=0.9 to use the 90% percentile.
#'   Default is 0, i.e., all data are used.
#' @param maxX numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param tol numeric giving the desired accuracy of the optimization
#'   process.
#' @param shapeIni,scaleIni initial values for the shape and scale parameters,
#'   i.e., values where the optimization starts.
#' @param shapeMin,scaleMin lower bound on the shape and scale parameters.
#'   If != -Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @param shapeMax,scaleMax upper bound on the shape and scale parameters.
#'   If != Inf, \code{optimMethod} is set to "L-BFGS-B".
#' @references
#'   \insertRef{Song2012quantile}{permAprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal

.fit_gpd_nls2 <- function(x,
                     q = 0,
                     maxX = NULL,
                     maxXOrig = NULL,
                     twosteps = TRUE,
                     optimMethod = "Nelder-Mead",
                     tol = 1e-8,
                     shapeIni = NULL,
                     scaleIni = NULL,
                     shapeMin = -Inf,
                     scaleMin = -Inf,
                     shapeMax = Inf,
                     scaleMax = Inf,
                     ...) {

  #This function returns 2 sets of estimated parameters
  #First set is from NLS-1, second set is from NLS-2.
  #Our simulation study shows that NLS-2 gives better results.

  #-----------------------------------------------------------------------------
  # SP
  if (!(shapeMin == -Inf && scaleMin == -Inf &&
       shapeMax == Inf && scaleMax == Inf)) {
    optimMethod <- "L-BFGS-B"
  }

  if (optimMethod == "L-BFGS-B") {
    optcontr <- list(factr = tol)
  } else {
    optcontr <- list(reltol = tol)
  }

  # Compute initial values via method of moments
  s2 <- var(x)
  m  <- mean(x)

  if (is.null(shapeIni)) {
    if (is.null(maxX)) {
      shapeIni <- -(0.5 * (m^2 / s2 - 1))
    } else {
      shapeIni <- 0.01
    }
  }

  if (is.null(scaleIni)) {
    scaleIni <- 0.5 * m * (m^2 / s2 + 1)
  }

  theta <- c(shapeIni, scaleIni)

  n <- length(x)

  if (q == 0) {
    x.u <- x

  } else {
    x.q <- quantile(x, q)
    x.u <- x[x > x.q]
  }

  x.u <- sort(x.u)

  # if (is.null(maxX)) {
  #   maxX <- x[1]
  # }

  # Actual maximum value (GPD density must be non-zero at this value)
  maxXact <- max(c(x, maxX))

  if (is.null(maxX)) {
    maxX <- maxXOrig <- maxXact
  }

  # Modified optim call to estimate initial values
  res.ini <- optim(c(scaleIni, shapeIni), fn = .NLS2_rss1,
                   method = optimMethod,
                   lower = c(scaleMin, shapeMin),
                   upper = c(scaleMax, shapeMax),
                   control = optcontr,
                   ...,
                   n = n, q = q, x.u = x.u, maxX = maxXact)$par

  scale <- res.ini[1]
  shape <- res.ini[2]

  res.fin <- NULL

  if (twosteps) {
    #res.fin <- optim(res.ini, .NLS2_rss2, x.u=x.u)$par

    # Modified optim call
    res.fin <- optim(res.ini, fn = .NLS2_rss2,
                     method = optimMethod,
                     lower = c(scaleMin, shapeMin),
                     upper = c(scaleMax, shapeMax),
                     control = optcontr,
                     ...,
                     n = n, q = q, x.u = x.u, maxX = maxXact)$par

    # Two-step output
    scale <- res.fin[1]
    shape <- res.fin[2]
  }

  bound <- - scale / shape

  densMax <- VGAM::dgpd(maxXOrig, scale = scale, shape = shape)

  out <- list(shape = shape,
              scale = scale,
              bound = bound,
              densMax = densMax,
              res.ini = res.ini,
              res.fin = res.fin)

  class(out) <- "GPDest"

  return(out)
}



#-------------------------------------------------------------------------------
# Helpers

# GPD distribution function F(x)
.NLS2_gpdf <- function(x, scale, shape) {

  # Edited by SP: term is set to 0 if < 0
  tmp <- pmax(1 + shape * x / scale, 0)

  #res <- 1 - (1+shape*x/scale)^(-1/shape)
  res <- 1 - tmp^(-1 / shape)
  return(res)
}

# Empirical distribution function
.NLS2_ecdf <- function(x, n, q) {
  #n is the original sample size
  res<-seq((q * n + 1) / n, 1, 1 / n)
  return(res)
}

# Theoretical distribution function (GPD)
.NLS2_gpdf2 <- function(x, scale, shape) {
  res <- -log(1 + shape * x / scale) / shape
  return(res)
}

# Empirical distribution function (EDF)
.NLS2_ecdf2<-function(x,n,q) {
  res <- log(1 - .NLS2_ecdf(x, n, q) + 1 / (10 * n))
  return(res)
}

# Function to optimize in the first step
# params: c(scale, shape)
.NLS2_rss1 <- function(params, n, q, x.u, maxX) {
  scale <- params[1]
  shape <- params[2]
  maxval <- max(c(x.u, maxX))

  cond1 <- scale <= 0
  cond2 <- (shape <= 0) && (maxval > (-scale/shape))

  if (cond1 || cond2) {
    temp1 <- 1e+6
  } else {
    # Squared sum of EDF - GPD
    temp1 <- sum((.NLS2_ecdf2(x.u, n, q) - .NLS2_gpdf2(x.u, scale, shape))^2)
  }

  return(temp1)
}

# Function to optimize in the second step
.NLS2_rss2 <- function(params, n, q, x.u, maxX) {
  if (params[2] < 0 & maxX >= (-params[1] / params[2])) {
    temp1 <- Inf

  } else {
    temp1<-sum((.NLS2_ecdf(x.u, n, q)-.NLS2_gpdf(x.u, params[1], params[2]))^2)
  }
  if (is.infinite(temp1)) temp1 <- 1e6
  return(temp1)
}


################################################


# # positive shape
# scale <- 1
# shape <- 0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                maxX = NULL, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                maxX = NULL, twosteps = TRUE,
#                shapeMin = 0, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                maxX = 8, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                maxX = NULL, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                maxX = NULL, twosteps = TRUE,
#                shapeMin = 0, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                maxX = 8, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# -gpdfit$scale / gpdfit$shape
# # -> works well

