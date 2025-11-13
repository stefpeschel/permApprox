################################################################################
# MOM
################################################################################

# Method of Moments estimator for the GPD

#' @title Method of Moments estimation for the GPD
#'
#' @description
#'   Method of Moments estimation for the two-parameter Generalized Pareto
#'   Distribution (GPD) proposed by \cite{Hosking & Wallis (1987)}.
#'
#' @param x data vector
#'
#' @references
#'   \insertRef{Hosking1987parameter}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal

.fit_gpd_mom <- function(x) {
  
  # x must be numeric
  x <- as.numeric(x)
  
  # Mean and variance
  meanx <- mean(x)
  varx <- var(x)
  
  shape <- -0.5 * (((meanx^2) / varx) - 1)
  
  scale <- 0.5 * meanx * (((meanx^2) / varx) + 1)
  
  return(list(shape = shape, scale = scale))
}

################################################################################
# PWM
################################################################################

# Probability Weighted Moments estimator for the GPD

#' @title Probability Weighted Moments estimation for the GPD
#'
#' @description
#'   Probability Weighted Moments (PWM) estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD) as described by
#'   \insertCite{Hosking1987parameter;textual}{permApprox}.
#'
#' @param x data vector
#'
#' @details
#'   The PWM estimators are computed using the plotting positions
#'   \eqn{p_i = (i - 0.35) / n} proposed by \insertCite{Hosking1987parameter;textual}{permApprox}.
#'   Let \eqn{a_0 = \bar{x}} and
#'   \eqn{a_1 = \frac{1}{n} \sum_{i=1}^{n} x_{(i)} (1 - p_i)}, where
#'   \eqn{x_{(i)}} are the order statistics of \eqn{x}.
#'   The estimators for the shape (\eqn{\xi}) and scale (\eqn{\sigma}) parameters are then
#'   \deqn{
#'     \hat{\xi} = 2 - \frac{a_0}{a_0 - 2 a_1}, \qquad
#'     \hat{\sigma} = \frac{2 a_0 a_1}{a_0 - 2 a_1}.
#'   }
#'
#' @return A list containing
#'   \item{shape}{Estimated shape parameter \eqn{\xi}.}
#'   \item{scale}{Estimated scale parameter \eqn{\sigma}.}
#'
#' @references
#'   \insertRef{Hosking1987parameter}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @export
#' @keywords internal
.fit_gpd_pwm <- function(x) {
  
  # x must be numeric
  x <- as.numeric(x)
  n <- length(x)
  
  if (n < 2L) {
    return(list(shape = NA_real_, scale = NA_real_))
  }
  
  # Sort and compute PWMs
  x_sorted <- sort(x)
  pvec <- (seq_len(n) - 0.35) / n
  a0 <- mean(x_sorted)
  a1 <- mean(x_sorted * (1 - pvec))
  
  denom <- a0 - 2 * a1
  if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
    return(list(shape = NA_real_, scale = NA_real_))
  }
  
  shape <- 2 - a0 / denom
  scale <- (2 * a0 * a1) / denom
  
  return(list(shape = shape, scale = scale))
}

################################################################################
# LME
################################################################################

#' @title Likelihood Moment Estimation for the GPD
#'
#' @description Likelihood Moment Estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD) proposed by \cite{Zhang (2007)}.
#'
#' @param x data vector
#' @param boundary numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param r numeric defining the r constant as part of the algorithm that must
#'   be smaller than 1/2 and different from zero. Default is -1/2.
#'   See details.
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
#'   \insertRef{Zhang2007lme}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @importFrom eva dgpd
#' @keywords internal

.fit_gpd_lme <- function(x, boundary = NULL, eval_point = NULL, r = -1/2, tol = 1e-8) {
  
  # Starting value for b
  b <- -1
  
  # Maximum value (GPD density must be non-zero at this value)
  xn <- max(c(x, boundary)) # Edited by SP
  
  if (is.null(boundary)) {
    boundary <- eval_point <- xn
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
  
  # GPD density at boundary
  densMax <- eva::dgpd(eval_point, scale = scale, shape = shape)
  
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
# gpdfit <- .fit_gpd_lme(x, boundary = NULL)
# gpdfit
#
# gpdfit <- .fit_gpd_lme(x, boundary = 8)
# gpdfit
#
# # -> fit is independent of boundary for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_lme(x, boundary = NULL, r = 0.25)
# gpdfit
#
# gpdfit <- .fit_gpd_lme(x, boundary = 8)
# gpdfit
#
# -gpdfit$scale / gpdfit$shape
#
# # -> works perfectly

################################################################################
# MLE1D
################################################################################

#' @title One-dimensional MLE for the GPD
#'
#' @description
#'   One-dimensional Maximum Likelihood Estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD) proposed by
#'   \cite{Castillo & Serra (2015)}.
#'
#' @param x data vector
#' @param boundary numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param shapePos logical. If \code{TRUE}, the optimization is done under the
#'   constraint of a positive shape parameter.
#' @param tol numeric giving the desired accuracy of the optimization
#'   process.
#' @references
#'   \insertRef{Castillo2015likelihood}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal

.fit_gpd_mle1d <- function(x,
                           boundary = NULL,
                           eval_point = NULL,
                           shapePos = FALSE,
                           tol = 1e-8) {
  
  # Positive-shape constraint
  if (shapePos) {
    int <- c(-100*max(x), 0)
  } else {
    int <- c(-100*max(x), 100*max(x))
  }
  
  # Actual maximum value (GPD density must be non-zero at this value)
  actual_boundary <- max(c(x, boundary))
  
  if (is.null(boundary)) {
    boundary <- eval_point <- actual_boundary
  }
  
  # Optimization
  sigma <- optimize(.MLE1D_fp, interval=int, boundary = actual_boundary, x = x,
                    maximum = FALSE, tol = tol)$minimum
  
  shape <- -.MLE1D_fk(sigma, x)
  scale <- -shape * sigma
  
  # GPD density at boundary
  densMax <- eva::dgpd(eval_point, scale = scale, shape = shape)
  
  out <- list(shape = shape,
              scale = scale,
              bound = sigma,
              densMax = densMax)
  
  class(out) <- "GPDest"
  
  return(out)
}


# Function for calculating shape from sigma
.MLE1D_fk <- function(sigma, x) {
  # edit by SP
  tmp <- 1 - x / sigma
  tmp[tmp < 0] <- 0
  tmp <- log(tmp)
  tmp[is.infinite(tmp)] <- NA
  -mean(tmp, na.rm = TRUE)
}


# Function to minimize
.MLE1D_fp <- function(sigma, boundary, x) {
  if ((sigma > 0) && (boundary > sigma)) {
    out <- 1e+6
  } else {
    out <- -length(x) * (-log(.MLE1D_fk(sigma, x) * sigma) +
                           .MLE1D_fk(sigma, x) - 1)
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
# gpdfit <- .fit_gpd_mle1d(x, boundary = NULL, shapePos = FALSE)
# gpdfit$shape
#
# gpdfit <- .fit_gpd_mle1d(x, boundary = 5, shapePos = FALSE)
# gpdfit$shape
#
# gpdfit <- .fit_gpd_mle1d(x, boundary = 5, shapePos = TRUE)
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
# gpdfit <- .fit_gpd_mle1d(x, boundary = NULL, shapePos = FALSE)
# gpdfit
#
# gpdfit <- .fit_gpd_mle1d(x, boundary = 5, shapePos = FALSE)
# gpdfit
# -gpdfit$scale / gpdfit$shape
#
# gpdfit <- .fit_gpd_mle1d(x, boundary = 5, shapePos = TRUE)
# gpdfit
#
# # -> works perfectly

################################################################################
# MLE2D
################################################################################

#' @title Two-dimensional MLE for the GPD
#'
#' @description
#'   Two-dimensional Maximum Likelihood Estimation for the two-parameter
#'   Generalized Pareto Distribution (GPD).
#'
#' @param x data vector
#' @param boundary numeric. Maximum value at which the GPD probability density
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
                           boundary = NULL,
                           eval_point = NULL,
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
  
  # --- Start values (PWM-based with simple fallback) --------------------------
  eps <- 1e-8
  
  if (is.null(shapeIni) || is.null(scaleIni)) {
    
    # 1) Try PWM 
    # (boundary is included to get more reliable start values in constrained case)
    initEst <- .fit_gpd_pwm(c(x, boundary))
    
    # 2) If PWM degenerate, fall back to Method of Moments
    if (anyNA(unlist(initEst)) || !is.finite(initEst$shape) ||
        !is.finite(initEst$scale) || initEst$scale <= 0) {
      initEst <- .fit_gpd_mom(c(x, boundary))
    }
    
    shapeIni <- initEst$shape
    scaleIni <- initEst$scale
  }
  
  # Final sanity + box constraints
  if (!is.finite(shapeIni)) shapeIni <- 0
  if (!is.finite(scaleIni) || scaleIni <= 0)
    scaleIni <- max(mean(pmax(x, 0)), eps)
  
  shapeIni <- min(max(shapeIni, shapeMin), shapeMax)
  scaleIni <- min(max(scaleIni, scaleMin, eps),
                  if (is.finite(scaleMax)) scaleMax else scaleIni)
  
  params <- c(shapeIni, scaleIni)
  
  # ----------------------------------------------------------------------------
  
  # Actual maximum value (GPD density must be non-zero at this value)
  actual_boundary <- max(c(x, boundary))
  
  if (is.null(boundary)) {
    boundary <- eval_point <- actual_boundary
  }
  
  fit <- optim(par = params, fn = .MLE2D_negloglik, method = optimMethod,
               lower = c(shapeMin, scaleMin), upper = c(shapeMax, scaleMax),
               control = optcontr,
               ...,
               x = x, boundary = actual_boundary)
  
  if (fit$convergence != 0) {
    warning("GPD fit may not have succeeded.")
  }
  
  shape <- fit$par[1]
  scale <- fit$par[2]
  bound <- - scale / shape
  
  densMax <- eva::dgpd(eval_point, scale = scale, shape = shape)
  
  out <- list(shape = shape,
              scale = scale,
              bound = bound,
              densMax = densMax,
              negLogLik = fit$value,
              optimRes = fit)
  
  class(out) <- "GPDest"
  
  return(out)
}


.MLE2D_negloglik <- function(params, x, boundary) {
  shape <- params[1]
  scale <- params[2]
  
  cond1 <- scale <= 0
  cond2 <- (shape <= 0) && (boundary > (-scale/shape))
  
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
# gpdfit <- .fit_gpd_mle2d(x, boundary = NULL)
# gpdfit$shape
#
# gpdfit <- .fit_gpd_mle2d(x, boundary = 5)
# gpdfit$shape
#
# gpdfit <- .fit_gpd_mle2d(x, boundary = 5, shapeMin = 0)
# gpdfit$shape
#
# # -> fit is nearly independent of boundary and shapePos constraint for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- eva::rgpd(n = 1000, scale = scale, shape = shape)
#
# gpdfit <- .fit_gpd_mle2d(x, boundary = NULL)
# gpdfit
#
# gpdfit <- .fit_gpd_mle2d(x, boundary = 5)
# gpdfit
# -gpdfit$scale / gpdfit$shape
#
# gpdfit <- .fit_gpd_mle2d(x, boundary = 5, shapeMin = 0)
# gpdfit
#
# # -> works perfectly

################################################################################
# NLS2
################################################################################

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
#' @param boundary numeric. Maximum value at which the GPD probability density
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
#'   \insertRef{Song2012quantile}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal

.fit_gpd_nls2 <- function(x,
                          q = 0,
                          boundary = NULL,
                          eval_point = NULL,
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
    if (is.null(boundary)) {
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
  
  # if (is.null(boundary)) {
  #   boundary <- x[1]
  # }
  
  # Actual maximum value (GPD density must be non-zero at this value)
  actual_boundary <- max(c(x, boundary))
  
  if (is.null(boundary)) {
    boundary <- eval_point <- actual_boundary
  }
  
  # Modified optim call to estimate initial values
  res.ini <- optim(c(scaleIni, shapeIni), fn = .NLS2_rss1,
                   method = optimMethod,
                   lower = c(scaleMin, shapeMin),
                   upper = c(scaleMax, shapeMax),
                   control = optcontr,
                   ...,
                   n = n, q = q, x.u = x.u, boundary = actual_boundary)$par
  
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
                     n = n, q = q, x.u = x.u, boundary = actual_boundary)$par
    
    # Two-step output
    scale <- res.fin[1]
    shape <- res.fin[2]
  }
  
  bound <- - scale / shape
  
  densMax <- eva::dgpd(eval_point, scale = scale, shape = shape)
  
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
.NLS2_rss1 <- function(params, n, q, x.u, boundary) {
  scale <- params[1]
  shape <- params[2]
  support_boundary <- max(c(x.u, boundary))
  
  cond1 <- scale <= 0
  cond2 <- (shape <= 0) && (support_boundary > (-scale/shape))
  
  if (cond1 || cond2) {
    temp1 <- 1e+6
  } else {
    # Squared sum of EDF - GPD
    temp1 <- sum((.NLS2_ecdf2(x.u, n, q) - .NLS2_gpdf2(x.u, scale, shape))^2)
  }
  
  return(temp1)
}

# Function to optimize in the second step
.NLS2_rss2 <- function(params, n, q, x.u, boundary) {
  if (params[2] < 0 & boundary >= (-params[1] / params[2])) {
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
#                boundary = NULL, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                boundary = NULL, twosteps = TRUE,
#                shapeMin = 0, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                boundary = 8, twosteps = TRUE,
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
#                boundary = NULL, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                boundary = NULL, twosteps = TRUE,
#                shapeMin = 0, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# gpdfit <- .fit_gpd_nls2(x = x, q = 0,
#                boundary = 8, twosteps = TRUE,
#                shapeMin = -Inf, scaleMin = -Inf,
#                shapeMax = Inf, scaleMax = Inf)
# gpdfit
#
# -gpdfit$scale / gpdfit$shape
# # -> works well

################################################################################
# WNLLSM
################################################################################

#' @title WNLSM and WNLLSM parameter estimation for the GPD
#'
#' @description
#'   Weighted nonlinear least squares moments (WNLSM) estimation and
#'   likelihood WNLSM (WNLLSM) estimation for the two-parameter Generalized
#'   Pareto Distribution (GPD) proposed by \cite{Zhao et al. (2019)}.
#'
#' @param x data vector
#' @param boundary numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @param tol numeric giving the desired accuracy of the optimization
#'   process.
#' @param shapeIni,scaleIni initial values for the shape and scale parameters,
#'   i.e., values where the optimization starts.
#'
#' @details A simplified version of the R code has been provided by
#'   \cite{Zhao et al. (2019)}. The function has been adapted to ensure that
#'   the GPD density is non-zero at boundary + eps or max(x) + eps.
#'
#' @references
#'   \insertRef{Zhao2019new}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal

.fit_gpd_wnllsm <- function(x,
                            boundary = NULL,
                            eval_point = NULL,
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
    
    lower <- -1 / max(c(x, boundary))
    
    if (is.null(boundary)) {
      boundary <- eval_point <- max(x)
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
    bMin <- -1 / max(c(x, boundary))
    shapeMin <- -Inf
    
    if (is.null(boundary)) {
      boundary <- eval_point <- max(x)
    }
    
    # increase lower limit because b must be larger than lower
    bMin <- bMin + sqrt(.Machine$double.eps)
    
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
    # Fit always with constraint to ensure a positive density at max(x) and boundary
    
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
  
  densMax <- eva::dgpd(eval_point, scale = scale, shape = shape)
  
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
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = NULL, method = "WNLLSM")
# gpdfit$shape
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = 8,  method = "WNLLSM")
# gpdfit$shape
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = NULL, method = "WNLSM")
# gpdfit$shape
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = 8,  method = "WNLSM")
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
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = NULL, method = "WNLLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = 8,  method = "WNLLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = 8, eps = 1e-5,  method = "WNLLSM")
# gpdfit$shape
# gpdfit$densMax
# # density depends on eps
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = NULL, method = "WNLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = 8, eps = 1e-8, method = "WNLSM")
# gpdfit$shape
# gpdfit$densMax
#
# gpdfit <- .fit_gpd_wnllsm(x = x, boundary = 8, eps = 1e-5, method = "WNLSM")
# gpdfit$shape
# gpdfit$densMax
#
# -gpdfit$scale / gpdfit$shape

################################################################################
# ZSE
################################################################################

#' @title Parameter estimation for the GPD proposed by Zhang & Stephens
#'
#' @description
#'   Parameter estimation for the two-parameter Generalized Pareto Distribution
#'   (GPD) proposed by \cite{Zhang & Stephens (2009)}.
#'
#' @param x data vector
#' @param boundary numeric. Maximum value at which the GPD probability density
#'   function must be positive. If \code{NULL}, the maximum of the data vector
#'   \code{x} is used.
#' @references
#'   \insertRef{Zhang2009new}{permApprox}
#'
#' @importFrom Rdpack reprompt
#' @keywords internal


.fit_gpd_zse <- function(x, boundary = NULL, eval_point = NULL, m = NULL) {
  
  # constr: "none", "shapePos", "boundary"
  # shapePos doesn't work
  
  n <- length(x)
  x <- sort(x)
  
  if (is.null(m)) {
    m <- 20 + floor(sqrt(n))
  } else {
    stopifnot(m > 20)
  }
  
  # Actual maximum value (GPD density must be non-zero at this value)
  actual_boundary <- max(c(x, boundary))
  
  if (is.null(boundary)) {
    boundary <- eval_point <- actual_boundary
  }
  
  b <- w <- L <-
    1 / actual_boundary + (1 - sqrt(m / (1:m - 0.5))) / 3 / x[floor(n / 4 + 0.5)]
  
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
  
  densMax <- eva::dgpd(eval_point, scale = scale, shape = shape)
  
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
# gpdfit <- ZSE(x, constr = "none", boundary = NULL)
# gpdfit
#
# gpdfit <- ZSE(x, constr = "boundary", boundary = 8)
# gpdfit
#
# # -> fit is independent of boundary for positive shape
#
#
# # negative shape
# shape <- -0.25
#
# set.seed(123456)
# x <- rgpd(n = 100, scale = scale, shape = shape)
#
# gpdfit <- ZSE(x, constr = "none", boundary = NULL)
# gpdfit
#
# gpdfit <- ZSE(x, constr = "shapePos", boundary = NULL)
# gpdfit
# # shapePos constraint doesn't work
#
# gpdfit <- ZSE(x, constr = "boundary", boundary = 8)
# gpdfit
#
# -gpdfit$scale / gpdfit$shape
# # -> much larger than boundary



