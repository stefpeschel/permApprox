#' @title Fit and test the Generalized Pareto Distribution (GPD)
#'
#' @description
#' Fits the Generalized Pareto Distribution (GPD) to exceedances over a given threshold
#' using one of several available estimation methods. Optionally applies a goodness-of-fit test
#' to assess the adequacy of the fitted model.
#'
#' @param data Numeric vector of observations.
#'
#' @param thresh Numeric scalar; exceedance threshold.
#'
#' @param fit_method Character string specifying the fitting method to be used.
#'   Must be one of \code{"lme"}, \code{"mle1d"}, \code{"mle2d"}, \code{"mom"},
#'   \code{"nls2"}, \code{"wnllsm"}, or \code{"zse"}.
#'
#' @param tol Numeric. Convergence tolerance (default 1e-8).
#'
#' @param constraint Character. Type of constraint to enforce during GPD fitting.
#'   Options are \code{"unconstrained"}, \code{"shape_nonneg"},
#'   \code{"support_at_obs"}, and \code{"support_at_max"}.
#'
#' @param support_boundary Numeric or NULL. Upper boundary of the GPD support,
#'   used in constrained fitting when the shape parameter is negative. If NULL,
#'   the maximum of the data vector is used.
#'
#' @param gof_test Character. Specifies the goodness-of-fit test to apply.
#'   Options are \code{"ad"} (Anderson-Darling), \code{"cvm"} (Cramér-von Mises),
#'   or \code{"none"}. Default is \code{"ad"}.
#'
#' @param ... Further arguments passed to low-level fitting helpers.
#'
#' @details
#' If the shape parameter is negative (indicating a bounded tail), a constraint
#' may be imposed to ensure that the support of the GPD does not exceed a specified
#' maximum value. The \code{support_boundary} argument defines this bound and is
#' only relevant for such constrained fitting. The default behavior uses the
#' maximum observed data value unless otherwise specified.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{shape}}{Estimated shape parameter of the GPD.}
#'     \item{\code{scale}}{Estimated scale parameter of the GPD.}
#'     \item{\code{p_value}}{P-value from the goodness-of-fit test, if applicable.}
#'     \item{\code{epsilon}}{Epsilon used for each test.}
#'   }
#'
#' @export


fit_gpd <- function(data,
                    thresh = NULL,
                    fit_method = "mle1d",
                    tol = 1e-8,
                    epsilon = 0,
                    constraint = "unconstrained",
                    support_boundary = NA_real_,
                    gof_test = "ad",
                    ...) {
  
  stopifnot(is.vector(data) & is.numeric(data))
  stopifnot(is.numeric(thresh) & length(thresh) == 1)
  
  fit_method <- match.arg(fit_method,
                          choices = c("lme",
                                      "mle1d",
                                      "mle2d",
                                      "mom",
                                      "nls2",
                                      "pwm",
                                      "wnllsm",
                                      "zse"))
  
  constraint <- match.arg(constraint,
                          choices = c("unconstrained",
                                      "shape_nonneg",
                                      "support_at_obs",
                                      "support_at_max"))
  
  gof_test <- match.arg(gof_test, choices = c("ad", "cvm", "none"))
  
  #-----------------------------------------------------------------------------
  # Exceedances
  #-----------------------------------------------------------------------------
  
  exceedances <- data[data > thresh]
  excess      <- exceedances - thresh
  
  # Evaluation point / support boundary shift
  eval_point      <- excess_boundary <- NULL
  if (!is.na(support_boundary)) {
    eval_point      <- support_boundary - thresh
    excess_boundary <- eval_point + epsilon
  }
  
  #-----------------------------------------------------------------------------
  # Shape-parameter constraints
  #-----------------------------------------------------------------------------
  
  if (constraint == "shape_nonneg" &&
      !fit_method %in% c("mle1d","mle2d","nls2"))
    stop("Constraint 'shape_nonneg' only available for mle1d, mle2d, nls2.")
  
  shapeMin <- if (constraint == "shape_nonneg") 0 else -Inf
  shapePos <- identical(shapeMin, 0)
  
  #-----------------------------------------------------------------------------
  # Fit switch
  #-----------------------------------------------------------------------------
  
  fit <- switch(fit_method,
                mle2d  = try(.fit_gpd_mle2d(excess, boundary = excess_boundary, 
                                            eval_point = eval_point,
                                            tol = tol, shapeMin = shapeMin), 
                             silent = TRUE),
                mle1d  = try(.fit_gpd_mle1d(excess, boundary = excess_boundary, 
                                            eval_point = eval_point,
                                            tol = tol, shapePos = shapePos), 
                             silent = TRUE),
                lme    = try(.fit_gpd_lme  (excess, boundary = excess_boundary, 
                                            eval_point = eval_point,
                                            tol = tol), 
                             silent = TRUE),
                nls2   = try(.fit_gpd_nls2 (excess, boundary = excess_boundary, 
                                            eval_point = eval_point,
                                            tol = tol, shapeMin = shapeMin), 
                             silent = TRUE),
                wnllsm = try(.fit_gpd_wnllsm(excess, boundary = excess_boundary, 
                                             eval_point = eval_point,
                                             tol = tol), 
                             silent = TRUE),
                zse    = try(.fit_gpd_zse  (excess, boundary = excess_boundary, 
                                            eval_point = eval_point),
                             silent = TRUE),
                mom    = try(.fit_gpd_mom  (excess), 
                             silent = TRUE),
                pwm    = try(.fit_gpd_pwm  (excess), 
                             silent = TRUE)
  )
  
  #-----------------------------------------------------------------------------
  # Collect results
  #-----------------------------------------------------------------------------
  
  if (inherits(fit, "try-error")) {
    list(shape = NA, scale = NA, p_value = 0, epsilon = epsilon)
  } else {
    shape <- fit$shape; scale <- fit$scale
    
    p_val <- if (gof_test == "none") {
      NULL
    } else {
      tfun <- if (gof_test == "ad") .gof_gpd_ad else .gof_gpd_cvm
      tmp  <- try(tfun(excess, scale = scale, shape = shape), silent = TRUE)
      if (inherits(tmp, "try-error")) 0 else tmp$p.value
    }
    
    list(shape = shape, scale = scale, p_value = p_val, epsilon = epsilon)
  }
}

