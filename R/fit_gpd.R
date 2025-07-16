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
#'   Must be one of \code{"LME"}, \code{"MLE1D"}, \code{"MLE2D"}, \code{"MOM"},
#'   \code{"NLS2"}, \code{"WNLLSM"}, or \code{"ZSE"}.
#'
#' @param tol Numeric. Convergence tolerance (default 1e-8).
#'
#' @param eps_fun Function that returns the ε to use. Possible options are 
#'   \code{\link{eps_fixed}}, \code{\link{eps_fixed}}, \code{\link{eps_fixed}},
#'   or a user-defined function. It is called as
#'        `eps_fun(data = data, support_boundary = support_boundary,
#'                 thresh = thresh, !!!eps_par)`.
#'                 
#' @param eps_par List of additional named arguments forwarded to `eps_fun`.
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
                    fit_method = "MLE1D",
                    tol = 1e-8,
                    eps_fun = eps_power,
                    eps_par = list(),
                    constraint = "unconstrained",
                    support_boundary = NULL,
                    gof_test = "ad",
                    ...) {
  stopifnot(is.vector(data) & is.numeric(data))
  stopifnot(is.numeric(thresh) & length(thresh) == 1)

  fit_method <- match.arg(fit_method,
                          choices = c("LME",
                                      "MLE1D",
                                      "MLE2D",
                                      "MOM",
                                      "NLS2",
                                      "WNLLSM",
                                      "ZSE"))

  constraint <- match.arg(constraint,
                          choices = c("unconstrained",
                                      "shape_nonneg",
                                      "support_at_obs",
                                      "support_at_max"))

  gof_test <- match.arg(gof_test, choices = c("ad", "cvm", "none"))
  
  stopifnot(is.function(eps_fun))
  
  #-----------------------------------------------------------------------------
  # Compute epsilon
  #-----------------------------------------------------------------------------
  
  eps <- do.call(
    eps_fun,
    c(list(data = data,
           support_boundary = support_boundary,
           thresh = thresh),
      eps_par)
  )
  if (!is.numeric(eps) || length(eps) != 1L || eps < 0)
    stop("`eps_fun` must return a single non-negative numeric value.")

  #-----------------------------------------------------------------------------
  # Exceedances
  #-----------------------------------------------------------------------------

  exceedances <- data[data > thresh]
  excess      <- exceedances - thresh
  
  # Evaluation point / support boundary shift
  eval_point      <- excess_boundary <- NULL
  if (!is.null(support_boundary)) {
    eval_point      <- support_boundary - thresh
    excess_boundary <- eval_point + eps
  }

  #-----------------------------------------------------------------------------
  # Shape-parameter constraints
  #-----------------------------------------------------------------------------
  
  if (constraint == "shape_nonneg" &&
      !fit_method %in% c("MLE1D","MLE2D","NLS2"))
    stop("Constraint 'shape_nonneg' only available for MLE1D, MLE2D, NLS2.")
  
  shapeMin <- if (constraint == "shape_nonneg") 0 else -Inf
  shapePos <- identical(shapeMin, 0)
  
  #-----------------------------------------------------------------------------
  # Fit switch
  #-----------------------------------------------------------------------------
  
  fit <- switch(fit_method,
                MLE2D  = try(.fit_gpd_mle2d(excess, excess_boundary, eval_point,
                                            tol = tol, shapeMin = shapeMin), 
                             silent = TRUE),
                MLE1D  = try(.fit_gpd_mle1d(excess, excess_boundary, eval_point,
                                            tol = tol, shapePos = shapePos), 
                             silent = TRUE),
                LME    = try(.fit_gpd_lme  (excess, excess_boundary, eval_point,
                                            tol = tol), 
                             silent = TRUE),
                NLS2   = try(.fit_gpd_nls2 (excess, excess_boundary, eval_point,
                                            tol = tol, shapeMin = shapeMin), 
                             silent = TRUE),
                WNLLSM = try(.fit_gpd_wnllsm(excess, excess_boundary, eval_point,
                                             tol = tol), 
                             silent = TRUE),
                ZSE    = try(.fit_gpd_zse  (excess, excess_boundary, eval_point),
                             silent = TRUE),
                MOM    = try(.fit_gpd_mom  (excess), 
                             silent = TRUE)
  )
  
  #-----------------------------------------------------------------------------
  # Collect results
  #-----------------------------------------------------------------------------
  
  if (inherits(fit, "try-error")) {
    list(shape = NA, scale = NA, p_value = 0, epsilon = eps)
  } else {
    shape <- fit$shape; scale <- fit$scale
    
    p_val <- if (gof_test == "none") {
      NULL
    } else {
      tfun <- if (gof_test == "ad") .gof_gpd_ad else .gof_gpd_cvm
      tmp  <- try(tfun(excess, scale = scale, shape = shape), silent = TRUE)
      if (inherits(tmp, "try-error")) 0 else tmp$p.value
    }
    
    list(shape = shape, scale = scale, p_value = p_val, epsilon = eps)
  }
}

