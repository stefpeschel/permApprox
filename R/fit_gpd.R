#' @title Fit and test the Generalized Pareto Distribution (GPD)
#'
#' @description
#' Fits the Generalized Pareto Distribution (GPD) to exceedances over a given threshold
#' using one of several available estimation methods. Optionally applies a goodness-of-fit test
#' to assess the adequacy of the fitted model.
#'
#' @param data Numeric vector. Input data to which the GPD is fitted.
#'
#' @param thresh Numeric scalar. Threshold above which the data are considered exceedances.
#'
#' @param fit_method Character string specifying the fitting method to be used.
#'   Must be one of \code{"LME"}, \code{"MLE1D"}, \code{"MLE2D"}, \code{"MOM"},
#'   \code{"NLS2"}, \code{"WNLLSM"}, or \code{"ZSE"}.
#'
#' @param tol Numeric. Convergence tolerance for iterative fitting procedures.
#'   Default is \code{1e-8}.
#'
#' @param eps Numeric. A small value used in support-related constraints (e.g., added to
#'   the upper support limit). Default is \code{0.8}.
#'
#' @param eps_type Character. Type of epsilon adjustment to apply. Can be
#'   \code{"quantile"} or \code{"fix"}.
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
#' @param ... Additional arguments passed to internal fitting functions.
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
#'   }
#'
#' @export


fit_gpd <- function(data,
                    thresh = NULL,
                    fit_method = "MLE1D",
                    tol = 1e-8,
                    eps = 0.8,
                    eps_type = "quantile",
                    constraint = "none",
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

  stopifnot(is.numeric(eps))
  stopifnot(eps_type %in% c("quantile", "fix"))

  gof_test <- match.arg(gof_test, choices = c("ad", "cvm", "none"))

  #-----------------------------------------------------------------------------
  # Fit GPD to the data
  #-----------------------------------------------------------------------------

  # Exceedances (test statistics above threshold)
  exceedances <- data[data > thresh]
  excess <- exceedances - thresh

  # Consider maximum value at which the GPD density must be positive
  if (!is.null(support_boundary)) {
    # Point at which the support is checked/evaluated
    eval_point <- support_boundary - thresh
    excess_boundary <- eval_point

    if (eps_type == "quantile") {
      eps <- quantile(data, eps)
    }

    excess_boundary <- eval_point + eps

  } else {
    eval_point <- excess_boundary <- NULL
  }

  # Contraint of a positive shape parameter
  if (constraint == "shape_nonneg") {
    if (!fit_method %in% c("MLE1D", "MLE2D", "NLS2")) {
      stop("Constraint \"shape_nonneg\" only available for methods ",
           "MLE1D, MLE2D, and NLS2.")
    }

    shapeMin <- 0
    shapePos <- TRUE
  } else {
    shapeMin <- -Inf
    shapePos <- FALSE
  }

  # GPD fit
  if (fit_method == "MLE2D") {
    fit <- try(.fit_gpd_mle2d(x = excess, boundary = excess_boundary,
                              eval_point = eval_point,
                             tol = tol, shapeMin = shapeMin),
                  silent = TRUE)

  } else if (fit_method == "MLE1D") {
    fit <- try(.fit_gpd_mle1d(x = excess, boundary = excess_boundary,
                              eval_point = eval_point,
                             tol = tol, shapePos = shapePos),
                  silent = TRUE)

  } else if (fit_method == "LME") {
    fit <- try(.fit_gpd_lme(x = excess, boundary = excess_boundary,
                            eval_point = eval_point,
                           tol = tol),
                  silent = TRUE)

  } else if (fit_method == "NLS2") {
    fit <- try(.fit_gpd_nls2(x = excess, boundary = excess_boundary,
                             eval_point = eval_point,
                            tol = tol, shapeMin = shapeMin),
                  silent = TRUE)

  } else if (fit_method == "WNLLSM") {
    fit <- try(.fit_gpd_wnllsm(x = excess, boundary = excess_boundary,
                               eval_point = eval_point,
                              tol = tol),
                  silent = TRUE)

  } else if (fit_method == "ZSE") {
    fit <- try(.fit_gpd_zse(x = excess, boundary = excess_boundary,
                            eval_point = eval_point),
                  silent = TRUE)

  } else if (fit_method == "MOM") {
    fit <- try(.fit_gpd_mom(x = excess),
                  silent = TRUE)
  }

  if ("try-error" %in% class(fit)) {
    shape <- scale <- negLogLik <- NA
    p_value <- 0

  } else {
    shape <- fit$shape
    scale <- fit$scale
    names(shape) <- names(scale) <- NULL

    #negLogLik <- fit$negLogLik

    if (gof_test == "none") {
      p_value <- NULL

    } else {
      if (gof_test == "ad") {
        testres <- try(.gof_gpd_ad(excess,
                                   scale = scale, shape = shape), silent = TRUE)
      } else if (gof_test == "cvm") {
        testres <- try(.gof_gpd_cvm(excess,
                                    scale = scale, shape = shape), silent = TRUE)
      }

      p_value <- ifelse("try-error" %in% class(testres), 0, testres$p.value)
    }

  }

  return(list(shape = shape, scale = scale, p_value = p_value))#, negLogLik = negLogLik))
}

