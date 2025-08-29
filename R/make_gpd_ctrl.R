#' Create a control list for GPD approximation
#'
#' Constructs a list of control parameters for the Generalized Pareto
#' Distribution (GPD) tail approximation method used in permApprox.
#'
#' @section Constraints and epsilon:
#' Epsilon-related settings (\code{eps_rule}, \code{eps_par}, \code{sample_size},
#' and \code{eps_retry}/\code{zero_guard}) are only relevant when
#' \code{constraint != "unconstrained"}. If \code{constraint == "unconstrained"},
#' epsilon is not used and these arguments are ignored (\code{eps_rule} is set to
#' \code{"constant"} with \code{eps_par = 0}, and \code{zero_guard} is disabled).
#'
#' @param fit_method Character. Method for GPD fitting. Options: "LME", "MLE1D",
#'   "MLE2D", "MOM", "NLS2", "WNLLSM", "ZSE". Default: "LME".
#'
#' @param include_obs Logical scalar. If \code{TRUE}, the observed test statistic is
#'   included into the tail approximation (considered as permutation test statistic).
#'   Default: \code{FALSE}.
#'
#' @param constraint Character. Constraint for the fitting process. Default:
#'   \code{"support_at_obs"}. Options:
#'   \describe{
#'     \item{\code{"unconstrained"}}{No constraint (epsilon is not used).}
#'     \item{\code{"shape_nonneg"}}{Shape parameter must be non-negative.}
#'     \item{\code{"support_at_obs"}}{Positive density at the observed test statistic.}
#'     \item{\code{"support_at_max"}}{Positive density at the maximum observed
#'       test statistic (multiple testing case).}
#'   }
#'
#' @param eps_rule Character. Rule for computing \eqn{\varepsilon}. Used only if
#'   \code{constraint != "unconstrained"}. One of:
#'   \itemize{
#'     \item \code{"constant"}: fixed \eqn{\varepsilon}, set by \code{eps_par}.
#'     \item \code{"factor"}: \eqn{\varepsilon = c\,|t_{\mathrm{obs}}|} with \code{c = eps_par}.
#'     \item \code{"slls"}: Standardized Lifted Log-Saturation (Z-scale with
#'           permutation-based cap); \code{eps_par} acts as \emph{target_factor}.
#'   }
#'
#' @param eps_par Numeric. Parameter for \code{eps_rule}. If \code{eps_rule = "constant"},
#'   this may be a scalar (same \eqn{\varepsilon} for all tests) or a vector (per-test).
#'   For \code{"factor"} and \code{"slls"} this must be a non-negative scalar.
#'
#' @param zero_guard Logical. If \code{TRUE}, enable adaptive epsilon refinement to
#'   eliminate machine-underflow p-values during GPD fitting (guarding at
#'   \code{.Machine$double.xmin}). Only used if \code{constraint != "unconstrained"}.
#'
#' @param eps_retry List of tuning parameters for the adaptive epsilon refinement.
#'   Only used when \code{zero_guard = TRUE} and \code{constraint != "unconstrained"}.
#'   Recognized fields (defaults in parentheses):
#'   \describe{
#'     \item{\code{step_init}}{Initial additive step for the target factor (10).}
#'     \item{\code{grow}}{Geometric growth factor for the step (>1; default
#'                        \eqn{(1+\sqrt{5})/2 \approx 1.618}).}
#'     \item{\code{max_expand_iter}}{Max expansion iterations (20).}
#'     \item{\code{bisect_iter_max}}{Max bisection iterations (30).}
#'     \item{\code{bisect_tol}}{Bisection tolerance (0.1).}
#'   }
#'
#' @param sample_size Optional numeric. Required for \code{eps_rule = "slls"} when
#'   \code{constraint != "unconstrained"}. For groups with different sizes, use the
#'   smaller one.
#'
#' @param tol Numeric > 0. Convergence tolerance for fitting GPD parameters.
#'   Default: 1e-8.
#'
#' @param thresh_method Character. Method for threshold detection. Default:
#'   \code{"ftr_min5"}. Options: \code{"fix"}, \code{"ftr"}, \code{"ftr_min5"},
#'   \code{"pr_below_alpha"}, \code{"fwd_stop"}, \code{"gof_cp"}.
#'
#' @param thresh0 Numeric or NULL. Initial threshold (fixed when
#'   \code{thresh_method = "fix"}). Must be NULL if \code{exceed0} is provided.
#'
#' @param thresh_step Integer >= 1. Step size for adaptive threshold search.
#'
#' @param exceed0 Numeric or NULL. Initial number (or proportion < 1) of exceedances.
#'   Must be NULL if \code{thresh0} is provided. Default: 0.25.
#'
#' @param exceed_min Integer >= 0. Minimum exceedances required for fitting.
#'   Default: 10.
#'
#' @param gof_test Character. Goodness-of-fit test for GPD: \code{"ad"}, \code{"cvm"},
#'   or \code{"none"}. Default: \code{"ad"}.
#'
#' @param gof_alpha Numeric in (0,1). Significance level for GOF test. Default: 0.05.
#'
#' @return A named list of class \code{"gpd_ctrl"} containing GPD settings.
#'
#' @export
make_gpd_ctrl <- function(
    fit_method = "LME",
    include_obs = FALSE,
    constraint = "support_at_obs",
    eps_rule = "slls",
    eps_par = 0.25,
    zero_guard = TRUE,
    eps_retry = list(
      step_init       = 10,
      grow            = (1 + sqrt(5)) / 2,   # ~1.618
      max_expand_iter = 20L,
      bisect_iter_max = 30L,
      bisect_tol      = 0.1
    ),
    sample_size = NULL,
    tol = 1e-8,
    thresh_method = "ftr_min5",
    thresh0 = NULL,
    thresh_step = 10,
    exceed0 = 0.25,
    exceed_min = 10,
    gof_test = "ad",
    gof_alpha = 0.05
) {
  ## normalize / validate simple scalars
  fit_method <- match.arg(fit_method, 
                          c("LME","MLE1D","MLE2D","MOM","NLS2","WNLLSM","ZSE"))
  
  if (!is.logical(include_obs) || length(include_obs) != 1L || is.na(include_obs))
    stop("'include_obs' must be a single logical.")
  
  constraint <- match.arg(constraint, 
                          c("unconstrained", "shape_nonneg", "support_at_obs", 
                            "support_at_max"))
  
  if (!(is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0))
    stop("'tol' must be a single positive finite numeric.")
  
  if (!(is.numeric(thresh_step) && length(thresh_step) == 1L && 
        is.finite(thresh_step) && thresh_step >= 1))
    stop("'thresh_step' must be a single integer >= 1.")
  
  if (!(is.numeric(exceed_min) && length(exceed_min) == 1L && 
        is.finite(exceed_min) && exceed_min >= 0))
    stop("'exceed_min' must be a single integer >= 0.")
  
  ## constraint-specific check
  if (constraint == "shape_nonneg" && !fit_method %in% c("MLE1D", "MLE2D", "NLS2")) {
    stop("Constraint 'shape_nonneg' only available for methods MLE1D, MLE2D, and NLS2.")
  }
  
  ## threshold method
  thresh_method <- match.arg(thresh_method, 
                             c("fix", "ftr", "ftr_min5", "pr_below_alpha",
                               "fwd_stop", "gof_cp"))
  if (is.null(thresh0) && is.null(exceed0)) {
    stop("You must specify either 'thresh0' or 'exceed0'. Both are currently NULL.")
  }
  
  if (!is.null(thresh0)) {
    if (!is.null(exceed0)) stop("Specify either 'thresh0' or 'exceed0', not both.")
    if (!(is.numeric(thresh0) && length(thresh0) == 1L && is.finite(thresh0)))
      stop("'thresh0' must be a single numeric or NULL.")
  }
  
  if (!is.null(exceed0)) {
    if (!is.null(thresh0)) stop("Specify either 'thresh0' or 'exceed0', not both.")
    if (!(is.numeric(exceed0) && length(exceed0) == 1L && is.finite(exceed0)))
      stop("'exceed0' must be a single numeric or NULL.")
    if (exceed0 < 0) stop("'exceed0' must be non-negative (count) or proportion < 1.")
  }
  
  ## GOF
  gof_test <- match.arg(gof_test, c("ad","cvm","none"))
  if (!(is.numeric(gof_alpha) && length(gof_alpha) == 1L && 
        is.finite(gof_alpha) && gof_alpha > 0 && gof_alpha < 1))
    stop("'gof_alpha' must be a single numeric in (0,1).")
  
  ## epsilon controls (only if constrained)
  if (constraint == "unconstrained") {
    # epsilon not used in unconstrained mode
    eps_rule   <- "constant"
    eps_par    <- 0
    zero_guard <- FALSE
    
  } else {
    if (!eps_rule %in% c("constant","factor","slls"))
      stop("Invalid 'eps_rule'. Must be one of: 'constant', 'factor', 'slls'.")
    
    if (identical(eps_rule, "constant")) {
      if (!(is.numeric(eps_par) && all(is.finite(eps_par)) && 
            all(eps_par >= 0)))
        stop("For eps_rule = 'constant', ", 
             "'eps_par' must be non-negative numeric (scalar or vector).")
    }
    if (identical(eps_rule, "factor")) {
      if (!(is.numeric(eps_par) && length(eps_par) == 1L && 
            is.finite(eps_par) && eps_par >= 0))
        stop("For eps_rule = 'factor', ", 
             "'eps_par' must be a single non-negative numeric.")
    }
    if (identical(eps_rule, "slls")) {
      if (!(is.numeric(eps_par) && length(eps_par) == 1L && 
            is.finite(eps_par) && eps_par >= 0))
        stop("For eps_rule = 'slls', ", 
             "'eps_par' must be a single non-negative numeric (target_factor).")
      if (is.null(sample_size))
        stop("'sample_size' must be provided for eps_rule = 'slls' ", 
             "(use the smaller group size).")
      if (!(is.numeric(sample_size) && length(sample_size) == 1L && 
            is.finite(sample_size) && sample_size > 0))
        stop("'sample_size' must be a single positive numeric for eps_rule = 'slls'.")
    }
    
    ## merge/validate eps_retry if zero_guard is enabled
    eps_retry_defaults <- list(
      step_init       = 10,
      grow            = (1 + sqrt(5)) / 2,
      max_expand_iter = 20L,
      bisect_iter_max = 30L,
      bisect_tol      = 0.01
    )
    
    if (is.null(eps_retry)) eps_retry <- list()
    
    eps_retry <- utils::modifyList(eps_retry_defaults, eps_retry)
    
    if (isTRUE(zero_guard)) {
      if (!(is.numeric(eps_retry$step_init) && 
            length(eps_retry$step_init) == 1L &&
            is.finite(eps_retry$step_init) && eps_retry$step_init > 0))
        stop("'eps_retry$step_init' must be a single positive numeric.")
      if (!(is.numeric(eps_retry$grow) && length(eps_retry$grow) == 1L &&
            is.finite(eps_retry$grow) && eps_retry$grow > 1))
        stop("'eps_retry$grow' must be a single numeric > 1.")
      if (!(is.numeric(eps_retry$max_expand_iter) && 
            length(eps_retry$max_expand_iter) == 1L &&
            is.finite(eps_retry$max_expand_iter) && 
            eps_retry$max_expand_iter >= 0))
        stop("'eps_retry$max_expand_iter' must be a single integer >= 0.")
      eps_retry$max_expand_iter <- as.integer(eps_retry$max_expand_iter)
      if (!(is.numeric(eps_retry$bisect_iter_max) && 
            length(eps_retry$bisect_iter_max) == 1L &&
            is.finite(eps_retry$bisect_iter_max) && 
            eps_retry$bisect_iter_max >= 0))
        stop("'eps_retry$bisect_iter_max' must be a single integer >= 0.")
      eps_retry$bisect_iter_max <- as.integer(eps_retry$bisect_iter_max)
      if (!(is.numeric(eps_retry$bisect_tol) && 
            length(eps_retry$bisect_tol) == 1L &&
            is.finite(eps_retry$bisect_tol) && eps_retry$bisect_tol > 0))
        stop("'eps_retry$bisect_tol' must be a single positive numeric.")
    }
  }
  
  control <- list(
    fit_method    = fit_method,
    include_obs   = include_obs,
    constraint    = constraint,
    eps_rule      = eps_rule,
    eps_par       = eps_par,
    zero_guard    = zero_guard,
    eps_retry     = eps_retry,
    sample_size   = sample_size,
    tol           = tol,
    thresh_method = thresh_method,
    thresh0       = thresh0,
    thresh_step   = as.integer(thresh_step),
    exceed0       = exceed0,
    exceed_min    = as.integer(exceed_min),
    gof_test      = gof_test,
    gof_alpha     = gof_alpha,
    cores         = cores,
    verbose       = verbose
  )
  
  class(control) <- "gpd_ctrl"
  control
}

