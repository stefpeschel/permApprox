#' Create a control list for GPD approximation
#'
#' Constructs a list of control parameters for the Generalized Pareto
#' Distribution (GPD) tail approximation method used in \pkg{permApprox}.
#'
#' @section Constraints and epsilon:
#' Epsilon-related settings (\code{eps_fun}, \code{eps_tune}, \code{eps_args},
#' \code{sample_size}, and \code{eps_retry}/\code{zero_guard}) are only used
#' when a \emph{support constraint} is active, that is, for
#' \code{constraint = "support_at_obs"} or \code{"support_at_max"}.
#' In these cases, a small positive constant \eqn{\varepsilon} is added to the
#' evaluation point to ensure that it lies strictly within the GPD support.
#'
#' If \code{constraint} is \code{"unconstrained"} or \code{"shape_nonneg"},
#' epsilon is not used: \code{eps_fun} is ignored, \code{eps_tune} is set to
#' \code{0}, \code{eps_args} is cleared, and \code{zero_guard} is disabled.
#'
#' @param fit_method Character. Method for GPD fitting.
#'   Options: \code{"lme"}, \code{"mle1d"}, \code{"mle2d"}, \code{"mom"},
#'   \code{"nls2"}, \code{"wnllsm"}, \code{"zse"}. Default: \code{"lme"}.
#'
#' @param include_obs Logical scalar. If \code{TRUE}, the observed test
#'   statistic is included into the tail approximation (treated as an additional
#'   permutation statistic). Default: \code{FALSE}.
#'   
##' @param discrete_screen Logical scalar. If \code{TRUE} (default), performs a
#'   discreteness screening step before GPD fitting and falls back to empirical
#'   p-values for tests flagged as too discrete for reliable tail modeling.
#'   If \code{FALSE}, this screening step is skipped.
#'
#' @param constraint Character. Constraint for the fitting process.
#'   Default: \code{"support_at_obs"}. Options:
#'   \describe{
#'     \item{\code{"none"} or \code{"unconstrained"}}{No constraint.}
#'     \item{\code{"shape_nonneg"}}{Shape parameter must be non-negative.}
#'     \item{\code{"support_at_obs"}}{Positive density at the observed test
#'       statistic for each test.}
#'     \item{\code{"support_at_max"}}{Positive density at the maximum observed
#'       test statistic across tests (multiple testing case).}
#'   }
#'
#' @param eps_fun Epsilon function used when a support constraint is active.
#'   Default: \code{\link{eps_slls}}. The function must accept at least the
#'   arguments \code{obs_stats}, \code{perm_stats}, \code{sample_size},
#'   and \code{tune}, and must return a numeric vector of the same length as
#'   \code{obs_stats}. For user-defined epsilon rules, \code{sample_size} may
#'   be ignored or treated flexibly. If \code{constraint} is not a support
#'   constraint, \code{eps_fun} is ignored.
#'
#' @param eps_tune Numeric scalar. Global tuning parameter passed as
#'   \code{tune} to \code{eps_fun}. This is the scalar knob that may be adapted
#'   by the zero-guard mechanism to avoid underflowing p-values.
#'   Default: \code{0.25}.
#'
#' @param eps_args List of additional arguments passed on to \code{eps_fun}.
#'   For the default \code{\link{eps_slls}}, this can be used to change
#'   \code{cap_base}, \code{alpha}, \code{q_ref}, \code{rho_lift}, etc.
#'   Ignored if \code{constraint} is not a support constraint.
#'
#' @param zero_guard Logical. If \code{TRUE}, enable adaptive epsilon refinement
#'   (via \code{eps_tune}) to eliminate machine-underflow p-values during GPD
#'   fitting (guarding at \code{.Machine$double.xmin}). Only used when a support
#'   constraint is active.
#'
#' @param eps_retry List of tuning parameters for the adaptive epsilon
#'   refinement (zero-guard). Only used when \code{zero_guard = TRUE} and a
#'   support constraint is active. Recognized fields (defaults in parentheses):
#'   \describe{
#'     \item{\code{step_init}}{Initial additive step for the tuning parameter
#'       (\code{10}).}
#'     \item{\code{grow}}{Geometric growth factor for the step (> 1;
#'       default \eqn{(1+\sqrt{5})/2 \approx 1.618}).}
#'     \item{\code{max_expand_iter}}{Maximum expansion iterations (\code{20}).}
#'     \item{\code{bisect_iter_max}}{Maximum bisection iterations (\code{30}).}
#'     \item{\code{bisect_tol}}{Bisection tolerance for the tuning parameter
#'       (\code{0.1}).}
#'   }
#'
#' @param sample_size Optional numeric. Per-test sample size proxy passed to
#'   \code{eps_fun}. For the default \code{\link{eps_slls}} under a support
#'   constraint, this must be provided and strictly positive (for groups with
#'   unequal sizes, use the smaller one). For custom \code{eps_fun}, it may be
#'   left \code{NULL} or used as needed. Ignored if \code{constraint} is not a
#'   support constraint.
#'
#' @param tol Numeric > 0. Convergence tolerance for fitting GPD parameters.
#'   Default: \code{1e-7}.
#'
#' @param thresh_method Character. Method for threshold detection.
#'   Default: \code{"rob_ftr"}. Options: \code{"fix"}, \code{"ftr"},
#'   \code{"rob_ftr"}, \code{"rpc"}, \code{"fwd_stop"}, \code{"gof_cp"}.
#'
#' @param thresh0 Numeric or \code{NULL}. Initial threshold (fixed when
#'   \code{thresh_method = "fix"}). Must be \code{NULL} if \code{exceed0} is
#'   provided.
#'
#' @param thresh_step Integer >= 1. Step size for adaptive threshold search.
#'
#' @param exceed0 Numeric or \code{NULL}. Initial number (or proportion < 1) of
#'   exceedances. Must be \code{NULL} if \code{thresh0} is provided.
#'   Default: \code{0.25}.
#'
#' @param exceed0_min Integer or \code{NULL}. Minimum initial number of
#'   exceedances. The actual initial number is taken as
#'   \code{max(exceed0, exceed0_min)} after translating proportions to counts.
#'   If \code{NULL}, no minimum is enforced. Default: \code{250}.
#'
#' @param exceed_min Integer >= 0. Minimum number of exceedances required for
#'   fitting. Default: \code{10}.
#'   
#' @param gof_test_thresh Character. Goodness-of-fit test used for threshold 
#'   search. Possible values are \code{"ad"} (Anderson-Darling) and 
#'   \code{"cvm"} (Cramér-von-Mises).
#'
#' @param gof_test Character. Goodness-of-fit test for the final GPD parameters.
#'   Possible values: \code{"ad"} (Anderson-Darling), \code{"cvm"} 
#'   (Cramér-von-Mises), or \code{"none"} (no test; default).
#'
#' @param gof_alpha Numeric in (0, 1). Significance level for the GOF test.
#'   Default: \code{0.05}.
#'
#' @return A named list of class \code{"gpd_ctrl"} containing GPD settings.
#'
#' @export
make_gpd_ctrl <- function(
    fit_method      = "lme",
    include_obs     = FALSE,
    discrete_screen = TRUE,
    constraint      = "support_at_obs",
    eps_fun         = eps_slls,
    eps_tune        = 0.25,
    eps_args        = list(),
    zero_guard      = TRUE,
    eps_retry       = list(
      step_init       = 10,
      grow            = (1 + sqrt(5)) / 2,   # ~1.618
      max_expand_iter = 20L,
      bisect_iter_max = 30L,
      bisect_tol      = 0.1
    ),
    sample_size     = NULL,
    tol             = 1e-7,
    thresh_method   = "rob_ftr",
    thresh0         = NULL,
    thresh_step     = 10,
    exceed0         = 0.25,
    exceed0_min     = 250,
    exceed_min      = 10,
    gof_test_thresh = "ad",
    gof_test        = "none",
    gof_alpha       = 0.05
) {
  ## ------------------------------------------------------------------------
  ## Basic validation
  ## ------------------------------------------------------------------------
  fit_method <- match.arg(
    fit_method,
    c("lme", "mle1d", "mle2d", "mom", "nls2", "wnllsm", "zse")
  )
  
  if (!is.logical(include_obs) || length(include_obs) != 1L || is.na(include_obs))
    stop("'include_obs' must be a single logical.")
  
  if (!is.logical(discrete_screen) ||
      length(discrete_screen) != 1L ||
      is.na(discrete_screen)) {
    stop("'discrete_screen' must be a single logical.")
  }
  
  constraint <- match.arg(
    constraint,
    c("none", "unconstrained", "shape_nonneg",
      "support_at_obs", "support_at_max")
  )
  if (constraint == "none") {
    constraint <- "unconstrained"
  }
  
  if (!(is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0))
    stop("'tol' must be a single positive finite numeric.")
  
  if (!(is.numeric(thresh_step) && length(thresh_step) == 1L &&
        is.finite(thresh_step) && thresh_step >= 1))
    stop("'thresh_step' must be a single numeric >= 1.")
  
  if (!(is.numeric(exceed_min) && length(exceed_min) == 1L &&
        is.finite(exceed_min) && exceed_min >= 0))
    stop("'exceed_min' must be a single numeric >= 0.")
  
  ## constraint-specific check
  if (constraint == "shape_nonneg" && !fit_method %in% c("mle1d", "mle2d", "nls2")) {
    stop("Constraint 'shape_nonneg' only available for methods ",
         "mle1d, mle2d, and NLS2.")
  }
  
  ## ------------------------------------------------------------------------
  ## Threshold method + initial thresholds
  ## ------------------------------------------------------------------------
  thresh_method <- match.arg(
    thresh_method,
    c("fix", "ftr", "rob_ftr", "rpc", "fwd_stop", "gof_cp")
  )
  
  if (is.null(thresh0) && is.null(exceed0)) {
    stop("You must specify either 'thresh0' or 'exceed0'. Both are currently NULL.")
  }
  
  if (!is.null(thresh0)) {
    if (!is.null(exceed0))
      stop("Specify either 'thresh0' or 'exceed0', not both.")
    if (!(is.numeric(thresh0) && length(thresh0) == 1L && is.finite(thresh0)))
      stop("'thresh0' must be a single numeric or NULL.")
  }
  
  if (!is.null(exceed0)) {
    if (!is.null(thresh0))
      stop("Specify either 'thresh0' or 'exceed0', not both.")
    if (!(is.numeric(exceed0) && length(exceed0) == 1L && is.finite(exceed0)))
      stop("'exceed0' must be a single numeric or NULL.")
    if (exceed0 < 0)
      stop("'exceed0' must be non-negative (count) or a proportion < 1.")
  }
  
  ## exceed0_min: may be NULL or a single non-negative numeric
  if (!is.null(exceed0_min)) {
    if (!(is.numeric(exceed0_min) &&
          length(exceed0_min) == 1L &&
          is.finite(exceed0_min) &&
          exceed0_min >= 0))
      stop("'exceed0_min' must be NULL or a single non-negative numeric.")
    
    exceed0_min <- as.integer(exceed0_min)
  }
  
  ## ------------------------------------------------------------------------
  ## GOF
  ## ------------------------------------------------------------------------
  gof_test <- match.arg(gof_test, c("ad", "cvm", "none"))
  gof_test_thresh <- match.arg(gof_test_thresh, c("ad", "cvm"))
  
  if (!(is.numeric(gof_alpha) && length(gof_alpha) == 1L &&
        is.finite(gof_alpha) && gof_alpha > 0 && gof_alpha < 1))
    stop("'gof_alpha' must be a single numeric in (0,1).")
  
  ## ------------------------------------------------------------------------
  ## Epsilon controls: function, tuning, args
  ## ------------------------------------------------------------------------
  if (!is.null(eps_fun) && !is.function(eps_fun))
    stop("'eps_fun' must be NULL or a function.")
  
  if (is.null(eps_args)) {
    eps_args <- list()
  } else if (!is.list(eps_args)) {
    stop("'eps_args' must be a list or NULL.")
  }
  
  if (!(is.numeric(eps_tune) && length(eps_tune) == 1L && is.finite(eps_tune)))
    stop("'eps_tune' must be a single finite numeric.")
  
  support_constraints <- c("support_at_obs", "support_at_max")
  
  if (!constraint %in% support_constraints) {
    ## No support constraint → epsilon not used
    eps_fun    <- NULL
    eps_tune   <- 0
    eps_args   <- list()
    zero_guard <- FALSE
    
  } else {
    ## Support constraints → epsilon may be used
    if (is.null(eps_fun)) {
      stop("For constraints 'support_at_obs' or 'support_at_max', ",
           "'eps_fun' must be provided (e.g., eps_slls).")
    }
    
    ## For default eps_slls we require a positive sample_size, so users see it early.
    if (identical(eps_fun, eps_slls)) {
      if (is.null(sample_size)) {
        stop("For the default eps_slls under support constraints, ",
             "'sample_size' must be provided (use the smaller group size).")
      }
      if (!(is.numeric(sample_size) && length(sample_size) == 1L &&
            is.finite(sample_size) && sample_size > 0)) {
        stop("'sample_size' must be a single positive numeric for eps_slls.")
      }
    } else {
      ## Custom eps_fun: sample_size is optional; if provided, check sanity.
      if (!is.null(sample_size)) {
        if (!(is.numeric(sample_size) && length(sample_size) == 1L &&
              is.finite(sample_size) && sample_size > 0))
          stop("'sample_size' must be a single positive numeric if supplied.")
      }
    }
    
    ## merge/validate eps_retry if zero_guard is enabled
    eps_retry_defaults <- list(
      step_init       = 10,
      grow            = (1 + sqrt(5)) / 2,
      max_expand_iter = 20L,
      bisect_iter_max = 30L,
      bisect_tol      = 0.1
    )
    
    if (is.null(eps_retry)) eps_retry <- list()
    eps_retry <- utils::modifyList(eps_retry_defaults, eps_retry)
    
    if (isTRUE(zero_guard)) {
      if (!(is.numeric(eps_retry$step_init) &&
            length(eps_retry$step_init) == 1L &&
            is.finite(eps_retry$step_init) && eps_retry$step_init > 0))
        stop("'eps_retry$step_init' must be a single positive numeric.")
      if (!(is.numeric(eps_retry$grow) &&
            length(eps_retry$grow) == 1L &&
            is.finite(eps_retry$grow) && eps_retry$grow > 1))
        stop("'eps_retry$grow' must be a single numeric > 1.")
      if (!(is.numeric(eps_retry$max_expand_iter) &&
            length(eps_retry$max_expand_iter) == 1L &&
            is.finite(eps_retry$max_expand_iter) &&
            eps_retry$max_expand_iter >= 0))
        stop("'eps_retry$max_expand_iter' must be a single numeric >= 0.")
      eps_retry$max_expand_iter <- as.integer(eps_retry$max_expand_iter)
      if (!(is.numeric(eps_retry$bisect_iter_max) &&
            length(eps_retry$bisect_iter_max) == 1L &&
            is.finite(eps_retry$bisect_iter_max) &&
            eps_retry$bisect_iter_max >= 0))
        stop("'eps_retry$bisect_iter_max' must be a single numeric >= 0.")
      eps_retry$bisect_iter_max <- as.integer(eps_retry$bisect_iter_max)
      if (!(is.numeric(eps_retry$bisect_tol) &&
            length(eps_retry$bisect_tol) == 1L &&
            is.finite(eps_retry$bisect_tol) && eps_retry$bisect_tol > 0))
        stop("'eps_retry$bisect_tol' must be a single positive numeric.")
    }
  }
  
  ## ------------------------------------------------------------------------
  ## Assemble control object
  ## ------------------------------------------------------------------------
  control <- list(
    fit_method       = fit_method,
    include_obs      = include_obs,
    discrete_screen  = discrete_screen,
    constraint       = constraint,
    eps_fun          = eps_fun,
    eps_tune         = eps_tune,
    eps_args         = eps_args,
    zero_guard       = zero_guard,
    eps_retry        = eps_retry,
    sample_size      = sample_size,
    tol              = tol,
    thresh_method    = thresh_method,
    thresh0          = thresh0,
    thresh_step      = as.integer(thresh_step),
    exceed0          = exceed0,
    exceed0_min      = exceed0_min,
    exceed_min       = as.integer(exceed_min),
    gof_test_thresh  = gof_test_thresh,
    gof_test         = gof_test,
    gof_alpha        = gof_alpha
  )
  
  class(control) <- "gpd_ctrl"
  control
}
