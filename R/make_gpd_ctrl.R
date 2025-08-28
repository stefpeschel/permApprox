#' Create a control list for GPD approximation
#'
#' Constructs a list of control parameters for the Generalized Pareto
#' Distribution (GPD) tail approximation method used in permApprox.
#'
#' @param fit_method Character. Method for GPD fitting. Options: "LME", "MLE1D",
#'   "MLE2D", "MOM", "NLS2", "WNLLSM", "ZSE". Default: "LME".
#'
#' @param include_obs Logical. If \code{TRUE}, the observed test statistic is
#'   included into the tail approximation (considered as permutation test
#'   statistic). Default: \code{FALSE}.
#'
#' @param constraint Character. Constraint for the fitting process.
#'   Default: "support_at_max". Options:
#'   \describe{
#'     \item{\code{"unconstrained"}}{No constraint.}
#'     \item{\code{"shape_nonneg"}}{Shape parameter must be non-negative.}
#'     \item{\code{"support_at_obs"}}{Positive density for observed test statistic.}
#'     \item{\code{"support_at_max"}}{Positive density for maximum of all
#'     observed test statistics (in the multiple testing case).}
#'   }
#'
#' @param eps_rule Choice of rule for computing \eqn{\varepsilon}. 
#'   Allowed character values are:
#'   \itemize{
#'     \item \code{"constant"}: Use a fixed \eqn{\varepsilon}. 
#'       The constant is defined via \code{eps_par}.
#'     \item \code{"factor"}: Set \eqn{\varepsilon = c \cdot |t_{\mathrm{obs}}|} 
#'       with \code{c} provided via \code{eps_par}.
#'     \item \code{"slls"}: Use the Standardized Lifted Log-Saturation (SLLS) 
#'       rule on the \emph{Z}-scale with permutation-based cap; \code{eps_par} 
#'       controls its conservativeness via \code{target_factor}. The sample 
#'       size must be provided via \code{sample_size} for this rule.
#'   }
#'                 
#' @param eps_par Adaptive parameterization for \code{eps_rule}. Its meaning
#'   depends on \code{eps_rule}:
#'   \itemize{
#'     \item \code{"constant"}: numeric. If length 1, the same constant 
#'       \eqn{\varepsilon} is used for all tests; if length equals the number of 
#'       tests, values are applied per test.
#'     \item \code{"factor"}: numeric scalar \code{c} used as 
#'       \eqn{\varepsilon = c \cdot |t_{\mathrm{obs}}|}.
#'     \item \code{"slls"}: Numeric scalar that is used as \code{target_factor} 
#'       of the epsilon function (higher \code{target_factor}
#'           \textrightarrow{} more conservative plateau).
#'   }
#' 
#' @param sample_size Optional numeric giving the sample size. Needed for the 
#'   "slls" rule. For groups with different sample sizes use the smaller one. 
#'   If the sample size is unknown, use a small value (e.g., 30) to be 
#'   conservative.
#'
#' @param tol Numeric. Convergence tolerance for fitting GPD parameters.
#'   Default: 1e-8.
#'
#' @param thresh_method Character. Method for threshold detection.
#'   Default: "ftr_min5". Options:
#'   \describe{
#'     \item{\code{"fix"}}{Fix threshold (defined via \code{thresh0}).}
#'     \item{\code{"ftr"}}{Failure to reject.}
#'     \item{\code{"ftr_min5"}}{Failure to reject (min 5 subsequent acceptances).}
#'     \item{\code{"pr_below_alpha"}}{Proportion of rejections below alpha.}
#'     \item{\code{"fwd_stop"}}{Forward Stop}
#'     \item{\code{"gof_cp"}}{GOF change point}
#'   }
#'
#' @param thresh0 Numeric or NULL. Initial threshold for defining exceedances
#'   in the GPD fitting. Is used as fixed threshold when
#'   \code{thresh_method = "fix"}. Set to \code{NULL} if \code{exceed0} is
#'   provided instead. Default is \code{NULL}.
#'   Either \code{thresh0} or \code{exceed0} must be specified (non-NULL).
#'
#' @param thresh_step Integer. Step size for adaptive threshold search.
#'   If \code{thresh_step} is 2, for instance, every second possible threshold
#'   is considered to save runtime.
#'
#' @param exceed0 Numeric or NULL. Initial number of exceedances.
#'   If less than 1, it is interpreted as a proportion of the number of
#'   permutations. Default is \code{0.25}. Either \code{exceed0} or
#'   \code{thresh0} must be specified (non-NULL).
#'
#' @param exceed_min Integer. Minimum exceedances needed for fitting. Default: 10.
#'
#' @param gof_test Character. Goodness-of-fit test for GPD. Options:
#'   "ad" (Anderson-Darling), "cvm" (Cramer-von-Mises), "none". Default: "ad".
#'
#' @param gof_alpha Numeric. Significance level for GOF test (0 < gof_alpha < 1).
#'   Default: 0.05.
#'
#' @param cores Integer. Number of CPU cores to use for parallel computations.
#'   Must be >=1. Default: \code{1}.
#'
#' @param verbose Logical. If \code{TRUE}, progress messages are printed.
#'   Default: \code{TRUE}.
#'
#' @return A named list of class "gpd_ctrl" containing GPD settings.
#'
#' @export
#' 
make_gpd_ctrl <- function(
    fit_method = "LME",
    include_obs = FALSE,
    constraint = "support_at_max",
    eps_rule = "slls",
    eps_par = 0.25,
    zero_guard = TRUE,
    eps_retry = list(
      step_init       = 10,
      grow            = (1 + sqrt(5)) / 2,   # ~1.618
      max_expand_iter = 20L,
      bisect_iter_max = 30L,
      bisect_tol      = 1e-6    
    ),
    sample_size = NULL,
    tol = 1e-8,
    thresh_method = "ftr_min5",
    thresh0 = NULL,
    thresh_step = 10,
    exceed0 = 0.25,
    exceed_min = 10,
    gof_test = "ad",
    gof_alpha = 0.05,
    cores = 1,
    verbose = TRUE
) {
  
  fit_method <- match.arg(fit_method,
                          choices = c("LME",
                                      "MLE1D",
                                      "MLE2D",
                                      "MOM",
                                      "NLS2",
                                      "WNLLSM",
                                      "ZSE"))
  
  stopifnot(is.logical(include_obs))
  
  constraint <- match.arg(constraint,
                          choices = c("unconstrained",
                                      "shape_nonneg",
                                      "support_at_obs",
                                      "support_at_max"))
  
  if (constraint == "shape_nonneg" && !fit_method %in% c("MLE1D", "MLE2D", "NLS2")) {
    stop("Constraint \"shape_nonneg\" only available for methods ",
         "MLE1D, MLE2D, and NLS2.")
  }
  
  # --- Error handling for eps_rule / eps_par ---
  if (!eps_rule %in% c("constant", "factor", "slls")) {
    stop("Invalid 'eps_rule'. Must be one of: 'constant', 'factor', 'slls'.")
  }
  
  if (eps_rule == "constant") {
    if (!is.numeric(eps_par) || eps_par < 0) {
      stop("'eps_par' must be a numeric, non-negative scalar or vector for eps_rule = 'constant'.")
    }
  }
  
  if (eps_rule == "factor") {
    if (!is.numeric(eps_par) || length(eps_par) != 1 || eps_par < 0) {
      stop("'eps_par' must be a numeric, non-negative scalar for eps_rule = 'factor'.")
    }
  }
  
  if (eps_rule == "slls") {
    if (!is.numeric(eps_par) || length(eps_par) != 1 || eps_par < 0) {
      stop("'eps_par' must be a numeric, non-negative scalar (target_factor) for eps_rule = 'slls'.")
    }
    if (is.null(sample_size)) {
      stop("'sample_size' must be provided for eps_rule = 'slls'. ", 
           "For groups with different sample sizes use the smaller one. ")
    }
  }
  # ---
  
  thresh_method <- match.arg(thresh_method,
                             choices = c("fix",
                                         "ftr",
                                         "ftr_min5",
                                         "pr_below_alpha",
                                         "fwd_stop",
                                         "gof_cp"))
  
  # Error handling for thresh0 and exceed0
  if (is.null(thresh0) && is.null(exceed0)) {
    stop("You must specify either 'thresh0' or 'exceed0'. Both are currently NULL.")
  }
  if (!is.null(thresh0)) {
    if (!is.null(exceed0)) {
      stop("Specify either 'thresh0' or 'exceed0', not both.")
    }
    if (!is.numeric(thresh0) || length(thresh0) != 1 || is.na(thresh0)) {
      stop("'thresh0' must be a single numeric value or NULL.")
    }
  }
  if (!is.null(exceed0)) {
    if (!is.null(thresh0)) {
      stop("Specify either 'thresh0' or 'exceed0', not both.")
    }
    if (!is.numeric(exceed0) || length(exceed0) != 1 || is.na(exceed0)) {
      stop("'exceed0' must be a single numeric value or NULL.")
    }
    if (exceed0 < 0) {
      stop("'exceed0' must be a non-negative number (or NULL).")
    }
  }
  
  stopifnot(is.numeric(exceed_min) & exceed_min >= 0)
  
  gof_test <- match.arg(gof_test, choices = c("ad", "cvm", "none"))
  
  stopifnot(is.numeric(gof_alpha) & gof_alpha > 0 & gof_alpha < 1)
  
  control <- list(
    fit_method = fit_method,
    include_obs = include_obs,
    constraint = constraint,
    eps_rule = eps_rule,
    eps_par = eps_par,
    zero_guard = zero_guard,
    eps_retry = eps_retry,
    sample_size = sample_size,
    tol = tol,
    thresh_method = thresh_method,
    thresh0 = thresh0,
    thresh_step = thresh_step,
    exceed0 = exceed0,
    exceed_min = exceed_min,
    gof_test = gof_test,
    gof_alpha = gof_alpha,
    cores = cores,
    verbose = verbose
  )
  class(control) <- "gpd_ctrl"
  control
}
