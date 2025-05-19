#' @title Compute empirical and approximated p-values from permutation tests
#'
#' @description
#' Computes empirical p-values for permutation tests.
#' When p-values are small, a Gamma or Generalized Pareto Distribution
#' (GPD) is fitted to the (tail of the) permutation distribution to improve
#' resolution if the number of permutations is small.
#'
#' @param obs_stats Numeric vector of observed test statistic(s).
#'
#' @param perm_stats Numeric vector or matrix of permutation test statistics.
#'   For a single test, supply a vector; for multiple tests, supply a matrix
#'   with one row per test.
#'
#' @param method Character. Method used to compute p-values. Default is \code{"gpd"}.
#'   Options are:
#'   \itemize{
#'     \item \code{"empirical"}: Empirical p-values are returned directly,
#'     based on the number of permutation test statistics that are as or
#'     more extreme than the observed statistic.
#'     \item \code{"gamma"}: A Gamma distribution is fitted to the permutation
#'       distribution if the empirical p-value falls below \code{fit_thresh}.
#'     \item \code{"gpd"}: A Generalized Pareto Distribution (GPD) is fitted
#'       to the tail of the permutation distribution if the empirical p-value
#'       falls below \code{fit_thresh}.
#'   }
#'
#' @param fit_thresh Numeric. Threshold on empirical p-values below which
#'   parametric fitting is applied. Default: \code{0.2} (parametric
#'   approximation if the empirical p-value is smaller than 0.2).
#'
#' @param alternative Character. One of \code{"greater"}, \code{"less"}, or
#'   \code{"two_sided"} (default), indicating the tail of the test.
#'
#' @param null_center Numeric or character. Specifies the value around which the
#'   null distribution is centered. If set to \code{"mean"} or \code{"median"},
#'   the per-row mean or median of \code{perm_stats} is used instead.
#'   This allows testing against a null hypothesis other than zero or centering
#'   based on the empirical distribution.

#' @param gpd_ctrl A control object created by \code{\link{make_gpd_ctrl}}.
#'   Contains settings for the GPD approximation, such as the fitting method,
#'   constraints, and thresholding strategy. Defaults to \code{make_gpd_ctrl()}.
#'
#' @param gamma_ctrl A control object created by \code{\link{make_gamma_ctrl}}.
#'   Contains settings for the Gamma approximation, including goodness-of-fit
#'   test and inclusion of the observed statistic. Defaults to
#'   \code{make_gamma_ctrl()}.
#'
#' @param adjust_ctrl A control object created by \code{\link{make_adjust_ctrl}}.
#'   Contains settings for multiple testing correction, such as the adjustment
#'   method and estimation of the proportion of true null hypotheses.
#'   Defaults to \code{make_adjust_ctrl()}.
#'
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing:
#' \describe{
#'   \item{p_values}{Numeric vector of approximated p-values. These are either
#'   unadjusted or adjusted for multiple testing, depending on the chosen
#'   adjustment method in the \code{control$adjust} settings.}
#'   \item{p_unadjusted}{Numeric vector of approximated but unadjusted p-values.}
#'   \item{p_empirical}{Numeric vector of raw empirical p-values.}
#'   \item{gpd_fit}{Details of GPD fit (if used), or \code{NULL}.}
#'   \item{gamma_fit}{Details of Gamma fit (if used), or \code{NULL}.}
#'   \item{method_used}{Character. Method used per test.}
#'   \item{adjust_result}{Full output of \code{mult_adjust()} if adjustment used.}
#'   \item{control}{List with the used control arguments.}
#' }
#'
#' @details
#' The function always computes empirical p-values for each test statistic by
#' comparing the observed statistic to the permutation distribution.
#' Specifically, the empirical p-value is calculated as
#' \deqn{p = \frac{r + 1}{B + 1}}
#' where \eqn{r} is the number of permutation statistics that are as or more
#' extreme than the observed statistic (according to the specified alternative),
#' and \eqn{B} is the number of permutations. This correction ensures that no
#' p-value is exactly zero.
#'
#' If \code{method} is set to \code{"gpd"} or \code{"gamma"}, a
#' parametric distribution is fitted to (the tail of) the permutation
#' distribution for cases where the empirical p-value falls below the threshold
#' specified by \code{fit_thresh}.
#' In such cases, the returned \code{p_values} vector contains the approximated
#' values (empirical for larger p-values, approximated for smaller ones).
#'
#' If \code{method = "none"}, only empirical p-values are returned.
#'
#' @examples
#' # Generate observed and permuted test statistics
#' set.seed(12345)
#' obs <- c(2.0, 3.0, 4.0, 5.0)
#' perm <- matrix(rnorm(4000), nrow = 4)
#'
#' # Empirical p-values
#' res_emp <- compute_p_values(obs_stats = obs,
#'                             perm_stats = perm,
#'                             method = "empirical")
#'
#' # Gamma approximation
#' gamma_ctrl <- make_gamma_ctrl(gof_test = "none")
#' res_gamma <- compute_p_values(obs_stats = obs,
#'                               perm_stats = perm,
#'                               method = "gamma",
#'                               gamma_ctrl = gamma_ctrl)
#'
#' # GPD approximation
#' res_gpd <- compute_p_values(obs_stats = obs,
#'                             perm_stats = perm,
#'                             method = "gpd")
#'
#' # GPD approximation with constraint
#' gpd_ctrl <- make_gpd_ctrl(constraint = "support_at_obs")
#'
#' res_gpd_constr <- compute_p_values(obs_stats = obs,
#'                                    perm_stats = perm,
#'                                    method = "gpd",
#'                                    gpd_ctrl = gpd_ctrl)
#'
#' # Data frame with (unadjusted) p-values
#' p_values <- data.frame(empirical = res_emp$p_unadjusted,
#'                        gamma = res_gamma$p_unadjusted,
#'                        gpd = res_gpd$p_unadjusted,
#'                        gpd_constr = res_gpd_constr$p_unadjusted)
#'
#' p_values
#'
#' @export

compute_p_values <- function(
    obs_stats,
    perm_stats,
    method = "gpd",
    fit_thresh = 0.2,
    alternative = "two_sided",
    null_center = 0,
    gpd_ctrl = make_gpd_ctrl(),
    gamma_ctrl = make_gamma_ctrl(),
    adjust_ctrl = make_adjust_ctrl(),
    ...
) {

  alternative <- match.arg(alternative,
                           choices = c("greater", "less", "two_sided"))

  method <- match.arg(method,
                      choices = c("gpd", "gamma", "empirical"))

  # Validate control arguments
  if (!inherits(gpd_ctrl, "gpd_ctrl")) {
    stop("'gpd_ctrl' must be created with make_gpd_ctrl().")
  }

  if (!inherits(gamma_ctrl, "gamma_ctrl")) {
    stop("'gamma_ctrl' must be created with make_gamma_ctrl().")
  }

  if (!inherits(adjust_ctrl, "adjust_ctrl")) {
    stop("'adjust_ctrl' must be created with make_adjust_ctrl().")
  }

  # Number of permutations and tests
  n_perm <- ncol(perm_stats)
  n_test <- length(obs_stats)

  # Ensure matrix format for multiple tests
  if (length(obs_stats) == 1) {
    perm_stats <- matrix(perm_stats, nrow = 1)
  }

  # Determine centering vector
  if (is.numeric(null_center)) {
    center_vec <- null_center
  } else {
    center_vec <- switch(
      null_center,
      mean = if (n_test == 1) mean(perm_stats) else rowMeans(perm_stats),
      median = if (n_test == 1) median(perm_stats) else apply(perm_stats, 1, median),
      stop("Invalid 'null_center': must be numeric, 'mean', or 'median'")
    )
  }

  # Center permutation statistics
  if (n_test == 1) {
    perm_stats <- perm_stats - center_vec
  } else {
    perm_stats <- sweep(perm_stats, 1, center_vec, FUN = "-")
  }

  # Center observed statistic(s)
  obs_stats <- obs_stats - center_vec

  # Empirical p-values
  pvals_emp_list <- .compute_pvals_emp(obs_stats = obs_stats,
                                       perm_stats = perm_stats,
                                       n_test = n_test,
                                       n_perm = n_perm,
                                       alternative = alternative)

  p_empirical <- pvals_emp_list$pvals
  n_perm_exceeding <- pvals_emp_list$n_perm_exceeding

  # Initialize gamme_fit and gpd_fit
  gamma_fit <- gpd_fit <- NULL

  #-----------------------------------------------------------------------------
  # Compute p-values

  if (method == "empirical" | all(p_empirical > fit_thresh)) { # Empirical p-values

    p_values <- p_empirical

    method_used <- rep("empirical", n_test)

  } else if (method == "gamma") { # Gamma approximation
    gamma_fit <- .compute_pvals_gamma(p_empirical = p_empirical,
                                      perm_stats = perm_stats,
                                      obs_stats = obs_stats,
                                      n_test = n_test,
                                      fit_thresh = fit_thresh,
                                      alternative = alternative,
                                      control = gamma_ctrl)

    p_values <- gamma_fit$p_values
    method_used <- gamma_fit$method_used

  } else if (method == "gpd") { # Tail approximation using the GPD

    gpd_fit <- .compute_pvals_gpd(p_empirical = p_empirical,
                                  perm_stats = perm_stats,
                                  obs_stats = obs_stats,
                                  n_test = n_test,
                                  fit_thresh = fit_thresh,
                                  alternative = alternative,
                                  control = gpd_ctrl)

    p_values <- gpd_fit$p_values
    method_used <- gpd_fit$method_used
  }

  #-----------------------------------------------------------------------------
  # Multiple testing adjustment

  # Store unadjusted p-values
  p_unadjusted <- p_values

  adjust_method <- adjust_ctrl$method

  if (adjust_method == "none") {
    adjust_result <- NULL

  } else {

    adjust_result <- mult_adjust(p_values = p_values,
                                 method = adjust_method,
                                 true_null_method = adjust_ctrl$true_null_method,
                                 p_true_null = adjust_ctrl$p_true_null,
                                 seq_length = adjust_ctrl$seq_length,
                                 perm_stats = adjust_ctrl$perm_stats,
                                 cores = adjust_ctrl$cores,
                                 verbose = adjust_ctrl$verbose)

    p_values <- adjust_result$p_adjusted
  }

  #-----------------------------------------------------------------------------
  callArgs <- mget(names(formals()),sys.frame(sys.nframe()))
  callArgs$gpdEstimate <- NULL

  output <- list(p_values = p_values,
                 p_unadjusted = p_unadjusted,
                 p_empirical = p_empirical,
                 gpd_fit = gpd_fit,
                 gamma_fit = gamma_fit,
                 method_used = method_used,
                 adjust_result = adjust_result,
                 control = list(
                   gpd = gpd_ctrl,
                   gamma = gamma_ctrl,
                   adjust = adjust_ctrl
                 )
  )

  return(output)
}


