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
#' @param perm_stats Permutation test statistics. Should be provided as a 
#'   numeric matrix or data frame with permutations in rows and tests in 
#'   columns.
#'
#' @param method Character. Method used to compute p-values. Default is
#'   \code{"gpd"}. Options are:
#'   \itemize{
#'     \item \code{"empirical"}: Empirical p-values are returned directly,
#'     based on the number of permutation test statistics that are as or
#'     more extreme than the observed statistic.
#'     \item \code{"gamma"}: A Gamma distribution is fitted to the permutation
#'       distribution if the empirical p-value falls below \code{approx_thresh}.
#'     \item \code{"gpd"}: A Generalized Pareto Distribution (GPD) is fitted
#'       to the tail of the permutation distribution if the empirical p-value
#'       falls below \code{approx_thresh}.
#'   }
#'
#' @param approx_thresh Numeric. Threshold on empirical p-values below which
#'   parametric fitting is applied. Default: \code{0.1} (parametric
#'   approximation if the empirical p-value is smaller than 0.1).
#'
#' @param alternative Character. One of \code{"greater"}, \code{"less"}, or
#'   \code{"two_sided"} (default), indicating the tail of the test.
#'   
#' @param null_center Numeric or character. Specifies the value around which the
#'   null distribution is centered. If set to \code{"mean"} or \code{"median"},
#'   the per-row mean or median of \code{perm_stats} is used instead.
#'   This allows testing against a null hypothesis other than zero or centering
#'   based on the empirical distribution.
#'   
#' @param power Numeric scalar. Optional monotone power transformation applied
#'   to both permutation statistics and the observed statistic before threshold
#'   detection and GPD fitting. Must be positive. 
#'   Default: 1 (no transformation).
#'
#' @param adjust_method Character. Method for multiple testing adjustment.
#'   Options include:
#'   \itemize{
#'     \item \code{"none"}: No adjustment.
#'     \item \code{"lfdr"}: Local false discovery rates via \code{fdrtool}.
#'     \item \code{"adapt_BH"}: Adaptive Benjamini-Hochberg (requires estimation
#'     of the proportion of true nulls).
#'     \item Any method supported by \code{stats::p.adjust} (e.g., "holm", "BH",
#'     "BY").
#'   }
#'
#' @param cores Integer >= 1. Number of CPU cores for parallel computations.
#'   Default: 1.
#'   
#' @param parallel_min Integer. Minimum number of tests in a given
#'   computation step (e.g., threshold detection, GPD/Gamma fitting, epsilon
#'   refinement) for which parallel computation is used when \code{cores > 1}.
#'   If the number of tests is smaller than \code{parallel_min}, the
#'   corresponding step is run sequentially to avoid parallel overhead.
#'   Default: \code{10}.
#'
#' @param verbose Logical scalar. If \code{TRUE}, progress messages are printed.
#'   Default: TRUE.
#'
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
#' @return
#' A list with the following components:
#'
#' \describe{
#'   \item{\code{p_values}}{
#'     Numeric vector of final p-values. These are empirical p-values or,
#'     for tests where parametric tail approximation was successful,
#'     Gamma- or GPD-approximated p-values. If multiple testing correction
#'     was requested, \code{p_values} contains the adjusted p-values.
#'   }
#'   \item{\code{p_unadjusted}}{
#'     Numeric vector of unadjusted p-values, combining empirical and
#'     approximated values (where applicable).
#'   }
#'   \item{\code{p_empirical}}{
#'     Numeric vector of raw empirical p-values computed directly from the
#'     permutation distribution (before any approximation or adjustment).
#'   }
#'   \item{\code{n_perm_exceed}}{
#'   Integer vector of length \code{n_test} giving, for each test,
#'   the number of permutation test statistics that are as or more
#'   extreme than the observed statistic, according to the specified
#'   test \code{alternative}. 
#'   This corresponds to the numerator \eqn{r} in the empirical
#'   permutation p-value formula \eqn{p = (r + 1) / (B + 1)}.
#'   }
#'   \item{\code{fit_method}}{
#'     Character scalar indicating the parametric tail-approximation method
#'     used: one of \code{"gpd"}, \code{"gamma"}, or \code{"empirical"}.
#'     When \code{"empirical"}, no tail approximation was performed.
#'   }
#'   \item{\code{fit_result}}{
#'     A list containing the detailed output of the chosen tail-approximation
#'     method.  
#'
#'     If \code{fit_method = "gpd"}, this is the result returned by
#'     \code{.compute_pvals_gpd()} (including GPD parameters, thresholds,
#'     exceedances, epsilon values, status flags, and fitted exceedances).
#'
#'     If \code{fit_method = "gamma"}, this is the result returned by
#'     \code{.compute_pvals_gamma()} (including Gamma parameters,
#'     goodness-of-fit results, status flags, and permutation values used
#'     for the fit).  
#'
#'     If \code{fit_method = "empirical"}, this entry is \code{NULL}.
#'   }
#'   \item{\code{method_used}}{
#'     Character vector of length \code{n_test}. For each hypothesis test,
#'     indicates whether the final p-value was empirical (\code{"empirical"})
#'     or came from Gamma or GPD approximation (\code{"gamma"} or \code{"gpd"}).
#'   }
#'   \item{\code{adjust_result}}{
#'     The full output of \code{mult_adjust()} if multiple testing correction
#'     was applied, otherwise \code{NULL}.
#'   }
#'   \item{\code{control}}{
#'     A list containing the control settings actually used in the computation,
#'     including the selected approximation control (Gamma or GPD) and the
#'     multiple testing adjustment control.
#'   }
#' }
#'
#' @details
#' For each test, \code{perm_approx()} first constructs a centered and,
#' optionally, transformed permutation distribution and then computes
#' empirical p-values. Parametric tail approximation (via Gamma or GPD) is
#' applied only for tests whose empirical p-value falls below
#' \code{approx_thresh}; all other tests remain purely empirical.
#'
#' \subsection{Centering and transformation}{
#' The input permutation statistics \code{perm_stats} are organized with
#' permutations in rows and tests in columns. The observed statistics
#' \code{obs_stats} are aligned with the columns of \code{perm_stats}.
#'
#' The null distribution is centered according to \code{null_center}:
#' \itemize{
#'   \item If \code{null_center} is numeric, the same value is subtracted
#'         from both \code{obs_stats} and \code{perm_stats}.
#'   \item If \code{null_center = "mean"} or \code{"median"}, the per-test
#'         mean or median of \code{perm_stats} is used as the centering
#'         value.
#' }
#' Centering is applied to both observed and permuted statistics.
#'
#' Optionally, a monotone power transformation with exponent \code{power > 0}
#' is applied to all centered statistics. This can increase sensitivity in
#' the tail when the scale of the test statistic is highly skewed.
#'
#' Finally, test statistics are transformed to a one-sided scale according
#' to \code{alternative} (\code{"greater"}, \code{"less"}, or
#' \code{"two_sided"}) via an internal helper function
#' \code{.transform_stats()}, so that all subsequent computations operate on
#' a common tail direction.
#' }
#'
#' \subsection{Empirical p-values}{
#' Empirical p-values are computed for each test by comparing the observed
#' statistic to its permutation distribution on the transformed scale.
#' Specifically, if \eqn{r} denotes the number of permutation statistics
#' that are at least as extreme as the observed statistic (in the direction
#' defined by \code{alternative}) and \eqn{B} is the number of permutations,
#' the empirical p-value is
#' \deqn{
#'   \hat{p}_{\mathrm{emp}} = \frac{r + 1}{B + 1} \, .
#' }
#' This correction ensures that empirical p-values are never exactly zero.
#'
#' These empirical p-values are always computed and returned in
#' \code{p_empirical}, regardless of the chosen \code{method}.
#' }
#'
#' \subsection{Parametric tail approximation}{
#' If \code{method = "gpd"} or \code{"gamma"}, a parametric approximation
#' is applied only to tests with empirical p-values below
#' \code{approx_thresh}. The indices of these tests are passed to the
#' corresponding internal workhorse:
#' \itemize{
#'   \item \code{.compute_pvals_gpd()} for Generalized Pareto (GPD) tail
#'         approximation of the permutation distribution,
#'   \item \code{.compute_pvals_gamma()} for Gamma approximation of the full
#'         permutation distribution.
#' }
#'
#' Both functions return a list with one entry per test (see their
#' respective documentation), including approximated tail p-values and
#' diagnostic information such as fitted parameters, goodness-of-fit
#' p-values, and per-test status flags. This list is returned as
#' \code{fit_result}, and the chosen approximation is indicated by
#' \code{fit_method} (\code{"gpd"} or \code{"gamma"}). If
#' \code{method = "empirical"} or no tests satisfy the threshold
#' \code{approx_thresh}, \code{fit_method = "empirical"} and
#' \code{fit_result} is \code{NULL}.
#'
#' For tests where the parametric fit is deemed successful
#' (typically those with \code{status == "success"} in \code{fit_result}),
#' the empirical p-value is replaced by the approximated tail p-value.
#' For all other tests, the empirical p-value is retained. The resulting
#' mixture of empirical and approximated p-values is returned in
#' \code{p_unadjusted}. The per-test origin of the final p-value is tracked
#' in \code{method_used}, indicating \code{"empirical"}, \code{"gamma"}, or
#' \code{"gpd"} for each test.
#' }
#'
#' \subsection{Multiple testing adjustment}{
#' If \code{adjust_method = "none"}, no multiple testing correction is
#' applied, and \code{p_values} is identical to \code{p_unadjusted}.
#' Otherwise, \code{p_values} contains multiplicity-adjusted p-values
#' obtained by calling \code{mult_adjust()} with the specified adjustment
#' method and control parameters in \code{adjust_ctrl}. The full output of
#' \code{mult_adjust()} is returned in \code{adjust_result}.
#'
#' The \code{control} element in the output collects the control settings
#' actually used for tail approximation (GPD or Gamma) and multiple testing
#' adjustment, which facilitates reproducibility and downstream inspection.
#' }
#'
#' @examples
#'
#' ## ---------------------------------------------------------------------
#' ## 10 tests, one truly non-null (two-sample mean differences)
#' ## ---------------------------------------------------------------------
#'
#' set.seed(42)
#'
#' n_per_group <- 50
#' m_tests     <- 10
#'
#' # Group labels: 0 = control, 1 = treatment
#' group <- rep(c(0, 1), each = n_per_group)
#'
#' # Data matrix: rows = samples, cols = tests/features
#' X <- matrix(
#'   rnorm(2 * n_per_group * m_tests, mean = 0, sd = 1),
#'   ncol = m_tests
#' )
#'
#' # Introduce a true effect in the first test (column 1)
#' X[group == 1, 1] <- X[group == 1, 1] + 0.8
#'
#' # Observed test statistics: difference in means (treated - control)
#' obs_stats <- colMeans(X[group == 1, , drop = FALSE]) -
#'              colMeans(X[group == 0, , drop = FALSE])
#'
#' # Permutation distribution: shuffle group labels and recompute all 10 stats
#' B <- 500
#' perm_mat <- matrix(NA_real_, nrow = B, ncol = m_tests)
#'
#' for (b in seq_len(B)) {
#'   grp_perm <- sample(group)
#'   perm_mat[b, ] <- colMeans(X[grp_perm == 1, , drop = FALSE]) -
#'                    colMeans(X[grp_perm == 0, , drop = FALSE])
#' }
#'
#' ## ---------------------------------------------------------------------
#' ## Empirical p-values only (BH-adjusted by default)
#' ## ---------------------------------------------------------------------
#'
#' res_emp <- perm_approx(
#'   obs_stats  = obs_stats,
#'   perm_stats = perm_mat,
#'   method     = "empirical",
#'   verbose    = FALSE
#' )
#'
#' # Print and summary methods
#' res_emp
#' summary(res_emp)
#'
#' # Adjusted p-values (BH by default)
#' res_emp$p_values
#'
#' # Unadjusted p-values
#' res_emp$p_unadjusted
#'
#'
#' ## ---------------------------------------------------------------------
#' ## GPD approximation for the same 10 tests
#' ## ---------------------------------------------------------------------
#'
#' # GPD-related arguments are changed using make_gpd_ctrl()
#' gpd_ctrl <- make_gpd_ctrl(
#'   constraint  = "support_at_obs",
#'   sample_size = n_per_group
#' )
#'
#' # Run p-value approximation with GPD
#' res_gpd <- perm_approx(
#'   obs_stats  = obs_stats,
#'   perm_stats = perm_mat,
#'   method     = "gpd",
#'   gpd_ctrl   = gpd_ctrl,
#'   verbose    = FALSE
#' )
#'
#' # Print and summary for GPD-based approximation
#' res_gpd
#' summary(res_gpd)
#'
#' # Tail fit details
#' fit <- res_gpd$fit_result
#' fit$status      # success / discrete / no_threshold / gof_reject / fit_failed
#' fit$shape       # GPD shape estimates per test
#' fit$thresh      # thresholds
#' fit$n_exceed    # number of exceedances
#'
#' # Compare empirical vs approximated (unadjusted) p-values
#' emp_gpd_df <- data.frame(
#'   empirical = res_emp$p_unadjusted,
#'   GPD       = res_gpd$p_unadjusted
#' )
#' emp_gpd_df
#'
#'
#' ## ---------------------------------------------------------------------
#' ## Gamma approximation for the same 10 tests
#' ## ---------------------------------------------------------------------
#'
#' # Gamma-related arguments are changed using make_gamma_ctrl()
#' gamma_ctrl <- make_gamma_ctrl(gof_test = "none")
#'
#' res_gamma <- perm_approx(
#'   obs_stats   = obs_stats,
#'   perm_stats  = perm_mat,
#'   method      = "gamma",
#'   gamma_ctrl  = gamma_ctrl,
#'   verbose     = FALSE
#' )
#'
#' # Print and summary for Gamma approximation
#' res_gamma
#' summary(res_gamma)
#'
#' # Gamma fit parameters (successful fits)
#' gamma_fit <- res_gamma$fit_result
#' gamma_fit$shape
#' gamma_fit$rate
#'
#'
#' ## ---------------------------------------------------------------------
#' ## Multiple testing adjustment with adapt_BH for GPD
#' ## ---------------------------------------------------------------------
#'
#' # Adaptive BH method for multiple testing adjustment
#' adjust_ctrl <- make_adjust_ctrl(true_null_method = "lfdr")
#'
#' res_adapt_BH <- perm_approx(
#'   obs_stats     = obs_stats,
#'   perm_stats    = perm_mat,
#'   method        = "gpd",
#'   gpd_ctrl      = gpd_ctrl,
#'   adjust_method = "adapt_BH",
#'   adjust_ctrl   = adjust_ctrl,
#'   verbose       = FALSE
#' )
#'
#' # Adjusted vs. unadjusted GPD-based p-values
#' adj_unadj_df <- data.frame(
#' GPD_unadj    = res_gpd$p_unadjusted,
#' GPD_BH       = res_gpd$p_values, 
#' GPD_adapt_BH = res_adapt_BH$p_values
#' )
#' adj_unadj_df
#'
#' @export

perm_approx <- function(
    obs_stats,
    perm_stats,
    method = "gpd",
    approx_thresh = 0.1,
    alternative = "two_sided",
    null_center = 0,
    power = 1,
    adjust_method = "BH",
    cores = 1,
    parallel_min = 10L,
    verbose = TRUE,
    gpd_ctrl = make_gpd_ctrl(),
    gamma_ctrl = make_gamma_ctrl(),
    adjust_ctrl = make_adjust_ctrl(),
    ...
) {
  
  ## ---------------------------------------------------------------------------
  ## Basic argument validation
  ## ---------------------------------------------------------------------------
  
  # Validate 'alternative'
  alternative <- match.arg(alternative,
                           choices = c("greater", "less", "two_sided"))
  
  # Validate 'method'
  method <- match.arg(method,
                      choices = c("gpd", "gamma", "empirical"))
  
  # Validate 'power'
  stopifnot(is.numeric(power), length(power) == 1L, power > 0)
  
  # Validate 'cores'
  if (!(is.numeric(cores) && length(cores) == 1L && 
        is.finite(cores) && cores >= 1))
    stop("'cores' must be a single integer >= 1.")
  cores <- as.integer(cores)
  
  # Validate 'parallel_min'
  parallel_min <- as.integer(parallel_min)
  if (!is.finite(parallel_min) || parallel_min < 1L) {
    stop("'parallel_min' must be a single integer >= 1.")
  }
  
  # Validate 'verbose'
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("'verbose' must be a single logical.")
  
  # Validate control arguments
  if (method == "gpd") {
    if (is.null(gpd_ctrl)) {
      gpd_ctrl <- make_gpd_ctrl()
    } else if (!inherits(gpd_ctrl, "gpd_ctrl")) {
      stop("'gpd_ctrl' must be created with make_gpd_ctrl().")
    }
  }
  if (method == "gamma") {
    if (is.null(gamma_ctrl)) {
      gamma_ctrl <- make_gamma_ctrl()
    } else if (!inherits(gamma_ctrl, "gamma_ctrl")) {
      stop("'gamma_ctrl' must be created with make_gamma_ctrl().")
    }
  }
  if (!inherits(adjust_ctrl, "adjust_ctrl")) {
    stop("'adjust_ctrl' must be created with make_adjust_ctrl().")
  }
  
  # Validate multiple testing adjustment method
  adjust_method <- match.arg(adjust_method, c("none", p.adjust.methods,
                                              "lfdr", "adapt_BH"))
  
  ## ---------------------------------------------------------------------------
  ## Ensure matrix format: permutations in rows, tests in columns
  ## ---------------------------------------------------------------------------
  
  if (is.null(dim(perm_stats))) {
    # Single test: vector of permutation statistics
    perm_stats <- matrix(perm_stats, ncol = 1L)
  } else {
    perm_stats <- as.matrix(perm_stats)
  }
  
  n_perm <- nrow(perm_stats)    # number of permutations
  n_test <- ncol(perm_stats)    # number of tests
  
  # Match observed statistics length
  if (length(obs_stats) != n_test) {
    stop("Length of 'obs_stats' must match number of columns in 'perm_stats'.")
  }
  
  ## ---------------------------------------------------------------------------
  ## Center test statistics around H0
  ## ---------------------------------------------------------------------------
  
  # Determine centering vector
  if (is.numeric(null_center)) {
    center_vec <- null_center
  } else {
    center_vec <- switch(
      null_center,
      mean = if (n_test == 1L) mean(perm_stats) else colMeans(perm_stats),
      median = if (n_test == 1L) median(perm_stats) else apply(perm_stats, 2L, median),
      stop("Invalid 'null_center': must be numeric, 'mean', or 'median'")
    )
  }
  
  # Center permutation statistics
  if (n_test == 1L) {
    perm_stats <- perm_stats - center_vec
  } else {
    perm_stats <- sweep(perm_stats, 2L, center_vec, FUN = "-")
  }
  
  # Center observed statistic(s)
  obs_stats <- obs_stats - center_vec
  
  ## ---------------------------------------------------------------------------
  ## Power transformation
  ## ---------------------------------------------------------------------------
  if (power != 1) {
    perm_stats <- perm_stats^power
    obs_stats  <- obs_stats^power
  }
  
  ## ---------------------------------------------------------------------------
  ## Transform according to the alternative (one-sided scale)
  ## ---------------------------------------------------------------------------
  transformed <- lapply(seq_len(n_test), function(i) {
    .transform_stats(
      perm_stats  = perm_stats[, i],
      obs_stat    = obs_stats[i],
      alternative = alternative
    )
  })
  t_obs       <- vapply(transformed, `[[`, numeric(1), "obs_stat")
  t_perm_list <- lapply(transformed, `[[`, "perm_stats")
  t_perm      <- do.call(cbind, t_perm_list)
  dimnames(t_perm) <- dimnames(perm_stats)
  
  ## ---------------------------------------------------------------------------
  ## Empirical p-values
  ## ---------------------------------------------------------------------------
  pvals_emp_list <- .compute_pvals_emp(
    obs_stats  = t_obs,
    perm_stats = t_perm
  )
  
  p_empirical   <- pvals_emp_list$pvals
  n_perm_exceed <- pvals_emp_list$n_perm_exceed # exceeding permutation statistics
  method_used   <- rep("empirical", n_test)
  
  ## ---------------------------------------------------------------------------
  ## Decide which tests are candidates for parametric approximation
  ## Only consider tests with:
  ## - finite empirical p-value
  ## - p_empirical < approx_thresh
  ## - positive transformed observed statistic (t_obs > 0)
  ##   => ensures we're in the upper tail on the working scale
  ## ---------------------------------------------------------------------------
  idx_fit <- which(
    !is.na(p_empirical) &
      (p_empirical < approx_thresh) &
      (t_obs > 0)
  )
  
  if (length(idx_fit) == 0L && method != "empirical" && isTRUE(verbose)) {
    message(
      "No tests selected for parametric approximation ",
      "(either p_empirical >= approx_thresh or transformed obs_stat <= 0); ",
      "returning empirical p-values."
    )
  }
  
  ## ---------------------------------------------------------------------------
  ## Approximate p-values
  ## ---------------------------------------------------------------------------
  fit_result <- NULL
  fit_method <- method   # store which one was actually used
  
  # Control list for output
  control_out <- list(
    adjust        = adjust_ctrl,
    approx_thresh = approx_thresh,
    adjust_method = adjust_method
  )
  
  p_values <- p_empirical
  
  if (method == "empirical") {
    
    fit_result <- NULL
    
  } else if (method == "gpd") {
    
    control_out$gpd <- gpd_ctrl
    
    # Call GPD fitter even if idx_fit is empty: skeleton list if so
    fit_result <- .compute_pvals_gpd(
      obs_stats    = t_obs,
      perm_stats   = t_perm,
      n_perm       = n_perm,
      idx_fit      = idx_fit,
      control      = gpd_ctrl,
      cores        = cores,
      parallel_min = parallel_min,
      verbose      = verbose,
      ...
    )
    
    success <- fit_result$status == "success"
    p_values[success]    <- fit_result$p_value[success]
    method_used[success] <- "gpd"
    
  } else if (method == "gamma") {
    
    control_out$gamma <- gamma_ctrl
    
    # Call gamma fitter even if idx_fit is empty: it will return a skeleton list
    fit_result <- .compute_pvals_gamma(
      obs_stats    = t_obs,
      perm_stats   = t_perm,
      idx_fit      = idx_fit,
      control      = gamma_ctrl,
      cores        = cores,
      parallel_min = parallel_min,
      verbose      = verbose
    )
    
    success <- fit_result$status == "success"
    p_values[success]    <- fit_result$p_value[success]
    method_used[success] <- "gamma"
    
  }
  
  ## ---------------------------------------------------------------------------
  ## Multiple testing adjustment
  ## ---------------------------------------------------------------------------
  
  p_unadjusted <- p_values
  adjust_result <- NULL
  
  if (adjust_method != "none") {
    adjust_result <- mult_adjust(
      p_values       = p_values,
      method         = adjust_method,
      true_null_method = adjust_ctrl$true_null_method,
      p_true_null    = adjust_ctrl$p_true_null,
      verbose        = verbose
    )
    p_values <- adjust_result$p_adjusted
  }
  
  ## ---------------------------------------------------------------------------
  ## Final output
  ## ---------------------------------------------------------------------------
  
  output <- list(
    p_values      = p_values,
    p_unadjusted  = p_unadjusted,
    p_empirical   = p_empirical,
    n_perm_exceed = n_perm_exceed,
    fit_method    = fit_method,
    fit_result    = fit_result,
    method_used   = method_used,
    adjust_result = adjust_result,
    control       = control_out
  )
  
  class(output) <- "perm_approx"
  output
}
