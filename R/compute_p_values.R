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
#'       distribution if the empirical p-value falls below \code{fit_thresh}.
#'     \item \code{"gpd"}: A Generalized Pareto Distribution (GPD) is fitted
#'       to the tail of the permutation distribution if the empirical p-value
#'       falls below \code{fit_thresh}.
#'   }
#'
#' @param fit_thresh Numeric. Threshold on empirical p-values below which
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
#' # GPD approximation without constraint
#' gpd_ctrl <- make_gpd_ctrl(constraint = "unconstrained")
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
#' # GPD approximation with constraint and fixed epsilon
#' gpd_ctrl <- make_gpd_ctrl(constraint = "support_at_max",
#'                           eps_fun = eps_fixed, 
#'                           eps_par = list(value = 0.1))
#' 
#' res_gpd_constr_eps0.1 <- compute_p_values(obs_stats = obs,
#'                                    perm_stats = perm,
#'                                    method = "gpd",
#'                                    gpd_ctrl = gpd_ctrl)
#' 
#' # Data frame with (unadjusted) p-values
#' p_values <- data.frame(empirical = res_emp$p_unadjusted,
#'                        gamma = res_gamma$p_unadjusted,
#'                        gpd = res_gpd$p_unadjusted,
#'                        gpd_constr = res_gpd_constr$p_unadjusted,
#'                        gpd_constr_eps0.1 = res_gpd_constr_eps0.1$p_unadjusted)
#' 
#' p_values
#'
#' @export

compute_p_values <- function(
    obs_stats,
    perm_stats,
    method = "gpd",
    fit_thresh = 0.1,
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
  
  # Validate 'alternative'
  alternative <- match.arg(alternative,
                           choices = c("greater", "less", "two_sided"))
  
  # Validate 'method'
  method <- match.arg(method,
                      choices = c("gpd", "gamma", "empirical"))
  
  # Validate 'power'
  stopifnot(is.numeric(power))
  
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
  } else if (method == "gamma") {
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
    perm_stats <- matrix(perm_stats, ncol = 1)
  } else {
    perm_stats <- as.matrix(perm_stats)
  }
  
  n_perm <- nrow(perm_stats)    # permutations
  n_test <- ncol(perm_stats)    # tests
  
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
      mean = if (n_test == 1) mean(perm_stats) else colMeans(perm_stats),
      median = if (n_test == 1) median(perm_stats) else apply(perm_stats, 2, median),
      stop("Invalid 'null_center': must be numeric, 'mean', or 'median'")
    )
  }
  
  # Center permutation statistics
  if (n_test == 1) {
    perm_stats <- perm_stats - center_vec
  } else {
    perm_stats <- sweep(perm_stats, 2, center_vec, FUN = "-")
  }
  
  # Center observed statistic(s)
  obs_stats <- obs_stats - center_vec
  
  ## ---------------------------------------------------------------------------
  ## Power transformation
  ## ---------------------------------------------------------------------------
  if (power != 1) {
    perm_stats <- perm_stats^power
    obs_stats <- obs_stats^power
  }
  
  ## ---------------------------------------------------------------------------
  ## Transform according to the alternative (one-sided scale)
  ## ---------------------------------------------------------------------------
  transformed <- lapply(seq_len(n_test), function(i) {
    .transform_stats(
      perm_stats  = perm_stats[, i],
      obs_stats   = obs_stats[i],
      alternative = alternative
    )
  })
  t_obs       <- vapply(transformed, `[[`, numeric(1), "obs_stats")
  t_perm_list <- lapply(transformed, `[[`, "perm_stats")
  t_perm      <- do.call(cbind, t_perm_list)
  dimnames(t_perm) <- dimnames(perm_stats)
  
  ## ---------------------------------------------------------------------------
  ## Empirical p-values
  ## ---------------------------------------------------------------------------
  pvals_emp_list <- .compute_pvals_emp(
    obs_stats = t_obs,
    perm_stats = t_perm
  )
  
  p_empirical      <- pvals_emp_list$pvals
  n_perm_exceeding <- pvals_emp_list$n_perm_exceeding
  method_used      <- rep("empirical", n_test)
  
  # Decide which tests are candidates for parametric approximation
  idx_fit <- which(!is.na(p_empirical) & (p_empirical < fit_thresh))
  
  if (length(idx_fit) == 0L && method != "empirical" && verbose) {
    message("No empirical p-values below 'fit_thresh'; returning empirical p-values.")
  }
  
  ## ---------------------------------------------------------------------------
  ## Approximate p-values
  ## ---------------------------------------------------------------------------
  
  # Initialize gamma_fit and gpd_fit
  gamma_fit <- gpd_fit <- NULL
  
  # Initialize control output list
  control <- list()
  
  if (method == "empirical") { # Empirical p-values
    
    p_values <- p_empirical
    
  } else if (method == "gamma") { # Gamma approximation
    
    if (length(idx_fit) == 0L) {
      # nothing to approximate
      p_values <- p_empirical
      
    } else {
      control$gamma <- gamma_ctrl
      
      gamma_fit <- .compute_pvals_gamma(
        obs_stats    = t_obs,
        perm_stats   = t_perm,
        idx_fit      = idx_fit,
        control      = control$gamma,
        cores        = cores,
        parallel_min = parallel_min,
        verbose      = verbose
      )
      
      # Replace p-values with successful Gamma fit and reset method
      p_values <- p_empirical
      success <- gamma_fit$status == "success"
      p_values[success]    <- gamma_fit$p_value[success]
      method_used[success] <- "gpd"
    }
    
  } else if (method == "gpd") { # Tail approximation using the GPD
    
    if (length(idx_fit) == 0L) {
      p_values <- p_empirical
      
    } else {
      control$gpd <- gpd_ctrl
      
      gpd_fit <- .compute_pvals_gpd(
        obs_stats    = t_obs,
        perm_stats   = t_perm,
        idx_fit      = idx_fit,
        p_empirical  = p_empirical,
        control      = control$gpd,
        cores        = cores,
        parallel_min = parallel_min,
        verbose      = verbose,
        ...
      )
      
      # Replace p-values with successful GPD fit and reset method
      p_values <- p_empirical
      success <- gpd_fit$status == "success"
      p_values[success]    <- gpd_fit$p_value[success]
      method_used[success] <- "gpd"
    }
  }
  
  ## ---------------------------------------------------------------------------
  ## Multiple testing adjustment
  ## ---------------------------------------------------------------------------
  
  # Store unadjusted p-values
  p_unadjusted <- p_values
  
  if (adjust_method == "none") {
    adjust_result <- NULL
    
  } else {
    
    adjust_result <- mult_adjust(
      p_values = p_values,
      method = adjust_method,
      true_null_method = adjust_ctrl$true_null_method,
      p_true_null = adjust_ctrl$p_true_null,
      verbose = verbose
    )
    
    p_values <- adjust_result$p_adjusted
  }
  
  #-----------------------------------------------------------------------------
  
  output <- list(p_values = p_values,
                 p_unadjusted = p_unadjusted,
                 p_empirical = p_empirical,
                 gpd_fit = gpd_fit,
                 gamma_fit = gamma_fit,
                 method_used = method_used,
                 adjust_result = adjust_result,
                 control = control
  )
  
  return(output)
}


