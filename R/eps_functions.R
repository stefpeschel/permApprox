################################################################################
### Factor of the observed test statistic
################################################################################

#' Simple multiplicative epsilon rule: eps = tune * |t_obs|
#'
#' @description
#' Computes \eqn{\varepsilon}-values as a simple multiplicative factor of the
#' absolute observed test statistic:
#' \deqn{\varepsilon_j = \mathrm{tune} \cdot |T_{\mathrm{obs}, j}|.}
#'
#' This is the simplest possible epsilon rule and can be used as an alternative
#' to \code{\link{eps_slls}}. It is monotone in \code{tune}, making it fully
#' compatible with the zero-guard refinement mechanism of
#' \code{\link{make_gpd_ctrl}}.
#'
#' @param obs_stats Numeric vector of observed test statistics.
#' @param perm_stats Numeric matrix of permutation statistics (ignored by this
#'   epsilon rule; included only for compatibility with the \code{eps_fun}
#'   interface).
#' @param sample_size Numeric scalar or \code{NULL}. Ignored by this epsilon
#'   rule; included only for compatibility.
#' @param tune Numeric scalar. Multiplicative factor controlling the overall
#'   size of \eqn{\varepsilon}. Larger values yield more conservative
#'   (larger) epsilon values.
#' @param ... Further arguments ignored.
#'
#' @return
#' Numeric vector of \eqn{\varepsilon}-values, one per test.
#'
#' @details
#' This epsilon rule is useful when a simple, interpretable, and fast epsilon
#' generator is desired. It does not depend on permutation values, sample size,
#' or shape of the distribution––only on the magnitude of each observed test
#' statistic.
#'
#' For tests where \eqn{T_{\mathrm{obs}} = 0}, a small fallback value
#' \eqn{\varepsilon = \mathrm{tune} \cdot 1e-12} is returned to avoid
#' zero-support issues.
#'
#' @examples
#' obs  <- c(0.5, 2.0, 3.5)
#' perm <- matrix(rnorm(1000 * 3), nrow = 1000)
#'
#' # epsilon = 0.1 * |T_obs|
#' eps_factor(obs_stats = obs, perm_stats = perm, sample_size = NULL, tune = 0.1)
#'
#' # use in permApprox control
#' ctrl <- make_gpd_ctrl(
#'   constraint = "support_at_obs",
#'   eps_fun    = eps_factor,
#'   eps_tune   = 0.05
#' )
#' ctrl
#'
#' @export
eps_factor <- function(obs_stats,
                       perm_stats,
                       sample_size,
                       tune = 0.1,
                       ...) {
  
  if (!is.numeric(obs_stats))
    stop("'obs_stats' must be numeric")
  
  if (!(is.numeric(tune) && length(tune) == 1L && is.finite(tune)))
    stop("'tune' must be a single numeric scalar")
  
  # basic rule: epsilon = tune * |t_obs|
  eps <- tune * abs(obs_stats)
  
  eps
}

################################################################################
### SLLS
################################################################################

#' Standardized Lifted Log-Saturation (SLLS) epsilon rule
#'
#' @description
#' Computes \eqn{\varepsilon}-values using the
#' \emph{Standardized Lifted Log-Saturation (SLLS)} rule. The construction
#' operates on a standardized \emph{Z}-scale, combines a logarithmic
#' saturation curve with a sample-size–dependent plateau, and adds a smooth
#' Wendland-\eqn{C^2} lift term to stabilize \eqn{\varepsilon} for small and
#' moderate test statistics. The saturation anchor is defined by a high
#' quantile of the (standardized) permutation distribution, inflated by a
#' simple sample-size dependent factor.
#'
#' This function is the default epsilon generator (\code{eps_fun}) in
#' \code{\link{make_gpd_ctrl}} when a support constraint is active.
#'
#' @param obs_stats Numeric vector of observed test statistics (one per test).
#' @param perm_stats Numeric matrix of permutation test statistics with one
#'   column per test and one row per permutation replicate. A numeric vector
#'   is allowed for a single test and will be internally reshaped to a
#'   matrix with one column.
#' @param sample_size Positive numeric scalar giving an effective sample size
#'   per test (typically the smaller group size). Controls the curvature and
#'   plateau height via simple \eqn{n}-dependent scaling.
#' @param tune Numeric scalar. Global tuning parameter controlling the target
#'   plateau height of the saturation curve. Larger values yield more
#'   conservative (larger) \eqn{\varepsilon}-values. This corresponds to the
#'   \emph{target factor} and is the quantity adapted by the zero-guard
#'   procedure in \code{\link{make_gpd_ctrl}}.
#' @param cap_base Character string specifying how the saturation cap is
#'   constructed on the standardized Z-scale. One of:
#'   \itemize{
#'     \item \code{"perm"} (default): cap based on a high quantile of
#'       \eqn{|Z_{\text{perm}}|} inflated by a sample-size dependent factor.
#'     \item \code{"tmax"}: cap based on the maximum of \eqn{|Z_{\text{obs}}|}
#'       across tests (common cap).
#'   }
#' @param alpha Numeric scalar controlling the inflation factor
#'   \eqn{g(n) = \max(1, \alpha \sqrt{n})} applied to the permutation-based
#'   cap. Larger values produce a higher cap and therefore a slower
#'   saturation. Default is \code{0.5}.
#' @param q_ref Numeric scalar between 0 and 1. Quantile of the absolute
#'   standardized permutation statistics \eqn{|Z_{\text{perm}}|} used as the
#'   baseline cap when \code{cap_base = "perm"}. Default is \code{0.99}.
#' @param rho_lift Numeric scalar. Multiplier for the Wendland-\eqn{C^2} lift
#'   term, which prevents underestimation of \eqn{\varepsilon} for small and
#'   moderate test statistics. Default is \code{3}.
#' @param k_factor Numeric scalar. Base curvature parameter of the saturation
#'   curve. The effective curvature scales as \code{k_factor * (500 / n)},
#'   where \code{n = sample_size}. Larger values bend the curve earlier.
#'   Default is \code{1000}.
#' @param floor_epsZ Numeric scalar. Minimum allowed value on the Z-scale to
#'   avoid vanishingly small \eqn{\varepsilon}. Default is \code{1e-6}.
#' @param ... Further arguments passed to internal helpers (currently unused).
#'
#' @return
#' Numeric vector of \eqn{\varepsilon}-values on the scale of the original
#' test statistics, with the same length as \code{obs_stats}.
#'
#' @details
#' The SLLS rule proceeds in four steps:
#' \enumerate{
#'   \item Standardize each test statistic using the mean and standard
#'     deviation of its permutation distribution, yielding
#'     \eqn{Z_{\text{obs}}} and \eqn{Z_{\text{perm}}}.
#'   \item Construct a cap \eqn{T_{\text{cap}}} on the Z-scale via
#'     \code{cap_base}. For \code{"perm"}, this uses a high quantile
#'     \eqn{Q_{q_{\mathrm{ref}}}(|Z_{\text{perm}}|)} inflated by
#'     \eqn{g(n) = \max(1, \alpha \sqrt{n})}; for \code{"tmax"}, it uses
#'     \eqn{\max |Z_{\text{obs}}|}.
#'   \item Calibrate a log-saturation curve such that its plateau height is
#'     proportional to \code{tune * (500 / sample_size)}, and add a
#'     Wendland-\eqn{C^2} lift term for small and moderate standardized
#'     statistics. The resulting \eqn{\varepsilon_Z} is bounded below by
#'     \code{floor_epsZ}.
#'   \item Map \eqn{\varepsilon_Z} back to the original test statistic scale
#'     via the per-test standard deviation of the permutation distribution.
#' }
#'
#' This construction yields a smooth, monotone epsilon that grows with the
#' magnitude of the test statistic and can be tuned globally via \code{tune}.
#' It is particularly suited for enforcing support constraints in GPD-based
#' tail approximation under limited numbers of permutations.
#'
#' @examples
#' set.seed(123)
#' perm <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' obs  <- rnorm(5, mean = 2)
#' eps  <- eps_slls(obs_stats = obs,
#'                  perm_stats = perm,
#'                  sample_size = 50,
#'                  tune = 0.25)
#'
#' @export

eps_slls <- function(obs_stats,
                     perm_stats,
                     sample_size,
                     tune          = 0.25, # targer_factor
                     cap_base      = c("perm", "tmax"),  # Z-cap choice
                     alpha         = 0.5, # g(n) = max(1, alpha*sqrt(n))
                     q_ref         = 0.99, # quantile of |Z_perm|
                     rho_lift      = 3,
                     k_factor      = 1000,
                     floor_epsZ    = 1e-6,
                     ...) {
  
  cap_base <- if (is.character(cap_base)) match.arg(cap_base) else cap_base
  n <- sample_size
  stopifnot(n > 0)
  
  # Allow vector perm_stats for single test
  if (is.null(dim(perm_stats))) perm_stats <- matrix(perm_stats, ncol = 1)
  Tobs <- as.numeric(obs_stats)
  
  # Z-standardization per test
  mu  <- colMeans(perm_stats, na.rm = TRUE)
  sig <- apply(perm_stats, 2, sd, na.rm = TRUE)
  
  # Guard against non-finite or zero SDs
  valid_sig <- sig[is.finite(sig) & sig > 0]
  sig_fallback <- if (length(valid_sig)) max(min(valid_sig), 1e-12) else 1e-12
  sig_ok <- sig
  sig_ok[!is.finite(sig_ok) | sig_ok <= 0] <- sig_fallback
  
  Zobs  <- abs((Tobs - mu) / sig_ok)
  Zperm <- sweep(sweep(perm_stats, 2, mu, "-"), 2, sig_ok, "/")
  
  # Build cap on Z-scale
  T_cap_Z <- .build_Tcap(
    cap_base = cap_base,
    t_obs    = Zobs,
    perm_stats = Zperm,
    n        = n,
    q_ref    = q_ref,
    alpha    = alpha
  )
  
  # LLS core on Z
  k        <- k_factor      * (500 / n)
  E_target <- tune * (500 / n) # tune = target_factor
  
  if (!any(is.finite(T_cap_Z)) || all(T_cap_Z <= 0)) {
    # Conservative constant if the cap failed
    epsZ <- rep(max(E_target + rho_lift, floor_epsZ), length(Zobs))
  } else {
    s   <- pmin(Zobs / T_cap_Z, 1)        # scale-free index in [0,1]
    psi <- (1 - s)^4 * (1 + 4*s)          # Wendland C^2 lift
    epsZ_core <- E_target * log1p(k * s) / log1p(k) + rho_lift * psi
    epsZ <- pmax(epsZ_core, floor_epsZ)
  }
  
  # Map back to T-scale
  epsT <- epsZ * sig_ok
  epsT
}

.build_Tcap <- function(cap_base, t_obs, perm_stats, n, 
                        q_ref = 0.99, alpha = 0.5) {
  if (identical(cap_base, "tmax")) {
    # Observed cap (same for all tests)
    T_cap <- rep(max(abs(t_obs), na.rm = TRUE), length(t_obs))
    
  } else if (identical(cap_base, "perm")) {
    # Permutation cap: high quantile of |perm|, inflated by g(n)
    if (is.null(dim(perm_stats))) perm_stats <- matrix(perm_stats, ncol = 1)
    Tabs <- abs(perm_stats)
    Tq <- apply(Tabs, 2, stats::quantile, 
                probs = q_ref, na.rm = TRUE, names = FALSE)
    
    # Inflation g(n) = max(1, alpha * sqrt(n))
    infl <- pmax(1, alpha * sqrt(n))
    T_cap <- Tq * infl
    
    # Broadcast if needed
    if (length(T_cap) == 1L) T_cap <- rep(T_cap, length(t_obs))
    if (length(T_cap) != length(t_obs)) T_cap <- rep(T_cap[1L], length(t_obs))
    
  } else {
    stop("cap_base must be 'perm' or 'tmax'.")
  }
  
  # Guardrails: enforce positive finite caps
  bad <- !is.finite(T_cap) | T_cap <= 0
  if (any(bad)) {
    fallback <- max(T_cap[!bad], 1)
    T_cap[bad] <- fallback
  }
  T_cap
}

#-------------------------------------------------------------------------------
#' @title Define epsilon values based on the user-defined rule
#'
#' @keywords internal
.define_eps <- function(perm_stats, obs_stats, sample_size,
                        constraint, eps_fun, eps_tune, eps_args) {
  
  # Epsilon only used for support constraints
  if (!constraint %in% c("support_at_obs", "support_at_max")) {
    return(rep(0, length(obs_stats)))
  }
  
  if (is.null(eps_fun))
    stop("For constrained fitting, 'eps_fun' must be non-NULL.")
  
  if (is.null(eps_args)) eps_args <- list()
  
  base_args <- list(
    obs_stats   = obs_stats,
    perm_stats  = perm_stats,
    sample_size = sample_size,
    tune        = eps_tune
  )
  
  eps <- do.call(eps_fun, c(base_args, eps_args))
  
  if (!is.numeric(eps) || length(eps) != length(obs_stats))
    stop("The epsilon vector returned by 'eps_fun' must be numeric ",
         "and of length(obs_stats).")
  
  eps
}



