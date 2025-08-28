#' @title Standardized Lifted Log-Saturation (SLLS) epsilon rule
#'
#' @description
#' Computes \eqn{\varepsilon}-values using the 
#' \emph{Standardized Lifted Log-Saturation (SLLS)} 
#' rule. The rule operates on the standardized \emph{Z}-scale, combines 
#' logarithmic saturation with a sample-size–dependent plateau, and adds a 
#' smooth Wendland-\eqn{C^2} lift term to stabilize small and moderate test 
#' statistics. The saturation anchor is defined by a high quantile of the 
#' permutation distribution, inflated by a simple sample-size dependent scaling.
#'
#' @param obs_stats Numeric vector of observed test statistics.
#' @param perm_stats Matrix of permutation test statistics with one column per test 
#'   and one row per permutation replicate.
#' @param sample_size Per-group sample size of the test (used to scale curvature and plateau).
#' @param k_factor Numeric scalar. Controls the curvature of the saturation curve.
#'   Larger values bend the curve earlier. Default is \code{1000}.
#' @param target_factor Numeric scalar. Controls the plateau height at the anchor 
#'   point. Larger values yield more conservative \eqn{\varepsilon}. Default is \code{0.25}.
#' @param rho_lift Numeric scalar. Multiplier for the Wendland-\eqn{C^2} lift term, 
#'   which prevents underestimation of \eqn{\varepsilon} for small and moderate test 
#'   statistics. Default is \code{3}.
#' @param floor_epsZ Numeric scalar. Minimum allowed value on the \emph{Z}-scale 
#'   to avoid vanishingly small \eqn{\varepsilon}. Default is \code{1e-6}.
#' @param q_ref Numeric scalar between 0 and 1. Quantile of the absolute standardized 
#'   permutation statistics used as the baseline cap. Default is \code{0.99}.
#' @param n_ref Reference sample size used to anchor the inflation factor. 
#'   Default is \code{100}.
#' @param anchor_mult Numeric scalar. Multiplier at the reference sample size 
#'   (\code{n_ref}). Default is \code{5}.
#' @param beta Numeric scalar. Growth exponent for the inflation factor. 
#'   \code{beta = 0.5} corresponds to \eqn{\sqrt{n}} scaling. Default is \code{0.5}.
#' @param ... Further arguments passed to internal methods (currently unused).
#'
#' @return Numeric vector of \eqn{\varepsilon}-values on the scale of the original 
#'   test statistics (same length as \code{obs_stats}).
#'
#' @details
#' The SLLS rule first standardizes test statistics using the mean and standard 
#' deviation of their permutation distribution. It then defines a reference point 
#' \eqn{Z_\mathrm{cap}} as a high quantile of the standardized permutation values, 
#' inflated by a sample-size dependent factor. A log-saturation curve is calibrated 
#' to this anchor, and a Wendland-\eqn{C^2} kernel is added as a lift term for 
#' small and moderate values. The result is mapped back to the original test 
#' statistic scale.
#'
#' @examples
#' set.seed(123)
#' perm <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' obs  <- rnorm(5, mean = 2)
#' eps_slls(obs, perm, sample_size = 50)
#'
#' @export
eps_slls <- function(obs_stats,
                     perm_stats,
                     sample_size,
                     k_factor = 1000,       # curvature
                     target_factor = 0.25,  # target plateau height
                     rho_lift = 3,          # lift multiplier
                     floor_epsZ = 1e-6, 
                     # permutation-based cap controls:
                     q_ref         = 0.99,  # quantile of |Z_perm|
                     n_ref         = 100,   # anchor n for inflation
                     anchor_mult   = 5,     # multiplier at n_ref
                     beta          = 0.5,   # growth exponent (~sqrt(n))
                     ...) {
  # Standardize to Z-scale
  mu  <- colMeans(perm_stats)
  sig <- apply(perm_stats, 2, sd)
  Zobs <- abs((obs_stats - mu) / sig)
  
  # Permutation-based Zcap
  Zperm_abs <- abs(sweep(sweep(perm_stats, 2, mu, "-"), 2, sig, "/"))
  Zcap_base <- apply(Zperm_abs, 2, stats::quantile, 
                     probs = q_ref, na.rm = TRUE, names = FALSE)
  
  inflation <- 1 + (anchor_mult - 1) * (sample_size / n_ref)^beta
  Zcap <- Zcap_base * inflation
  
  # Curvature and target plateau height in Z-scale
  k        <- k_factor * (500 / sample_size)
  E_target <- target_factor * (500 / sample_size)
  
  # Calibrate to Zcap
  B <- k / Zcap
  A <- E_target / log1p(B * Zcap)
  
  # Relative position and lift term
  s         <- pmin(Zobs / Zcap, 1)
  psi_wendland_c2 <- (1 - s)^4 * (1 + 4*s)
  lift_term <- rho_lift * (1 - s)^4 * (1 + 4*s)
  
  # Epsilon in Z-space
  epsZ <- pmax(A * log1p(B * Zobs) + lift_term, floor_epsZ)
  
  # Convert back to T-scale
  epsT <- epsZ * sig
  epsT
}

#-------------------------------------------------------------------------------
#' @title Define epsilon values based on the user-defined rule
#'
#' @keywords internal
.define_eps <- function(perm_stats, obs_stats, sample_size,
                        constraint, eps_rule, eps_par) {
  
  if (constraint != "unconstrained") {
    
    if (eps_rule == "constant") {
      if (length(eps_par) == 1) {
        eps <- rep(eps_par, length(obs_stats))
      } else {
        eps <- eps_par
      }
      
    } else if (eps_rule == "factor") {
      eps <- eps_par * obs_stats
      
    } else if (eps_rule == "slls") {
      eps <- eps_slls(obs_stats = obs_stats,
                      perm_stats = perm_stats,
                      sample_size = sample_size,
                      target_factor = eps_par)
      
    } else {
      stop("'eps_rule' not supported.")
    }
    
    if (!is.numeric(eps) || length(eps) != length(obs_stats)) {
      stop("The epsilon vector has not the expected format.")
    }
    
  } else {
    eps <- rep(0, length(obs_stats))
  }
  
  eps
}



