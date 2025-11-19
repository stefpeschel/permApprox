#' Create a control list for Gamma approximation in permApprox
#'
#' Constructs a list of control parameters for the Gamma distribution
#' approximation method used in permApprox.
#'
#' @param include_obs Logical. If \code{TRUE}, the observed test statistic is
#'   included into the tail approximation (considered as permutation test
#'   statistic). Default: \code{FALSE}.
#'
#' @param gof_test Character. Goodness-of-fit test for the Gamma distribution.
#'   Options: "cvm" (Cramer-von-Mises), "none". Default: "cvm".
#'
#' @param gof_alpha Numeric. Significance level for GOF test (0 < gof_alpha < 1).
#'   Default: 0.05.
#'   
#' @param fit_args Optional named list of additional arguments passed to
#'   \code{fitdistrplus::fitdist()} when fitting the Gamma distribution.
#'   Entries in \code{fit_args} must not include \code{data} or \code{distr},
#'   which are always set internally. Typical uses are to change the
#'   \code{method} (e.g., "mle", "mme") or provide \code{start} values.
#'
#' @return A named list of class "gamma_ctrl" containing Gamma settings.
#'
#' @export
make_gamma_ctrl <- function(
    include_obs = FALSE,
    gof_test    = c("none", "cvm"),
    gof_alpha   = 0.05,
    fit_args    = NULL
) {
  
  stopifnot(is.logical(include_obs))
  
  gof_test <- match.arg(gof_test)
  
  stopifnot(is.numeric(gof_alpha) & gof_alpha > 0 & gof_alpha < 1)
  
  if (!is.null(fit_args) && !is.list(fit_args)) {
    stop("'fit_args' must be NULL or a list of arguments for ", 
         "fitdistrplus::fitdist().")
  }
  
  control <- list(
    include_obs = include_obs,
    gof_test   = gof_test,
    gof_alpha  = gof_alpha,
    fit_args   = fit_args
  )
  
  class(control) <- "gamma_ctrl"
  
  control
}
