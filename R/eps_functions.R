#' @title Epsilon rule: fixed constant
#'
#' @description
#' Returns the same non-negative constant for every call.  Useful when you want to
#' supply a *single, non–data-dependent* value of ε.
#'
#' @param value Numeric scalar \eqn{>= 0}. The constant \eqn{\varepsilon}
#'   to return.  Default is `0.05`.
#' @param ... Unused, included for forward compatibility.
#'
#' @return A numeric scalar equal to `value`.
#'
#' @examples
#' eps_fixed()          # 0.05
#' eps_fixed(value = 1) # 1
#'
#' @export

eps_fixed <- function(value = 0.05, ...) {
  if (!is.numeric(value) || length(value) != 1L || value < 0)
    stop("`value` must be a single non-negative number.")
  value
}

# ---------------------------------------------------------------------------

#' @title Epsilon rule: factor of support boundary
#'
#' @description
#' Computes \eqn{\varepsilon} = \code{factor} x \code{support_boundary}.
#' It is intended for constrained GPD fits where the support boundary is known
#' (or fixed to the maximum observed statistic).
#'
#' @param support_boundary Numeric scalar \eqn{>= 0}.
#'   Upper boundary of the GPD support.
#' @param factor Numeric scalar \eqn{>= 0}.  Multiplier applied to
#'   `support_boundary`.  Default is `0.05`.
#' @param ... Unused, included for forward compatibility.
#'
#' @details
#' If `support_boundary` is `NULL`, the function stops with an error so the user
#' realises a value is required.
#'
#' @return A numeric scalar \eqn{\varepsilon}.
#'
#' @examples
#' eps_factor(support_boundary = 3.2)               # 0.16
#' eps_factor(support_boundary = 3.2, factor = 0.1) # 0.32
#'
#' @export

eps_factor <- function(support_boundary = NULL,
                       factor = 0.05,
                       ...
                       ) {
  
  if (is.null(support_boundary))
    stop("`support_boundary` must be supplied for `eps_factor()`.")
  
  if (!is.numeric(factor) || length(factor) != 1L || factor < 0)
    stop("`factor` must be a single non-negative number.")
  
  factor * support_boundary
}

# ---------------------------------------------------------------------------

#' @title Epsilon rule: inverse power of sample size
#'
#' @description
#' Implements the general power-law
#' \deqn{\varepsilon = \frac{A}{n^{B}} \times x_{bound},}
#' which captures the empirically observed decrease of the optimal ε with
#' growing sample size.
#'
#' @param support_boundary Numeric scalar \eqn{>= 0}.
#'   Upper boundary of the GPD support.
#' @param n Integer (length-one). Sample size.
#' @param A Numeric scalar \eqn{> 0}.  Numerator (default `7.5e3`).
#' @param B Numeric scalar \eqn{> 0}.  Exponent (default `2`).
#' @param ... Unused, included for forward compatibility.
#'
#' @return A numeric scalar \eqn{\varepsilon}.
#'
#' @examples
#' eps_power(n = 100)                  # 10000/100^2 = 1
#' eps_power(n = 250, A = 5000, B = 2) # 0.08
#'
#' @export

eps_power <- function(support_boundary = NULL, n = NULL, A = 1e4, B = 2, ...) {
  if (is.null(support_boundary))
    stop("`support_boundary` must be supplied for `eps_factor()`.")
  
  if (is.null(n))
    stop("Sample size `n` must be given for this epsilon rule.")
  
  if (!is.numeric(n) || length(n) != 1L || n <= 0)
    stop("`n` must be a single positive number.")
  
  if (!is.numeric(A) || length(A) != 1L || A <= 0)
    stop("`A` must be a single positive number.")
  
  if (!is.numeric(B) || length(B) != 1L || B <= 0)
    stop("`B` must be a single positive number.")
  
  support_boundary * (A / n^B)
  
}



