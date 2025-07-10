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
#' eps_fixed(n = 100)            # 0.05
#' eps_fixed(n = 250, value = 1) # 1
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
#'   Upper boundary of the GPD support (on the original test-statistic scale).
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

eps_factor <- function(support_boundary,
                       factor = 0.05,
                       ...
                       ) {
  
  if (is.null(support_boundary))
    stop("`support_boundary` must be supplied for `eps_factor()`.")
  
  if (!is.numeric(factor) || length(factor) != 1L || factor <= 0)
    stop("`factor` must be a single positive number.")
  
  factor * support_boundary
}

# ---------------------------------------------------------------------------

#' @title Epsilon rule: inverse power of sample size
#'
#' @description
#' Implements the general power-law
#' \deqn{\varepsilon = \frac{A}{n^{B}}}{epsilon = A / n^B}
#' which captures the empirically observed decrease of the optimal ε with
#' growing sample size.
#'
#' @param n Integer (length-one). Sample size.
#' @param A Numeric scalar \eqn{> 0}.  Numerator (default `7.5e3`).
#' @param B Numeric scalar \eqn{> 0}.  Exponent (default `2`).
#' @param ... Unused, included for forward compatibility.
#'
#' @return A numeric scalar \eqn{\varepsilon}.
#'
#' @examples
#' eps_power(n = 100)                 # 7500 / 100^2 = 0.75
#' eps_power(n = 250, A = 5000, B = 2) # 0.08
#'
#' @export

eps_power <- function(n, A = 7.5e3, B = 2, ...) {
  if (!is.numeric(n) || length(n) != 1L || n <= 0)
    stop("`n` must be a single positive number.")
  
  if (!is.numeric(A) || length(A) != 1L || A <= 0)
    stop("`A` must be a single positive number.")
  
  if (!is.numeric(B) || length(B) != 1L || B <= 0)
    stop("`B` must be a single positive number.")
  
  A / n^B
}
