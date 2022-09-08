# Method of Moments estimator for the GPD

#' @title Method of Moments estimation for the GPD
#'
#' @description
#'   Method of Moments estimation for the two-parameter Generalized Pareto
#'   Distribution (GPD) proposed by \cite{Hosking & Wallis (1987)}.
#'
#' @param x data vector
#'
#' @references
#'   \insertRef{Hosking1987parameter}{permpap}
#'
#' @importFrom Rdpack reprompt
#' @export

gpd_MOM <- function(x) {

  # x must be numeric
  x <- as.numeric(x)

  # Mean and variance
  meanx <- mean(x)
  varx <- var(x)

  shape <- -0.5 * (((meanx^2) / varx) - 1)

  scale <- 0.5 * meanx * (((meanx^2) / varx) + 1)

  return(list(shape = shape, scale = scale))
}
