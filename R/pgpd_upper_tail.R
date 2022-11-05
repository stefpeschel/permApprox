#' @title Quantile function for the upper tail of the GPD
#'
#' @description
#' Enables the computation of very small p-values.
#'
#' @param q quantile
#' @param location location parameter of the GPD
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @export

pgpd_upper_tail <- function(q, location = 0, shape, scale){
  zedd <- (q - location)/scale
  use.zedd <- pmax(zedd, 0)
  pmax(1 + shape * use.zedd, 0)^(-1/shape)
}
