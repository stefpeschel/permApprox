#' Create a control object for multiple testing adjustment in permApprox
#'
#' Constructs a list of parameters controlling multiple testing correction
#' methods, so that these options are decoupled from the main
#' \code{compute_p_values()} function.
#'
#' @inheritParams mult_adjust
#'
#' @return A named list (class \code{"adjust_ctrl"}) with the specified settings.
#' @export
#'
make_adjust_ctrl <- function(
    true_null_method = "convest",
    p_true_null = NULL
) {
  
  true_null_method <- match.arg(true_null_method, c("farco", "lfdr", "mean",
                                                    "hist", "convest"))
  
  if (!is.null(p_true_null)) {
    if (!is.numeric(p_true_null) || length(p_true_null) != 1 || 
        is.na(p_true_null) || p_true_null < 0 || p_true_null > 1) {
      stop("'p_true_null' must be a single numeric value between 0 and 1, or NULL.")
    }
  }
  
  control <- list(
    true_null_method = true_null_method,
    p_true_null = p_true_null
  )
  
  class(control) <- "adjust_ctrl"
  return(control)
}
