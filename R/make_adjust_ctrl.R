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
    p_true_null = NULL,
    seq_length = 100,
    cores = 1,
    verbose = FALSE
) {

  true_null_method <- match.arg(true_null_method, c("farco", "lfdr", "mean",
                                                    "hist", "convest"))

  if (!is.null(p_true_null)) {
    if (!is.numeric(p_true_null) || length(p_true_null) != 1 || is.na(p_true_null) ||
        p_true_null < 0 || p_true_null > 1) {
      stop("'p_true_null' must be a single numeric value between 0 and 1, or NULL.")
    }
  }

  if (!is.numeric(seq_length) || length(seq_length) != 1 || seq_length < 1) {
    stop("'seq_length' must be a positive integer.")
  }
  seq_length <- as.integer(seq_length)

  if (!is.numeric(cores) || length(cores) != 1 || cores < 1) {
    stop("'cores' must be an integer >= 1.")
  }
  cores <- as.integer(cores)

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value (TRUE or FALSE)." )
  }

  control <- list(
    true_null_method = true_null_method,
    p_true_null = p_true_null,
    seq_length = seq_length,
    cores = cores,
    verbose = verbose
  )
  class(control) <- "adjust_ctrl"
  return(control)
}
