#' Create a control object for multiple testing adjustment in permAprox
#'
#' Constructs a list of parameters controlling multiple testing correction
#' methods, so that these options are decoupled from the main \code{permaprox()}
#' function.
#'
#' @inheritParams mult_adjust
#'
#' @return A named list (class \code{"controlMultAdjust"}) with the specified settings.
#' @export
#'
make_ctrl_adjust <- function(
    method = "adaptBH",
    trueNullMethod = "convest",
    pTrueNull = NULL,
    nseq = 100,
    cores = 1,
    verbose = FALSE
) {
  # --- Validation ---
  method <- match.arg(method)
  trueNullMethod <- match.arg(trueNullMethod)

  if (!is.null(pTrueNull)) {
    if (!is.numeric(pTrueNull) || length(pTrueNull) != 1 || is.na(pTrueNull) ||
        pTrueNull < 0 || pTrueNull > 1) {
      stop("'pTrueNull' must be a single numeric value between 0 and 1, or NULL.")
    }
  }

  if (!is.numeric(nseq) || length(nseq) != 1 || nseq < 1) {
    stop("'nseq' must be a positive integer.")
  }
  nseq <- as.integer(nseq)

  if (!is.numeric(cores) || length(cores) != 1 || cores < 1) {
    stop("'cores' must be an integer >= 1.")
  }
  cores <- as.integer(cores)

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value (TRUE or FALSE)." )
  }

  # --- Construct control object ---
  ctrl <- list(
    method = method,
    trueNullMethod = trueNullMethod,
    pTrueNull = pTrueNull,
    nseq = nseq,
    cores = cores,
    verbose = verbose
  )
  class(ctrl) <- "controlMultAdjust"
  return(ctrl)
}
