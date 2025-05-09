#' @title Adapted version of gpdAd from eva package
#'
#' @keywords internal

.gof_gpd_ad <- function(data, scale, shape, bootstrap = FALSE, bootnum = NULL,
                        allowParallel = FALSE, numCores = 1){

  if (bootstrap == TRUE & is.null(bootnum))
    stop("Must specify some number of boostrap samples")
  n <- length(data)
  # fit <- tryCatch(eva::gpdFit(data, nextremes = n, method = "mle"),
  #                 error = function(w) {
  #                   return(NULL)
  #                 }, warning = function(w) {
  #                   return(NULL)
  #                 })
  # if (is.null(fit))
  #   stop("Maximum likelihood failed to converge at initial step")
  # scale <- fit$par.ests[1]
  # shape <- fit$par.ests[2]
  theta <- c(scale, shape)
  if (bootstrap == FALSE & shape > 1)
    stop("Estimated parameters are outside the table range, please use the bootstrap version")
  thresh <- eva:::findthresh(data, n)
  newdata <- eva::pgpd(data, loc = thresh, scale = scale, shape = shape)
  newdata <- sort(newdata)

  # edit by SP:
  n <- length(newdata)

  i <- seq(1, n, 1)
  stat <- -n - (1/n) * sum((2 * i - 1) * (log(newdata) + log1p(-rev(newdata))))
  if (bootstrap == TRUE) {
    if (allowParallel == TRUE) {
      cl <- makeCluster(numCores)
      fun <- function(cl) {
        parSapply(cl, 1:bootnum, function(i, ...) {
          gpdAdGen(n, theta)
        })
      }
      teststat <- fun(cl)
      stopCluster(cl)
    }
    else {
      teststat <- replicate(bootnum, gpdAdGen(n, theta))
    }
    teststat <- teststat[!is.na(teststat)]
    eff <- length(teststat)
    p <- (sum(teststat > stat) + 1)/(eff + 2)
  }
  else {
    row <- which(rownames(eva:::ADQuantiles) == max(round(shape,
                                                    2), -0.5))
    if (stat > eva:::ADQuantiles[row, 999]) {
      pvals <- -log(as.numeric(colnames(eva:::ADQuantiles[950:999])))
      x <- as.numeric(eva:::ADQuantiles[row, 950:999])
      y <- lm(pvals ~ x)
      stat <- as.data.frame(stat)
      colnames(stat) <- c("x")
      p <- as.numeric(exp(-predict(y, stat)))
    }
    else {
      bound <- as.numeric(colnames(eva:::ADQuantiles)[which.max(stat <
                                                            eva:::ADQuantiles[row, ])])
      if (bound == 0.999) {
        p <- 0.999
      }
      else {
        lower <- eva:::ADQuantiles[row, which(colnames(eva:::ADQuantiles) ==
                                          bound + 0.001)]
        upper <- eva:::ADQuantiles[row, which(colnames(eva:::ADQuantiles) ==
                                          bound)]
        dif <- (upper - stat)/(upper - lower)
        val <- (dif * (-log(bound) - -log(bound + 0.001))) +
          log(bound)
        p <- exp(val)
      }
    }
  }
  names(theta) <- c("Scale", "Shape")
  if (!bootstrap) {
    out <- list(as.numeric(stat), as.numeric(p), theta)
    names(out) <- c("statistic", "p.value", "theta")
  }
  else {
    out <- list(as.numeric(stat), as.numeric(p), theta,
                eff)
    names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  }
  out
}
