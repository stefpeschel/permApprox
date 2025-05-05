get_thresh_idx <- function(thresh_method, gofPvalVec, gof_alpha){

  if(all(gofPvalVec <= gof_alpha)){

    idxUse <- NA

    return(list(idxUse = idxUse))
  }

  idxH0accept <- which(gofPvalVec > gof_alpha)
  idxH0reject <- which(gofPvalVec <= gof_alpha)

  if(thresh_method == "PRbelowAlpha"){

    # Actual number of iterations
    n <- length(gofPvalVec)

    # Proportion of rejected GOF tests for all thresholds
    propReject <- sapply(1:n, function(i){
      sum(gofPvalVec[i:n] <= gof_alpha) / (n - i + 1)
    })

    # Ensure that H0 is accepted at the chosen threshold
    propReject2 <- propReject
    propReject2[idxH0reject] <- 1

    # Select the first threshold with a PR below alpha
    idxUse <- which(propReject2 <= gof_alpha)[1]

    if (is.na(idxUse)) {
      # Use the minimum if no threshold leads to a PR below alpha
      idxUse <- which.min(propReject2)
    }

    out <- list(idxUse = idxUse, propReject = propReject)

  } else if(thresh_method == "fwdStop"){

    # Log-transformed p-values
    y <- -log(1 - gofPvalVec)

    # Transform log-p-values as described in Barder et. al 2018
    ysum <- sapply(1:length(gofPvalVec), function(k) {
      1 / k * sum(y[1:k])
    })

    if (all(ysum <= gof_alpha)) {
      idxUse <- length(gofPvalVec)

    } else if (all(ysum > gof_alpha)) {
      idxUse <- 1

    } else{
      # Select the first value above alpha
      idxUse <- which(ysum > gof_alpha)[1]
    }

    out <- list(idxUse = idxUse, ysum = ysum)

  } else if (thresh_method == "gofCP") {

    # Add 100 fake p-values (sampled from U(0, 0.01)) to ensure a correct
    # estimate if (nearly) all hypotheses are true
    gofPvalTmp <- c(runif(100, min = 0, max = 0.01), gofPvalVec)

    # Changepoint detection
    cp <- changepoint::cpt.meanvar(gofPvalTmp)
    idxSel <- cp@cpts[1] - 100

    # Find the next larger index for which the AD test is accepted
    idxUse <- idxH0accept[idxH0accept >= idxSel][1]

    out <- list(idxUse = idxUse)

  } else {
    stop("Threshold method not supported.")
  }

  return(out)
}
