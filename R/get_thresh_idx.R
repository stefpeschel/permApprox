get_thresh_idx <- function(threshMethod, idxVec, shapeVec,
                           gofPvalVec, gofAlpha, exceedMin,
                           gofTailRMMeth, gofTailRMPar){

  if(all(gofPvalVec <= gofAlpha)){

    idxUse <- length(gofPvalVec)

    return(list(idxUse = idxUse))
  }

  idxH0accept <- which(gofPvalVec > gofAlpha)
  idxH0reject <- which(gofPvalVec <= gofAlpha)

  if(threshMethod == "PRbelowAlpha"){

    # Actual number of iterations
    n <- length(gofPvalVec)

    # Proportion of rejected GOF tests for all thresholds
    propReject <- sapply(1:n, function(i){
      sum(gofPvalVec[i:n] <= gofAlpha) / (n - i + 1)
    })

    # Ensure that H0 is accepted at the chosen threshold
    propReject2 <- propReject
    propReject2[idxH0reject] <- 1

    # Select the first threshold with a PR below alpha
    idxUse <- which(propReject2 <= gofAlpha)[1]

    if (is.na(idxUse)) {
      # Use the minimum if no threshold leads to a PR below alpha
      idxUse <- which.min(propReject2)
    }

    out <- list(idxUse = idxUse, propReject = propReject)

  } else if(threshMethod == "fwdStop"){

    # Log-transformed p-values
    y <- -log(1 - gofPvalVec)

    # Transform log-p-values as described in Barder et. al 2018
    ysum <- sapply(1:length(gofPvalVec), function(k) {
      1 / k * sum(y[1:k])
    })

    if (all(ysum <= gofAlpha)) {
      idxUse <- length(gofPvalVec)

    } else if (all(ysum > gofAlpha)) {
      idxUse <- 1

    } else{
      # Select the first value above alpha
      idxUse <- which(ysum > gofAlpha)[1]
    }

    out <- list(idxUse = idxUse, ysum = ysum)

  } else if (threshMethod == "gofCP") {

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
