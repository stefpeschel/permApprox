get_thresh_idx <- function(threshMethod, idxVec, shapeVec,
                         gofPvalVec, gofAlpha, nExceedMin,
                         gofTailRMMeth, gofTailRMPar){

  if(all(gofPvalVec <= gofAlpha)){
    idxUse <- NULL
    out <- list(idxUse = idxUse)

  } else{
    idxH0accept <- which(gofPvalVec > gofAlpha)
    idxH0reject <- which(gofPvalVec <= gofAlpha)

    if(threshMethod %in% c("minPR", "PRbelowAlpha")){

      # Actual number of iterations
      n <- length(idxVec)

      if(gofTailRMMeth != "none"){

        if(gofTailRMMeth == "num"){
          torem <- gofTailRMPar

        } else if(gofTailRMMeth == "prop"){
          torem <- floor(length(gofPvalVec) * gofTailRMPar)

        } else if(gofTailRMMeth == "allrej"){
          idxSplit <- split(idxH0reject, cumsum(c(1, diff(idxH0reject) != 1)))
          torem <- length(idxSplit[[length(idxSplit)]])

        } else if(gofTailRMMeth == "changepoint"){
          # replace NAs by value before
          if(is.na(shapeVec[1])) shapeVec[1] <- 0
          for(i in 2:length(shapeVec)){
            if(is.na(shapeVec[i])) shapeVec[i] <- shapeVec[i-1]
          }
          cp.est <- changepoint::cpt.meanvar(shapeVec, penalty="SIC", method="AMOC")
          torem <- n - cp.est@cpts[1] + 1
        }

        gofPvalVec.sel <- gofPvalVec[1:(length(gofPvalVec) - torem)]
        rmTailIdx <- length(gofPvalVec) - torem + 1

      } else{
        torem <- 0
        startidx <- NA
        gofPvalVec.sel <- gofPvalVec
        rmTailIdx <- NA
      }

      # Proportion of removed tail
      propTailRem <- torem / length(gofPvalVec)

      # Proportion of rejected GOF tests
      n.sel <- length(gofPvalVec.sel)

      propReject <- sapply(1:n.sel, function(i){
        sum(gofPvalVec.sel[i:n.sel] <= gofAlpha) / (n.sel - i + 1)
      })

      # Ensure that H0 is accepted at the chosen threshold
      propReject2 <- propReject
      propReject2[idxH0reject] <- 1

      # Choose index for threshold
      if(threshMethod == "PRbelowAlpha"){
        if(!any(propReject2 <= gofAlpha)){
          idxUse <- NULL
        } else{
          idxUse <- which(propReject2 <= gofAlpha)[1]
        }

      } else{
        idxUse <- which.min(propReject2)
      }

      out <- list(idxUse = idxUse, propReject = propReject, tailRem = torem,
                  propTailRem = propTailRem, rmTailIdx = rmTailIdx)

    } else if(threshMethod == "ftr"){
      # Failure-to-reject method
      idxUse <- idxH0accept[1]

      out <- list(idxUse = idxUse)

    } else if(threshMethod == "fwdStop"){
      pvalssort <- sort(gofPvalVec)
      logpvals <- log(1-pvalssort)

      transSum <- sapply(1:length(pvalssort), function(k){
        -1/k * sum(logpvals[1:k])
      })

      if(all(transSum <= gofAlpha)){
        idxUse <- NULL
      } else if(all(transSum > gofAlpha)){
        idxUse <- 1
      } else{
        belowAlpha <- which(transSum <= gofAlpha)
        idxUseSort <- belowAlpha[length(belowAlpha)] + 1
        idxUse <- which(gofPvalVec == pvalssort[idxUseSort])
      }

      out <- list(idxUse = idxUse)

    }
  }

  return(out)
}
