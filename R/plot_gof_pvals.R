plot_gof_pvals <- function(gofPvalVec, gofAlpha, thresh, 
                           threshMethod, propReject = NULL,
                           xvar = "thresh",
                           idxVec = NULL, threshVec = NULL,
                           title = NULL,
                           mar = NULL, cexaxis = NULL,
                           cexlab = NULL, cexmain = NULL, cexlegend = NULL,
                           plotPropReject = TRUE){
  oldpar <- par()
  
  if(is.null(mar)){
    mar <- c(5, 5, 3, 2)
  }
  par(mar = mar)
  
  if(is.null(title)){
    title <- "GOF p-values"
  }
  
  if(xvar == "thresh"){
    xvar.plot <- threshVec
  } else{
    xvar.plot <- idxVec
  }
  
  plot(gofPvalVec ~ xvar.plot, type = "l", ylim = c(0,1), main = title,
       cex.axis = cexaxis, 
       cex.main = cexmain, cex.lab = cexlab)
  
  abline(h = gofAlpha, col = "red")
  
  if(xvar == "thresh"){
    abline(v = thresh, col = "green")
  } else{
    abline(v = which(threshVec == thresh), col = "green")
  }
  
  if(plotPropReject && threshMethod %in% c("minPR", "PRbelowAlpha")){
    lines(propReject ~ xvar.plot[1:length(propReject)], col = "blue")
  }
  
  lines(gofPvalVec ~ xvar.plot)
  
  if(is.null(cexlegend)) cexlegend <- 1
  
  if(threshMethod %in% c("minPR", "PRbelowAlpha")){
    legend("topleft", 
           legend = c("gofAlpha", "threshold", "prop. of rejections"),
           col = c("red", "green", "blue"), lty = 1, cex = cexlegend)
  } else{
    legend("topleft", 
           legend = c("gofAlpha", "threshold"),
           col = c("red", "green"), lty = 1, cex = cexlegend)
  }
  
  par(mar = oldpar$mar)
}