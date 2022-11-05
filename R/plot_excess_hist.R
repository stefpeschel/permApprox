plot_excess_hist <- function(excessPerm, excessObs,
                             nExceed, shape, scale, gofPval, pval, idx = NULL,
                             breaks = 100,
                             maintext = "", main = NULL,
                             mar = NULL, xlim = NULL, ylim = NULL,
                             col = NULL, cexaxis = NULL,
                             cexlab = NULL, cexmain = NULL, cexlegend = NULL){
  
  if(!is.null(idx)){
    excessPerm <- excessPerm[idx]
    excessObs <- excessObs[idx]
    nExceed <- nExceed[idx]
    shape <- shape[idx]
    scale <- scale[idx]
    gofPval <- gofPval[idx]
    pval <- pval[idx]
  }
  
  # Set graphic parameters
  oldpar <- par()
  
  if(is.null(mar)){
    mar <- c(5, 5, 10, 4)
  } 
  par(mar = mar)
  
  if(is.null(main)){
    main <- paste("\n number of exceedances:", nExceed,
                   "\n GPD shape param:", round(shape, 7),
                   "\n GPD scale param:", round(scale, 7),
                   "\n GOF test p-value:", round(gofPval, 7), 
                   "\n permutation p-value:", pval,
                   maintext)
  }
  
  if(is.null(col)) col <- "lightgray"
  if(is.null(cexlegend)) cexlegend <- 1
  
  upperlim <- max(excessPerm, excessObs)*1.1
  xseq <- seq(0, upperlim, 0.01)
  
  if(is.null(xlim)) xlim <- c(0, upperlim)
  
  # Plot histogram
  h <- hist(excessPerm, breaks = breaks,
            probability = TRUE, xlim = xlim, ylim = ylim,
            main = main, col = col, cex.axis = cexaxis, 
            cex.main = cexmain, cex.lab = cexlab)
  
  lines(xseq, VGAM::dgpd(xseq, scale = scale, shape = shape), col = "blue")
  
  abline(v = excessObs, col = "red")
  
  legend(excessObs*0.7, max(h$density)*0.9, 
         legend = c("fitted GPD", "obs.teststat"),
         col = c("blue", "red"), lty = 1, cex = cexlegend)
  
  par(mar = oldpar$mar)
}