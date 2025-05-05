# Histogram of Gamma distribution


# if(plotHist){
#   par(mar = c(5, 5, 7, 4))
#   xseq <- seq(0, obs_stats+0.1, 0.01)
#
#   if(is.null(histXlim)) histXlim <- c(0, obs_stats+0.1)
#
#   h <- hist(perm_stats, breaks = histBreaks,
#             probability = TRUE, xlim = histXlim, ylim = histYlim,
#             main = paste("\n Gamma shape param:", round(shape, 7),
#                          "\n Gamma rate param:", round(rate, 7),
#                          "\n CvM test p-value:", round(gofPval, 7),
#                          "\n permutation p-value:", pval,
#                          maintext),
#             col = histCol,
#             cex.axis = 1.7,
#             cex.main = 2, cex.lab = 2)
#
#   lines(xseq, dgamma(xseq, shape = shape, rate = rate), col = "blue")
#   abline(v = obs_stats, col = "red")
#   legend(obs_stats*0.7, max(h$density)*0.9,
#          legend = c("fitted Gamma distr.", "obs.teststat"),
#          col = c("blue", "red"), lty = 1, cex = 2)
#
# }
