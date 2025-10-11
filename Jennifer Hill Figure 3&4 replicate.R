load("~/Desktop/phd/HBB project/BART Replication/code/.RData")
library(ggplot2)
library(reshape2)
library(knitr)

ls()
ls(pattern = "tau")
table_ATE <- data.frame(
  Model = c("BART-A", "BART-B", "BART-C"),
  Mean_ATE = c(mean(tauAs), mean(tauBs), mean(tauCs)),
  L95 = c(quantile(tauAs, 0.025), quantile(tauBs, 0.025), quantile(tauCs, 0.025)),
  U95 = c(quantile(tauAs, 0.975), quantile(tauBs, 0.975), quantile(tauCs, 0.975))
)
table_ATE


## --------------------------------------------------------------
## Hill (2011) Figure 3/4 replication
##  Scatter of deviation from truth vs. interval length
## --------------------------------------------------------------

plot_hill_centered <- function(resA, resB = NULL,
                               title_text = "Hill (2011) Real-world Simulation Replication") {
  
  draw_row <- function(res, row_label = "Complete Overlap") {
    methods <- list(
      "BART"              = list(te = "b.te",  cil = "b.cil"),
      "Linear Regression" = list(te = "r.te",  cil = "r.cil"),
      "Propensity Score Matching" = list(te = "ps.te", cil = "ps.cil"),
      "Propensity Weighted"       = list(te = "ipw.te", cil = "ipw.cil")
    )
    

    tau_true <- mean(res[, "tau.est"])
    tau_centered <- res[, "tau.est"] - tau_true
    
    for (m in names(methods)) {
      te  <- res[, methods[[m]]$te]
      cil <- res[, methods[[m]]$cil]
      deviation <- te - tau_centered
      
      rmse   <- sqrt(mean((deviation)^2))
      cover  <- mean(res[, gsub("\\.te", ".cov", methods[[m]]$te)], na.rm = TRUE)
      
      plot(cil, deviation,
           pch = 20, cex = 0.6,
           xlab = "Interval Length", ylab = "Deviation from truth",
           main = m,
           ylim = c(-0.4, 0.4),
           xlim = range(cil) * c(0.95, 1.05))
      abline(h = 0, lty = 2, lwd = 1.5)           # true = 0
      abline(a = 0, b = 1, lty = 3, col = "gray60")
      abline(a = 0, b = -1, lty = 3, col = "gray60")
      
      legend("bottomleft",
             legend = c(sprintf("RMSE %.2f", rmse),
                        sprintf("Coverage %.2f", cover)),
             bty = "n", cex = 0.8)
    }
    mtext(row_label, outer = FALSE, side = 2, line = 1.5, cex = 1)
  }
  
  if (!is.null(resB)) {
    par(mfrow = c(2, 4), mar = c(4.5, 4.5, 3, 1), oma = c(1, 0, 3, 0))
  } else {
    par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3, 1), oma = c(1, 0, 3, 0))
  }
  
  draw_row(resA, "Complete Overlap")
  if (!is.null(resB)) draw_row(resB, "Incomplete Overlap")
  
  mtext(title_text, outer = TRUE, side = 3, line = 0.5, cex = 1.2)
}


plot_hill_centered(results.a, results.b)

png("Hill2011_Figure3_4_centered.png", width = 1600, height = 900, res = 150)
plot_hill_centered(results.a, results.b)
dev.off()
