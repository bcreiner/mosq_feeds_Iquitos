tmp_data <- data.frame(viremia = rep(mean(df_D2$viremia), 5), DoI = 1:5, bin_EIP = rep(2,5))
pred <- predict(mod_gam_linear_bin, newdata = tmp_data, se = TRUE)

A <- exp(pred$fit - pred$fit[1])
B <- exp(pred$fit - 1.96 * pred$se - pred$fit[1])
C <- exp(pred$fit + 1.96 * pred$se - pred$fit[1])

YSEQ <- pretty(range(log10(c(B, C))))
options(scipen = 999)
YLAB <- head(log10(sort(as.vector(outer(10^YSEQ , c(1,5), '*')))), -2)[-c(2,4)]


pdf(file="output/DoI_OR.pdf", height = 6, width = 6)
plot(log10(A), ylim = range(YLAB), type = 'n', xlab = "Day of illness",
     ylab = "Odds ratio relative to Day of illness of 1", axes = FALSE)
box()
abline(h = 0, lty = 3)
axis(1,1:5)
axis(2, YLAB, 10^YLAB)
for (i in 1:5){
  lines(c(i-.1, i+.1), rep(log10(A[i]), 2), lwd = 3)
  if (i > 1){
    lines(rep(i, 2), log10(c(A[i], B[i])), lwd = 2)
    lines(c(i-.05, i+.05), log10(rep(B[i], 2)), lwd = 2)
    lines(rep(i, 2), log10(c(A[i], C[i])), lwd = 2)
    lines(c(i-.05, i+.05), log10(rep(C[i], 2)), lwd = 2)
  }
}
dev.off()