# =====================================================================
# 30_figure2_table2.r
# Figure 2 (12 X-COP clusters, gas-only HMG) and Table 2
# (fitted r_nei, chi^2, p-value, significance).
# Uses Einasto mass profiles (EIN3_mass.fits) from Eckert et al. (2022).
# =====================================================================

sdir <- tryCatch({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else getwd()
}, error = function(e) getwd())
source(file.path(sdir, "00_hmg_common.r"))
suppressMessages(library(FITSio))
P <- hmg_paths()

dir_ein <- file.path(P$data, "einasto_mass_profiles")

## Cluster display order matching the paper (A85, A644, A1644, ...)
clu_order <- c(10, 9, 1:8, 11, 12)
clusters  <- list.files(dir_ein)[clu_order]
R500      <- c(1214, 1398, 1031, 1160, 1340, 1453, 1202, 1424, 1146, 1381, 1109, 1346)

# Read the Einasto mass profile for one cluster and derive the stellar
# mass column as total minus missing (dark) minus gas.
read_cluster <- function(clust) {
  f  <- file.path(dir_ein, clust, list.files(file.path(dir_ein, clust), pattern = ".fits"))
  d  <- as.data.frame(readFITS(f[1])$col)
  names(d) <- readFITS(f[1])$colNames
  d$M_STAR <- d$MASS - d$M_DM - d$MGAS
  d
}

## ---- Fit r_nei by scanning R0 and minimising chi^2 -------------------
R0s <- seq(200, 2000, 10)
R0_best <- R1_best <- R2_best <- rep(NA, length(clusters))
chis1 <- chis2 <- chis1.p <- chis2.p <- rep(NA, length(clusters))
chimond <- chimond.p <- rep(NA, length(clusters))
nbins_R1000 <- rep(NA, length(clusters))
nbins_R500  <- rep(NA, length(clusters))

for (iclu in seq_along(clusters)) {
  d <- read_cluster(clusters[iclu])
  M_gas_kg <- Msol * d$MGAS

  ## observed inward gravity = missing mass + gas (Eckert dynamical)
  obs    <- GN*Msol*(d$M_DM    + d$MGAS)    / (d$R_REF*kpc)^2
  obs_HI <- GN*Msol*(d$M_DM_HI + d$MGAS_HI) / (d$R_IN *kpc)^2
  obs_LO <- GN*Msol*(d$M_DM_LO + d$MGAS_LO) / (d$R_OUT*kpc)^2
  obs_err <- (obs_HI - obs_LO) / 2

  # Three outer-radius windows used to fit r_nei; R>=1000 (lgR2) defines Table 2
  lgR0 <- d$R_REF >= 1500
  lgR1 <- d$R_REF >= 500
  lgR2 <- d$R_REF >= 1000
  nbins_R1000[iclu] <- sum(lgR2)
  nbins_R500[iclu]  <- sum(lgR1)

  chi00 <- chi01 <- chi02 <- rep(NA, length(R0s))
  aMOND2_last <- NULL
  for (iir in seq_along(R0s)) {
    pr <- hmg_predict(M_gas_kg, d$R_REF, R0s[iir], r_grav = 50)
    ap <- pr$a_pred1
    chi00[iir] <- sum((ap-obs)[lgR0]^2/obs_err[lgR0]^2, na.rm = TRUE)
    chi01[iir] <- sum((ap-obs)[lgR1]^2/obs_err[lgR1]^2, na.rm = TRUE)
    chi02[iir] <- sum((ap-obs)[lgR2]^2/obs_err[lgR2]^2, na.rm = TRUE)
    aMOND2_last <- pr$a_MOND2
  }
  chimond[iclu]   <- sum((aMOND2_last-obs)[lgR1]^2/obs_err[lgR1]^2, na.rm = TRUE)
  chimond.p[iclu] <- round(100*pchisq(chimond[iclu], df = sum(lgR1)-1))

  R0_best[iclu] <- R0s[order(chi00)[1]]
  R1_best[iclu] <- R0s[order(chi01)[1]]
  R2_best[iclu] <- R0s[order(chi02)[1]]

  chis1[iclu]   <- min(chi01)
  chis2[iclu]   <- min(chi02)
  chis1.p[iclu] <- round(100*pchisq(chis1[iclu], df = sum(lgR1)-1))
  chis2.p[iclu] <- round(100*pchisq(chis2[iclu], df = sum(lgR2)-1))
}

## Paper convention: for clusters whose chi^2 CDF rounds to 100% (a perfect or
## over-fit), the reported r_nei is the fitted value + 200 kpc; chi^2 is not
## recomputed.
R2_best[chis2.p == 100] <- R2_best[chis2.p == 100] + 200
## r_nei uncertainty (empirical window-sensitivity spread, not a formal
## confidence interval): sqrt(10^2 + sample standard deviation, over the three
## outer-radius windows (R>=500, >=1000, >=1500 kpc), of the best-fit r_nei^2),
## with the 10 kpc term set by the R0 scan grid step. This variant matches the
## paper's tabulated uncertainties. RXC1825 is an outlier relative to the
## paper's quoted +/-110.
R_err <- round(sqrt(10^2 + apply(cbind(R2_best, R1_best, R0_best)^2, 1, sd)))

## Degrees of freedom in Table 2:
##  - the tabulated chi^2 (chis2) is computed over the R>=1000 kpc window
##    (per-cluster ~22-31 bins);
##  - the chi2_cdf_pct column is the lower-tail chi^2 CDF percentile (pchisq)
##    with per-cluster dof, df = n_i - 1 (sum(lgR2)-1). This is the paper's
##    convention (high value = rejection), not a standard upper-tail
##    goodness-of-fit p-value;
##  - the 95% significance flag compares chi^2 to the fixed threshold
##    chi^2(95%, df=42) = 58.12, matching the paper's n=43 (R>=500 kpc)
##    convention.
signif95 <- ifelse(chis2 < CHI2_THRESHOLD_95, "*", "-")
pfmt <- function(p) ifelse(p < 1, "<1", ifelse(p > 99, ">99", as.character(p)))

table2 <- data.frame(
  Cluster      = clusters,
  r_nei        = R2_best,
  r_nei_err    = R_err,
  chi2         = round(chis2, 1),
  chi2_cdf_pct = pfmt(chis2.p),
  signif95     = signif95,
  stringsAsFactors = FALSE
)

cat("\n==== TABLE 2 ====\n")
print(table2, row.names = FALSE)
cat(sprintf("\nchi^2 threshold (95%%, df=42) = %.2f\n", CHI2_THRESHOLD_95))
## MOND (MLS interpolation, no EFE) evaluated per cluster over the R>=500 kpc
## window. The smallest per-cluster chi^2 is MOND's best case across the sample;
## even this minimum lies far above the 95% threshold.
imin <- which.min(chimond)
cat(sprintf("MOND best-case (minimum) per-cluster chi^2 = %.0f (%s, n=%d bins, R>=500 kpc)\n",
            chimond[imin], clusters[imin], nbins_R500[imin]))

dir.create(file.path(P$root, "tables"), showWarnings = FALSE)
write.csv(table2, file.path(P$root, "tables", "table2_hmg_fit.csv"),
          row.names = FALSE)

## ---- Figure 2 -------------------------------------------------------
plot_fig2 <- function(dev_open) {
  dev_open()
  par(mfrow = c(3,4), mar = c(0.1,0.1,0.1,0.1), oma = c(7.9,3.1,3.1,3.1))
  for (iclu in seq_along(clusters)) {
    d <- read_cluster(clusters[iclu])
    M_gas_kg <- Msol * d$MGAS
    pr <- hmg_predict(M_gas_kg, d$R_REF, R2_best[iclu], r_grav = 50)

    plot(d$R_REF, log10(GN*Msol*d$MASS/(d$R_REF*kpc)^2),
         xlim = c(200,3000), ylim = c(-11.5,-9), col = "white", log = "x",
         axes = FALSE, xlab = "", ylab = "")
    rect(0.1, -14, 10000, -7, col = rgb(1,1,0.95,0.5)); box()
    if (iclu %in% 1:4)   axis(3, c(1,2,5,10,20,50,100,200,300,500,750,1000,1500,2000,2500), padj = 0.8)
    if (iclu %in% 9:12)  axis(1, c(1,2,5,10,20,50,100,200,300,500,750,1000,1500,2000,2500), padj = -0.8)
    if (iclu %in% c(1,5,9))  axis(2, seq(-13,-9,1), padj = 0.8)
    if (iclu %in% c(4,8,12)) axis(4, seq(-13,-9,1), padj = -0.8)
    if (iclu %in% 1:4)   mtext(side=3,"Radius [kpc]",line=1.5)
    if (iclu %in% 9:12)  mtext(side=1,"Radius [kpc]",line=1.7)
    if (iclu %in% c(1,5,9))  mtext(side=2,expression(paste(log[10],"(",a," / ",ms^{-2},")")),line=1.5)
    if (iclu %in% c(4,8,12)) mtext(side=4,expression(paste(log[10],"(",a," / ",ms^{-2},")")),line=1.9)
    abline(v=c(1,2,5,10,20,50,100,200,400,500,750,1000,1500), h=seq(-13,-8,0.5), lty=2, col="gray95")
    abline(v = R2_best[iclu], lty = 2, col = "gray80")

    lines(d$R_REF, log10(GN*Msol*d$MGAS/(d$R_REF*kpc)^2), col="indianred", lwd=2)
    lines(d$R_REF, log10(GN*Msol*d$M_DM/(d$R_REF*kpc)^2), col="gray90", lwd=5)
    lines(d$R_REF, log10(GN*Msol*(d$M_DM+d$MGAS)/(d$R_REF*kpc)^2), col="gray75", lwd=3)
    lines(d$R_REF, log10(pr$a_MOND),  col="green3",     lwd=1.5, lty=2)
    lines(d$R_REF, log10(pr$a_pred_s),col="darkviolet", lwd=1.5, lty=2)
    lines(d$R_REF, log10(pr$a_pred0), col="darkviolet", lwd=0.8, lty=3)
    lines(d$R_REF, log10(pr$a_pred1), col="darkviolet", lwd=3)
    text(1720, -9.2, clusters[iclu], cex = 1.15)
  }
  par(fig=c(0,1,0,0.2), mar=c(0.3,1.2,0.6,0.3), oma=c(0,0,0,0), new=TRUE)
  plot.new()
  legend("bottomleft", ncol=1, bty="n", cex=0.95,
         legend=c("Hot-gas baryonic gravity (Eckert+2022)",
                  "Estimated missing-mass gravity (Eckert+2022)",
                  "Estimated missing-mass + baryonic gravity (Eckert+2022)"),
         lty=1, lwd=c(2,5,3), col=c("indianred","gray90","gray75"))
  legend("bottomright", ncol=1, bty="n", cex=0.95,
         legend=c("MOND model for gas + MLS interpolation function (no EFE)",
                  "HMG model for gas with only spatial contribution",
                  "HMG model for gas with also time-like contribution",
                  expression(paste("HMG model with Hubble-Newton equilibium (",gamma[sys]==pi/3,")"))),
         lty=c(2,2,1,3), lwd=c(1.6,1.6,3,1),
         col=c("green3","darkviolet","darkviolet","darkviolet"))
  dev.off()
}

plot_fig2(function() pdf(file.path(P$fig, "Figure2.pdf"), width=7, height=7))
plot_fig2(function() png(file.path(P$fig, "Figure2.png"), width=7, height=7, units="in", res=150))
cat("Figure2 written to figures/.\n")
