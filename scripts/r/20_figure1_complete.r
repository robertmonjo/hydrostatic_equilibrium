# =====================================================================
# 20_figure1_complete.r
# Figure 1: acceleration profiles for the five X-COP clusters with
# complete stellar-mass information (A1795, A2029, A2142, A2319, A644).
# Baryonic (star / gas / total), missing mass, total mass, MOND, HMG.
# Data: EIN3_bar_mass.fits + gbar_gobs_*.fits (Eckert et al. 2022).
# =====================================================================

sdir <- tryCatch({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else getwd()
}, error = function(e) getwd())
source(file.path(sdir, "00_hmg_common.r"))
suppressMessages(library(FITSio))
P <- hmg_paths()

dir_e22  <- file.path(P$data, "data_E22")
clusters <- list.files(dir_e22)     # A1795 A2029 A2142 A2319 A644

# Read the two FITS tables for one cluster:
#   d1 = EIN3_bar_mass (baryonic + total mass profiles)
#   d2 = gbar_gobs     (observed inward acceleration with error band)
read_bar <- function(clust) {
  f  <- file.path(dir_e22, clust, list.files(file.path(dir_e22, clust), pattern = ".fits"))
  d1 <- as.data.frame(readFITS(f[1])$col); names(d1) <- readFITS(f[1])$colNames  # bar_mass
  d2 <- as.data.frame(readFITS(f[2])$col); names(d2) <- readFITS(f[2])$colNames  # gbar_gobs
  list(d1 = d1, d2 = d2)
}

## ---- Fit r_nei: scan R0 and keep the value minimising chi^2 (R>=1500)
R0s <- seq(50, 2000, 10)
R22_best <- rep(NA, length(clusters))
chiss2   <- chiss2.p <- rep(NA, length(clusters))
for (iclu in seq_along(clusters)) {
  dat <- read_bar(clusters[iclu]); d1 <- dat$d1
  M_kg <- Msol*(d1$M_STAR + d1$MGAS)               # enclosed baryonic mass [kg]
  # Observed inward acceleration and its half-width from the mass band
  obs    <- GN*Msol*d1$MASS    / (d1$R_REF*kpc)^2
  obs_HI <- GN*Msol*d1$MASS_HI / (d1$R_IN *kpc)^2
  obs_LO <- GN*Msol*d1$MASS_LO / (d1$R_OUT*kpc)^2
  obs_err <- (obs_HI - obs_LO)/2
  lgR2 <- d1$R_REF >= 1500
  chi002 <- sapply(R0s, function(R0) {
    ap <- hmg_predict(M_kg, d1$R_REF, R0, r_grav = 40)$a_pred1
    sum((ap-obs)[lgR2]^2/obs_err[lgR2]^2, na.rm = TRUE)
  })
  R22_best[iclu] <- R0s[order(chi002)[1]]
  chiss2[iclu]   <- min(chi002)
  chiss2.p[iclu] <- round(100*pchisq(min(chi002), df = sum(lgR2)-1))
}

## ---- Figure 1 -------------------------------------------------------
plot_fig1 <- function(dev_open) {
  dev_open()
  par(mfrow = c(3,2), mar = c(0.1,0.1,0.1,0.1), oma = c(3.1,3.2,3.1,3.2))
  for (iclu in seq_along(clusters)) {
    dat <- read_bar(clusters[iclu]); d1 <- dat$d1; d2 <- dat$d2
    M_kg <- Msol*(d1$M_STAR + d1$MGAS)
    R0 <- R22_best[iclu]; if (R0 < 400) R0 <- R0 + 200
    pr <- hmg_predict(M_kg, d1$R_REF, R0, r_grav = 40)

    plot(d1$R_REF, log10(GN*Msol*d1$MASS/(d1$R_REF*kpc)^2), ylim = c(-11.5,-8.5),
         col = "white", log = "x", axes = FALSE, xlim = c(1,3000), xlab = "", ylab = "")
    rect(0.1, -14, 10000, -7, col = rgb(1,1,0.95,0.5)); box()
    if (iclu %in% 1:2)     axis(3, c(1,2,5,10,20,50,100,200,500,1000,2000), padj = 0.8)
    if (iclu %in% c(4,5))  axis(1, c(1,2,5,10,20,50,100,200,500,1000,2000), padj = -0.8)
    if (iclu %in% c(1,3,5))axis(2, seq(-13,-9,1), padj = 0.8)
    if (iclu %in% c(2,4))  axis(4, seq(-13,-9,1), padj = -0.8)
    if (iclu %in% 1:2)     mtext(side=3,"Radius [kpc]",line=1.5)
    if (iclu %in% c(4,5))  mtext(side=1,"Radius [kpc]",line=1.7)
    if (iclu %in% c(1,3,5))mtext(side=2,expression(paste(log[10],"(",a," / ",ms^{-2},")")),line=1.5)
    if (iclu %in% c(2,4))  mtext(side=4,expression(paste(log[10],"(",a," / ",ms^{-2},")")),line=1.9)
    abline(v=c(1,2,5,10,20,50,100,200,500,1000), h=seq(-13,-8,0.5), lty=2, col="gray95")
    abline(v = R0, lty = 2, col = "gray82")

    lines(d1$R_REF, log10(GN*Msol*d1$M_DM/(d1$R_REF*kpc)^2), col="gray88", lwd=5)
    lines(d1$R_REF, log10(GN*Msol*d1$MGAS/(d1$R_REF*kpc)^2), col="indianred", lwd=3)
    lines(d1$R_REF, log10(GN*Msol*d1$M_STAR/(d1$R_REF*kpc)^2), col="cyan3", lwd=2)
    lines(d1$R_REF, log10(GN*Msol*(d1$M_STAR+d1$MGAS)/(d1$R_REF*kpc)^2), col="blue2", lwd=1)
    lines(d1$R_REF, log10(GN*Msol*d1$MASS/(d1$R_REF*kpc)^2), col="gray60", lwd=2)
    lines(d1$R_REF, log10(pr$a_MOND),   col="green3",     lwd=1.5, lty=2)
    lines(d1$R_REF, log10(pr$a_pred_s), col="darkviolet", lwd=1.5, lty=2)
    lines(d1$R_REF, log10(pr$a_pred0),  col="darkviolet", lwd=0.8, lty=3)
    lines(d1$R_REF, log10(pr$a_pred1),  col="darkviolet", lwd=3)
    points(d2$RADIUS, log10(d2$GOBS))
    segments(d2$RADIUS, log10(d2$GOBS-d2$GOBS_ERR_LO),
             d2$RADIUS, log10(d2$GOBS+d2$GOBS_ERR_HI))
    text(1500, -8.7, clusters[iclu], cex = 1.15)
  }
  par(mar=c(0.3,1.2,0.6,0.3), oma=c(0,0,0,0), new=TRUE); plot.new()
  legend("left", ncol=1, bty="n", pch=1,
    legend=c("Inward gravity for hydrostatic equilibium  ",
             "Hot-gas baryonic gravity",
             "Stars and dust baryonic gravity (accounted to date)",
             "Total baryonic (Newtonian) gravity",
             "Estimated missing mass",
             "Total estimated mass: baryonic + missing mass",
             "MOND model + MLS interpolation function (no EFE)",
             "HMG model with only spatial contribution",
             "HMG model with also time-like contribution",
             expression(paste("HMG model with Hubble-Newton equilibium (",gamma[sys]==pi/3,")"))),
    lty=c(1,1,1,1,1,1,2,2,1,3), pt.cex=c(1.5,0,0,0,0,0,0,0,0,0),
    lwd=c(1,3,2,1,5,2,1.6,1.6,3,1),
    col=c("gray10","indianred","cyan3","blue2","gray85","gray60","green3","darkviolet","darkviolet","darkviolet"))
  dev.off()
}

plot_fig1(function() pdf(file.path(P$fig, "Figure1.pdf"), width=7, height=7))
plot_fig1(function() png(file.path(P$fig, "Figure1.png"), width=7, height=7, units="in", res=150))

cat("\n==== Figure 1 fitted r_nei (5 complete clusters) ====\n")
print(data.frame(Cluster=clusters, r_nei=R22_best,
                 chi2_R1500=round(chiss2,1), p=chiss2.p), row.names = FALSE)
cat("Figure1 written to figures/.\n")
