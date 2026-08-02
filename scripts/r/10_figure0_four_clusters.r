# =====================================================================
# 10_figure0_four_clusters.r
# Supplementary four-cluster figures: acceleration profiles for four
# clusters/groups (NGC-533, NGC-5044, Abell-2717, Abell-2029) using the
# mass profiles from Angus et al. (2008), file data/four_clusters.txt.
# The HMG 0-parameter (Hubble-Newton equilibrium) model is compared to
# the inward gravity required for dynamical/hydrostatic equilibrium.
# Outputs supplementary_fig_angus_four_clusters[_scheme_b].{pdf,png}.
# =====================================================================

sdir <- tryCatch({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else getwd()
}, error = function(e) getwd())
source(file.path(sdir, "00_hmg_common.r"))
P <- hmg_paths()

angus <- read.table(file.path(P$data, "four_clusters.txt"), quote = "#", sep = "\t")
colnames(angus) <- c("clust","radkpc","total","thn_bnk","thk_bnk",
                     "thn_blu","thk_blu","green","red")
clusts <- unique(angus$clust)                       # N533 N5044 A2717 A2029
clusters <- c("NGC-533 group","NGC-5044 group",
              "Abell-2717 cluster","Abell-2029 cluster")

acc_obs <- GN*Msol*angus$total/(angus$radkpc*kpc)^2

## ---- Recover a Newtonian baryonic acceleration + its error ----------
## Invert a_C^2 = a_N^2 + 2 a_N (c/t)/g0  with g0 ~ 7.9 to get a_N from the
## thick/thin baryonic-mass curves (banked mass), then model stars/gas.
a  <- 1
b  <- 2*((c0/T0)/7.9)
d1 <- -(GN*Msol*angus$thk_bnk/(angus$radkpc*kpc)^2)^2
d2 <- -(GN*Msol*angus$thn_bnk/(angus$radkpc*kpc)^2)^2
acc_N    <- (-b + sqrt(b^2 - 4*a*d1))/2
acc_Nerr <- abs(((-b + sqrt(b^2 - 4*a*d2))/2 - acc_N)/acc_N)

## Fill gas (red) NA gaps by a per-cluster log-quadratic fit
gas2 <- angus$red
for (ic in 1:4) {
  lg <- angus$clust == clusts[ic]
  lm2 <- lm(log(angus$red[lg]) ~ log(angus$radkpc[lg]) + I(log(angus$radkpc[lg])^2))
  gas2[lg] <- exp(predict(lm2, newdata = data.frame(angus$radkpc[lg])))
  gas2[lg] <- exp(lm2$coefficients[1] +
                  lm2$coefficients[2]*log(angus$radkpc[lg]) +
                  lm2$coefficients[3]*log(angus$radkpc[lg])^2)
}
angus$red[is.na(angus$red) & !is.na(angus$total)] <-
  gas2[is.na(angus$red) & !is.na(angus$total)]

## Stellar acceleration = total baryonic Newtonian - gas
acc_stars <- acc_N - (GN*Msol*angus$red/(angus$radkpc*kpc)^2)
acc_stars[acc_stars <= 0] <- NA
star1 <- acc_stars
for (ic in 1:4) {
  lg <- angus$clust == clusts[ic]
  lm1 <- lm(log(acc_stars[lg]) ~ I(log(angus$radkpc[lg])^1) +
              I(log(angus$radkpc[lg])^2) + I(log(angus$radkpc[lg])^4))
  star1[lg] <- exp(lm1$coefficients[1] +
                   lm1$coefficients[2]*log(angus$radkpc[lg]^1) +
                   lm1$coefficients[3]*log(angus$radkpc[lg])^2 +
                   lm1$coefficients[4]*log(angus$radkpc[lg])^4)
}
acc_stars[is.na(acc_stars) & !is.na(angus$total)] <-
  star1[is.na(acc_stars) & !is.na(angus$total)]
acc_newton <- acc_stars + (GN*Msol*angus$red/(angus$radkpc*kpc)^2)

lg3 <- !is.na(acc_Nerr)
lm_err <- lm(log(acc_Nerr[lg3]) ~ log(angus$radkpc[lg3]))
acc_obs_err <- acc_obs*(0.005 + exp(lm_err$coefficients[1] +
                                    lm_err$coefficients[2]*log(angus$radkpc)))

## ---- HMG 0-parameter model (eps0 without the 1/6 term) --------------
ve_clus2  <- 2*acc_newton*angus$radkpc*kpc
vh_clus2  <- (angus$radkpc*kpc)^2/T0^2
M_clus    <- acc_newton/GN*(angus$radkpc*kpc)^2
dens_clus <- M_clus/(4/3*pi*(angus$radkpc*kpc)^3)
eps0      <- sqrt(dens_clus/rho_vac)
q0  <- abs(ve_clus2 - vh_clus2*eps0^2)/(vh_clus2*eps0^2 + ve_clus2)
g0  <- hmg_gamma0(q0)
a_pred0 <- sqrt(acc_newton^2 + 2*acc_newton*((c0/T0)/g0))
acc_stars[86] <- NA           # cosmetic mask of a single outlier point (A2029 outermost)

## ---- Plot helper (shared by Fig0 and Fig0b) -------------------------
draw_panels <- function(gas_col, star_col, newt_col, newt_lty,
                        model_col, band_thn, band_thk, band_lwd) {
  xlim2 <- c(250,300,800,1000)
  par(mfrow = c(2,2), mar = c(3.1,3.1,0.3,0.3), oma = c(4.6,0,0,0))
  for (ic in 1:4) {
    lg <- angus$clust == clusts[ic]
    plot(angus$radkpc[lg], acc_obs[lg], col="gray85", type="l", axes=FALSE,
         xlab="", ylab="", lty=1, lwd=6, log="xy",
         xlim=c(20,xlim2[ic]), ylim=c(5e-13,1.18e-9))
    rect(1,1e-14,2000,1e-8, col=rgb(1,1,0.95,0.5)); box()
    axis(1, c(20,50,100,200,500,1000), padj=-0.8)
    axis(2, 10^(seq(-13,-9,1)), seq(-13,-9,1), padj=0.8)
    mtext(side=1,"Radius [kpc]",line=1.7)
    mtext(side=2,expression(paste(log[10],"(",a," / ",ms^{-2},")")),line=1.5)
    abline(v=c(20,50,100,200,500,1000),
           h=c(10^(seq(-13,-9,1)),5*10^(seq(-13,-9,1))), lty=2, col="gray95")
    lines(angus$radkpc[lg], acc_obs[lg], col="gray95", lwd=6)
    points(angus$radkpc[lg], acc_obs[lg], col="gray42", pch=20, cex=1.5)
    lines(angus$radkpc[lg], (GN*Msol*angus$thn_bnk/(angus$radkpc*kpc)^2)[lg], col=band_thn, lwd=band_lwd[1])
    lines(angus$radkpc[lg], (GN*Msol*angus$thk_bnk/(angus$radkpc*kpc)^2)[lg], col=band_thk, lwd=band_lwd[2])
    lines(angus$radkpc[lg], (GN*Msol*angus$red/(angus$radkpc*kpc)^2)[lg], col=gas_col, lwd=2)
    lines(angus$radkpc[lg], acc_stars[lg], col=star_col, lwd=2)
    lines(angus$radkpc[lg], acc_newton[lg], col=newt_col, lty=newt_lty, lwd=2)
    lines(angus$radkpc[lg], a_pred0[lg], col=model_col, lwd=2, lty=1)
    lines(angus$radkpc[lg], (a_pred0+0.5*acc_obs_err)[lg], col=model_col, lwd=0.8, lty=1)
    segments(angus$radkpc[lg], (acc_obs-0.5*acc_obs_err)[lg],
             angus$radkpc[lg], (acc_obs+acc_obs_err)[lg], col="gray30")
  }
  ## panel titles
  coords <- list(c(0,0.5,0.5,1), c(0.5,1,0.5,1), c(0,0.5,0,0.5), c(0.5,1,0,0.5))
  for (i in 1:4) {
    par(fig=coords[[i]], mar=c(3.1,3.1,0.3,0.3), oma=c(4.6,0,0,0), new=TRUE)
    plot.new(); text(0.5,0.98, clusters[i])
  }
}

## ---- Supplementary figure (default colour scheme) -------------------
draw_fig0 <- function(dev_open) {
  dev_open()
  draw_panels(gas_col="indianred", star_col="cyan3", newt_col="blue2", newt_lty=2,
              model_col="darkviolet", band_thn="green3", band_thk="green3",
              band_lwd=c(0.8,2))
  par(fig=c(0,1,0,1), mar=c(0,3.1,0.3,0.3), oma=c(0,0,0,0), new=TRUE); plot.new()
  legend("bottom", ncol=2, bty="n", pch=20, pt.cex=c(1.5,0,0,0,0,0,0,0),
    legend=c("Inward gravity for dyn. and hydros. equilibium  ",
             "Hot-gas baryonic gravity","Stars and dust baryonic gravity",
             "Total Newtonian baryonic gravity","MOND model (Angus et al. 2008)",
             expression(paste("MOND model +1",sigma," (Angus et al. 2008)")),
             "HMG model (0-parameter, this work)",
             expression(paste("HMG model +1",sigma," (0-parameter, this work)"))),
    lty=c(1,1,1,2,1,1,1,1), lwd=c(1,2,2,2,2,1,2,1),
    col=c("gray40","indianred","cyan3","blue2","green3","green3","darkviolet","darkviolet"))
  dev.off()
}

## ---- Supplementary figure (alternate colour scheme, scheme b) -------
draw_fig0b <- function(dev_open) {
  dev_open()
  draw_panels(gas_col="yellow3", star_col="cyan3", newt_col="darkred", newt_lty=2,
              model_col="indianred",
              band_thn=rgb(1,0.5,0.5,0.2), band_thk=rgb(1,0.5,0.5,0.2),
              band_lwd=c(1,4))
  par(fig=c(0,1,0,1), mar=c(0,3.1,0.3,0.3), oma=c(0,0,0,0), new=TRUE); plot.new()
  legend("bottom", ncol=2, bty="n", pch=20, pt.cex=c(1.5,0,0,0,0,0,0,0),
    legend=c("Inward gravity for dyn. and hydros. equilibium  ",
             "Hot-gas baryonic gravity","Stars and dust baryonic gravity",
             "Total Newtonian baryonic gravity","MOND model (Angus et al. 2008)",
             expression(paste("MOND model +1",sigma," (Angus et al. 2008)")),
             "HMG model (0-parameter, this work)",
             expression(paste("HMG model +1",sigma," (0-parameter, this work)"))),
    lty=c(1,1,1,2,1,1,1,1), lwd=c(1,2,2,2,3,1.2,2,1),
    col=c("gray40","yellow3","cyan3","darkred",
          rgb(1,0.5,0.5,0.2),rgb(1,0.5,0.5,0.2),"indianred","indianred"))
  dev.off()
}

draw_fig0(function()  pdf(file.path(P$fig,"supplementary_fig_angus_four_clusters.pdf"),  width=7, height=7))
draw_fig0(function()  png(file.path(P$fig,"supplementary_fig_angus_four_clusters.png"),  width=7, height=7, units="in", res=150))
draw_fig0b(function() pdf(file.path(P$fig,"supplementary_fig_angus_four_clusters_scheme_b.pdf"), width=7, height=7))
draw_fig0b(function() png(file.path(P$fig,"supplementary_fig_angus_four_clusters_scheme_b.png"), width=7, height=7, units="in", res=150))
cat("Supplementary four-cluster figures written to figures/.\n")
