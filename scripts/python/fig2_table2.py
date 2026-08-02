"""
fig2_table2.py
Figure 2 (main paper): gas-only HMG acceleration profiles for the twelve
X-COP clusters, plus Table 2 (fitted r_nei, chi^2, p-value, significance).
Reads the Einasto mass profiles (Eckert et al. 2022), fits r_nei per
cluster over three outer-radius windows, and writes the figure and the
Table 2 CSV.
Outputs Figure2_python.{pdf,png} and tables/table2_python.csv.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist
import hmg_common as h

SUBDIR = "einasto_mass_profiles"
clusters_alpha = sorted(os.listdir(os.path.join(h.DATA, SUBDIR)))  # alphabetical
# Display order matching the paper (A85, A644, A1644, ...); 1-based indices
clu_order = [10, 9, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12]
clusters = [clusters_alpha[i-1] for i in clu_order]


def read_cluster(clust):
    """Read one Einasto mass profile; derive stellar mass = total - DM - gas."""
    f = h.cluster_fits(SUBDIR, clust)
    d = h.read_fits_table(f[0])
    d["M_STAR"] = d["MASS"] - d["M_DM"] - d["MGAS"]
    return d


R0s = np.arange(200, 2001, 10)
n = len(clusters)
R0_best = np.full(n, np.nan)
R1_best = np.full(n, np.nan)
R2_best = np.full(n, np.nan)
chis1 = np.full(n, np.nan); chis2 = np.full(n, np.nan)
chis1p = np.full(n, np.nan); chis2p = np.full(n, np.nan)
chimond = np.full(n, np.nan)
nbins_R500 = np.full(n, 0, dtype=int)

for i, clust in enumerate(clusters):
    d = read_cluster(clust)
    M_gas = h.Msol * d["MGAS"]
    Rref = d["R_REF"]

    obs    = h.GN*h.Msol*(d["M_DM"]    + d["MGAS"])    / (Rref     *h.kpc)**2
    # R_IN can be 0 for the innermost bin; suppress the divide warning (this bin
    # lies inside the R>=1000 kpc fit window's lower edge and does not affect it).
    with np.errstate(divide="ignore", invalid="ignore"):
        obs_HI = h.GN*h.Msol*(d["M_DM_HI"] + d["MGAS_HI"]) / (d["R_IN"] *h.kpc)**2
    obs_LO = h.GN*h.Msol*(d["M_DM_LO"] + d["MGAS_LO"]) / (d["R_OUT"]*h.kpc)**2
    obs_err = (obs_HI - obs_LO) / 2

    lgR0 = Rref >= 1500
    lgR1 = Rref >= 500
    lgR2 = Rref >= 1000

    chi00 = np.empty_like(R0s, float)
    chi01 = np.empty_like(R0s, float)
    chi02 = np.empty_like(R0s, float)
    aMOND2 = None
    for j, R0 in enumerate(R0s):
        pr = h.hmg_predict(M_gas, Rref, R0, r_grav=50)
        ap = pr["a_pred1"]
        chi00[j] = np.nansum(((ap-obs)[lgR0])**2 / obs_err[lgR0]**2)
        chi01[j] = np.nansum(((ap-obs)[lgR1])**2 / obs_err[lgR1]**2)
        chi02[j] = np.nansum(((ap-obs)[lgR2])**2 / obs_err[lgR2]**2)
        aMOND2 = pr["a_MOND2"]
    chimond[i] = np.nansum(((aMOND2-obs)[lgR1])**2 / obs_err[lgR1]**2)
    nbins_R500[i] = int(lgR1.sum())

    R0_best[i] = R0s[np.argmin(chi00)]
    R1_best[i] = R0s[np.argmin(chi01)]
    R2_best[i] = R0s[np.argmin(chi02)]
    chis1[i] = chi01.min(); chis2[i] = chi02.min()
    chis1p[i] = round(100*chi2_dist.cdf(chis1[i], df=int(lgR1.sum())-1))
    chis2p[i] = round(100*chi2_dist.cdf(chis2[i], df=int(lgR2.sum())-1))

# Paper convention: for clusters whose chi^2 CDF rounds to 100% (a perfect or
# over-fit), the reported r_nei is the fitted value + 200 kpc; chi^2 is not
# recomputed.
R2_best[chis2p == 100] += 200
# r_nei uncertainty (empirical window-sensitivity spread, not a formal
# confidence interval): sqrt(10^2 + sample standard deviation, over the three
# outer-radius windows (R>=500, >=1000, >=1500 kpc), of the best-fit r_nei^2),
# with the 10 kpc term set by the R0 scan grid step. This variant matches the
# paper's tabulated uncertainties. RXC1825 is an outlier relative to the
# paper's quoted +/-110.
R_err = np.round(np.sqrt(10**2 +
        np.std(np.vstack([R2_best**2, R1_best**2, R0_best**2]), axis=0, ddof=1)))

# Degrees of freedom in Table 2:
#  - the tabulated chi^2 (chis2) is computed over the R>=1000 kpc window
#    (per-cluster ~22-31 bins);
#  - the chi2_cdf_pct column is the lower-tail chi^2 CDF percentile (chi2.cdf)
#    with per-cluster dof, df = n_i - 1 (lgR2.sum()-1). This is the paper's
#    convention (high value = rejection), not a standard upper-tail
#    goodness-of-fit p-value;
#  - the 95% significance flag compares chi^2 to the fixed threshold
#    chi^2(95%, df=42) = 58.12, matching the paper's n=43 (R>=500 kpc)
#    convention.
signif = np.where(chis2 < h.CHI2_THRESHOLD_95, "*", "-")


def pfmt(p):
    return "<1" if p < 1 else (">99" if p > 99 else str(int(p)))


print("\n==== TABLE 2 (Python) ====")
print(f"{'Cluster':<9}{'r_nei':>7}{'+/-':>6}{'chi2':>8}{'chi2_cdf_pct':>14}{'sig':>5}")
for i, c in enumerate(clusters):
    print(f"{c:<9}{int(R2_best[i]):>7}{int(R_err[i]):>6}"
          f"{chis2[i]:>8.1f}{pfmt(chis2p[i]):>14}{signif[i]:>5}")
print(f"\nchi^2 threshold (95%, df=42) = {h.CHI2_THRESHOLD_95:.2f}")
# MOND (MLS interpolation, no EFE) evaluated per cluster over the R>=500 kpc
# window. The smallest per-cluster chi^2 is MOND's best case across the sample;
# even this minimum lies far above the 95% threshold.
imin = int(np.nanargmin(chimond))
print(f"MOND best-case (minimum) per-cluster chi^2 = {chimond[imin]:.0f} "
      f"({clusters[imin]}, n={nbins_R500[imin]} bins, R>=500 kpc)")

# CSV
import csv
with open(os.path.join(h.TAB, "table2_python.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Cluster", "r_nei", "r_nei_err", "chi2", "chi2_cdf_pct", "signif95"])
    for i, c in enumerate(clusters):
        w.writerow([c, int(R2_best[i]), int(R_err[i]),
                    round(chis2[i], 1), pfmt(chis2p[i]), signif[i]])

# ---- Figure 2 -------------------------------------------------------------
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
XT = [200, 500, 1000, 2000]
YT = [-11, -10, -9]
fig, axes = plt.subplots(3, 4, figsize=(7, 7.7))
plt.subplots_adjust(wspace=0, hspace=0, left=0.085, right=0.915, top=0.905, bottom=0.185)
for i, clust in enumerate(clusters):
    ax = axes.flat[i]
    row, col = divmod(i, 4)
    d = read_cluster(clust)
    Rref = d["R_REF"]
    pr = h.hmg_predict(h.Msol*d["MGAS"], Rref, R2_best[i], r_grav=50)

    def L(y): return np.log10(y)
    ax.set_xscale("log")
    ax.set_xlim(200, 3000); ax.set_ylim(-11.5, -9)
    h.panel_style(ax)
    ax.axvline(R2_best[i], ls="--", color="0.80", lw=0.8)
    ax.plot(Rref, L(h.GN*h.Msol*d["MGAS"]/(Rref*h.kpc)**2), color="indianred", lw=2)
    ax.plot(Rref, L(h.GN*h.Msol*d["M_DM"]/(Rref*h.kpc)**2), color="0.90", lw=5)
    ax.plot(Rref, L(h.GN*h.Msol*(d["M_DM"]+d["MGAS"])/(Rref*h.kpc)**2), color="0.75", lw=3)
    ax.plot(Rref, L(pr["a_MOND"]),   color="green",      lw=1.5, ls="--")
    ax.plot(Rref, L(pr["a_pred_s"]), color="darkviolet", lw=1.5, ls="--")
    ax.plot(Rref, L(pr["a_pred0"]),  color="darkviolet", lw=0.8, ls=":")
    ax.plot(Rref, L(pr["a_pred1"]),  color="darkviolet", lw=3)
    ax.text(0.95, 0.88, clust, transform=ax.transAxes, fontsize=8.5, ha="right")

    # Shared outer axes: radius on top+bottom, log a on left+right.
    # Keep "200" on every panel but drop the "2000" that would collide with
    # the next panel's "200"; keep "2000" only on the right-most column.
    lab = [("" if (v == 2000 and col != 3) else str(v)) for v in XT]
    ax.xaxis.set_major_locator(FixedLocator(XT))
    ax.xaxis.set_major_formatter(FixedFormatter(lab))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_locator(FixedLocator(YT))
    ax.yaxis.set_major_formatter(FixedFormatter([str(v) for v in YT]))
    ax.tick_params(axis="x", which="both", direction="out", labelsize=7,
                   top=(row == 0), bottom=(row == 2),
                   labeltop=(row == 0), labelbottom=(row == 2))
    ax.tick_params(axis="y", which="both", direction="out", labelsize=7,
                   left=(col == 0), right=(col == 3),
                   labelleft=(col == 0), labelright=(col == 3))
    if row == 0:
        ax.set_title("Radius [kpc]", fontsize=8.5, pad=16)
    if row == 2:
        ax.set_xlabel("Radius [kpc]", fontsize=8.5)
    if col == 0:
        ax.set_ylabel(h.YLAB, fontsize=8.5)
    if col == 3:
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(h.YLAB, fontsize=8.5, rotation=90, labelpad=16)

# Two-block legend describing the baryonic/missing-mass curves and the models
left_handles = [
    plt.Line2D([], [], color="indianred", lw=2, label="Hot-gas baryonic gravity (Eckert+2022)"),
    plt.Line2D([], [], color="0.90", lw=5, label="Estimated missing-mass gravity (Eckert+2022)"),
    plt.Line2D([], [], color="0.75", lw=3, label="Estimated missing-mass + baryonic gravity (Eckert+2022)"),
]
right_handles = [
    plt.Line2D([], [], color="green", lw=1.5, ls="--", label="MOND model for gas + MLS interpolation function (no EFE)"),
    plt.Line2D([], [], color="darkviolet", lw=1.5, ls="--", label="HMG model for gas with only spatial contribution"),
    plt.Line2D([], [], color="darkviolet", lw=3, label="HMG model for gas with also time-like contribution"),
    plt.Line2D([], [], color="darkviolet", lw=1, ls=":", label=r"HMG model with Hubble-Newton equilibrium ($\gamma_{sys}=\pi/3$)"),
]
fig.legend(handles=left_handles, loc="lower left", bbox_to_anchor=(0.055, 0.012),
           fontsize=6.9, frameon=False, handlelength=2.4, labelspacing=0.4)
fig.legend(handles=right_handles, loc="lower left", bbox_to_anchor=(0.505, 0.012),
           fontsize=6.9, frameon=False, handlelength=2.4, labelspacing=0.4)
fig.savefig(os.path.join(h.FIG, "Figure2_python.pdf"))
fig.savefig(os.path.join(h.FIG, "Figure2_python.png"), dpi=150)
print("Figure2_python written.")
