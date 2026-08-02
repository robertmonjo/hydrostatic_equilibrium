"""
fig1_complete.py
Figure 1 (main paper): acceleration profiles for the five X-COP clusters
with complete stellar-mass information (A1795, A2029, A2142, A2319, A644).
Reads the Eckert et al. (2022) baryonic-mass and gbar/gobs FITS tables,
fits the HMG r_nei parameter per cluster, and draws baryonic, missing and
total mass together with the MOND and HMG predictions.
Outputs Figure1_python.{pdf,png}.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist
import hmg_common as h

SUBDIR = "data_E22"
clusters = sorted(os.listdir(os.path.join(h.DATA, SUBDIR)))  # A1795 A2029 A2142 A2319 A644


def read_bar(clust):
    """Return (bar_mass, gbar_gobs) FITS tables for one cluster."""
    f = h.cluster_fits(SUBDIR, clust)
    fbar = [x for x in f if "bar_mass" in os.path.basename(x)][0]
    fgg  = [x for x in f if "gbar_gobs" in os.path.basename(x)][0]
    d1 = h.read_fits_table(fbar)   # EIN3_bar_mass: baryonic + total mass profiles
    d2 = h.read_fits_table(fgg)    # gbar_gobs: observed inward acceleration
    return d1, d2


# ---- Fit r_nei (min chi^2 for R>=1500) ------------------------------------
R0s = np.arange(50, 2001, 10)
R22 = np.full(len(clusters), np.nan)
chi2v = np.full(len(clusters), np.nan)
chi2p = np.full(len(clusters), np.nan)
for i, clust in enumerate(clusters):
    d1, _ = read_bar(clust)
    M_kg = h.Msol*(d1["M_STAR"] + d1["MGAS"])
    Rref = d1["R_REF"]
    obs    = h.GN*h.Msol*d1["MASS"]    / (Rref*h.kpc)**2
    obs_HI = h.GN*h.Msol*d1["MASS_HI"] / (d1["R_IN"]*h.kpc)**2
    obs_LO = h.GN*h.Msol*d1["MASS_LO"] / (d1["R_OUT"]*h.kpc)**2
    obs_err = (obs_HI - obs_LO)/2
    lgR2 = Rref >= 1500
    chi = np.array([np.nansum(((h.hmg_predict(M_kg, Rref, R0, r_grav=40)["a_pred1"]-obs)[lgR2])**2
                              / obs_err[lgR2]**2) for R0 in R0s])
    R22[i] = R0s[np.argmin(chi)]
    chi2v[i] = chi.min()
    chi2p[i] = round(100*chi2_dist.cdf(chi.min(), df=int(lgR2.sum())-1))

# ---- Figure 1 -------------------------------------------------------------
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
fig, axes = plt.subplots(3, 2, figsize=(7, 7))
plt.subplots_adjust(wspace=0, hspace=0, left=0.085, right=0.915, top=0.915, bottom=0.065)
XT = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
YT = [-11, -10, -9]
for i, clust in enumerate(clusters):
    ax = axes.flat[i]
    row, col = divmod(i, 2)
    d1, d2 = read_bar(clust)
    Rref = d1["R_REF"]
    R0 = R22[i] + (200 if R22[i] < 400 else 0)
    pr = h.hmg_predict(h.Msol*(d1["M_STAR"]+d1["MGAS"]), Rref, R0, r_grav=40)

    def L(y):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log10(y)
    ax.set_xscale("log"); ax.set_xlim(1, 3000); ax.set_ylim(-11.5, -8.5)
    h.panel_style(ax)
    ax.axvline(R0, ls="--", color="0.82", lw=0.8)
    ax.plot(Rref, L(h.GN*h.Msol*d1["M_DM"]/(Rref*h.kpc)**2), color="0.88", lw=5)
    ax.plot(Rref, L(h.GN*h.Msol*d1["MGAS"]/(Rref*h.kpc)**2), color="indianred", lw=3)
    ax.plot(Rref, L(h.GN*h.Msol*d1["M_STAR"]/(Rref*h.kpc)**2), color="c", lw=2)
    ax.plot(Rref, L(h.GN*h.Msol*(d1["M_STAR"]+d1["MGAS"])/(Rref*h.kpc)**2), color="blue", lw=1)
    ax.plot(Rref, L(h.GN*h.Msol*d1["MASS"]/(Rref*h.kpc)**2), color="0.60", lw=2)
    ax.plot(Rref, L(pr["a_MOND"]),   color="green",      lw=1.5, ls="--")
    ax.plot(Rref, L(pr["a_pred_s"]), color="darkviolet", lw=1.5, ls="--")
    ax.plot(Rref, L(pr["a_pred0"]),  color="darkviolet", lw=0.8, ls=":")
    ax.plot(Rref, L(pr["a_pred1"]),  color="darkviolet", lw=3)
    ax.errorbar(d2["RADIUS"], L(d2["GOBS"]),
                yerr=[L(d2["GOBS"])-L(d2["GOBS"]-d2["GOBS_ERR_LO"]),
                      L(d2["GOBS"]+d2["GOBS_ERR_HI"])-L(d2["GOBS"])],
                fmt="o", mfc="none", mec="k", ecolor="k", ms=4, lw=0.8)
    ax.text(0.74, 0.88, clust, transform=ax.transAxes, fontsize=9)

    # Radius labels sit on the top row (panels 1,2) and under panels 4,5;
    # log a labels on the left (1,3,5) and right (2,4). Trim the shared
    # boundary label so "2000" (left panel) and "1" (right panel) do not collide.
    top = i in (0, 1)
    bot = i in (3, 4)
    lft = col == 0
    rgt = i in (1, 3)
    lab = list(XT)
    if col == 1:                       # right column: drop leftmost "1"
        lab = [("" if v == 1 else str(v)) for v in XT]
    else:                              # left column: drop rightmost "2000"
        lab = [("" if v == 2000 else str(v)) for v in XT]
    ax.xaxis.set_major_locator(FixedLocator(XT))
    ax.xaxis.set_major_formatter(FixedFormatter(lab))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_locator(FixedLocator(YT))
    ax.yaxis.set_major_formatter(FixedFormatter([str(v) for v in YT]))
    ax.tick_params(axis="x", which="both", direction="out", labelsize=6,
                   top=top, bottom=bot, labeltop=top, labelbottom=bot)
    ax.tick_params(axis="y", which="both", direction="out", labelsize=7,
                   left=lft, right=rgt, labelleft=lft, labelright=rgt)
    if top:
        ax.set_title("Radius [kpc]", fontsize=8.5, pad=15)
    if bot:
        ax.set_xlabel("Radius [kpc]", fontsize=8.5)
    if lft:
        ax.set_ylabel(h.YLAB, fontsize=8.5)
    if rgt:
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(h.YLAB, fontsize=8.5, rotation=270, labelpad=14)
axes.flat[5].axis("off")
axes.flat[5].legend(handles=[
    plt.Line2D([], [], marker="o", mfc="none", mec="k", ls="", label="Inward gravity for hydrostatic equilibrium"),
    plt.Line2D([], [], color="indianred", lw=3, label="Hot-gas baryonic gravity"),
    plt.Line2D([], [], color="c", lw=2, label="Stars and dust baryonic gravity (accounted to date)"),
    plt.Line2D([], [], color="blue", lw=1, label="Total baryonic (Newtonian) gravity"),
    plt.Line2D([], [], color="0.85", lw=5, label="Estimated missing mass"),
    plt.Line2D([], [], color="0.60", lw=2, label="Total estimated mass: baryonic + missing mass"),
    plt.Line2D([], [], color="green", lw=1.5, ls="--", label="MOND model + MLS interpolation function (no EFE)"),
    plt.Line2D([], [], color="darkviolet", lw=1.5, ls="--", label="HMG model with only spatial contribution"),
    plt.Line2D([], [], color="darkviolet", lw=3, label="HMG model with also time-like contribution"),
    plt.Line2D([], [], color="darkviolet", lw=1, ls=":", label=r"HMG model with Hubble-Newton equilibrium ($\gamma_{sys}=\pi/3$)"),
], loc="center", bbox_to_anchor=(0.5, 0.42), fontsize=6.6, frameon=False,
   labelspacing=0.5, handlelength=2.4)
fig.savefig(os.path.join(h.FIG, "Figure1_python.pdf"))
fig.savefig(os.path.join(h.FIG, "Figure1_python.png"), dpi=150)

print("\n==== Figure 1 fitted r_nei (Python) ====")
for i, c in enumerate(clusters):
    print(f"{c:<7} r_nei={int(R22[i]):>5}  chi2(R>=1500)={chi2v[i]:6.1f}  p={int(chi2p[i])}")
print("Figure1_python written.")
