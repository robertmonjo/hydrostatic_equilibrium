"""
fig0_four_clusters.py
Supplementary four-cluster figures: acceleration profiles for four
clusters/groups (NGC-533, NGC-5044, Abell-2717, Abell-2029) from the
Angus et al. (2008) mass profiles (data/four_clusters.txt), comparing the
observed inward gravity with the HMG 0-parameter model.
Outputs supplementary_fig_angus_four_clusters[_scheme_b]_python.{pdf,png}.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import hmg_common as h

a = h.load_angus()
clusts = ["N533", "N5044", "A2717", "A2029"]
titles = ["NGC-533 group", "NGC-5044 group", "Abell-2717 cluster", "Abell-2029 cluster"]
rad = a["radkpc"]
kpc, GN, Msol, c0, T0 = h.kpc, h.GN, h.Msol, h.c0, h.T0

acc_obs = GN*Msol*a["total"]/(rad*kpc)**2


def logfit(x, y, powers):
    """Least-squares fit of log(y) ~ sum_k c_k (log x)^p_k ; return fitted y for
    every x. Non-finite points are excluded from the fit, but predictions are
    returned for all x, so NaN gaps in y get filled by the fitted curve."""
    lx, ly = np.log(x), np.log(y)
    m = np.isfinite(lx) & np.isfinite(ly)
    Afit = np.column_stack([np.ones(int(m.sum()))] + [lx[m]**p for p in powers])
    coef, *_ = np.linalg.lstsq(Afit, ly[m], rcond=None)
    A = np.column_stack([np.ones_like(lx)] + [lx**p for p in powers])
    return np.exp(A @ coef)


# ---- Recover Newtonian baryonic acceleration + error ----------------------
b = 2*((c0/T0)/7.9)
d1 = -(GN*Msol*a["thk_bnk"]/(rad*kpc)**2)**2
d2 = -(GN*Msol*a["thn_bnk"]/(rad*kpc)**2)**2
acc_N = (-b + np.sqrt(b**2 - 4*d1))/2
acc_Nerr = np.abs(((-b + np.sqrt(b**2 - 4*d2))/2 - acc_N)/acc_N)

# Fill gas (red) NA gaps by per-cluster log-quadratic fit
gas2 = a["red"].copy()
for c in clusts:
    lg = a["clust"] == c
    gas2[lg] = logfit(rad[lg], a["red"][lg], [1, 2])
fill = np.isnan(a["red"]) & ~np.isnan(a["total"])
a["red"] = a["red"].copy()
a["red"][fill] = gas2[fill]

# Stellar acceleration
acc_stars = acc_N - (GN*Msol*a["red"]/(rad*kpc)**2)
acc_stars = np.where(acc_stars <= 0, np.nan, acc_stars)
star1 = acc_stars.copy()
for c in clusts:
    lg = (a["clust"] == c)
    valid = lg & ~np.isnan(acc_stars)
    lx = np.log(rad[valid])
    A = np.column_stack([np.ones_like(lx), lx**1, lx**2, lx**4])
    coef, *_ = np.linalg.lstsq(A, np.log(acc_stars[valid]), rcond=None)
    lxl = np.log(rad[lg])
    star1[lg] = np.exp(coef[0] + coef[1]*np.log(rad[lg]**1) +
                       coef[2]*np.log(rad[lg])**2 + coef[3]*np.log(rad[lg])**4)
fill2 = np.isnan(acc_stars) & ~np.isnan(a["total"])
acc_stars[fill2] = star1[fill2]
acc_newton = acc_stars + (GN*Msol*a["red"]/(rad*kpc)**2)

lg3 = ~np.isnan(acc_Nerr)
coef_err = np.polyfit(np.log(rad[lg3]), np.log(acc_Nerr[lg3]), 1)  # slope, intercept
acc_obs_err = acc_obs*(0.005 + np.exp(coef_err[1] + coef_err[0]*np.log(rad)))

# ---- HMG 0-parameter model (eps0 without +1/6) ----------------------------
ve2 = 2*acc_newton*rad*kpc
vh2 = (rad*kpc)**2/T0**2
M_clus = acc_newton/GN*(rad*kpc)**2
dens = M_clus/(4/3*np.pi*(rad*kpc)**3)
eps0 = np.sqrt(dens/h.rho_vac)
q0 = np.abs(ve2 - vh2*eps0**2)/(vh2*eps0**2 + ve2)
g0 = h.hmg_gamma0(q0)
a_pred0 = np.sqrt(acc_newton**2 + 2*acc_newton*((c0/T0)/g0))
acc_stars[(a["clust"] == "A2029") & (rad == 1000)] = np.nan   # cosmetic mask: A2029 outermost point

XLIM2 = [250, 300, 800, 1000]


from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
YE = list(range(-13, -8))          # log10(a) tick exponents on the y axis


def draw(fname, gas_col, star_col, newt_col, model_col, band_col, band_lw):
    fig, axes = plt.subplots(2, 2, figsize=(7, 7))
    plt.subplots_adjust(hspace=0.28, wspace=0.28, bottom=0.20, top=0.95)
    for ic, c in enumerate(clusts):
        ax = axes.flat[ic]
        lg = a["clust"] == c
        order = np.arange(len(rad))[lg]
        r = rad[lg]
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(20, XLIM2[ic]); ax.set_ylim(5e-13, 1.18e-9)
        h.panel_style(ax)
        ax.plot(r, acc_obs[lg], color="0.90", lw=6)
        ax.plot(r, acc_obs[lg], color="0.42", marker="o", ls="", ms=4)
        ax.plot(r, (GN*Msol*a["thn_bnk"]/(rad*kpc)**2)[lg], color=band_col, lw=band_lw[0])
        ax.plot(r, (GN*Msol*a["thk_bnk"]/(rad*kpc)**2)[lg], color=band_col, lw=band_lw[1])
        ax.plot(r, (GN*Msol*a["red"]/(rad*kpc)**2)[lg], color=gas_col, lw=2)
        ax.plot(r, acc_stars[lg], color=star_col, lw=2)
        ax.plot(r, acc_newton[lg], color=newt_col, lw=2, ls="--")
        ax.plot(r, a_pred0[lg], color=model_col, lw=2)
        ax.plot(r, (a_pred0+0.5*acc_obs_err)[lg], color=model_col, lw=0.8)
        ax.errorbar(r, acc_obs[lg],
                    yerr=[0.5*acc_obs_err[lg], acc_obs_err[lg]],
                    fmt="none", ecolor="0.30", lw=0.7)
        ax.set_title(titles[ic], fontsize=10)
        ax.set_xlabel("Radius [kpc]")
        ax.set_ylabel(h.YLAB)
        # x radius labels as plain integers
        xt = [20, 50, 100, 200, 500, 1000]
        ax.xaxis.set_major_locator(FixedLocator(xt))
        ax.xaxis.set_major_formatter(FixedFormatter([str(v) for v in xt]))
        ax.xaxis.set_minor_locator(NullLocator())
        # y ticks labelled as log10 exponents (-13..-9)
        ax.yaxis.set_major_locator(FixedLocator([10.0**e for e in YE]))
        ax.yaxis.set_major_formatter(FixedFormatter([str(e) for e in YE]))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.tick_params(axis="both", which="both", direction="out", labelsize=8)
    fig.legend(handles=[
        plt.Line2D([], [], color="0.42", marker="o", ls="", label="Inward gravity for dyn. and hydros. equilibrium"),
        plt.Line2D([], [], color=gas_col, lw=2, label="Hot-gas baryonic gravity"),
        plt.Line2D([], [], color=star_col, lw=2, label="Stars and dust baryonic gravity"),
        plt.Line2D([], [], color=newt_col, lw=2, ls="--", label="Total Newtonian baryonic gravity"),
        plt.Line2D([], [], color=band_col, lw=band_lw[1], label="MOND model (Angus et al. 2008)"),
        plt.Line2D([], [], color=band_col, lw=band_lw[0], label=r"MOND model +1$\sigma$ (Angus et al. 2008)"),
        plt.Line2D([], [], color=model_col, lw=2, label="HMG model (0-parameter, this work)"),
        plt.Line2D([], [], color=model_col, lw=0.8, label=r"HMG model +1$\sigma$ (0-parameter, this work)"),
    ], loc="lower center", ncol=2, fontsize=7.5, frameon=False, labelspacing=0.5)
    fig.savefig(os.path.join(h.FIG, fname + ".pdf"))
    fig.savefig(os.path.join(h.FIG, fname + ".png"), dpi=150)
    plt.close(fig)


draw("supplementary_fig_angus_four_clusters_python",
     gas_col="indianred", star_col="c", newt_col="blue",
     model_col="darkviolet", band_col="green", band_lw=(0.8, 2))
draw("supplementary_fig_angus_four_clusters_scheme_b_python",
     gas_col="y", star_col="c", newt_col="darkred",
     model_col="indianred", band_col=(1, 0.5, 0.5, 0.4), band_lw=(1, 4))
print("Supplementary four-cluster figures written.")
