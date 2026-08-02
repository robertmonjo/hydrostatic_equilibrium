"""
hmg_common.py
Shared constants, the Hyperconical Modified Gravity (HMG) acceleration
model and data loaders for the figure and table scripts of

  Monjo, R. (2025), ApJ 981, 195. DOI 10.3847/1538-4357/adb723
  "Hydrostatic Equilibrium Constraints in X-COP Galaxy Clusters"

Reads only files under ../../data.
"""
import os
import numpy as np
from astropy.io import fits
from scipy.stats import chi2 as chi2_dist

# ---- Shared plotting style ------------------------------------------------
CREAM = (1.0, 1.0, 0.975)          # panel background fill
GRID_COL = "0.90"                  # faint gridline colour
YLAB = r"$\log_{10}(a\,/\,{\rm ms^{-2}})$"


def panel_style(ax, cream=True, grid=True):
    """Apply the cream background, faint gridlines and a full box to a panel."""
    if cream:
        ax.set_facecolor(CREAM)
    if grid:
        ax.grid(True, which="major", color=GRID_COL, ls="-", lw=0.5, zorder=0)
    for s in ax.spines.values():
        s.set_linewidth(0.8)
        s.set_zorder(3)

# ---- Physical constants (SI unless noted) ---------------------------------
c0   = 3e8                                    # speed of light, m/s (paper value)
Msol = 1.9891e30                              # solar mass, kg
kpc  = 3261.8116478174 * 365*24*3600 * c0     # kpc in metres
T0   = 13.7e9 * 365*24*3600                   # age of the universe, s (13.7 Gyr)
GN   = 6.674e-11                              # gravitational constant, SI

rho_vac   = 3.0 / (8*np.pi*GN*T0**2)          # HMG vacuum/critical density
GA_EMPTY  = np.pi/3                            # gamma_U universal angle
GA_CENTER = np.pi/2                            # gamma_center (1-parameter model)

CHI2_THRESHOLD_95 = chi2_dist.ppf(0.95, df=42)  # ~58.12, reported in Table 2

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data")
FIG  = os.path.join(ROOT, "figures")
TAB  = os.path.join(ROOT, "tables")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)


def hmg_gamma0(quotient, gaempty=GA_EMPTY, gacenter=GA_CENTER):
    """HMG projection factor gamma_0 = gamma_sys / cos(gamma_sys)."""
    g_sys = np.arcsin(np.sqrt(np.sin(gaempty)**2 +
                              (np.sin(gacenter)**2 - np.sin(gaempty)**2) * quotient))
    return g_sys / np.cos(g_sys)


def hmg_predict(M_enc_kg, R_ref_kpc, R0_kpc, r_grav=40.0, add_e02=True):
    """
    Full HMG prediction on a radial grid (all accelerations in m/s^2).

    M_enc_kg   : enclosed baryonic mass [kg], array over the grid
    R_ref_kpc  : reference radius [kpc], same length
    R0_kpc     : r_nei fitting parameter [kpc]
    r_grav     : inner gravitational-domination scale [kpc]
    add_e02    : include the +1/6 term in eps^2 (True for X-COP figures)
    """
    M = np.asarray(M_enc_kg, float)
    R_ref_m = np.asarray(R_ref_kpc, float) * kpc
    acc_newton = GN * M / R_ref_m**2
    ve2  = 2 * acc_newton * R_ref_m           # v_N^2 = 2 G M / r
    vh2  = R_ref_m**2 / T0**2                  # v_H^2 = (r/t)^2
    dens = M / (4/3 * np.pi * R_ref_m**3)

    e02  = 1/6 if add_e02 else 0.0
    eps0 = np.sqrt(e02 + dens / rho_vac)

    win = (np.asarray(R_ref_kpc) > r_grav) & (np.asarray(R_ref_kpc) < R0_kpc)
    eps0_low = np.nanmax(eps0[win])
    eps0_hig = np.nanmin(eps0[win])
    eps1 = eps0.copy()
    eps1[np.asarray(R_ref_kpc) < r_grav] = eps0_low
    eps1[np.asarray(R_ref_kpc) > R0_kpc] = eps0_hig

    q1 = np.abs(ve2 - vh2*eps1**2) / (vh2*eps1**2 + ve2)
    g1 = hmg_gamma0(q1)
    a_pred1  = np.sqrt(acc_newton**2 + 2*acc_newton*((c0/T0)/g1))   # total
    a_pred_s = acc_newton + (c0/T0)/g1                             # spatial only

    q0 = np.abs(ve2 - vh2*eps0**2) / (vh2*eps0**2 + ve2)
    g0 = hmg_gamma0(q0)
    a_pred0 = np.sqrt(acc_newton**2 + 2*acc_newton*((c0/T0)/g0))    # balanced

    a_MOND  = np.sqrt(acc_newton**2 + 2*acc_newton*((c0/T0)/7.5))
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        a_MOND2 = acc_newton / (1 - np.exp(-np.sqrt(acc_newton/((c0/T0)/5.5))))

    return dict(acc_newton=acc_newton, a_pred1=a_pred1, a_pred_s=a_pred_s,
                a_pred0=a_pred0, a_MOND=a_MOND, a_MOND2=a_MOND2,
                eps0=eps0, eps1=eps1, g1=g1, g0=g0)


# ---- FITS loaders ---------------------------------------------------------
def read_fits_table(path):
    """Read a binary FITS table into a dict of float64 numpy columns."""
    with fits.open(path) as h:
        d = h[1].data
        cols = {name: np.asarray(d[name], dtype=float) for name in d.columns.names}
    return cols


def cluster_fits(subdir, clust):
    """List .fits files inside data/<subdir>/<clust>, sorted by name."""
    p = os.path.join(DATA, subdir, clust)
    return [os.path.join(p, f) for f in sorted(os.listdir(p)) if f.endswith(".fits")]


# ---- four_clusters.txt (Angus 2008) ---------------------------------------
COLS_ANGUS = ["clust", "radkpc", "total", "thn_bnk", "thk_bnk",
              "thn_blu", "thk_blu", "green", "red"]


def _parse_blocks(path):
    """Read four_clusters.txt into per-cluster ordered blocks (drop # lines)."""
    blocks = {}
    order = []
    with open(path, "r", encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            row = [name] + [np.nan if (v.strip() in ("NA", "")) else float(v)
                            for v in parts[1:9]]
            if name not in blocks:
                blocks[name] = []
                order.append(name)
            blocks[name].append(row)
    return blocks, order


def load_angus(path=None):
    """
    Build the 87-row Angus (2008) table used by the four-cluster figure.

    The row layout is a 15-row prefix (N533 rows with radius >= 40, then the
    first 4 N5044 rows) followed by the full per-cluster blocks in the order
    N533, N5044, A2717, A2029. This ordering fixes which points enter each
    per-cluster fit. Returns a dict of numpy arrays keyed by COLS_ANGUS.
    """
    if path is None:
        path = os.path.join(DATA, "four_clusters.txt")
    blocks, _ = _parse_blocks(path)

    n533  = blocks["N533"]
    n5044 = blocks["N5044"]
    a2717 = blocks["A2717"]
    a2029 = blocks["A2029"]

    prefix = [r for r in n533 if not np.isnan(r[1]) and r[1] >= 40] + n5044[:4]
    rows = prefix + n533 + n5044 + a2717 + a2029

    out = {c: [] for c in COLS_ANGUS}
    for r in rows:
        for c, v in zip(COLS_ANGUS, r):
            out[c].append(v)
    out["clust"] = np.array(out["clust"], dtype=object)
    for c in COLS_ANGUS[1:]:
        out[c] = np.array(out[c], dtype=float)
    return out
