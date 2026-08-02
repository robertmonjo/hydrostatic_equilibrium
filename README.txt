Hydrostatic equilibrium in X-COP galaxy clusters (HMG)
======================================================

Code and data to reproduce the figures and table of:

  R. Monjo, "Hydrostatic Equilibrium Constraints in X-COP Galaxy Clusters",
  The Astrophysical Journal 981, 195 (2025).
  https://doi.org/10.3847/1538-4357/adb723

Hyperconical modified gravity (HMG) adds a projection-induced radial
acceleration a0 = 2c/(gamma0 * t). This code tests it against the hydrostatic
mass profiles of the X-COP galaxy clusters and compares it with the Newtonian
and MOND expectations.


Contents
--------
  data/      X-COP mass profiles (FITS) and digitised cluster data
  scripts/   python/ and r/ pipelines; run_all regenerates everything
  figures/   output figures (PNG + PDF); R files and their *_python counterparts
  tables/    Table 2 (per-cluster r_nei, chi^2, p-value, significance) as CSV


Run (from this folder)
----------------------
  python scripts/python/run_all.py
  Rscript scripts/r/run_all.r

Both write to figures/ and tables/. No network access is required.

Requirements:
  Python 3 with numpy, scipy, matplotlib, astropy (reads the FITS profiles).
  R with base packages + FITSio.


Outputs
-------
  Main paper figures
    figures/Figure1   five complete X-COP clusters: baryonic, missing and
                      total mass, with MOND and HMG
    figures/Figure2   twelve X-COP clusters, gas-only HMG fit
    tables/table2_hmg_fit         per-cluster r_nei, chi^2, p-value and
                      significance (Table 2)

  Supplementary figures
    figures/supplementary_fig_angus_four_clusters            Angus et al. (2008)
    figures/supplementary_fig_angus_four_clusters_scheme_b   four-cluster
                      0-parameter HMG test (two colour schemes)

  Each figure is written as PNG + PDF; the R pipeline writes the plain names
  and the Python pipeline writes the matching *_python names.


Data sources (bundled in data/)
-------------------------------
  data_E22/, einasto_mass_profiles/   X-COP mass profiles (FITS).
      Eckert et al. (2022), A&A 662, A123 (arXiv:2205.01110).
  four_clusters.txt   Angus, Famaey & Buote (2008), MNRAS 387, 1470.
  Chandra.txt         Abell 2029 profile (Horne 2006).
  table1_xcop_properties.csv   Table 1 (literature values), from Monjo (2025),
      adapted from Eckert et al. (2022).

  See data/PROVENANCE.md for the full manifest: per-file source, citation
  (with DOI/URL), epoch, reformatting notes, SHA-256 hashes, and which files
  the pipeline consumes versus which are bundled only for reference.

Datasets are the property of their original authors and are redistributed here
for reproducibility with attribution.


Methods notes
-------------
  r_nei uncertainty (Table 2, r_nei_err). This is a heuristic sensitivity
  measure, not a formal confidence interval: the sample standard deviation of
  the best-fit r_nei across the three outer-radius window definitions
  (R>=500, >=1000, >=1500 kpc), combined in quadrature with a 10 kpc floor set
  by the R0 scan grid step. It quantifies how r_nei depends on the chosen
  fitting window.

  Degrees of freedom (Table 2). The tabulated chi^2 is computed over the
  R>=1000 kpc window (per-cluster ~22-31 bins). The p-value column uses
  per-cluster degrees of freedom, df = n_i - 1. The 95% significance flag
  compares chi^2 to the fixed threshold chi^2(95%, df=42) = 58.12, which
  corresponds to the paper's n=43 (R>=500 kpc) convention.
