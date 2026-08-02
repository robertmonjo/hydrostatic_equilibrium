# Data provenance

This file documents the origin, citation, formatting and integrity of every
file under `data/`. All datasets remain the property of their original authors
and are redistributed here for reproducibility with attribution.

For each file we give: source and citation (with DOI/URL where known), version
or epoch if known, whether it was reformatted from the original, whether the
pipeline consumes it, and a SHA-256 hash.

## Primary sources

- **Eckert et al. (2022)**, "Hydrostatic mass profiles of X-COP galaxy
  clusters", Astronomy & Astrophysics 662, A123.
  arXiv:2205.01110. DOI: 10.1051/0004-6361/202142905.
  X-COP project: https://www.astro.unige.ch/xcop/
- **Angus, Famaey & Buote (2008)**, "X-ray group and cluster mass profiles in
  MOND: unexplained mass on the group scale", MNRAS 387, 1470.
  DOI: 10.1111/j.1365-2966.2008.13332.x.
- **Monjo (2025)**, "Hydrostatic Equilibrium Constraints in X-COP Galaxy
  Clusters", The Astrophysical Journal 981, 195.
  DOI: 10.3847/1538-4357/adb723.

## Pipeline usage summary

Files consumed by the scripts:

- `einasto_mass_profiles/<cluster>/EIN3_mass.fits` — Figure 2 and Table 2.
- `data_E22/<cluster>/EIN3_bar_mass.fits` — Figure 1.
- `data_E22/<cluster>/gbar_gobs_<cluster>.fits` — Figure 1.
- `four_clusters.txt` — supplementary four-cluster figures.

Files bundled for reference / completeness but **not consumed by the
pipeline**:

- `data_E22/<cluster>/GP_thermo.fits`
- `data_E22/<cluster>/<cluster>_BCG.txt`
- `Chandra.txt`
- `table1_xcop_properties.csv`

---

## einasto_mass_profiles/

**Source:** Eckert et al. (2022), X-COP. Einasto (EIN3) hydrostatic mass
profiles per cluster: total, dark ("missing") and gas mass columns on a radial
grid, with high/low error variants. FITS binary tables. Epoch of the FITS
files: 2022-04-25.
**Reformatted:** repackaged into one `EIN3_mass.fits` per cluster directory;
column contents are as produced by the X-COP analysis.
**Consumed by:** `fig2_table2.py`, `30_figure2_table2.r` (Figure 2, Table 2).

| File | SHA-256 |
|------|---------|
| einasto_mass_profiles/A1644/EIN3_mass.fits | a63b3169f4712355c26dfa61a63b9abd95fdde64af94626de9dfc392b36b5cea |
| einasto_mass_profiles/A1795/EIN3_mass.fits | 86b0d9bc8853dcfc24e2b58cfb8a8f012c134fc932eb679216f5c828d3296497 |
| einasto_mass_profiles/A2029/EIN3_mass.fits | c7126585e447c0828392cf8b4baac5aee053585ae60b492b0da6bfba2eee26e9 |
| einasto_mass_profiles/A2142/EIN3_mass.fits | a6bdba4d8ce3ae63cf08f5afc411094162bd7c2110122db86bfa044c0aa98574 |
| einasto_mass_profiles/A2255/EIN3_mass.fits | c42ea08b8fae913f13f513012f9389d44dbdf19c0e78a8c5d93a58a663704a7a |
| einasto_mass_profiles/A2319/EIN3_mass.fits | bbb5a688173b95eeaa47d42a7bde0e65bdf57dd90072b8af304ef10999fa8877 |
| einasto_mass_profiles/A3158/EIN3_mass.fits | f8cc2192e90d82a77ca982b2521d4dae5a8d774deb118374b7cab9cd36371637 |
| einasto_mass_profiles/A3266/EIN3_mass.fits | 37682784ae5519ba4875b9052056fb866bc80ce8385cb5ab76ab6018d56a5d17 |
| einasto_mass_profiles/A644/EIN3_mass.fits  | aa4179309f903379fe3e31acead0e6e2858c979874bd89a7e4796253bb85ea58 |
| einasto_mass_profiles/A85/EIN3_mass.fits   | 6958fad8a2e2282603dd8b5886c314d40431c075fbd7acbf120e3151a8c2ca31 |
| einasto_mass_profiles/RXC1825/EIN3_mass.fits | 919ca21d1ea6440838d6612c8ef6d141f21424097e03a08c62682b147be7b34b |
| einasto_mass_profiles/ZW1215/EIN3_mass.fits | c4b6ab1810d28ee98b6f692e949312a4a4bb62060e6e993a2f57b9d6b7794ac2 |

## data_E22/

**Source:** Eckert et al. (2022), X-COP. Per-cluster subdirectories for the five
"complete" clusters used in Figure 1. FITS epochs: `EIN3_bar_mass.fits` and
`GP_thermo.fits` 2022 (2022-05-31 / 2022-04-25), `gbar_gobs_*.fits` 2023-04-20.
**Reformatted:** repackaged X-COP products.

### EIN3_bar_mass.fits — baryonic and total mass profiles
**Consumed by:** `fig1_complete.py`, `20_figure1_complete.r` (Figure 1).

| File | SHA-256 |
|------|---------|
| data_E22/A1795/EIN3_bar_mass.fits | 1095e52ca61f9e83718e4237fe8c8d8cfe5cd17f39c8a692cd1f2a007e69c281 |
| data_E22/A2029/EIN3_bar_mass.fits | f76d5ea2bc105aea218ada34dd06a6daaac6be42309c7fd560f7105fb28bee23 |
| data_E22/A2142/EIN3_bar_mass.fits | aaeaf4978a17ec970e5c34cf575e0abc68f5675a3aa7f2d53e829b929c085a85 |
| data_E22/A2319/EIN3_bar_mass.fits | bbcccd3d7a2173135ad9a0fb38b286eabd2f6cba61e63872925d4db0327e2e57 |
| data_E22/A644/EIN3_bar_mass.fits  | 34f8efc307989effee29d5e768d56dfb6fc9db963274878b795248e2763bb5b3 |

### gbar_gobs_<cluster>.fits — observed inward acceleration with error band
**Consumed by:** `fig1_complete.py`, `20_figure1_complete.r` (Figure 1).

| File | SHA-256 |
|------|---------|
| data_E22/A1795/gbar_gobs_A1795.fits | 8a6ff9f7ae360047aa0f3be5c51442d0e217f989e81c15e75e7fa54766c4982c |
| data_E22/A2029/gbar_gobs_A2029.fits | 8ed65afd32061e822e45730338ce55ae34b9098eecfddf9fba4cde5f88d83ed4 |
| data_E22/A2142/gbar_gobs_A2142.fits | f700885754486cc6bddca599cb71b7e59431f9dba48135f88f58434bf0b7cb87 |
| data_E22/A2319/gbar_gobs_A2319.fits | cb89b6b3a0fc3c4fa361a671713b33b47dc35dca9b837c17ad16cac9cc6c5cd9 |
| data_E22/A644/gbar_gobs_A644.fits   | 6994602130138d7f67ee1bdc542143efaa8e3788c0673968fb90e911a9da4b07 |

### GP_thermo.fits — gas thermodynamic profiles
Bundled for reference / completeness; **not consumed by the pipeline**.

| File | SHA-256 |
|------|---------|
| data_E22/A1795/GP_thermo.fits | 1de7b2d36c4dcef3505ad5df672788409c46190a8699d5c3ecdeafeaa3ef65e8 |
| data_E22/A2029/GP_thermo.fits | 0aa85a3b7f9eb3b86ad12aaef2a04f16bc5e40e3fd5c12d8c6f4db5a62e77a2b |
| data_E22/A2142/GP_thermo.fits | fa648295427d80fbe4da39f421255c5dc87a22f75409241a6950f659f6912396 |
| data_E22/A2319/GP_thermo.fits | 1022d575889912e5ac997ec8a4f7494e705a2918c71a9cc035e075aaac82be2e |
| data_E22/A644/GP_thermo.fits  | 16a7f339b6dfd48171fa8d8da551341faab7463e41ce7e421318cf20097a92e6 |

### <cluster>_BCG.txt — brightest cluster galaxy profile
Whitespace-separated text, no header. Bundled for reference / completeness;
**not consumed by the pipeline**.

| File | SHA-256 |
|------|---------|
| data_E22/A1795/A1795_BCG.txt | 9e7c85d9476422450e1182fe0b418e12989f61bfa32984b361a1c16171189137 |
| data_E22/A2029/A2029_BCG.txt | 6e5d39efd46408ba7dfcf3667e8d454b2d5b5c52d6a6c9ae370fb0d722e59d38 |
| data_E22/A2142/A2142_BCG.txt | dfacfb90b28145f01dd187b866b6fdca09dc0f919bda4099aab3c507ebf9bc39 |
| data_E22/A2319/A2319_BCG.txt | e596f88fb5de64db070596ebbbaec477a5065b57360ce8be7b4e3dd4d573a667 |
| data_E22/A644/A644_BCG.txt   | 5268448f6888ef04754a5d47b25642c495fd787202fd224745c95470ca8ea3fc |

## four_clusters.txt

**Source:** Angus, Famaey & Buote (2008), MNRAS 387, 1470. Digitised mass
profiles for four systems (N533, N5044, A2717, A2029) read from the paper's
figures. Tab-separated; column names appear only in a leading `#`-comment line
that the loaders skip, so the columns are read positionally (no consumed header
row).
**Reformatted:** digitised from published figures into a plain-text table.
**Consumed by:** `fig0_four_clusters.py`, `10_figure0_four_clusters.r`
(supplementary four-cluster figures).

| File | SHA-256 |
|------|---------|
| four_clusters.txt | 4480d1d46046d34a03809f28141a3a6e0cc3c2045c3408665c34c5b5020f2990 |

## Chandra.txt

**Source:** Abell 2029 acceleration/mass profile (Horne 2006). Tab-separated
with a header row. Bundled for reference / completeness; **not consumed by the
pipeline**.

| File | SHA-256 |
|------|---------|
| Chandra.txt | 85acdeb0d466693c21ab7d4f530231ba611aa7ea99414573cdab96370d1fcc3c |

## table1_xcop_properties.csv

**Source:** Table 1 of Monjo (2025), adapted from Eckert et al. (2022):
redshift, column density, R500, M500 (total and baryonic) and completeness flag
per cluster. CSV with a header row.
**Reformatted:** transcribed from the published Table 1. Bundled for reference /
completeness; **not consumed by the pipeline** (the manuscript Table 1 is
typeset directly, not generated from this file).

| File | SHA-256 |
|------|---------|
| table1_xcop_properties.csv | 2dffa634f7af2bed403ba0d9448335a94fdf1fd39d178ac3dee766b31bbc84b6 |
