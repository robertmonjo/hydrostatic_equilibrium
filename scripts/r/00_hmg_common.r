# =====================================================================
# 00_hmg_common.r
# Shared physical constants and Hyperconical Modified Gravity (HMG)
# model functions used by the figure and table scripts of
#   Monjo, R. (2025), ApJ 981, 195. DOI 10.3847/1538-4357/adb723
#   "Hydrostatic Equilibrium Constraints in X-COP Galaxy Clusters"
#
# Defines the HMG acceleration model, physical constants and path
# helpers. Reads data from ../../data (relative to this script).
# =====================================================================

## ---- Physical constants (SI unless noted) ----------------------------
c0   <- 3e8                                   # speed of light, m/s (paper value)
Msol <- 1.9891e30                             # solar mass, kg
kpc  <- 3261.8116478174 * 365*24*3600 * c0    # kpc in metres (light-yr based)
pc   <- kpc / 1000
kms  <- 1000
T0   <- 13.7e9 * 365*24*3600                  # age of universe, s (13.7 Gyr)
GN   <- 6.674e-11                             # gravitational constant, SI

## HMG vacuum / critical density for H = 1/t:  rho_vac = 3/(8 pi G t^2)
rho_vac <- 3 / (8*pi*GN*T0^2)

## Universal HMG angles
GA_EMPTY  <- pi/3   # gamma_U  (universal constant of the model)
GA_CENTER <- pi/2   # gamma_center fixed to pi/2 -> 1-parameter model

## ---- HMG projection factor gamma_0(r) --------------------------------
## quotient q(r) = | v_N^2 - eps^2 v_H^2 | / ( eps^2 v_H^2 + v_N^2 )
## sin^2 gamma_sys = sin^2(gaU) + (sin^2(gac) - sin^2(gaU)) * q
## gamma_0 = gamma_sys / cos(gamma_sys)
hmg_gamma0 <- function(quotient, gaempty = GA_EMPTY, gacenter = GA_CENTER) {
  g_sys <- asin(sqrt(sin(gaempty)^2 +
                     (sin(gacenter)^2 - sin(gaempty)^2) * quotient))
  g_sys / cos(g_sys)
}

## ---- Full HMG prediction on a radial grid ----------------------------
## M_enc_kg : enclosed (baryonic) mass in kg, vector over the grid
## R_ref_kpc: reference radius in kpc, same length
## R0_kpc   : r_nei fitting parameter (kpc)
## r_grav   : inner gravitational-domination scale (kpc)
## Returns list of acceleration vectors (m/s^2).
hmg_predict <- function(M_enc_kg, R_ref_kpc, R0_kpc, r_grav = 40,
                        add_e02 = TRUE) {
  R_ref_m    <- R_ref_kpc * kpc
  acc_newton <- GN * M_enc_kg / R_ref_m^2
  ve2        <- 2 * acc_newton * R_ref_m          # v_N^2 = 2 G M / r
  vh2        <- R_ref_m^2 / T0^2                   # v_H^2 = (r/t)^2
  dens       <- M_enc_kg / (4/3 * pi * R_ref_m^3)

  e02  <- if (add_e02) 1/6 else 0
  eps0 <- sqrt(e02 + dens / rho_vac)              # unclamped eps (balanced)

  ## clamped eps for the fitted 1-parameter model (r_nei = R0)
  win <- R_ref_kpc > r_grav & R_ref_kpc < R0_kpc
  eps0_low <- max(eps0[win])
  eps0_hig <- min(eps0[win])
  eps1 <- eps0
  eps1[R_ref_kpc < r_grav] <- eps0_low
  eps1[R_ref_kpc > R0_kpc] <- eps0_hig

  ## fitted model (clamped)
  q1 <- abs(ve2 - vh2*eps1^2) / (vh2*eps1^2 + ve2)
  g1 <- hmg_gamma0(q1)
  a_pred1  <- sqrt(acc_newton^2 + 2*acc_newton*((c0/T0)/g1))   # total (centrifugal)
  a_pred_s <- acc_newton + (c0/T0)/g1                          # spatial only

  ## balanced model gamma_sys = pi/3 (unclamped eps0)
  q0 <- abs(ve2 - vh2*eps0^2) / (vh2*eps0^2 + ve2)
  g0 <- hmg_gamma0(q0)
  a_pred0 <- sqrt(acc_newton^2 + 2*acc_newton*((c0/T0)/g0))    # dotted balanced line

  ## MOND
  a_MOND  <- sqrt(acc_newton^2 + 2*acc_newton*((c0/T0)/7.5))          # theoretical
  a_MOND2 <- acc_newton / (1 - exp(-sqrt(acc_newton/((c0/T0)/5.5))))  # empirical MLS

  list(acc_newton = acc_newton,
       a_pred1 = a_pred1, a_pred_s = a_pred_s, a_pred0 = a_pred0,
       a_MOND = a_MOND, a_MOND2 = a_MOND2,
       eps0 = eps0, eps1 = eps1, g1 = g1, g0 = g0)
}

## ---- Paths -----------------------------------------------------------
hmg_paths <- function() {
  ## resolve data / figures dirs relative to this script location
  args <- commandArgs(trailingOnly = FALSE)
  fa   <- grep("--file=", args, value = TRUE)
  if (length(fa)) {
    here <- dirname(normalizePath(sub("--file=", "", fa[1])))
  } else {
    here <- normalizePath(".")
  }
  root <- normalizePath(file.path(here, "..", ".."))
  list(root = root,
       data = file.path(root, "data"),
       fig  = file.path(root, "figures"),
       tab  = file.path(root, "tables"))
}

## qchisq threshold reported in the paper (n = 43 bins)
CHI2_THRESHOLD_95 <- qchisq(0.95, df = 42)   # ~ 58.12
