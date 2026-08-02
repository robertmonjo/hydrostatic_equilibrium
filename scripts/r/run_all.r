# =====================================================================
# run_all.r  --  regenerate every figure and table of
#   Monjo (2025), ApJ 981, 195, from data/ using the R pipeline.
# Usage:  Rscript run_all.r
# =====================================================================
sdir <- tryCatch(dirname(normalizePath(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)[1]))),
        error = function(e) ".")
message("== Supplementary four-cluster figures (Angus 2008) ==")
source(file.path(sdir, "10_figure0_four_clusters.r"))
message("== Figure 1 (five complete X-COP clusters) ==")
source(file.path(sdir, "20_figure1_complete.r"))
message("== Figure 2 + Table 2 (twelve X-COP clusters) ==")
source(file.path(sdir, "30_figure2_table2.r"))
message("== DONE. Outputs in figures/. ==")
