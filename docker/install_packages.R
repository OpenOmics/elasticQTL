#!/usr/bin/env Rscript
# install_packages.R
# Auto-generated from sessionInfo.txt
# Installs all required packages at the versions listed in the session info.
#
# R version required: 4.5.0
# Platform: x86_64-pc-linux-gnu (Rocky Linux 8.7)

# --- Install remotes (needed for install_version) ---
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# --- Helper: install a specific version if not already installed at that version ---
install_if_needed <- function(pkg, ver, repos = "https://cran.r-project.org") {
  current <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA)
  if (is.na(current) || current != ver) {
    message(sprintf("Installing %s %s ...", pkg, ver))
    remotes::install_version(pkg, version = ver, repos = repos, upgrade = "never")
  } else {
    message(sprintf("  %s %s already installed — skipping.", pkg, ver))
  }
}

# ============================================================
# Attached packages (non-base)
# ============================================================
install_if_needed("Matrix",     "1.7-3")
install_if_needed("glmnet",     "4.1-10")
install_if_needed("data.table", "1.18.2.1")
install_if_needed("dplyr",      "1.1.4")
install_if_needed("tidyr",      "1.3.2")
install_if_needed("ggplot2",    "4.0.2")
install_if_needed("optparse",   "1.7.5")

# ============================================================
# Namespace-loaded (dependency) packages
# ============================================================
install_if_needed("vctrs",        "0.6.5")
install_if_needed("cli",          "3.6.5")
install_if_needed("rlang",        "1.1.7")
install_if_needed("purrr",        "1.2.1")
install_if_needed("generics",     "0.1.4")
install_if_needed("S7",           "0.2.1")
install_if_needed("glue",         "1.8.0")
install_if_needed("scales",       "1.4.0")
install_if_needed("tibble",       "3.3.1")
install_if_needed("foreach",      "1.5.2")
install_if_needed("lifecycle",    "1.0.5")
install_if_needed("getopt",       "1.20.4")
install_if_needed("codetools",    "0.2-20")
install_if_needed("RColorBrewer", "1.1-3")
install_if_needed("Rcpp",         "1.1.1")
install_if_needed("pkgconfig",    "2.0.3")
install_if_needed("farver",       "2.1.2")
install_if_needed("lattice",      "0.22-6")
install_if_needed("R6",           "2.6.1")
install_if_needed("dichromat",    "2.0-0.1")
install_if_needed("tidyselect",   "1.2.1")
install_if_needed("pillar",       "1.11.1")
install_if_needed("shape",        "1.4.6.1")
install_if_needed("magrittr",     "2.0.4")
install_if_needed("withr",        "3.0.2")
install_if_needed("gtable",       "0.3.6")
install_if_needed("iterators",    "1.0.14")
install_if_needed("survival",     "3.8-3")

# ============================================================
# Verify installation
# ============================================================
message("\n--- Verification ---")
pkgs <- c(
  "Matrix", "glmnet", "data.table", "dplyr", "tidyr", "ggplot2", "optparse",
  "vctrs", "cli", "rlang", "purrr", "generics", "S7", "glue", "scales",
  "tibble", "foreach", "lifecycle", "getopt", "codetools", "RColorBrewer",
  "Rcpp", "pkgconfig", "farver", "lattice", "R6", "dichromat", "tidyselect",
  "pillar", "shape", "magrittr", "withr", "gtable", "iterators", "survival"
)
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  message(sprintf("  %-15s %s", p, v))
}
message("\nDone.")
