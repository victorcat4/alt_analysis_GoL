# =============================================================================
# Package Setup
# =============================================================================

required_packages <- c(
  "tidyverse",
  "sf",
  "readxl",
  "lubridate",
  "dbscan",
  # Parallel processing
  "future",
  "furrr",
  "progressr"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

theme_set(theme_bw())
message("  Packages loaded successfully")