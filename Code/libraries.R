load_or_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of packages to be loaded/installed
packages <- c("sf", "sp", "spdep", # Handling spatial data
              "dplyr", "tidyr",
              "PerformanceAnalytics", "randomForest", "pdp",
              "ggplot2", "gridExtra", "grid", "biscale", "cowplot", "ggtern", "viridis", # plotting
              "mgcv", "spaMM", "vegan", "car", "mvpart", "caret") # statistics

load_or_install_packages(packages)


packages2 <- c("raster", "terra", "lme4", "effects", "PerformanceAnalytics", "randomForest", "pdp")
