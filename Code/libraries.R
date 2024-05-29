load_or_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

packages <- c("ggplot2", "viridis", "sf", "sp", "spdep",
              "PerformanceAnalytics", "gridExtra", "biscale", 
              "cowplot", "randomForest", "pdp")

# List of packages to be loaded/installed
packages2 <- c("raster", "terra", "vegan", "car", "spaMM", "lme4", "effects")

packages_notused <- c(
  "lmerTest", "regclass", 
  "data.table", "scales", "visreg", 
  "spgwr", "rgl", 
  "Matrix", "regclass", "Rcpp", "RcppEigen",
  "reshape2"
)

load_or_install_packages(packages)