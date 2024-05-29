source("functions.R")

city = "São Paulo"

# Set individual predictors
predictors <- c("SVI", 
                "SGVI", 
                "SUVI", 
                "prop_veg")

# Set all models with predictor combinations
pred_1km <- c("1", predictors, 
              "SVI + prop_veg", 
              "SVI * prop_veg",
              "SVI + SGVI",
              "SVI * SGVI",
              "SVI + SUVI",
              "SVI * SUVI",
              "SVI * SEVI",
              "SGVI * SEVI",
              "SUVI * SEVI")

# set response variable
rich <- "richness"

# Load data
data_rich <- edit_shp(st_read("../Cities/SP/1km/richness_effort_1km.shp"))

## Set random effect to effort
data_rich$sp_codigo <- data_rich$effort

#plot map of samples coloured by prop_veg
plot_map(data_rich)

# see correlations
cor_mat <- as.data.frame(data_rich)[c(rich, predictors)]
suppressWarnings(chart.Correlation(cor_mat, histogram=TRUE, pch=16, method = "pearson"))

# model with fitme and spatial correlation matrix
richn <- spatial_models(data_rich, rich, pred_1km, "poisson")
# select best models (AIC < 2)
best_rich <- select_best_models(richn)

# AIC table - cam be saved if needed with write csv
AIC_eff_tab(richn, "Richness")

# Plot linear regression by SEVI category
plot_categories(data_rich, rich)

# plot linear regression by vegetation amount category
plot_PV_categories(data_rich, rich)

# plot richest samples
EVI <- "../Cities/SP/EVI/EVI_SP.tif"
plot_top_3(data_rich, rich, EVI, "id")

