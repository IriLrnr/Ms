source("functions2.R")
# set city context
shp <- "../Cities/SP/1km/T_two_ES.shp"
data <- edit_shp(st_read(shp))

# testing
library(spgwr)

library(vegan)
library(car)
library(mgcv)
library(e1071)
library(nortest)
library(plotly)

data <- as(data, "Spatial")
GWRbandwidth <- gwr.sel(TMEAN ~ prop_veg, data = data, adapt = T)



######### INLA

# set city context
shp <- "../Cities/SP/1km/T_two_ES.shp"
out <- "DT"
city <- "São Paulo"

# Load data
data <- edit_shp(st_read(shp))

data <- data[!is.na(data$DT_wait),]
data$RM_fac <- factor(data$RM) 
data <- data[data$RM < 15,]
data_sf <- st_transform(data, crs = 31983)
d<-data_sf
plot_map(d)
d$index <- 1:nrow(d)


#density plot

ggplot(cdata, aes(x = T_mean, fill = factor(RM))) +
  geom_density(alpha = 0.5) +
  labs(x = "Mean Temperature", y = "Density", fill = "RM") +
  custom_theme


# varpart

data_df <- as.data.frame(cdata)
part.lm = varpart(data_df$LTM, data_df[, c("EVI_mean")], data_df[, c("prop_veg")], data_df[, c("SVI", "SGVI", "SUVI")])
# Plot with custom settings
plot(part.lm, digits = 2,
     bg = c(rgb(48, 225, 210, 80, maxColorValue = 255), 
            rgb(80, 48, 225, 210, maxColorValue = 255),
            rgb(210, 80, 48, 225, maxColorValue = 255)))

out1 <- rda(data_df$DT, data_df[, c("EVI_mean")], data_df[, c("prop_veg")])
out2 <- rda(data_df$DT, data_df[, c("prop_veg")], data_df[, c("EVI_mean")])

anova(out1, step = 1000, perm.max = 1000)
anova(out2, step = 1000, perm.max = 1000)


optimize_lag_distance <- function(data, response, predictor, distances) {
  library(spdep)
  
  coords <- st_centroid(cdata)
  results <- data.frame(distance = distances, AIC = NA)
  
  for (i in seq_along(distances)) {
    dist <- distances[i]
    
    # Create neighbor list and weights matrix
    nb <- dnearneigh(coords, i, i+100) 
    listw <- nb2listw(nb, style="W", zero.policy=TRUE)
    
    MC <- moran.mc(cdata$EVI_mean, lw, nsim=99, zero.policy=T) 
    plot(MC, main="", las=1) 
    
    
    # Calculate spatial lag for the predictor
    data[[paste0(predictor, "_lag")]] <- lag.listw(listw, data[[predictor]])
    
    # Fit a linear model
    formula <- as.formula(paste(response, "~", predictor, "+", paste0(predictor, "_lag")))
    model <- lm(formula, data = data)
    
    # Store AIC value
    results$p[i] <- MC$p.value
    print(results$p[i])
  }
  
  print(results$distance[which.min(results$p)])
  
  optimal_distance <- results$distance[which.min(results$p)]
  
  return(optimal_distance)
}

# Example usage
distances <- seq(100, 1000, by = 100)

# Optimize lag distance for each predictor
optimal_distance <- optimize_lag_distance(cdata, "LTM", "EVI_mean", distances)
optimal_distance_SGVI <- optimize_lag_distance(cdata, "T_mean", "SGVI", distances)
optimal_distance_SUVI <- optimize_lag_distance(cdata, "T_mean", "SUVI", distances)

print(optimal_distance_SVI)
print(optimal_distance_SGVI)
print(optimal_distance_SUVI)


library(gstat)

variogram_model <- variogram(T_mean ~ 1, cdata)

# Fit a variogram model to determine the range             
vgm_model <- fit.variogram(variogram_model, vgm(psill = 1, model = "Exp", nugget = 0.1, range = 10000))
plot(variogram_model, model = vgm_model)
