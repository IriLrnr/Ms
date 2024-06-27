source("libraries.R")

# Add centroids to data
edit_shp <- function (data) {
  # Exclude polygons with no neighbours from your data
  clean_data <- data
  
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  clean_data$Category <- cut(clean_data$prop_veg, breaks = c(0, 0.3, 0.6, 1), 
                     labels = c("urban: pv < 30%", "mixed: 30% < pv < 60%", "green: pv > 60%"), include.lowest = TRUE)
  
  clean_data$SEVI <- cut(clean_data$prop_veg, breaks = c(0, 0.5, 1), 
                             labels = c("Urban unbalance", "Green unbalance"), include.lowest = TRUE)
  
  return(clean_data)
}

# fitme models
spatial_models <- function(data_df, response, predictors, family=gaussian) {
  # Create a list of formulas based on predictors adding random effects
  formulas <- lapply(predictors, function(x) {
    as.formula(paste(response, "~", x, "+ (1 | sp_codigo) + Matern(1 | centx + centy)"))
  })
  
  # Fit the models using fitme
  
  avail_thr <- parallel::detectCores(logical=FALSE) - 1L 

  models <- lapply(formulas, function(f) {
    fitme(formula = f, data = data_df, family = family,
          fixed=list(longdep=0.5,shape=0.5,rho=0.05),
          control.HLfit=list(NbThreads=max(avail_thr, 1L), 
                             algebra="spcorr"))
  })
  
  # Name the models for easy reference
  names(models) <- predictors
  
  # Return the list of models
  return(models)
}

non_spatial_models <- function(data_df, response, predictors, family=gaussian) {
  formulas <- lapply(predictors, function(x) {
    as.formula(paste(response, "~", x, "+ (1 | sp_codigo)"))
  })
  
  avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
  
  # Fit the models using fitme
  models <- lapply(formulas, function(f) {
    fitme(formula = f, data = data_df, family = family, 
          control.HLfit=list(NbThreads=max(avail_thr, 1L)))
  })
  
  # Name the models for easy reference
  names(models) <- predictors
  
  # Return the list of models
  return(models)
}

select_best_models <- function(models) {
  # Extract AIC values from each model
  AICs <- sapply(models, get_AIC)
  return(models[(AICs - min(AICs)) <= 2])
}

# AIC comparison
AIC_tab <- function (models) {
  
  AICs <- sapply(models, get_AIC, simplify = FALSE)
  # Create data frame with model names and AICs
  AIC_df <- data.frame(Model = names(AICs), AIC = unlist(AICs), row.names = NULL)
  # Calculate delta AICs
  AIC_df$deltaAIC <- AIC_df$AIC - min(AIC_df$AIC)
  # compute AIC weights
  AIC_df$w = exp(-0.5 * AIC_df$deltaAIC) / sum(exp(-0.5 * AIC_df$deltaAIC))  
  
  AIC_df$evidence <- min(AIC_df$AIC) / AIC_df$AIC
  
  # Add city as column
  #AIC_df$city <- rep(city, nrow(AIC_df))
  # Sort the data frame by delta AIC
  AIC_df <- AIC_df[order(AIC_df$deltaAIC),]
  
  return(AIC_df)
}

get_AIC <- function(model) {
  as.numeric((AIC.HLfit(model))[1])
}

# This only work with maximum 3 predictors
AIC_eff_tab <- function(models, name) {
  # Initialize AICs and effect sizes
  AICs <- sapply(models, get_AIC, simplify = FALSE)
  
  AIC_df <- data.frame(Class = name, Model = names(AICs), AIC = round(unlist(AICs), 4), row.names = NULL)
  AIC_df$deltaAIC <- round(AIC_df$AIC - min(AIC_df$AIC), 4)
  
  AIC_df$pred1 <- ""
  AIC_df$pred2 <- ""
  AIC_df$inter <- ""
  
  for (i in 1:length(models)) {
    fe <- fixef(models[[i]])
    se <- sqrt(diag(vcov(models[[i]])))
    tvals <- fe / se
    
    for (j in 2:min(length(fe), 4)) { # Ensure we don't exceed the number of available predictors
      sig <- "" # Initialize sig as an empty string
      
      # Check if tvals[j] is not NA before proceeding with the comparison
      if (!is.na(tvals[j])) {
        if (abs(tvals[j]) < 1.96) sig <- ""
        if (abs(tvals[j]) >= 1.96) sig <- "*"      # p < 0.05
        if (abs(tvals[j]) > 2.576) sig <- "**" # p < 0.01
        if (abs(tvals[j]) > 3.291) sig <- "***"# p < 0.001
      }
      
      # Round the effect size before converting to string and appending significance markers
      rounded_effect_size <- round(as.numeric(fe[j]), 4)
      val_with_sig <- if (!is.na(fe[j])) paste0(rounded_effect_size, sig) else ""
      
      if (j == 2) {
        AIC_df$pred1[i] <- val_with_sig
      } else if (j == 3) {
        AIC_df$pred2[i] <- val_with_sig
      } else if (j == 4) {
        AIC_df$inter[i] <- val_with_sig
      }
    }
  }
  
  AIC_df <- AIC_df[order(AIC_df$deltaAIC),]
  
  return(as.data.frame(AIC_df))
}

# plots linar
plot_PV_categories <- function(data, response) {
  data <- as.data.frame(data)
  data$response <- data[,c(response)]
  ggplot(data = data, 
    aes(x = SVI, y = response, color = Category)) +
    geom_point(shape = 4, size = 1) +  # Add points to the plot
    geom_smooth(method = "glm", formula = y ~ x, se = FALSE) +  # Add regression lines
    theme_minimal() +  # Use a minimal theme for a cleaner look
    labs(x = "SVI", y = response, color = "Category") +  # Label the axes and legend
    scale_color_viridis_d()
    #scale_color_brewer(palette = "Set1")  # Use a color palette for clarity
}

plot_categories <- function(data, response) {
  data <- as.data.frame(data)
  data$response <- data[,c(response)]
  p1 <- ggplot(data = data, 
               aes(x = SVI, y = response, color = SEVI)) +
               geom_point(shape = 4, size = 1) +  # Add points to the plot
               geom_smooth(method = "glm", formula = y ~ x, se = FALSE) +  # Add regression lines
               theme_minimal() +  # Use a minimal theme for a cleaner look
               labs(x = "SVI", y = response, color = "Category") +  # Label the axes and legend
               scale_color_viridis_d() +
               theme(legend.position = "none")  # Hide the legend
  
  p2 <- ggplot(data = data, 
               aes(x = SGVI, y = response, color = SEVI)) +
               geom_point(shape = 4, size = 1) +  # Add points to the plot
               geom_smooth(method = "glm", formula = y ~ x, se = FALSE) +  # Add regression lines
               theme_minimal() +  # Use a minimal theme for a cleaner look
               labs(x = "SGVI", y = response, color = "Category") +  # Label the axes and legend
               scale_color_viridis_d() +
               theme(axis.title.y = element_blank(),
                     legend.position = "none")  # Hide the y-axis title
  
  p3 <- ggplot(data = data, 
               aes(x = prop_veg, y = response, color = SEVI)) +
    geom_point(shape = 4, size = 1) +  # Add points to the plot
    geom_smooth(method = "glm", formula = y ~ x, se = FALSE) +  # Add regression lines
    theme_minimal() +  # Use a minimal theme for a cleaner look
    labs(x = "prop_veg", y = response, color = "Category") +  # Label the axes and legend
    scale_color_viridis_d() +
    theme(legend.position = "none")  # Hide the y-axis title
  
  # Extract just the legend
  # Create a plot that forces a legend to appear, but minimize plot elements
  p_legend <- ggplot(data, aes(x = SVI, y = response, color = SEVI)) +
    geom_point(alpha = 1) +  # Make points transparent
    geom_smooth(method = "glm", formula = y ~ x, se = FALSE) +
    scale_color_viridis_d() +
    theme_minimal() +  # Remove most plot elements
    theme()
  
  # Extract just the legend
  legend <- cowplot::get_legend(p_legend)
  # Draw the legend
  p_legend <- cowplot::plot_grid(legend, align = 'v', axis = 'tb', nrow = 1)

  all_p <- grid.arrange(p1, p2, p3, p_legend, ncol = 2, nrow = 2)
  return (all_p)
}

plot_map <- function(data) {
  ggplot(data) +
    geom_sf(aes(fill = prop_veg)) +
    scale_fill_gradient(low = "red", high = "lightgreen") +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          legend.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank())
}

plot_region_EVI <- function (polygon, raster) {
  
  polygon_st <- st_geometry(polygon)
  clipped_raster <- crop(raster, vect(polygon_st))
  clipped_raster <- mask(clipped_raster, vect(polygon_st))
  plot(clipped_raster, main = paste(polygon$name))
}

plot_top_3 <- function (data, res, path_raster, names, n = 3) {
  data$response <- as.data.frame(data)[,c(res)]
  data$name <- as.data.frame(data)[,c(names)]
  
  top_lines_df <- data[order(-data$response)[1:n],]
  polygons <- st_cast(top_lines_df, "POLYGON")
  
  raster <- rast(path_raster)
  
  # reproject if in different crs
  if (st_crs(polygons) != crs(raster)) {
    raster <- project(raster, crs(polygons))
    polygons <- st_transform(polygons, crs(raster))  # Match CRS to the raster
  }
  
  par(mfrow = c(1, 3))
  
  p1 <- plot_region_EVI(polygons[1,], raster)
  
  p2 <- plot_region_EVI(polygons[2,], raster)
  
  p3 <- plot_region_EVI(polygons[3,], raster)
  
  
  # Reset par settings to defaults
  on.exit(par(mfrow = c(1, 1)))
}

# GAMM models
gam_spatial_models <- function(data_df, response, formulas) {
  formulas <- lapply(formulas, function(x) {
    as.formula(paste(response, "~", x))
  })
  
  models <- lapply(formulas, function(f) {
    
    gam(f, 
         random = list(sp_codigo = ~ 1), 
         data = data_df,
         correlation = corSpatial(form = ~ centx + centy, type = "exponential"),
         method = "ML")
  })
  
  # Name the models for easy reference
  names(models) <- predictors
  
  # Return the list of models
  return(models)
}

# not working
gamm_mrf_models <- function(data_df, response, predictors) {
  # Create neighborhood structure based on coordinates (assuming data_df already includes centx and centy)
  coordinates(data_df) <- ~centx + centy  # This line is just illustrative; use it if you are converting data_df to a spatial object
  nb <- poly2nb(data_df, queen = TRUE)
  
  # Ensure categorical variables are factors, particularly for random effects
  data_df$sp_codigo <- as.factor(data_df$sp_codigo)
  
  # Construct models
  models <- lapply(predictors, function(x) {
    formula_str = paste(response, "~", x, "+ s(centx, centy, bs = 'mrf', xt = list(nb = nb))")
    formula = as.formula(formula_str)
    gamm(formula, random = list(sp_codigo = ~ 1), data = data_df, method = "ML")
  })
  
  # Name the models for easy reference
  names(models) <- predictors
  
  return(models)
}

select_best_gamm_models <- function(models) {
  # Extract AIC values from each model
  AICs <- sapply(models, AIC)
  return(models[(AICs - min(AICs)) <= 2])
}

AIC_tab_gam <- function (models) {
  # Initialize lists to store results
  AICs <- sapply(models, function(x) AIC(x$lme), simplify = FALSE)  # Assuming use of lme part for AIC
  
  # Create data frame with model names, AICs, and delta AICs
  AIC_df <- data.frame(Model = names(AICs), 
                       AIC = unlist(AICs), 
                       row.names = NULL)
  
  # Calculate delta AICs and weights
  AIC_df$deltaAIC <- AIC_df$AIC - min(AIC_df$AIC)
  
  # Sort the data frame by delta AIC
  AIC_df <- AIC_df[order(AIC_df$deltaAIC),]
  
  return(AIC_df)
}


