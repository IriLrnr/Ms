source("libraries.R")

# Add centroids to data
edit_shp <- function (data) {
  # EVI Indexes on the sample
  data$EVI_nodata <- (data$EVI_NODATA)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SVI <- (data$B)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SUVI <- (data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SGVI <- (data$C)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SEVI <- ((data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)+1)/2
  data$SBVI <- abs(data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  
  clean_data <- data
  if (length(data$EVI_mean) > 0) {
    clean_data <- data[!is.na(data$EVI_mean),]
    clean_data <- clean_data[clean_data$EVI_nodata < 0.3,]
  }
  
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  #clean_data$BH[is.na(clean_data$BH)] <- 0
  #clean_data$VH[is.na(clean_data$VH)] <- 0
  
  #clean_data$Cat <- cut(clean_data$EVI_mean, breaks = c(0, 0.16, 0.56, 1), 
  #                 labels = c("SU", "SH", "SG"), include.lowest = TRUE)
  
  return(clean_data)
}

# not using anymore
analyse_model <- function(data, model, moran = T) {
  if (ad.test(residuals(model))$p.value < 0.05) {
    sk <- skewness(residuals(model))
    print(sk)
  }
  if (moran == T) {
    listw <- nb2listw(poly2nb(data), style = "W", zero.policy = T)
    print(moran.test(residuals(model), listw))
  }
  return(AIC(model))
} 

# not using much
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

plot_bi_map <- function(d, res, pred, pal = "PinkGrn") {
  
  d$response <- as.data.frame(d)[,c(res)]
  d$predictor <- as.data.frame(d)[,c(pred)]
  
  sub_data <- subset(d, select = c(response, predictor))
  sub_data <- bi_class(sub_data, x = response, y = predictor, style = "fisher", dim = 3)
  
  print(table(as.factor(sub_data$bi_class)))
  
  # Plotting using ggplot2 with sf support
  plot <- ggplot(sub_data) +
    geom_sf(aes(fill = bi_class), color = NA, size = 0.1, show.legend = FALSE) +  # Remove borders with colour = NA
    bi_scale_fill(pal = pal) +
    theme_minimal() +
    labs(x = "", y = "", title = "Bivariate Map of T_mean and prop_veg") +
    theme_minimal() +
    custom_theme

  # Create the legend separately
  legend <- bi_legend(pal = pal,
                      xlab = res,
                      ylab = pred,
                      size = 10,
                      arrows = FALSE)
  
  finalPlot <- ggdraw() +
    draw_plot(plot, 0, 0, 1, 1) +
    draw_plot(legend, 0.8, 0.2, 0.2, 0.4)
  
  finalPlot
  
  #return(sub_data)
  
  # Display the counts
  #print(table(as.factor(sdb$bi_class)))
}

# fitme models
spatial_models <- function(data_df, response, predictors, family=gaussian) {
  # Create a list of formulas based on predictors adding random effects
  formulas <- lapply(predictors, function(x) {
    as.formula(paste(response, "~", x, "+ Matern(1 | centx + centy)"))
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

get_AIC <- function(model) {
  as.numeric((AIC.HLfit(model))[1])
}

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

# GAMM models
gam_spatial_models <- function(data_df, response, predictors) {
  formulas <- lapply(predictors, function(x) {
    as.formula(paste(response, "~", x, "+ s(centx, centy)"))
  })
  
  models <- lapply(formulas, function(f) {
    gam(f, data = data_df, method = "ML")
  })
  
  # Name the models for easy reference
  names(models) <- predictors
  
  # Return the list of models
  return(models)
}

AIC_tab_gam <- function (models) {
  # Initialize lists to store results
  AICs <- sapply(models, function(x) AIC(x), simplify = FALSE)  # Assuming use of lme part for AIC
  
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

