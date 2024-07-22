source("libraries.R")

normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# edit data for analysis
edit_shp <- function (data) {
  data$prop_veg <- (data$PV_1 + data$PV_2 + data$PV_3)/(data$PV_1 + data$PV_2 + data$PV_3 + data$PV_NODATA)
  
  # EVI Indexes on the sample
  data$EVI_nodata <- (data$EVI_NODATA)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SVI <- (data$B)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SUVI <- (data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SGVI <- (data$C)/(data$C + data$B + data$A + data$EVI_NODATA)
  #data$SEVI <- ((data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)+1)/2
  #data$SBVI <- abs(data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  
  #data$pop_d <- data$pop/data$ocup_area
  
  data$SBmean <- data$SBsum/data$SBcount

  clean_data <- data
  #clean_data <- data[!is.na(data$EVI_mean),]
  clean_data <- data[data$EVI_nodata < 0.3,]
  
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  abdata <- clean_data[!is.na(clean_data$SBmean),]
  
  clean_data <- clean_data %>%
    mutate(CR = DT,
           AB = SBmean,
           pop = pop)
  
  clean_data <- clean_data %>%
    mutate(CRn = normalize(clean_data$DT),
           ABn = normalize(clean_data$SBmean),
           popn = normalize(clean_data$pop))
  
  clean_data <- clean_data[,c("id", "CR", "AB", "pop", "CRn", "ABn", "popn", "SVI", "SGVI", "SUVI", "prop_veg", "centx", "centy")]
  
  class(clean_data)
  
  return(clean_data)
}

plot_map <- function(data, res, var) {
  data$r <- as.data.frame(data)[,c(paste(res))]
  ggplot(data = data) +
    geom_sf(aes(fill = r), color = NA) +
    scale_fill_viridis_c() +
    custom_theme +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.text = element_blank()) +
    labs(title = var)
}

# Create a dummy plot to generate the continuous legend
plot_dummy_legend <- function() {
  legend_plot <- ggplot(data.frame(x = 0, y = 0, z = seq(0, 1, length.out = 100)), aes(x, y, fill = z)) +
                         geom_tile(NULL) +
                          scale_fill_viridis_c(
                        guide = guide_colorbar(
                          direction = "horizontal",
                          title = NULL,
                          label.position = "bottom",
                          barwidth = 20,
                          barheight = 1)) +
                          theme_void() +
                          theme(
                            legend.position = "bottom",
                            legend.title = element_blank(),
                            legend.text = element_text(size = 10))
  
  # Extract the legend as a grob
  shared_legend <- cowplot::get_plot_component(legend_plot, 'guide-box-bottom', return_all = TRUE)
  
  legend_grob <- cowplot::ggdraw(shared_legend)
  return(legend_grob)
}

# make a bi plot map to see overlaps - not using much
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
  return(as.numeric((extractAIC.HLfit(model))[2]))
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

add_centroids <- function (data) {
  # Exclude polygons with no neighbours from your data
  clean_data <- data
  
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  return(clean_data)
}

