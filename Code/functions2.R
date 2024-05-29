source("libraries.R")

# Add centroids to data
edit_shp <- function (data) {
  # EVI Indexes on the sample
  data$EVI_nodata <- (data$EVI_NODATA)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SVI <- (data$B)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SUVI <- (data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SGVI <- (data$C)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SEVI <- (data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SBVI <- abs(data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  
  clean_data <- data[!is.na(data$EVI_mean),]
  clean_data <- clean_data[clean_data$EVI_nodata < 0.3,]
  
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  clean_data$BH[is.na(clean_data$BH)] <- 0
  clean_data$VH[is.na(clean_data)] <- 0
  
  clean_data$Cat <- cut(clean_data$EVI_mean, breaks = c(0, 0.16, 0.56, 1), 
                   labels = c("SU", "SH", "SG"), include.lowest = TRUE)
  
  return(clean_data)
}

# Add centroids to data
edit_shp2 <- function (data) {
  centroids <- st_centroid(data)
  # Add centroid coordinates as new columns
  data$centx <- st_coordinates(centroids)[, 1]
  data$centy <- st_coordinates(centroids)[, 2]
  
  return(data)
}


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