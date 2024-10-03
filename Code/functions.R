source("libraries.R")

# edit data for analysis
edit_shp <- function (data) {
  
  # EVI Indexes on the sample
  data$r100_total <-(data$r100_C + data$r100_B + data$r100_A + data$r100_NODAT)
  data$r100_NODAT <- (data$r100_NODAT)/data$r100_total
  data$r100_SVI <- (data$r100_B)/data$r100_total
  data$r100_SUVI <- (data$r100_A)/data$r100_total
  data$r100_SGVI <- (data$r100_C)/data$r100_total
  data$r100_SEVI <- (data$r100_C - data$r100_A)/(data$r100_C + data$r100_A)
  data$r100_SBVI <- (1 + (data$r100_C - data$r100_A)/data$r100_total)/2
  
  # EVI Indexes on the sample
  data$r250_total <-(data$r250_C + data$r250_B + data$r250_A + data$r250_NODAT)
  data$r250_NA <- (data$r250_NODAT)/data$r250_total
  data$r250_SVI <- (data$r250_B)/data$r250_total
  data$r250_SUVI <- (data$r250_A)/data$r250_total
  data$r250_SGVI <- (data$r250_C)/data$r250_total
  data$r250_SEVI <- (data$r250_C - data$r250_A)/(data$r250_C + data$r250_A)
  data$r250_SBVI <- (1 + (data$r250_C - data$r250_A)/data$r250_total)/2
  
  # EVI Indexes on the sample
  data$r500_total <-(data$r500_C + data$r500_B + data$r500_A + data$r500_NODAT)
  data$r500_NA <- (data$r500_NODAT)/data$r500_total
  data$r500_SVI <- (data$r500_B)/data$r500_total
  data$r500_SUVI <- (data$r500_A)/data$r500_total
  data$r500_SGVI <- (data$r500_C)/data$r500_total
  data$r500_SEVI <- (data$r500_C - data$r500_A)/(data$r500_C + data$r500_A)
  data$r500_SBVI <- (1 + (data$r500_C - data$r500_A)/data$r500_total)/2
  
  data$EVI_class <- cut(data$EVI_mean, breaks = c(0, 0.16, 0.56, 1), labels = c("A", "B", "C"))
  
  # Calculate SVIproximal
  data <- data %>%
    mutate(
      SVIproximal = if_else(r100_total - 1 == 0, NA_real_, r100_B / (r100_total - 1)),
      SVIproximal = if_else(EVI_class == "B" & r100_total - 1 != 0, (r100_B - 1) / (r100_total - 1), SVIproximal),
      SGVIproximal = if_else(r100_total - 1 == 0, NA_real_, r100_C / (r100_total - 1)),
      SGVIproximal = if_else(EVI_class == "C" & r100_total - 1 != 0, (r100_C - 1) / (r100_total - 1), SGVIproximal),
      SUVIproximal = if_else(r100_total - 1 == 0, NA_real_, r100_A / (r100_total - 1)),
      SUVIproximal = if_else(EVI_class == "A" & r100_total - 1 != 0, (r100_A - 1) / (r100_total - 1), SUVIproximal),
      #SBVIproximal = if_else(r100_total - 1 == 0, NA_real_, (1 + ((r100_C - r100_A) / (r100_total - 1)))/2),
      #SBVIproximal = if_else(EVI_class == "C" & r100_total - 1 != 0, (1 + ((r100_C - r100_A - 1) / (r100_total - 1)))/2, SBVIproximal),
      #SBVIproximal = if_else(EVI_class == "A" & r100_total - 1 != 0, (1 + ((r100_C - r100_A + 1) / (r100_total - 1)))/2, SBVIproximal)
      SBVIproximal = if_else(r100_total - 1 == 0, NA_real_, (r100_C - r100_A) / (r100_total - 1)),
      SBVIproximal = if_else(EVI_class == "C" & r100_total - 1 != 0, (r100_C - r100_A - 1) / (r100_total - 1), SBVIproximal),
      SBVIproximal = if_else(EVI_class == "A" & r100_total - 1 != 0, (r100_C - r100_A + 1) / (r100_total - 1), SBVIproximal)
    )
  
  # Calculate SEVIproximal
  data <- data %>%
    mutate(
      SEVIproximal = if_else((r100_C + r100_A) == 0, NA_real_, (r100_C - r100_A) / (r100_C + r100_A)),
      SEVIproximal = if_else(EVI_class == "A" & (r100_C + r100_A - 1) != 0, (r100_C - r100_A + 1) / (r100_C + r100_A - 1), SEVIproximal),
      SEVIproximal = if_else(EVI_class == "C" & (r100_C + r100_A - 1) != 0, (r100_C - r100_A - 1) / (r100_C + r100_A - 1), SEVIproximal)
    )
  
  # Calculate SVIdistal
  data <- data %>%
    mutate(
      SVIdistal = if_else((r500_total - r100_total) == 0, NA_real_, (r500_B - r100_B) / (r500_total - r100_total)),
      SGVIdistal = if_else((r500_total - r100_total) == 0, NA_real_, (r500_C - r100_C) / (r500_total - r100_total)),
      SUVIdistal = if_else((r500_total - r100_total) == 0, NA_real_, (r500_A - r100_A) / (r500_total - r100_total)),
      #SBVIdistal = if_else((r500_total - r100_total) == 0, NA_real_, (1 + r500_SGVI - r500_SUVI)/2)
      SBVIdistal = if_else((r500_total - r100_total) == 0, NA_real_, r500_SGVI - r500_SUVI)
    )
  
  # Calculate SEVIdistal
  data <- data %>%
    mutate(
      SEVIdistal = if_else(((r500_C - r100_C) + (r500_A - r100_A)) == 0, NA_real_, ((r500_C - r100_C) - (r500_A - r100_A)) / ((r500_C - r100_C) + (r500_A - r100_A)))
    )
  
  clean_data <- data
  clean_data$pop <- clean_data$pop * 0.09
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  clean_data$EVI_class <- cut(clean_data$EVI_mean, breaks = c(0, 0.16, 0.56, 1), labels = c("A", "B", "C"))
  #data$SBVI <- abs(data$C - data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  
  return(clean_data)
}

z_scale <- function(x) {
  2 * ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))) - 1
}

normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

plot_map <- function(data, res, var) {
  data$r <- as.data.frame(data)[,c(paste(res))]
  ggplot(data = data) +
    geom_sf(aes(fill = r), color = NA) +
    scale_fill_viridis_c(option = "turbo", limits = c(-1, 1)) +
    custom_theme +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.text = element_blank()) +
    labs(title = "")
}

plot_tern <- function(data, var, color= "turbo") {
  color = "turbo"
  data$var <- (as.data.frame(data)[,c(paste(var))])
  ggtern(data, aes(x = SGVIproximal, y = SVIproximal, z = SUVIproximal, color = (var))) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_viridis_c(option = color, limits = c(-1,1)) +
    theme_hidetitles() +
    theme_hidelabels() +
    theme_hideticks() +
    theme(
      legend.position = "none",
    ) 
    # + labs(x = "SGVI", y = "SVI", z = "SUVI")
}

plot_tern2 <- function(data, var, color = "turbo") {
  color = "turbo"
  data$var <- (as.data.frame(data)[,c(paste(var))])
    ggtern(data, aes(x = SGVIdistal, y = SVIdistal, z = SUVIdistal, color = (var))) +
      geom_point(size = 1) +
      scale_color_viridis_c(option = color, limits = c(-1,1)) +
      theme_hidetitles() +
      theme_hidelabels() +
      theme_hideticks() +
      theme(
        legend.position = "none",
      ) 
      #+ labs(x = "SGVI", y = "SVI", z = "SUVI") 
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
