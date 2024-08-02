source("libraries.R")

# edit data for analysis
edit_shp2 <- function (data) {
  
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
      SBVIproximal = if_else(r100_total - 1 == 0, NA_real_, (1 + (r100_C - r100_A) / (r100_total - 1))/2),
      SBVIproximal = if_else(EVI_class == "C" & r100_total - 1 != 0, (1 + (r100_C - r100_A - 1) / (r100_total - 1))/2, SBVIproximal),
      SBVIproximal = if_else(EVI_class == "A" & r100_total - 1 != 0, (1 + (r100_C - r100_A + 1) / (r100_total - 1))/2, SBVIproximal)
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
      SBVIdistal = if_else((r500_total - r100_total) == 0, NA_real_, (1 + (r500_C - r100_C) - (r500_A - r100_A)/ (r500_total - r100_total))/2)
    )
  
  # Calculate SEVIdistal
  data <- data %>%
    mutate(
      SEVIdistal = if_else(((r500_C - r100_C) + (r500_A - r100_A)) == 0, NA_real_, ((r500_C - r100_C) - (r500_A - r100_A)) / ((r500_C - r100_C) + (r500_A - r100_A)))
    )
  
  #data$SEVIproximal <- (1+data$SEVIproximal)/2
  #data$SEVIproximal <- (1+data$SEVIproximal)/2
  
  clean_data <- data
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

calculate_surrounding_old <- function(nb_list, data, var_name, l) {
  data <- as.data.frame(data)
  index <- sapply(1:length(nb_list), function(i) {
    neighbors <- nb_list[[i]]
    sum_var <- sum(data[neighbors, var_name], na.rm = TRUE)
    sum_all <- sum(data[neighbors, "A"], na.rm = TRUE) +
      sum(data[neighbors, "B"], na.rm = TRUE) +
      sum(data[neighbors, "C"], na.rm = TRUE) +
      sum(data[neighbors, "EVI_nodata"])
    if (sum_all == 0) {
      return(NA)
    } else {
      return(sum_var / sum_all)
    }
  })
  return(index)
}

calculate_surrounding <- function(nb_list, data, var_name) {
  data <- as.data.frame(data)
  index <- sapply(1:length(nb_list), function(i) {
    neighbors <- nb_list[[i]]
    total_neighbors <- length(neighbors)
    
    if (total_neighbors == 0) {
      return(NA)
    } else {
      count_var <- sum(data[neighbors, "EVI_class"] == var_name, na.rm = TRUE)
      proportion_var <- count_var / l*l
      return(proportion_var)
    }
  })
  return(index)
}

calculate_surrounding_dif_old <- function(nb_list, data) {
  data <- as.data.frame(data)
  index <- sapply(1:length(nb_list), function(i) {
    neighbors <- nb_list[[i]]
    sum_C <- sum(data[neighbors, "C"], na.rm = TRUE)
    sum_A <- sum(data[neighbors, "A"], na.rm = TRUE)
    if (sum_A + sum_C == 0) {
      return(NA)
    } else {
      return((sum_C-sum_A)/(sum_C+sum_A))
    }
  })
  return(index)
}

calculate_surrounding_dif <- function(nb_list, data, var_name) {
  data <- as.data.frame(data)
  index <- sapply(1:length(nb_list), function(i) {
    neighbors <- nb_list[[i]]
    total_neighbors <- length(neighbors)
    
    if (total_neighbors == 0) {
      return(NA)
    } else {
      count_A <- sum(data[neighbors, "EVI_class"] == "A", na.rm = TRUE)
      count_C <- sum(data[neighbors, "EVI_class"] == "C", na.rm = TRUE)
      return((count_C-count_A)/(count_C+count_A))
    }
  })
  return(index)
}

# edit data for analysis
edit_shp <- function (data) {
  data$prop_veg <- (data$PV_1 + data$PV_2 + data$PV_3)/(data$PV_1 + data$PV_2 + data$PV_3 + data$PV_NODATA)
  
  # EVI Indexes on the sample
  data$EVI_nodata <- (data$EVI_NODATA)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SVI <- (data$B)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SUVI <- (data$A)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SGVI <- (data$C)/(data$C + data$B + data$A + data$EVI_NODATA)
  data$SEVI <- (data$C - data$A)/(data$C + data$A)
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
  
  #clean_data <- clean_data[,c("id", "CR", "AB", "pop", "CRn", "ABn", "popn", "SVI", "SGVI", "SUVI", "SEVI", "prop_veg", "centx", "centy")]
  
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

plot_tern <- function(data, var, color= "turbo") {
  color = "turbo"
  data$var <- (as.data.frame(data)[,c(paste(var))])
  ggtern(data, aes(x = SGVIproximal, y = SVIproximal, z = SUVIproximal, color = (var))) +
    geom_point(size = 1) +
    scale_color_viridis_c(option = color, limits = c(-1,1)) +
    theme_hidelabels() +
    theme_hideticks() +
    theme(
      legend.position = "none",
    ) +
    labs(x = "SGVI", y = "SVI", z = "SUVI")
}

plot_tern2 <- function(data, var, color = "turbo") {
  color = "turbo"
  data$var <- (as.data.frame(data)[,c(paste(var))])
    ggtern(data, aes(x = SGVIdistal, y = SVIdistal, z = SUVIdistal, color = (var))) +
      geom_point(size = 1) +
      scale_color_viridis_c(option = color, limits = c(-1,1)) +
      theme_hidelabels() +
      theme_hideticks() +
      theme(
        legend.position = "none",
      ) +
      labs(x = "SGVI", y = "SVI", z = "SUVI")
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

add_centroids <- function (data) {
  # Exclude polygons with no neighbours from your data
  clean_data <- data
  
  centroids <- st_centroid(clean_data)
  # Add centroid coordinates as new columns
  clean_data$centx <- st_coordinates(centroids)[, 1]
  clean_data$centy <- st_coordinates(centroids)[, 2]
  
  return(clean_data)
}

