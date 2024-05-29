source("functions.R")
generate_landscape <- function(prop_veg) {
  pv <- prop_veg*100
  urb <- 100 - pv
  cls_a <- flsgen_create_class_targets(
    "A", 
    NP = c(5, 20),    # Number of patches for vegetation
    AREA = c(100, 500), # Area range for vegetation patches
    PLAND = c(pv-5, pv+5)  # Landscape percentage covered by vegetation
  )
  cls_b <- flsgen_create_class_targets(
    "B", 
    NP = c(5, 15),     # Number of patches for urban areas
    AREA = c(100, 500), # Area range for urban patches
    PLAND = c(urb-10, urb-5)  # Landscape percentage covered by urban areas
  )
  
  # Create landscape targets with 1000x1000 pixels grid
  ls_targets <- flsgen_create_landscape_targets(
    nb_rows = 100, nb_cols = 100, classes = list(cls_a, cls_b)
  )
  
  # Generate landscape structure and the landscape itself
  structure <- flsgen_structure(ls_targets)
  landscape <- flsgen_generate(structure_str = structure, verbose = F)
  
  return(landscape)
}

# Generate the landscape
#limite prop_veg = (0.05, 0.95)
landscape <- generate_landscape(0.70)

landscape <- flsgen_terrain(100, 100, roughness = 0.2)

# Plot the generated landscape
plot(landscape)

structure <- flsgen_extract_structure_from_raster(landscape, c(0,1))



