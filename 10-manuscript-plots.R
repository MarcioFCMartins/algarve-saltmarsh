##        Plots for the manuscript               ##
##          Márcio Martins    2024-03-25         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##


# Setup  ------------------------------------------------------------------

libs <- c(
    "ggplot2",   # Plot framework
    "ggthemr",   # Pre-made themes for ggplot2
    "ggmosaic",  # Expands ggplot2 with an area mosaic plot
    "patchwork", # Arrange several ggplot2 plots
    "tidyverse", # Data wrangling
    "scales",    # Scale transformation and "prettyfying"
    "sf",        # Methods for vector spatial data
    "here",      # Path normalization
    "workflows", # Workflow for machine learning
    "janitor",   # Data cleaning
    "tidymodels",# Wrapper for several machine learning packages
    "terra"      # Spatial raster methods
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

# Load helper functions
source(here("code/helpers-plots.R"))

ggthemr("fresh")

FIGURES_DIR <- here("analysis/manuscript/saltmarsh/figures")

# FIG 3 -  Vegetation communities composition -------------------------------------
veg_quads <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-vegetation-quadrats.gpkg")
) |>
    mutate(
        water_system = factor(water_system, c("Alvor", "Arade", "Ria Formosa", "Guadiana")),
        cluster = as.factor(cluster),
        marsh_type = factor(marsh_type, levels = c("low", "middle", "high"))
    )

## Fig 2
plot_composition(
    data    = veg_quads, # Matrices
    sp_cols = c(9:32), # Columns with species presences
    group   = marsh_type
) # X axis shows saltmarsh community

ggsave(
    here(FIGURES_DIR, "Fig3-community-composition.png"),
    height = 10,
    width = 20,
    units = "cm"
)

# FIG 4 - Estimated submersion time vs DEM --------------------------------

elev <- rast("./outputs/rasters/submersion-time/ria-formosa-updated-submersion.tif")


lidar <- rast("~/Work/200-data-and-resources/203 - Spatial rasters/203.01 - Ria Formosa bathymetry/dem_slope_aspect_flowdir.tif")

# Project LIDAR to the same CRS as elev, and only keep elevation layer
lidar <- subset(lidar, 1) |>
    project("epsg:32629", threads = T)
# Make the raster geometry of LIDAR equal to elev
lidar <- crop(lidar, ext(elev))
lidar <- resample(lidar, elev)

compare <- data.frame(
    "elev" = unname(values(elev)[,1]),
    "lidar"= unname(values(lidar)[,1])
)

compare_subset <- na.omit(compare) |>
    filter(lidar < 2.5 & lidar > -2.5)


comp_plot <- compare_subset |>
    mutate(elev = round(elev, 1),
           lidar = round(lidar, 2)
    ) |>
    group_by(elev, lidar) |>
    summarise(z = n()) |>
    ungroup() |>
    complete(elev, lidar, fill = list(z = 0))

cor_test <- cor.test(compare_subset$elev, compare_subset$lidar)
cor_report <- paste0(
    "Pearson r: ",
    round(cor_test$estimate, 2),
    "\n",
    "p < 0.001"
)

ggplot(comp_plot) +
    geom_raster(aes(x = lidar, y = elev, fill = z)) +
    scale_y_continuous(
        expand = expansion(0,0),
        breaks = seq(0, 1, 0.1),
        labels = scales::percent) + 
    scale_x_continuous(expand = expansion(0,0)) +
    coord_fixed(ratio = 2) +
    scale_fill_viridis_c() +
    annotate(
        "text",
        x = -2.4, y = 0,
        label = cor_report,
        color = "white",
        hjust = 0,
        vjust = 0,
        size = 3.5
    ) +
    labs(
        x = "Elevation (meters vs Datum Cascais Helmert 38)",
        y = "Submersion time (%)",
        fill = "Pixel count")

ggsave(
    here(FIGURES_DIR, "Fig4b - subtime correlation plot.png"),
    width = 8,
    height = 4.5,
    dpi = 300
)



# FIG 4 - Candidate model performances ------------------------------------------
candidate_mod_perf <- readRDS(
    here("outputs/models/all_model_results.RDS")
) |> 
    # Predicting clusters makes no sense with them being different per system
    filter(response == "marsh_type")

model_name_dict <- c(
    "rf"      = "Random Forest",
    "knn"     = "K-Nearest Neighbors",
    "svm_pol" = "SVM polynomial kernel",
    "svm_rad" = "SVM radial kernel",
    "glm"     = "Generalized Linear Model",
    "xgb"     = "Extreme Gradient Boosting",
    "nn"      = "Neural Network",
    "lasso"   = "Lasso Regression"
)
## Fig 3
plot_candidate_mod_performances(
    data = candidate_mod_perf,
    model_names = model_name_dict,
    max_x = 10 # How many models to show in the X-axis
) 

ggsave(
    here(FIGURES_DIR, "Fig5-model-performances.png"),
    width = 10,
    height = 15,
    units = "cm"
)

