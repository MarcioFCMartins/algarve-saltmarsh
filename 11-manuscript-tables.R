##        Tables for the manuscript               ##
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

# Base directory to store tables
TABLES_DIR <- here("analysis/manuscript/saltmarsh/tables")

# TABLE 1 - Sampled data summary - transect and quadrat number ------------

veg_quads <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-vegetation-quadrats.gpkg")
) |>
    mutate(
        cluster = as.factor(cluster),
        marsh_type = factor(marsh_type, levels = c("low", "middle", "high"))
    )

sampling_summary <- veg_quads |>
    st_drop_geometry() |>
    group_by(water_system, transect_set) |>
    summarise(n_transect = length(unique(transect)), n_quads = n())

write.csv(
    sampling_summary,
    here(TABLES_DIR, "Tab1 - Sampling summary.csv"),
    row.names = FALSE
)

2# TABLE 2 - Final model performance  -----------------------------------------

final_mod <- readRDS(
    here("outputs/models/final-model.RDS")
)

test_data <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-training-data.gpkg"),
    quiet = TRUE
) |>
    filter(!train) |>
    st_drop_geometry()


test_data <- cbind(
    test_data,
    predict(final_mod, new_data = test_data)
)|>
    mutate(
        across(
            c(marsh_type, .pred_class),
            \(x) factor(x, levels = c("low", "middle", "high"))
        )
    )

# Performance indices for use in text
test_data |>
    group_by(water_system) |>
    nest() |>
    mutate(acc = map_dbl(data, \(x) accuracy(x, marsh_type, .pred_class)$.estimate))

accuracy(test_data, marsh_type, .pred_class)

# Kappa
kap(test_data, marsh_type, .pred_class)

# Confusion matrix - percentages
table_conf_mat <- janitor::tabyl(test_data,.pred_class, marsh_type) |>
    adorn_totals(c("col", "row")) |>
    adorn_title(
        row_name = "Predicted",
        col_name = "Truth",
        placement = "top"
    )

write.csv(
    table_conf_mat,
    here(TABLES_DIR, "Tab2 - Confusion Matrix.csv"),
    row.names = FALSE
)

# TABLE 3 - Saltmarsh community composition -----------------------------------------
study_systems <- c(
    "alvor",
    "arade",
    "ria-formosa",
    "guadiana"
)

roi_masks <- list.files(
    here("outputs/system-masks/"),
    full.names = TRUE
)

marsh_map <- rast(
    here("outputs/rasters/predicted-map.tif")
)

marsh_areas <- list()
for(system_i in study_systems){
    system_roi <- roi_masks[str_detect(roi_masks, system_i)] |>
        st_read() |>
        st_transform(crs(marsh_map))
    
    tmp <- tempfile(fileext = "_.tif")
    system_marsh <- crop(
        marsh_map,
        system_roi
    )
    
    # Divide the rasters into tiles to overcome memory issues
    raster_tiles <- makeTiles(
        system_marsh,
        rast(ncols = 2, nrows = 2, extent = ext(system_marsh)),
        filename = tmp)
    
    rm(system_marsh, system_roi)
    gc()
    
    # Process tiles one by one
    system_areas <- list()
    for(tile_i in raster_tiles){
        tile <- rast(tile_i)
        system_areas[[tile_i]] <- expanse(
            tile,
            unit = "ha",
            byValue = TRUE)
    }
    
    marsh_areas[[system_i]] <- system_areas
}

marsh_areas <- map(marsh_areas, \(x) do.call(rbind, x)) |> 
    bind_rows(.id = "system") |>
    group_by(system, value) |>
    summarise(area = sum(area, na.rm = T)) |>
    na.omit() |>
    mutate(
        value = case_when(
            value == 1 ~ "low",
            value == 2 ~ "middle",
            value == 3 ~ "high"))


table_marsh_areas <- marsh_areas |>
    mutate(factor(value, levels = c("low", "middle", "high"))) |>
    pivot_wider(
        names_from = system,
        values_from = area)  |>
    arrange(value) |>
    select(community = value, alvor, arade, `ria-formosa`, guadiana) |>
    adorn_totals("row") |>
    adorn_percentages("col") |>
    adorn_pct_formatting(digits = 1, rounding = "half up")  |>
    adorn_ns(position = "front", format_func = \(x) format(x, digits = 1))

# Round percentages to zero digits manually - use the following method:
#1. Rounding all values down to the nearest integer value;
#2. Determining the difference between the sum of the rounded values and total value;
#3. Distributing the difference between the rounded values in decreasing order of their decimal parts.
write.csv(
    table_marsh_areas,
    here(TABLES_DIR, "Tab3 - Saltmarsh area per system.csv"),
    row.names = FALSE
)


# TABLE ?? - Saltmarsh legally relevant areas ----------------------------------------

legal_areas <- list.files(
    here("data/clean/mapping/all-systems/legally-relevant-areas/"),
    pattern = "\\.gpkg",
    full.names = T
)

legal_areas <- data.frame(
    "agreement" = str_extract(basename(legal_areas), "^[a-z]+"),
    "file"      = legal_areas
)


marsh_map <- rast(
    here("outputs/rasters/predicted-map.tif")
)


agreement_areas <- list()
for(i in 3:nrow(legal_areas)){
    agreement          <- legal_areas$agreement[i]
    agreement_polygons <- st_read(legal_areas$file[i]) |>
        st_transform(32629)
    
    # Files use different variable names to identify sites, create a standard one
    names(agreement_polygons)[names(agreement_polygons) %in% c("officialna", "nome_ap", "site_name")] <- "site"
    
    marsh_map <- crop(
        marsh_map,
        agreement_polygons
    )
    
    # Go over systems in the legal agreement one by one
    system_areas <- list()
    for(system_i in agreement_polygons$site) {
        
        # Extract the polygon for the site we are analyzing
        system_polygon <- agreement_polygons[agreement_polygons$site == system_i, ]
        # Then crop the marsh mask and maskout all marsh not in the legal agreement area
        system_zones <- crop(marsh_map, system_polygon) |>
            mask(system_polygon)
        
        # Divide the rasters into tiles to overcome memory issues
        tmp <- tempfile(fileext = "_.tif")
        
        # Split map into tiles to avoid OOM issues
        raster_tiles <- makeTiles(
            system_zones,
            rast(ncols = 2, nrows = 2, extent = ext(system_zones)),
            filename = tmp)
        
        tile_areas <- list()
        for(tile_i in raster_tiles) {
            tile <- rast(tile_i)
            tile_areas[[tile_i]] <- expanse(
                tile,
                unit = "ha",
                byValue = TRUE)
        }
        
        system_areas[[system_i]] <- tile_areas
    }
    
    system_areas <- lapply(
        system_areas,
        function(x) {
            y <- do.call(rbind, x) 
            rownames(y) <- NULL
            return(y)
        }
    )
    
    system_areas <- bind_rows(system_areas, .id = "system")    |>
        group_by(system, value) |>
        summarise(area = sum(area))
    
    agreement_areas[[agreement]] <- system_areas
}

final_agreement_areas <- bind_rows(agreement_areas, .id = "agreement")

write.csv(
    final_agreement_areas, 
    here("analysis/tables/Tab3 - Legal relevant areas overlap - raw data.csv"))


