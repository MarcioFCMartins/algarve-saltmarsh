###################################################
##        Create prediction maps                 ##
##          Márcio Martins    2024-02-08         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##
###################################################


# Setup -------------------------------------------------------------------

libs <- c(
    "tidymodels", # Wrapper for several machine learning packages
    "sf", # Simple-features for vector spatial data
    "here", # Path normalization
    "tidyverse", # Data wrangling
    "terra" # Raster methods
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))


source(here("code/helpers-rasters.R"))

# Path to save the raster with the prediction variables, masked and cropped
PRED_RASTER_PATH <- here("outputs/rasters/prediction-variables-composite.tif")


# Read final model --------------------------------------------------------
pred_mod <- readRDS(here("outputs/models/final-model.RDS"))


# Move all outputs to data folder -----------------------------------------
# To ensure all data files are up to date, copy them to the clean data files
study_systems <- c(
    "alvor",
    "arade",
    "ria-formosa",
    "guadiana"
)

for (system_i in study_systems) {
    # List all output rasters
    files <- list.files(
        here("outputs/rasters"),
        recursive = TRUE,
        pattern = "tif",
        full.names = TRUE
    )
    
    # Filter for current system
    files <- files[str_detect(files, system_i)]
    
    # Select rasters actually needed for the analysis
    target_files <- c(
        "updated-composite-overlap-count",
        "updated-composite-q20",
        "updated-elev",
        "updated-marsh-mask-50",
        "updated-ndvi-q20",
        "updated-ndwi-high-tide",
        "updated-ndwi-low-tide",
        "updated-submersion"
    )
    
    files <- files[str_detect(files, paste(target_files, collapse = "|"))]
    
    target_path <- here("data/clean/mapping", system_i)
    
    file.copy(files, to = target_path, overwrite = TRUE, copy.date = TRUE)
}


# # Prepare prediction rasters ----------------------------------------------

message("\n\nCreating raster with selected predictors, masked to only contain marsh areas")
message("Path:", PRED_RASTER_PATH, "\n")

# Extract predictors from the preprocessor
predictors <- extract_preprocessor(pred_mod)$var_info$variable
# Remove the response variable
predictors <- predictors[-length(predictors)]

study_systems <- c(
    "alvor",
    "arade",
    "ria-formosa",
    "guadiana"
)

map_rasters <- list()

for (system_i in study_systems) {
    # Read rasters of interest for current water system
    files <- list.files(
        here("data/clean/mapping", system_i),
        pattern = "tif",
        full.names = TRUE
    )

    rasters <- lapply(
        files,
        terra::rast
    )

    ext <- get_full_extent(rasters)

    rasters <- lapply(
        rasters,
        \(x) extend(x, ext) |>
            crop(ext)
    )
    
    # Only keep layers which are used as predictors
    rasters <- subset(rast(rasters), predictors)
    
    # Exclude areas based on automatic masks
    system_mask <- rast(here("data/clean/mapping", system_i, paste0(system_i, "-updated-marsh-mask-50.tif")))
    rasters <- mask(rasters, system_mask)
    
    # Exclude areas based on manually curated masks
    mask_manual <- st_read(
        here("data/clean/mapping/all-systems/manual-clean-up-masks.gpkg"),
        quiet = TRUE
    ) |>
        filter(system == system_i) |>
        st_transform(crs = crs(rasters)) 
    # Masking of large files was causing OOM crash...using a temp file seems to help
    tmp <- tempfile(fileext = ".tif")
    rasters2 <- terra::mask(
        x = rasters, 
        mask = mask_manual, 
        inverse = TRUE,
        filename = tmp)

    map_rasters[[system_i]] <- rasters2
}

# Merge to a single raster and save to disk
map_rasters <- sprc(map_rasters)
map_rasters <- terra::merge(
    map_rasters,
    filename = here("outputs/rasters/prediction-variables-composite.tif"),
    overwrite = TRUE)

message("Predictor raster saved.", "\n\n")


# Now apply prediction model ----------------------------------------------
terraOptions(memfrac = 0.5)
pred_raster <- rast(PRED_RASTER_PATH)

message("Creating predicted map.", "\n\n")

response_map <- predict(
    pred_raster,
    pred_mod,
    fun = function(model, data) {
        # If all cells in this chunk are NA, return all preds as
        if (all(is.na(data))) {
            preds <- rep(NA, nrow(data))
            return(preds)
        }
        # If there are cells to actually predict:
        # Remove NAs so that the workflows pre-processing and prediction can be applied
        preds <- predict(model, na.omit(data))
        preds <- preds$.pred_class
        # Convert the response variable to numeric - seems to be a requirement of the geotif format
        preds <- case_when(
            preds == "low" ~ 1,
            preds == "middle" ~ 2,
            preds == "high" ~ 3
        )

        # The returned vector must have the same length as the original data
        # Initialize empty vector with same size as the input data
        preds_full <- as.double(rep(NA, nrow(data)))
        # Which cells from input data we want to predict for
        cells_with_data <- apply(data, 1, \(x) all(!is.na(x)))

        preds_full <- replace( # Replace the empty prediction vec with actual predictions
            preds_full,
            cells_with_data, # But only replace the cells with data in the original raster
            preds # Replace them with the predicitons
        )
        
        return(preds_full)
    }
)


writeRaster(response_map, here("outputs/rasters/predicted-map.tif"), overwrite = TRUE)

