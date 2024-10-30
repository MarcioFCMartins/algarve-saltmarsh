###################################################
## Extract training variables for training data  ##
##          Márcio Martins    2023/07/28         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##
###################################################


# Setup -------------------------------------------------------------------
libs <- c(
    "sf",
    "terra",
    "here",
    "dplyr",
    "tidyr",
    "vegan",
    "here",
    "ggplot2",
    "forcats",
    "scales",
    "stringr",
    "purrr",
    "readxl"
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

here::i_am("code/05-process-final-dataset.R")
source(here("code/helpers.R"))



# Load data ---------------------------------------------------------------

quadrats <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-vegetation-quadrats.gpkg")
) |>
    # Only keep information relevant for the model training stage
    select(water_system, transect, quad_id, cluster, marsh_type)

## Load system rasters and extract points one by one ----------------------
study_systems <- tolower(
    gsub(" ", "-", unique(quadrats$water_system), fixed = TRUE)
)
study_systems <- paste0(study_systems, "-updated")

points_data_full <- data.frame()

for (system in study_systems) {
    # Read rasters of interest for current water system
    files <- list.files(
        here("outputs/rasters"),
        pattern = system,
        recursive = TRUE,
        full.names = TRUE
    )
    files <- files[grepl("q20|submersion|texture|ndvi|ndwi", files)]
    files <- files[!grepl("clusters", files)]

    rasters <- lapply(
        files,
        terra::rast
    )

    # Only keep training points in current system and convert to CRS of rasters
    system_train_points <- quadrats |>
        filter(tolower(gsub(" ", "-", water_system, fixed = TRUE)) == gsub("-updated", "", system)) |>
        st_transform(crs(rasters[[1]]))

    points_data <- lapply(
        rasters,
        \(x) terra::extract(x, vect(system_train_points))
    )

    points_data <- do.call(cbind, points_data)
    points_data <- points_data[!colnames(points_data) == "ID"]

    points_data <- cbind(system_train_points, points_data)


    points_data_full <- rbind(points_data_full, points_data)
}


# Split data into test and train sets -------------------------------------
# All the new transects are part of the test data, also added some of the old
# transects (needed in the Ria, in the other systems to have more representative sets)
points_data_full$train <- !points_data_full$transect %in% c(
    "alvor_B1", "alvor_B1_extra", "alvor_C1", "?", "arade_A1", "arade_B1",
    "arade_C1", "arade_C2", "Transect 1", "transect 2", "Transect 3", "sal_alvor_WB1_s1_t1",
    "sal_arade_WB1_s1_t1", "sal_guad_WB1_s1_t3", "sal_ria_WB2_s1_t1", "sal_ria_WB4_s1_t1",
    "sal_ria_WB3_s2_t1", "sal_ria_WB1_s1_t2"
)


points_data_full <- relocate(points_data_full, train, .after = marsh_type)


st_write(
    points_data_full,
    here("./data/clean/mapping/all-systems/saltmarsh-training-data.gpkg"),
    delete_dsn = TRUE
)
