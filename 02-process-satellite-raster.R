##  Average satellite images over a year         ##
##        Márcio Martins    2023-06-19           ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##


# Setup -------------------------------------------------------------------
libs <- c(
    "PTtidaltools",
    "lubridate",
    "sf",
    "terra",
    "here",
    "stringr"
)

# Install missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

# Limit terra memory usage to prevent issues
terraOptions(memfrac = 0.4)

# Change me depending on desired quantile for composites
composite_quantile <- 0.2

region_names <- c(
    "alvor-updated",
    "arade-updated",
    "ria-formosa-updated",
    "guadiana-updated"
)


# Available images for 2022 -----------------------------------------------

# Make a list of available rasters for every water system
images <- list()
for (region in region_names) {
    available_queries <- list.files(
        "/mnt/10274c4b-4f18-41e0-a518-ff86b71a055f/planet_labs_imagery",
        pattern = paste0(region, ".+"),
        full.names = TRUE
    )

    # Only keep queries in the lower tidal interval
    # This regex will match _, digit, optional dot and optional digit - digit, optional dot and optional digit
    valid_queries <- available_queries[
        str_detect(available_queries, "_0\\.[3|6]-\\d\\.?\\d?")
    ]

    rasters <- list.files(
        path = paste0(valid_queries, "/", "PSScene"),
        pattern = "_AnalyticMS_SR_8b.+tif$",
        full.names = TRUE
    )

    images[[region]] <- rasters
}

rm(available_queries, valid_queries, rasters)



# Create composites of 2022 -----------------------------------------------

for (region in region_names) {
    message(paste0("Creating composites for ", region))

    files <- images[[region]]

    # Read system mask to delimit ROI
    roi <- list.files(
        here("data/clean/mapping/all-systems/system-masks/"),
        pattern = region,
        full.names = TRUE
    )
    # And re-project to the same projection as satellite images
    roi <- st_read(roi) |> st_transform(32629)

    # Calculate and save composite raster
    composite <- create_raster_composite(files, prob = composite_quantile, roi_mask = roi)
    file_name <- paste0(
        "./outputs/rasters/system-composites/",
        region,
        "-composite-q",
        composite_quantile * 100,
        ".tif"
    )
    writeRaster(composite, file_name, names = names(composite))
    rm(composite)
    gc()

    # Calculate and number of overlapping images per pixel
    overlap_counts <- calculate_overlap_per_pixel(files)
    file_name <- paste0(
        "./outputs/rasters/system-composites/",
        region,
        "-composite-overlap-count.tif"
    )
    writeRaster(overlap_counts, file_name, overwrite = TRUE)
    rm(overlap_counts)
    gc()
}


# NDVI based on 2022 composites -------------------------------------------

for (region in region_names) {
    composite <- paste0(
        "./outputs/rasters/system-composites/",
        region,
        "-composite-q",
        composite_quantile * 100,
        ".tif"
    ) |>
        rast()

    ndvi <- (composite$nir - composite$red) / (composite$nir + composite$red)

    file_name <- paste0(
        "./outputs/rasters/ndvi/",
        region,
        "-ndvi-q",
        composite_quantile * 100,
        ".tif"
    )

    writeRaster(ndvi, file_name, overwrite = TRUE, names = "ndvi")
}
