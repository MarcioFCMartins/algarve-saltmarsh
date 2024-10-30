## Determine elevation of intertidal sediment    ##
## using high frequency satellite imagery        ##
##          Márcio Martins    2023/06/15         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##



# Setup -------------------------------------------------------------------
libs <- c(
    "here",
    "PTtidaltools",
    "sf"
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

source(here("code/helpers-rasters.R"))
here::i_am("code/02-elevation-estimation.R")

# Available images --------------------------------------------------------

region_names <- c(
    "alvor-updated", "arade-updated",
    "ria-formosa-updated", "guadiana-updated"
)

# Port 18 is in Lagos, 19 in Faro-Olhão, 21 in Vila Real de Santo António
region_port <- c(
    "alvor-updated" = 18, "arade-updated" = 18,
    "ria-formosa-updated" = 19, "guadiana-updated" = 21
)

# Make a list of available rasters for every water system
images <- list()
for (region in region_names) {
    available_queries <- list.files(
        "/mnt/10274c4b-4f18-41e0-a518-ff86b71a055f/planet_labs_imagery",
        pattern = paste0(region, ".+"),
        full.names = TRUE
    )

    # Only keep queries that have tidal ranges in the name
    # This regex will match _, digit, optional dot and optional digit - digit, optional dot and optional digit
    valid_queries <- available_queries[
        str_detect(available_queries, "_\\d\\.?\\d?-\\d\\.?\\d?")
    ]

    rasters <- list.files(
        path = paste0(valid_queries, "/", "PSScene"),
        pattern = "_AnalyticMS_SR_8b.+tif$",
        full.names = TRUE
    )

    images[[region]] <- rasters
}

# This specific raster has a huge cloud shadow in it and is problematic. Remove it
images[[4]] <- images[[4]][!images[[4]] == "/mnt/10274c4b-4f18-41e0-a518-ff86b71a055f/planet_labs_imagery/guadiana-updated_20221001_20221231_1-1.3/PSScene/20221210_110352_14_2414_3B_AnalyticMS_SR_8b.tif"]



rm(available_queries, valid_queries, rasters)


# Estimate site elevation -------------------------------------------------

for (i in 1:length(images)) {
    raster_vec <- images[[i]]
    system_name <- names(images)[i]
    port_id <- region_port[system_name]

    system_mask <- st_read(
        here(
            paste0("data/clean/mapping/all-systems/system-masks/", system_name, ".geojson")
        )
    ) |>
        st_transform(32629)

    est_elev <- classify_tidal_regime(raster_vec, port_id, roi_mask = system_mask)

    file_names <- paste0(
        here("./outputs/rasters/elevation-estimates/"),
        system_name, "-", names(est_elev), ".tif"
    )


    for (j in 1:length(est_elev)) {
        writeRaster(est_elev[[j]], file_names[j], overwrite = TRUE)
    }
}

# Convert estimated elevations to submersion times ------------------------
tide_heights <- data.frame()
for (region in region_names) {
    system_name <- region
    port_id <- region_port[system_name]

    # Retrieve intertidal heights for 2022, at 15 min interval
    tides <- interpolate_tides(
        date_times = seq(ISOdate(2022, 1, 1), ISOdate(2022, 12, 31), by = 900),
        port_id    = port_id
    )

    # Read cell height estimates
    est_elev <- rast(
        here(paste0("outputs/rasters/elevation-estimates/", system_name, "-elev.tif"))
    )

    # Read estimated supratidal zones
    supratidal <- rast(
        here(paste0("outputs/rasters/elevation-estimates/", system_name, "-supratidal.tif"))
    )

    est_elev <- c(est_elev, supratidal)

    # Create new raster with submersion time for 2022 (as fraction)
    tide_nr <- length(tides)
    submersion_fraction <- app(
        est_elev,
        function(x) {
            if (x[2] == 1) { # If cell is supratidal return 0 submersion time
                return(0)
            } else { # Else, calculate fraction tides under cell height
                return(sum(x[1] < tides) / tide_nr)
            }
        }
    )

    writeRaster(
        submersion_fraction,
        here(paste0("outputs/rasters/submersion-time/", system_name, "-submersion.tif")),
        overwrite = TRUE,
        names = "subtime"
    )

    # Export tides used
    tides <- data.frame(
        "water_system" = gsub("-updated", "", region),
        "date_time"    = seq(ISOdate(2022, 1, 1), ISOdate(2022, 12, 31), by = 900),
        "height"       = tides
    )

    tide_heights <- rbind(tide_heights, tides)
}

write.csv2(
    tide_heights,
    here("outputs/tables/tides_2022.csv")
)

# Create tables with tidal information for each image ---------------------
tide_heights <- data.frame()

for (i in 1:length(images)) {
    raster_vec <- images[[i]]
    system_name <- names(images)[i]
    port_id <- region_port[system_name]

    system_mask <- st_read(
        here(
            paste0("data/clean/mapping/all-systems/system-masks/", system_name, ".geojson")
        )
    ) |>
        st_transform(32629)

    tides <- c()
    for (img in raster_vec) {
        tides <- c(
            tides,
            get_image_tide(img, port_id)
        )
    }

    tides <- data.frame(
        "water_system" = gsub("-updated", "", system_name),
        "height"       = tides
    )

    tide_heights <- rbind(tide_heights, tides)
}

write.csv2(
    tide_heights,
    here("outputs/tables/tide-heights-acquired-images.csv")
)


# Calculate NDWI at low and high tides ------------------------------------

for (region in region_names) {
    available_queries <- list.files(
        "/mnt/10274c4b-4f18-41e0-a518-ff86b71a055f/planet_labs_imagery",
        pattern = paste0(region, ".+"),
        full.names = TRUE
    )

    # Only keep queries that have tidal ranges in the name
    # This regex will match _, digit, optional dot and optional digit - digit, optional dot and optional digit
    valid_queries <- available_queries[
        str_detect(available_queries, "_\\d\\.?\\d?-\\d\\.?\\d?")
    ]

    query_range <- str_extract(valid_queries, "_\\d\\.?\\d?-\\d\\.?\\d?") |>
        str_replace("_", "")

    query_range <- str_split(query_range, "-")
    query_min <- sapply(query_range, \(x) as.numeric(x[[1]]))
    query_max <- sapply(query_range, \(x) as.numeric(x[[2]]))

    low_tide_queries <- valid_queries[query_max <= 1.0]
    high_tide_queries <- valid_queries[query_min >= 2.9]

    low_tide_rasters <- list.files(
        path = paste0(low_tide_queries, "/", "PSScene"),
        pattern = "_AnalyticMS_SR_8b.+tif$",
        full.names = TRUE
    )

    high_tide_rasters <- list.files(
        path = paste0(high_tide_queries, "/", "PSScene"),
        pattern = "_AnalyticMS_SR_8b.+tif$",
        full.names = TRUE
    )


    low_tide_rasters <- lapply(low_tide_rasters, terra::rast)
    high_tide_rasters <- lapply(high_tide_rasters, terra::rast)

    low_tide_ndwi <- lapply(
        low_tide_rasters,
        \(x) (x$green - x$nir) / (x$green + x$nir)
    )
    high_tide_ndwi <- lapply(
        high_tide_rasters,
        \(x) (x$green - x$nir) / (x$green + x$nir)
    )

    new_ext <- st_read(
        here(
            paste0("data/clean/mapping/all-systems/system-masks/", region, ".geojson")
        )
    ) |>
        st_transform(32629) |>
        ext()


    low_tide_ndwi <- lapply(low_tide_ndwi, \(x) terra::crop(x, new_ext, extend = TRUE))
    high_tide_ndwi <- lapply(high_tide_ndwi, \(x) terra::crop(x, new_ext, extend = TRUE))

    low_tide_composite <- rast(low_tide_ndwi) |>
        terra::median(na.rm = TRUE)

    high_tide_composite <- rast(high_tide_ndwi) |>
        terra::median(na.rm = TRUE)

    names(low_tide_composite) <- "ndwi_low"
    names(high_tide_composite) <- "ndwi_high"

    writeRaster(
        low_tide_composite,
        here(paste0("outputs/rasters/ndwi/", region, "-low-tide.tiff"))
    )

    writeRaster(
        high_tide_composite,
        here(paste0("outputs/rasters/ndwi/", region, "-high-tide.tiff"))
    )
}
