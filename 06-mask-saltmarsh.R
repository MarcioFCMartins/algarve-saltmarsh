####################################################
## Create saltmarsh masks used for classification ##
##          Márcio Martins    2023/10/11          ##
##  email     marciomartinsred@gmail.com          ##
## Github     MarcioFCMartins                     ##
####################################################


# Setup -------------------------------------------------------------------

libs <- c(
    "terra",
    "sf",
    "here",
    "readxl",
    "stringr",
    "dplyr",
    "rgrass"
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

# System names - used to find rasters for each system
study_systems <- c(
    "alvor-updated",
    "arade-updated",
    "ria-formosa-updated",
    "guadiana-updated"
)

# The final masks will be created based on the variables:
# blue,green,red,nir,ndwi_low,ndvi,ndwi_high
selected_cluster_method <- "k25-bgrnnnnxy"

# Submersion limit to be classified as water. Only areas with submersion time
# Below this will be considered as possible saltmarsh site
max_sub_time <- 0.9

here("code/06-saltmarsh-masks.R")

# Water masks -------------------------------------------------------------
# Load water masks which will be used in addition to the
# cluster based masks
water_masks <- list()

for (system in study_systems) {
    files <- list.files(
        here("data/clean/mapping"),
        pattern = paste0(system, "-submersion.tif"),
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(files) > 1) {
        stop(
            paste(
                "More than one submersion file found for",
                system
            )
        )
    }

    sub_time <- rast(files)

    water_masks[[system]] <- sub_time > max_sub_time

}


# Cluster masks -----------------------------------------------------------

# Selection of which clusters to keep
clusters_exclusion <- read_xlsx(
    here("data/raw/mapping/all-systems/marsh-mask-cluster-selection.xlsx"),
    sheet = "to-exclude"
) |>
    filter(cluster_variables == selected_cluster_method) |>
    mutate(
        water_system = str_replace(water_system, " ", "-") |>
            tolower() |>
            paste0("-updated")
    )

cluster_masks <- list()
for (system in study_systems) {
    files <- list.files(
        here("data/clean/mapping"),
        pattern = paste0(system, "-clusters.tif"),
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(files) > 1) {
        stop(
            paste(
                "More than one cluster file found for",
                system
            )
        )
    }

    clusters <- rast(files)

    clusters_exclusion_system <- clusters_exclusion |>
        filter(water_system == system) |>
        pull(marsh)

    # Exclude clusters determined as not marsh
    mask_system <- clusters %in% clusters_exclusion_system

    # Exclude area masked out as not part of system
    mask_system <- mask(
        mask_system,
        clusters,
        maskvalue = NA,
        updatevalue = FALSE
    )

    cluster_masks[[system]] <- mask_system
}



# Manual masks ------------------------------------------------------------
manual_masks <- st_read(
    here("data/clean/mapping/all-systems/manual-clean-up-masks.gpkg")
) |>
    st_transform(32629)


# Combine masks -----------------------------------------------------------
# First, remove the subtidal areas
combined_masks <- list()
for (system in study_systems) {
    final_mask <- mask(
        cluster_masks[[system]], # From the cluster mask
        water_masks[[system]], # remove the water mask
        maskvalue = TRUE, # where water mask = true
        updatevalue = FALSE # set as FALSE
    )

    final_mask <- mask(
        final_mask, # From the previous mask
        manual_masks, # Remove the areas manually masked out
        inverse = TRUE, # Polygons should be masked OUT
        updatevalue = FALSE # And cells set as false
    )


    combined_masks[[system]] <- final_mask
}

# Now, contract and expand the mask to deal with isolated pixels
# which are usually errors
for (system in study_systems) {
    # We need topographical algorithms to expand and grow cells from GRASS GIS
    # You need GRASS GIS installed and added to system path
    # run `sudo apt instal grass-gis` on debian based systems
    # see https://bookdown.org/robinlovelace/geocompr/gis.html#rgrass
    initGRASS(
        gisDbase = tempdir(),
        home = tempdir(),
        location = system,
        mapset = "PERMANENT",
        override = TRUE
    )

    # Set projection, extent and resolution for the files
    p4_string <- terra::crs(combined_masks[[system]], proj = TRUE)
    execGRASS("g.proj",
        flags = c("c", "quiet"),
        proj4 = p4_string
    )
    b_box <- st_bbox((combined_masks[[system]]))
    execGRASS("g.region",
        flags = c("quiet"),
        n = as.character(b_box["ymax"]), s = as.character(b_box["ymin"]),
        e = as.character(b_box["xmax"]), w = as.character(b_box["xmin"]),
        res = as.character(res(combined_masks[[system]])[[1]])
    )

    # Add the raster to the GRASS environment
    write_RAST(
        combined_masks[[system]],
        vname = "mask",
        flags = "overwrite"
    )

    # Set the non-marsh areas as NULL
    execGRASS(
        cmd = "r.null",
        map = "mask",
        setnull = "0"
    )

    # Apply contraction of 1 cell to mask
    execGRASS(
        cmd = "r.grow",
        input = "mask",
        output = "mask",
        radius = -1.01,
        flags = "overwrite"
    )

    execGRASS(
        cmd = "r.grow",
        input = "mask",
        output = "mask",
        radius = 1.01,
        flags = "overwrite"
    )

    system_mask <- read_RAST("mask")

    # Create file name that includes
    file_name <- paste0(
        system, # the system name
        "-marsh-mask-", # raster content
        as.integer(max_sub_time * 100), # submersion time limit
        ".tif"
    )

    writeRaster(
        system_mask,
        here("outputs/rasters/marsh-masks", file_name),
        overwrite = TRUE
    )
}
