libs <- c(
    "tmap", 
    "stringr",
    "jsonlite",
    "PTtidaltools",
    "terra",
    "viridis"
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

# Create a water mask based on a Normalized Water Index raster
mask_water <- function(planet_raster, ndwi_threshold = 0) {
    planet_raster <- terra::rast(planet_raster)

    ndwi <- (planet_raster$green - planet_raster$nir) /
        (planet_raster$green + planet_raster$nir)

    water_mask <- ndwi > ndwi_threshold

    return(water_mask)
}

# Calculate the tide for a planet raster, based on the acquisition time given
# in the metadata file (assumes standard folder and metadata structure from a planet request)
get_image_tide <- function(raster_vec, port_id) {
    # Extract components used to identify the raster item
    name_components <- stringr::str_split(basename(raster_vec), "_")
    # Combine components in way that identifies corresponding
    # metadata file
    metadata_file <- sapply(
        name_components,
        \(x) paste(
            x[1:4],
            collapse = "_"
        )
    )

    metadata_file <- paste0(
        dirname(raster_vec),
        "/",
        metadata_file,
        "_metadata.json"
    )

    capture_times <- vector("character")

    for (file in metadata_file) {
        metadata <- jsonlite::read_json(file)

        capture_time <- metadata$properties$acquired
        capture_time <- stringr::str_split(capture_time, "\\.")[[1]][1]
        capture_time <- stringr::str_replace(capture_time, "T", " ")

        capture_times <- c(capture_times, capture_time)
    }


    tides <- PTtidaltools::interpolate_tides(
        date_times = capture_times,
        port_id = port_id,
        timezone = "UTC"
    )

    return(tides)
}

# Given a list of rasters, calculate the extent necessary to cover all of them
get_full_extent <- function(rasters) {
    full_extent <- terra::ext(rasters[[1]])
    if (length(rasters) > 1) {
        for (i in 2:length(rasters)) {
            new_ext <- ext(rasters[[i]])

            if (new_ext[1] < full_extent[1]) {
                full_extent[1] <- new_ext[1]
            }

            if (new_ext[2] > full_extent[2]) {
                full_extent[2] <- new_ext[2]
            }

            if (new_ext[3] < full_extent[3]) {
                full_extent[3] <- new_ext[3]
            }

            if (new_ext[4] > full_extent[4]) {
                full_extent[4] <- new_ext[4]
            }
        }
    }

    return(full_extent)
}

classify_tidal_regime <- function(raster_vec, port_id, roi_mask = NULL) {
    tidal_regimes <- list()
    # Classify rasters as water / not water
    water_masks <- lapply(
        raster_vec,
        mask_water
    )

    if (!is.null(roi_mask)) {
        water_masks <- lapply(
            water_masks,
            \(x) crop(x, roi_mask)
        )
    }

    # Estimate tidal height at time of image capture
    tides <- get_image_tide(raster_vec, port_id)

    # Homogenize raster extents so rasters can be stacked
    # First get the extent that covers all target rasters
    total_extent <- get_full_extent(water_masks)

    # Now loop through all the water masks and classify water according to
    # the tidal height at time of image capture
    for (i in 1:length(water_masks)) {
        # Extend raster layer to target extent - done one by one to limit memory usage
        water_mask <- extend(water_masks[[i]], total_extent)
        # Estimated tidal height for this image
        tidal_height <- tides[i]

        # Reclassify water mask, so that water is assigned the tidal
        # height at time of capture, and dry area is assigned NA
        class_matrix <- matrix(
            c(
                TRUE, FALSE,
                tidal_height, NA
            ),
            ncol = 2
        )

        tidal_mask <- classify(
            water_mask,
            class_matrix
        )

        rm(water_mask)
        gc()

        # To establish the minimum detected elevation, rasters will be
        # Compared iteratively, instead of creating one large stack
        # to prevent memory overload, at the cost of speed
        # If we are looking at the first image, establish it as the new estimation
        if (i == 1) {
            estimated_elev <- tidal_mask
            # If we have more images, stack with previous estimate and only keep the minimum values
        } else {
            estimated_elev <- c(estimated_elev, tidal_mask)
            estimated_elev <- min(estimated_elev, na.rm = TRUE)
        }
    }

    tidal_regimes$elev <- estimated_elev

    tidal_regimes$supratidal <- is.na(estimated_elev)

    return(tidal_regimes)
}

create_raster_composite <- function(raster_vec, prob = 0.5, roi_mask = NULL) {
    rasters <- lapply(raster_vec, terra::rast)
    # Crop rasters if ROI was provided
    if (!is.null(roi_mask)) {
        roi_mask <- terra::vect(roi_mask)
        rasters <- lapply(rasters, \(x) crop(x, roi_mask))
    }

    # Ensure that all rasters have the same number of layers
    if (length(unique(sapply(rasters, dim)[3, ])) != 1) {
        stop("Not all rasters have the same number of layers,\n
             can not create composite.")
    }

    # Calculate the largest extent in our rasters
    total_extent <- get_full_extent(rasters)

    # Get the number of layers
    n_layers <- unique(sapply(rasters, dim)[3, ])
    # And the name of those layers
    layer_names <- names(rasters[[1]])

    # We now discard the full rasters, we'll import individual
    # layers to minimize memory usage
    rm(rasters)
    gc()

    # Loop over the number of layers in our rasters
    # and store each composite layer as a temporary file
    composite_raster <- c()

    for (n_layer in 1:n_layers) {
        message(paste("Creating composite for layer", n_layer))
        layers <- lapply(raster_vec, \(x) terra::rast(x, lyrs = n_layer))

        # And now extend AND crop all rasters to enforce equal extents
        layers <- lapply(layers, \(x) extend(x, total_extent))
        layers <- lapply(layers, \(x) crop(x, total_extent))

        tmp_file1 <- paste0(tempfile(), ".tif")
        # Now stack all the layers
        layers <- terra::rast(layers)
        writeRaster(layers, tmp_file1)
        layers <- terra::rast(tmp_file1)

        gc()

        # Save individual layers as temp files
        tmp_file2 <- paste0(tempfile(), ".tif")
        composite_raster[n_layer] <- tmp_file2


        layers <- terra::quantile(layers, probs = prob, na.rm = TRUE)
        writeRaster(layers, tmp_file2)
        rm(layers)
        gc()
    }

    composite <- terra::rast(composite_raster)
    names(composite) <- layer_names

    return(composite)
}

calculate_overlap_per_pixel <- function(raster_vec) {
    rasters <- lapply(
        raster_vec,
        rast
    )

    # Only one layer per raster is needed for this method
    layers <- lapply(
        rasters,
        \(x) subset(x, 1)
    )

    layers <- lapply(
        layers,
        \(x) !is.na(x)
    )

    # Homogenize raster extents so rasters can be stacked
    # First get the extent that covers all target rasters
    total_extent <- get_full_extent(layers)

    layers <- lapply(
        layers,
        \(x) extend(x, total_extent)
    )

    # Stack one layer per available raster
    layers <- rast(layers)

    # Now sum across layers to get the number of overlaps in each cell
    overlap_count <- sum(layers, na.rm = TRUE)

    return(overlap_count)
}

# @ rasters A single or list of terra::rast objects
# @ rgb_layers Vector indicating the name of the layers used to plot RGB images
# @ vars Vector with the name of the additional layers to plot along the RGB image
# @ roi A single or list of terra::vect objects from which extents will be taken to
#      plot spatial facets
plot_raster_matrix <- function(rasters,
                               rgb_layers = c("red", "green", "blue"),
                               vars = NULL,
                               roi,
                               scale = 1) {
    # Standardize extent of the rasters to the max extent of the ROIs
    if (!is.list(roi)) {
        roi <- list(roi)
    }

    max_ext <- get_full_extent(roi)

    if (is.list(rasters)) {
        rasters <- lapply(
            rasters,
            \(x) crop(aggregate(x, scale), max_ext, extend = TRUE)
        )
        rasters <- do.call(c, rasters)
    } else {
        rasters <- crop(aggregate(rasters, scale), max_ext, extend = TRUE)
    }


    # Calculate highest value from RGB layers for scaling values RGB plot
    max_rgb <- max(
        minmax(rasters[[rgb_layers[1]]]),
        minmax(rasters[[rgb_layers[2]]]),
        minmax(rasters[[rgb_layers[3]]])
    )

    # Create empty layer to allow RGB image to show in a facet
    rasters$RGB <- NA

    # Select variables to be plotted
    var_names <- c("RGB", names(rasters)[names(rasters) %in% vars])


    # Now loop over all ROIs and create plots for each one
    maps <- list()
    for (i in 1:length(roi)) {
        roi_raster <- crop(rasters, roi[[i]])

        roi_map <- tm_shape(roi_raster) +
            tm_rgb(
                r = which(names(roi_raster) == rgb_layers[1]),
                g = which(names(roi_raster) == rgb_layers[2]),
                b = which(names(roi_raster) == rgb_layers[3]),
                max.value = max_rgb
            ) +
            tm_shape(roi_raster[[var_names]]) +
            tm_raster(palette = viridis::viridis(10), title = "") +
            tm_facets(
                ncol = length(var_names),
                free.scales.raster = TRUE
            ) +
            tm_compass() +
            tm_scale_bar(bg.color = "white", bg.alpha = .8) +
            tm_layout(
                outer.margins = 0,
                legend.bg.color = "white",
                legend.bg.alpha = .8,
                panel.show = ifelse(i == 1, TRUE, FALSE) # Show panel titles for first map
            )

        maps[[i]] <- roi_map
    }

    # Convert the tmap objects to grob for usage with grid
    maps_grobs <- lapply(maps, tmap_grob)

    # Get the heights of the rows used by tmap to delimite a layout
    row_heights <- lapply(maps_grobs, \(x) x$vp[[2]]$layout$heights)
    # Get the row at which the rasters are drawn
    row_pos <- sapply(maps_grobs, \(x) unique(x$children$multiple_1$vp$layout.pos.row))

    # Sum the height of the multiple row with facet label row - gives us
    # the height of the actual drawn map
    map_heights <- mapply(
        FUN = \(height, row) height[row] + height[row - 1],
        height = row_heights,
        row = row_pos
    )

    # Height of all rows combined
    total_height <- sum(map_heights)

    # If height is over one, it won't fit the current device, apply scaling
    if (total_height <= 1) {
        vp_height <- 1 # height for view port
        free <- 1 - total_height # empty space in graph device
        occupied <- c(0, map_heights[-length(map_heights)]) # how much space is taken by previous rows
        occupied <- occupied + free / 2 # shift all rows to vertical center
    } else {
        vp_height <- 1 / total_height # scale viewport height to scale grob sized and ensure they all fit
        map_heights <- map_heights / sum(map_heights) # scale the height of all rows so the total is between 0 and 1
        occupied <- c(0, map_heights[-length(map_heights)]) # how much space is take by previous rows
    }

    # now we calculate the center point for all rows
    start_pos <- mapply(
        FUN = function(row_height, occupied) {
            free_height <- 1 - occupied
            free_center_offset <- free_height - 0.5
            grob_center <- 0.5 + (free_center_offset - row_height / 2)
            return(grob_center)
        },
        row_height = map_heights,
        occupied = occupied
    )


    pdf(NULL)
    dev.control(displaylist = "enable")
    grid::grid.newpage()
    grid::pushViewport(grid::viewport())
    for (i in 1:length(maps_grobs)) {
        row_grob <- maps_grobs[[i]]
        row_center <- start_pos[i]
        grid::pushViewport(grid::viewport(
            height = vp_height,
            y = row_center
        ))
        grid::grid.draw(row_grob)
        grid::popViewport()
    }
    final_map <- recordPlot()
    invisible(dev.off())

    return(final_map)
}
