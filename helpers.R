####################################################
## Helper functions for South Portugal BC mapping ##
####################################################

# Márcio Martins - marciomartinsred@gmail.com
# 2023-03-21



# 1. Community clustering ----------------------------------------------------
## 1.1 Reads excel files with percent cover matrices ---------------------------
# Appends all sheets, using sheet name to identify the transect
# And replaces cover percent with presence/absence
read_presence_matrix <- function(file) {
    presence_wide <- data.frame() # Store all quadrats

    for (sheet in excel_sheets(file)) {
        matrix <- read_xlsx(file, sheet = sheet)

        # Replace cover percentages with presence/absence
        matrix[, -1] <- replace(
            x = matrix[, -1],
            list = matrix[, -1] > 0,
            values = 1
        )

        # Fix any missing values
        matrix[is.na(matrix)] <- 0
        distance_from_water <- colnames(matrix)[-1]

        matrix <- t(matrix)
        colnames(matrix) <- matrix[1, ]
        matrix <- as.data.frame(matrix[-1, ]) %>%
            mutate_all(.funs = function(x) as.numeric(as.character(x)))

        matrix <- cbind(
            as.numeric(sheet),
            as.numeric(distance_from_water),
            matrix
        )
        names(matrix)[c(1, 2)] <- c("transect_id", "distance_from_water")


        presence_wide <- bind_rows(presence_wide, matrix)
        # Remove not identified species and walking platform records
        presence_wide <- select(presence_wide, -any_of(c("NID", "Walking platform")))
    }

    presence_wide[is.na(presence_wide)] <- 0

    # Add column which identifies if it is unvegetated
    is_vegetated <- rowSums(presence_wide[, -c(1, 2)])
    presence_wide$vegetated <- ifelse(is_vegetated > 0, 1, 0)

    # Reorder columns
    presence_wide <- select(
        presence_wide,
        transect_id, distance_from_water, vegetated, everything()
    )

    return(presence_wide)
}


# 2 Determine quadrat positions ---------------------------------------------

## 2.1 Read transect coordinates -----------------------------------------------
read_transects <- function(file, data_crs, region_crs) {
    # Start and end point coordinates of each transect
    points <- read_xlsx(file, sheet = "saltmarsh_data") |>
        sf::st_as_sf(coords = c("longitude", "latitude"), crs = data_crs) |>
        # Angle operations assume a planer projection, so we must convert
        sf::st_transform(crs = region_crs)

    points <- points |>
        dplyr::select(-notes) |>
        cbind(sf::st_coordinates(points)) |>
        sf::st_drop_geometry()

    # Transect lengths
    lengths <- read_xlsx(file, sheet = "saltmarsh_summary") |>
        dplyr::select(transect_id, name, length)

    # Identifying information for individual transects
    info <- read_xlsx(file, sheet = "saltmarsh_data") |>
        dplyr::filter(start_end == "Start") |>
        dplyr::select(1:6)

    crs <- c("data_crs" = data_crs, "region_crs" = region_crs)
    return(list(points, lengths, info, crs))
}


##  2.2 Calculate vector orientation -------------------------------------------
# https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r/24999820
# https://www.math.uh.edu/~jmorgan/Math6397/day13/LinearAlgebraR-Handout.pdf - page 5

# 1. translate the transect vector to the origin (i.e. calculate delta x and delta y)
# 2. create a second vector b, with a y of 1 and x of zero (aka North)
# 3. Returns the radians between those 2 vectors (orientation)
get_vector_angle <- function(delta_x, delta_y) {
    b <- c(delta_x, delta_y)
    a <- c(0, 1)
    # This calculates the angle to a vertical line (between 0 and pi)
    theta <- acos(sum(a * b) / (sqrt(sum(a * a)) * sqrt(sum(b * b))))
    # Correct to be between 0 and 2pi (conditions change depending on half of the unit circle it lands)
    if (delta_x >= 0) {
        theta
    } else {
        (pi - theta) + pi
    }
}


## 2.3 Break vector into known number of equal segments ------------------------
# Takes the start position, orientation and length of a vector
# Breaks it into n segments
# Returns a data.frame with points for center of each segment
break_vector <- function(name, start_x, start_y, theta, length, n_segments) {
    # How long should each segment be
    scale <- length / n_segments
    # Offset each point by half the length of the segment
    # So that we return the center point
    point_offset <- scale / 2

    points <- data.frame(
        "quad_index" = numeric(),
        "name" = character(),
        "x" = numeric(),
        "y" = numeric()
    )

    for (i in 0:(n_segments - 1)) {
        x <- start_x + (sin(theta) * scale * i)
        y <- start_y + (cos(theta) * scale * i)
        x <- x + sin(theta) * point_offset
        y <- y + cos(theta) * point_offset
        index <- i

        new_point <- data.frame(
            "quad_index" = index,
            "name" = name,
            "x" = x,
            "y" = y
        )
        points <- rbind(points, new_point)
    }
    points
}


## 2.4 Calculate quadrat central points ------------------------------------
# Takes the list returned by `read_transect` and calculates the central points
# for each quadrat in each transect
get_quatrat_centers <- function(transects) {
    points <- transects[[1]]
    lengths <- transects[[2]]
    info <- transects[[3]]
    crs_info <- transects[[4]]

    # Reshape data so the x and y of start and end points are individual columns
    orientations <- points |>
        # Combine lat and lon in a single column
        gather(key = variable, value = value, X, Y) |>
        # Combine identifier of lat/lon with combiner of transect start/end
        unite(col = temp, start_end, variable) %>%
        # Now spread coordinates with 4 coordinate columns: lat start, lat end, lon start, lon end
        spread(temp, value)

    orientations <- orientations |>
        # Translate vector to origin
        mutate(delta_x = End_X - Start_X, delta_y = End_Y - Start_Y) |>
        # Calculate angle vs north
        mutate(theta = map2_dbl(delta_x, delta_y, get_vector_angle)) |>
        select(name, theta)

    # Every quadrat was placed so that it began at the recorded meter
    #   i.e. quadrat 5 started at meter 5
    # The last quadrat also extended beyond the recorded GPS point
    #   (this was accounted for when calculating transect length by adding 1 meter
    #   to recorded value)

    # Convert points to a line
    lines <- points %>%
        st_as_sf(coords = c("X", "Y"), crs = crs_info["region_crs"]) |>
        group_by(name) %>%
        summarize(
            water_system = unique(water_system),
            water_body = unique(water_body),
            site_nr = unique(site_nr),
            transect_nr = unique(transect_nr),
            transect_id = unique(transect_id)
        ) %>%
        st_cast(to = "LINESTRING") %>%
        left_join(lengths)

    # Get the length and number of segments
    transect_geometry <- lines %>%
        # Get length of segment, so that transect is equally divided into k segments
        mutate(
            gps_length = map_dbl(geometry, st_length),
            # Each transect has a number of quadrats equal to length
            segment_number = length,
            segment_length = gps_length / segment_number
        ) %>%
        st_drop_geometry() %>%
        select(name, gps_length, segment_number, segment_length) |>
        left_join(orientations)

    starting_points <- points %>%
        filter(start_end == "Start") |>
        st_as_sf(coords = c("X", "Y"), crs = crs_info["region_crs"])

    # Join the starting point to the segment information
    transect_geometry <- transect_geometry %>%
        left_join(starting_points %>% select(name)) %>%
        st_as_sf(crs = crs_info["region_crs"])

    #  Move the coordinates to columns
    transect_geometry <- cbind(
        transect_geometry,
        st_coordinates(transect_geometry)
    ) %>%
        st_drop_geometry()

    # Break transects into segments (see helper function)
    quadrat_points <- lapply(
        split(transect_geometry, transect_geometry$name),
        FUN = function(x) {
            break_vector(
                name = x$name,
                start_x = x$X,
                start_y = x$Y,
                theta = x$theta,
                length = x$gps_length,
                n_segments = x$segment_number
            )
        }
    )

    quadrat_points <- do.call(rbind, quadrat_points) |>
        mutate(quad_index = as.numeric(quad_index)) |>
        st_as_sf(coords = c("x", "y"), crs = crs_info["region_crs"]) |>
        arrange(name, quad_index)

    quadrat_points <- left_join(quadrat_points, info, by = "name")
}
