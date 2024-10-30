##      Prepare vegetation data and clusters     ##
##          Márcio Martins    2023/07/28         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##

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
    "stringi",
    "purrr",
    "readxl",
    "lubridate"
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

here::i_am("code/01-calc-vegetation-communities.R")
source(here("code/helpers.R"))


# Create transect metadata table ------------------------------------------
## SET A ------------------------------------------------------------------
metadata_set_a <- read_xlsx(
    here("data/raw/mapping/all-systems/monipor_clean_data.xlsx"),
    sheet = "saltmarsh_data"
) |>
    filter(!duplicated(name)) |>
    group_by(water_system) |>
    mutate(
        date = dmy(date),
        # Create unique identifier by combining identifier of zone and site within zone
        site = paste0(water_body, site_nr),
        # Determine when unique site changes
        site_change = ifelse(site != lag(site), 1, 0)) 

metadata_set_a <- metadata_set_a |>
    mutate(
        # Add set identifier
        transect_set = "Set A",
        # Use the change points in unique sites to count number of unique sites
        site = cumsum(ifelse(is.na(site_change), 1, site_change))) |>
    select(water_system, transect_set, site, transect = name, date) |>
    ungroup()

# Combine both sets of notes per transect and add them to metadata
notes_set_a <- read_xlsx(
    here("data/raw/mapping/all-systems/monipor_clean_data.xlsx"),
    sheet = "saltmarsh_data"
) |>
    select(transect = name, notes) |>
    group_by(transect) |>
    summarise(notes = paste(unique(ifelse(is.na(notes), "", notes)), collapse = "; ")) 

metadata_set_a <- left_join(metadata_set_a, notes_set_a, by = "transect")

## SET B ---------------------------------------------------------
metadata_set_b <- read_xlsx(
    here("data/raw/mapping/all-systems/holdout-transects.xlsx"),
    sheet = "transects",
    skip = 3,
    na = "NA"
) |>
    mutate(
        water_system = str_to_title(system),
        transect_set = "Set B"
    ) |>
    select(
        water_system,
        transect_set,
        site = area,
        transect = trenasect,
        date
    ) |>
    filter(!duplicated(transect))
    
# Combine both sets of notes per transect and add them to metadata
notes_set_b <- read_xlsx(
    here("data/raw/mapping/all-systems/holdout-transects.xlsx"),
    sheet = "transects",
    skip = 3,
    na = "NA"
) |>
    select(transect =  trenasect, notes) |>
    group_by(transect) |>
    summarise(notes = paste(unique(ifelse(is.na(notes), "", notes)), collapse = "; ")) 

metadata_set_b <- left_join(metadata_set_b, notes_set_b, by = "transect")

### Correct missing site and transect name for alvor -------------------------
rows_to_correct <- metadata_set_b$water_system == "Alvor" & metadata_set_b$site == "?" & metadata_set_b$transect == "?"
metadata_set_b$site[rows_to_correct] <- "B"
metadata_set_b$transect[rows_to_correct] <- "alvor_B3"
metadata_set_b$notes[rows_to_correct] <- "This is not a normal transect, but rather points taken in high marsh along the water line"


## Combine and export data ---------------------------------------------------
metadata_transects <- rbind(metadata_set_a, metadata_set_b)
rm(metadata_set_a, metadata_set_b, rows_to_correct, notes_set_a, notes_set_b)

write.csv(
    metadata_transects,
    "./data/clean/mapping/all-systems/saltmarsh-transect-metadata.csv",
    row.names = F
)

# Read old MONIPOR transects and add spatial location ---------------------
presence_matrices <- list()

files <- list.files(
    here("data/raw/mapping/"),
    pattern = ".cover_matrices\\.xlsx",
    full.names = TRUE,
    recursive = TRUE
)

for (file in files) {
    system <- str_extract(file, "//(.+?)_cover")
    system <- sub("//", "", system)
    system <- strsplit(system, "/")[[1]][1]
    system <- sub("_(.+)", "", system)

    # Combines several excel sheets into a single data.frame
    # with transect IDs, quadrat distance from water
    # also removes unidentified species and walking platform records
    matrices <- read_presence_matrix(file)
    presence_matrices[[system]] <- matrices
    rm(matrices)
}

old_quadrats <- bind_rows(presence_matrices, .id = "water_system") |>
    mutate(water_system = case_when(
        str_detect(water_system, "alvor") ~ "Alvor",
        str_detect(water_system, "guadiana") ~ "Guadiana",
        str_detect(water_system, "arade") ~ "Arade",
        str_detect(water_system, "riaformosa") ~ "Ria Formosa"
    ))


# Add remaining information to these quadrats
info <- read_xlsx(
    here("data/raw/mapping/all-systems/monipor_clean_data.xlsx"),
    sheet = "saltmarsh_data"
) |>
    dplyr::filter(start_end == "Start") |>
    dplyr::select(1:6)

old_quadrats <- left_join(
    old_quadrats,
    info,
    by = c("water_system", "transect_id")
) |>
    select(water_system,
        site_name = site_nr, transect_id,
        distance_from_water, everything()
    )

old_quadrats[is.na(old_quadrats)] <- 0

# EPSG code for CRS PT-TM06/ETRS89
pttm06 <- 3763
# EPSG code for CRS WGS84
wgs84 <- 4326

transects <- read_transects(
    "./data/raw/mapping/all-systems/monipor_clean_data.xlsx",
    data_crs = wgs84,
    region_crs = pttm06
)

# Estimate the center point for each individual quadrat, based on the
# start point, end point and distance along the transect.
# Check source file for implementation details
quadrat_central_points <- get_quatrat_centers(transects) |>
    rename(distance_from_water = quad_index) |>
    select(water_system, name, distance_from_water)

old_quadrats <- left_join(
    old_quadrats,
    quadrat_central_points,
    by = c("water_system", "name", "distance_from_water")
) |>
    st_sf()

old_quadrats$transect_set <- "Set A"


# Read new transects ------------------------------------------------------

new_quadrats <- read_xlsx(
    here("data/raw/mapping/all-systems/holdout-transects.xlsx"),
    sheet = "transects",
    skip = 3,
    na = "NA"
) |>
    rename(
        name = trenasect, # Notice name typo
        distance_from_water = distance,
        water_system = system,
        site_name = area
    )

# Interpolate position of quadrats with missing points
new_quadrats <- new_quadrats |>
    group_by(name) |> # Only interpolate within the same transect
    mutate(
        longitude = approx(
            x = seq_along(longitude),
            y = longitude,
            xout = seq_along(longitude)
        )$y,
        latitude = approx(
            x = seq_along(latitude),
            y = latitude,
            xout = seq_along(latitude)
        )$y
    ) |>
    ungroup() |>
    mutate(
        water_system = case_when(
            str_detect(water_system, "alvor") ~ "Alvor",
            str_detect(water_system, "guadiana") ~ "Guadiana",
            str_detect(water_system, "arade") ~ "Arade",
            str_detect(water_system, "riaformosa") ~ "Ria Formosa"
        )
    )


## Convert new transects to presence/absence matrix ------------------------
# Start by creating a lookup table for corrected species names - the new
# names are based on the ones used for the training transects, even if the species
# name is outdated

# Names of species found in new data
sp_present <- stri_split(new_quadrats$species, regex = ",") |>
    unlist() |>
    stri_trim_both() |>
    unique() |>
    sort()

# Lookup names, based on the ones used for training
sp_dict <- list(
    "Unknown", # 1
    "Arthrocnemum macrostachyum", # 2
    "Arthrocnemum macrostachyum", # 3
    "Arthrocnemum macrostachyum", # 4
    "Arthrocnemum macrostachyum", # 5
    "Atriplex halimus", # 6
    "Atriplex halimus", # 7
    "Halimione portulacoides", # 8
    "Atriplex halimus", # 9
    "Unvegetated", # 10
    "Salsola vermiculata", # 11
    "Unknown", # 12
    "Halimione portulacoides", # 13
    "Limoniastrum monopetalum", # 14
    "Limonium ferulaceum", # 15
    "Limonium vulgare", # 16
    "Sarcocornia perennis", # 17
    "Salicornia sp.", # 18
    "Sarcocornia fruticosa", # 19
    "Sarcocornia fruticosa", # 20
    "Sarcocornia perennis", # 21
    "Sarcocornia perennis", # 22
    "Salicornia sp.", # 23
    c(
        "Sarcocornia perennis", # 24
        "Sarcocornia fruticosa"
    ),
    "Sarcocornia perennis", # 25
    "Sarcocornia fruticosa", # 26
    c(
        "Sarcocornia perennis", # 27
        "Sarcocornia fruticosa"
    ),
    "Sarcocornia perennis", # 28
    "Sarcocornia perennis", # 29
    "Sarcocornia fruticosa", # 30
    "Salsola vermiculata", # 31
    "Spartina maritima", # 33
    "Spartina densiflora", # 34
    "Spartina maritima", # 35
    "Suaeda vera", # 36
    "Suaeda vera", # 36
    "Unvegetated"
)

names(sp_dict) <- sp_present

# Create an empty presence / absence matrix
pres_matrix <- new_quadrats |>
    select(species)
# Names to be used in species columns
sp_correct <- unlist(sp_dict) |>
    unique()
pres_matrix_names <- c(names(pres_matrix), sp_correct)
# Add species columns
pres_matrix <- cbind(
    pres_matrix,
    matrix(0, nrow = nrow(pres_matrix), ncol = length(sp_correct))
)

names(pres_matrix) <- pres_matrix_names

# Now populate the matrix with presences
for (i in 1:nrow(pres_matrix)) {
    # Species present in this row
    row_sp <- pres_matrix[i, "species"] |>
        stringi::stri_split(regex = ",") |>
        unlist() |>
        stri_trim_both()

    # Corrected species names that match the ones found in the row
    sp_matches <- sp_dict[which(sp_present %in% row_sp)] |>
        unlist() |>
        unname()

    pres_matrix[i, sp_matches] <- 1
}

new_quadrats <- cbind(
    select(
        new_quadrats,
        water_system,
        name,
        distance_from_water,
        site_name,
        longitude,
        latitude
    ),
    select(pres_matrix, -species, -Unvegetated, -Unknown) # Remove original species column and "Unkown" species
)

new_quadrats <- st_as_sf(
    new_quadrats,
    coords = c("longitude", "latitude"),
    crs = 4326
) |>
    st_transform(crs = pttm06)

new_quadrats$transect_set <- "Set B"
# Combine old and new quadrats --------------------------------------------
quadrats <- bind_rows(
    mutate(old_quadrats, site_name = as.character(site_name)),
    mutate(new_quadrats, distance_from_water = as.numeric(distance_from_water))
) |>
    select(-water_body, -transect_nr, -site_name, -transect_id, transect = name) |>
    relocate(transect, transect_set, .after = water_system)

# Now fill in NAs on species with zero
quadrats <- mutate(quadrats, across(c(6:33), \(x) ifelse(is.na(x), 0, x))) |>
    # Add a vegetated bool
    rowwise() |>
    mutate(vegetated = as.numeric(sum(c_across(c(6:33))) >= 1)) |>
    ungroup() |>
    # add unique quadrats identifier
    mutate(quad_id = seq_along(water_system), .after = transect_set)


# Remove quadrats with unidentified Salicornia
quadrats <- filter(quadrats, `Salicornia sp.` != 1)
quadrats <- select(quadrats, -`Salicornia sp.`)

# Correct missing site and transect name for alvor -------------------------
rows_to_correct <- quadrats$water_system == "Alvor" & quadrats$transect == "?"
quadrats$transect[rows_to_correct] <- "alvor_B3"

# Correct species names to accepted ones ----------------------------------
accepted_names <- read.csv(
    here("data/clean/mapping/all-systems/saltmarsh-accepted-species.csv")
)

# Remove duplicated species
quadrats <- quadrats |>
    mutate(
        # Salsola vermiculata was misspelled at least once
        `Salsola vermiculata` = ifelse(
            `Salsola vermiculata` + `Salsola verniculata` >= 1, 1, 0
        ),
        # Most of the species identified as juncus maritimus were Spartina densiflora
        `Spartina densiflora` = ifelse(`Juncus maritimus` + `Spartina densiflora` >= 1, 1, 0),
        `Juncus effusus` = ifelse(`Juncus effusus` + `Juncus holoschoenos` >= 1, 1, 0),
    ) |>
    select(-c(`Salsola verniculata`, `Juncus holoschoenos`, `Juncus maritimus`))

# Indices for new names
name_index <- sapply(
    names(quadrats),
    \(x) which(accepted_names$old_name == x)
)

new_names <- c()
for (i in 1:length(name_index)) {
    if (length(name_index[[i]]) == 0) {
        new_names[i] <- names(quadrats)[i]
    } else {
        new_names[i] <- accepted_names$new_name[name_index[[i]]]
    }
}

names(quadrats) <- new_names

# Create vegetation clusters ----------------------------------------------
# Let's create a simple data.frame to avoid class conflicts
vegetated_quadrats <- quadrats |>
    filter(vegetated == 1) |>
    st_drop_geometry() |>
    as.data.frame()

# And now only keep the species, and split by water system
system_splits <- vegetated_quadrats$water_system
veg_matrices <- split(vegetated_quadrats, f = system_splits)

# Apply the clustering to the different systems
cluster_n_per_system <- c(
    "Alvor" = 3,
    "Arade" = 2,
    "Guadiana" = 6,
    "Ria Formosa" = 5
)


for (system in names(veg_matrices)) {
    # Get the species matrix for the system
    x <- veg_matrices[[system]]
    # And remove identifiers, only keep species presence/absence cols
    x <- as.data.frame(x[, -c(1:6)])

    cluster_n <- cluster_n_per_system[[system]]

    # Get pairwise Bray-Curtins binary distances (AKA Sorensen index)
    bray_distances <- vegdist(
        x,
        method = "bray",
        binary = TRUE
    )

    # Some references for usage of clustering:
    # Cluster. Ward's linkage was used based on
    # "Robustness of three hierarchical agglomerative clustering
    # techniques for ecological data" - See more references for what might be appropriate

    # See also:
    # "Use of the Bray-Curtis Similarity Measure in Cluster Analysis of Foraminiferal Data"

    # "Bacterial bioclusters relate to hydrochemistry in New Zealand groundwater"
    ## This paper even shows some more advanced clustering options

    # "Multivariate Analysis of Ecological Communities in R: vegan tutorial"; Chapter 6

    # SELECT CLUSTERS
    set.seed(42)
    dendrogram <- hclust(bray_distances, method = "ward.D")
    # plot(
    #     as.dendrogram(dendrogram),
    #     type = "triangle",
    #     center = TRUE)
    clusters <- cutree(dendrogram, k = cluster_n)

    veg_matrices[[system]]$cluster <- clusters
}

rm(x, bray_distances, cluster_n, dendrogram)

vegetated_quadrats <- do.call(rbind, veg_matrices)


# Assign cluster names to vegetation group --------------------------------

cluster_vegetation_type <- data.frame(
    water_system = c(
        rep("Alvor", cluster_n_per_system["Alvor"]),
        rep("Arade", cluster_n_per_system["Arade"]),
        rep("Guadiana", cluster_n_per_system["Guadiana"]),
        rep("Ria Formosa", cluster_n_per_system["Ria Formosa"])
    ),
    cluster = c(
        c(1:cluster_n_per_system["Alvor"]),
        c(1:cluster_n_per_system["Arade"]),
        c(1:cluster_n_per_system["Guadiana"]),
        c(1:cluster_n_per_system["Ria Formosa"])
    ),
    marsh_type = c(
        c("low", "middle", "high"),
        c("middle", "middle"),
        c("low", "low", "low", "low", "high", "middle"),
        c("low", "low", "middle", "middle", "high")
    )
)

write.csv(
    cluster_vegetation_type,
    here("outputs/cluster-community-dict.csv"),
    row.names = F
)

vegetated_quadrats <- left_join(
    vegetated_quadrats, cluster_vegetation_type
)


# And now add cluster information to the spatially referenced quad --------
vegetated_quadrats_sf <- quadrats |>
    # Only keep vegetated quadrats
    filter(vegetated == 1) |>
    # Add the cluster and marsh type column
    left_join(
        select(vegetated_quadrats, cluster, marsh_type, quad_id),
        by = "quad_id"
    ) |>
    # And relocate all information columns to the left of the matrix
    relocate(cluster, marsh_type, .after = vegetated)


# Save the data -----------------------------------------------------------
st_write(
    vegetated_quadrats_sf,
    here("data", "clean", "mapping", "all-systems", "saltmarsh-vegetation-quadrats.gpkg"),
    delete_dsn = TRUE
)
