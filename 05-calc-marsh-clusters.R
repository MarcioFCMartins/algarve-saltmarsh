###################################################
##  Creation of saltmarsh mask through k-means   ##
##          Márcio Martins    2023-08-28         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##
###################################################


# Setup -------------------------------------------------------------------
libs <- c(
    "ClusterR", # clustering algorithms - using rcpp implementation of k-means for speed
    "terra", # raster manipulation
    "sf",
    "foreach",
    "doParallel",
    "here"
)

# Install missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

here::i_am("./code/05-saltmarsh-mask.R")

# Define clustering to be ran ---------------------------------------------

# System names - used to find rasters for each system
study_systems <- c(
    "alvor-updated",
    "arade-updated",
    "ria-formosa-updated",
    "guadiana-updated"
)


### K-means -- Raw variables
# raw_kmeans <- F          # Should this be ran?
train_features <- c(
    "blue", "green", "red", "nir",
    "ndwi_low", "ndvi", "subtime"
)
include_coords <- T # Should coordinates be used for clusters
n_clusters <- 25 # Number of clusters to create

# Run clustering on raw data ----------------------------------------------
# Setup parallel processing
all_cores <- parallel::detectCores(logical = FALSE)
all_cores <- ifelse(all_cores > 12L, 12L, all_cores)
cl <- makePSOCKcluster(all_cores, outfile = "")
registerDoParallel(cl)


folder_name <- paste0(
    "k", n_clusters, "-",
    paste0(substr(train_features, 1, 1), collapse = ""),
    ifelse(include_coords, "xy", "")
)


foreach(
    i = 1:length(study_systems),
    .final = \(x) NULL,
    .packages = c("terra", "ClusterR", "here", "sf")
) %dopar% {
    system <- study_systems[i]

    target_path <- here(paste0(
        "outputs/rasters/system-clusters/",
        folder_name,
        "/",
        system,
        ".tif"
    ))

    # If directory does not exist, create it
    if (!file.exists(dirname(target_path))) {
        dir.create(dirname(target_path))
    }

    # Skip loop if target raster already exists
    if (file.exists(target_path)) {
        cat(paste(system, "raster already exists. Skipped\n"),
            file = stdout()
        )
        return(NULL)
    }

    if (!file.exists(here(dirname(target_path), "clustering-vars.txt"))) {
        write(
            train_features,
            here(dirname(target_path), "clustering-vars.txt")
        )
    }


    # Read rasters of interest for current water system
    files <- list.files(
        here("outputs/rasters"),
        pattern = system,
        recursive = TRUE,
        full.names = TRUE
    )

    files <- files[grepl("q20|submersion|texture|ndvi|ndwi", basename(files))]
    rasters <- lapply(
        files,
        terra::rast
    )

    # Stack rasters and only keep desired layers
    rasters <- do.call(c, rasters)
    rasters <- subset(rasters, train_features)

    roi_mask <- list.files(
        here("outputs/system-masks/"),
        pattern = system,
        full.names = TRUE
    ) |>
        st_read() |>
        st_transform(crs(rasters)) |>
        vect()

    rasters <- crop(rasters, roi_mask, mask = TRUE)

    # Prepare data for k-means
    if (include_coords) { # Convert raster to df
        data <- as.data.frame(rasters, xy = TRUE)
    } else {
        data <- as.data.frame(rasters)
    }
    data <- na.omit(data)
    data <- as.data.frame(scale(data)) # Center and scale variables

    # Cluster all pixels
    mask_clusters <- KMeans_rcpp(
        data,
        clusters = n_clusters,
        num_init = 10,
        seed = 99,
        verbose = TRUE
    )

    # Rebuild raster, with clustered pixels
    cluster_raster <- rast( # Initialize empty raster
        extent = ext(rasters),
        resolution = res(rasters)
    )
    crs(cluster_raster) <- crs(rasters) # Same crs
    empty_cell <- values(any(is.na(rasters))) # True for cells that have been masked out
    values(cluster_raster)[!empty_cell] <- mask_clusters$clusters # Set values for non-masked cells in new raster


    writeRaster(
        cluster_raster,
        target_path
    )
}

stopCluster(cl)
