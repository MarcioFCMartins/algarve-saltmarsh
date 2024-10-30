##   Manuscript supplementary material figures   ##
##          Márcio Martins    2024-03-27         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##


# Setup -------------------------------------------------------------------

libs <- c(
    "ggplot2",
    "ggthemr",
    "patchwork",
    "scales",
    "stringr",
    "sf",
    "dplyr",
    "tidyr",
    "here",
    "terra",
    "janitor",
    "tidymodels"
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

ggthemr("flat")

ZENODO_REPO_DIR <- here("analysis/manuscript/saltmarsh/zenodo repo")
FIGURES_DIR <- here("analysis/manuscript/saltmarsh/figures")
TABLES_DIR <- here("analysis/manuscript/saltmarsh/tables")

# Copy files for zenodo repository ----------------------------------------
files <- c(
    here("data/clean/mapping/all-systems/system-masks"),
    here("data/clean/mapping/all-systems/manual-clean-up-masks.gpkg"),
    here("data/clean/mapping/all-systems/saltmarsh-training-data.gpkg"),
    here("data/clean/mapping/all-systems/saltmarsh-vegetation-quadrats.gpkg"),
    here("data/clean/mapping/all-systems/saltmarsh-transect-metadata.csv"),
    here("outputs/rasters/predicted-map.tif")
)

sapply(
    files,
    function(x) {
        if(dir.exists(x)) {
            zip(
                here(ZENODO_REPO_DIR, paste0(basename(x), ".zip")),
                list.files(x, full.names = TRUE, recursive = TRUE),
                flags = "-rj9X"
            )
        } else {
            file.copy(
                x,
                here(ZENODO_REPO_DIR),
                recursive = TRUE,
                overwrite = TRUE
            )
        }
    }
)

# Figure S1 - image tides histogram ---------------------------------------

image_tides <- read.csv2(
    here("outputs/tables/tide-heights-acquired-images.csv")
)

year_tides <- read.csv2(
    here("outputs/tables/tides_2022.csv")
)

# Create labels with total number of image per system
facet_labels <- sapply(
    unique(image_tides$water_system),
    \(x) paste0(
        str_to_title(x),
        " (",
        table(image_tides$water_system)[names(table(image_tides$water_system)) == x], ")"
    )
)

ggplot() +
    geom_histogram(
        data = year_tides,
        aes(x = height)
    ) +
    geom_rug(
        data = image_tides,
        aes(x = height)
    ) +
    facet_grid(
        cols = vars(water_system),
        labeller = as_labeller(facet_labels)
    ) +
    labs(
        x = "Tide height (m vs HZ)",
        y = "Observed tidal height count"
    )

ggsave(
    here(FIGURES_DIR, "S1 - image tides histogram.png"),
    width = 8,
    height = 3,
    dpi = 300
)

# Table S1 - Sampling data summary --------------------------------------
transect_metadata <- read.csv("./data/clean/mapping/all-systems/saltmarsh-transect-metadata.csv")

veg_quads <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-vegetation-quadrats.gpkg")
) |>
    mutate(
        cluster = as.factor(cluster),
        marsh_type = factor(marsh_type, levels = c("low", "middle", "high"))
    )

quad_count <- veg_quads |>
    st_drop_geometry() |>
    group_by(transect) |>
    summarise(n_quads = n())

transect_summary <- transect_metadata |>
    left_join(quad_count, by = "transect") |>
    select(-notes)

write.csv(
    transect_summary,
    here(TABLES_DIR, "S1 - Transect summary.csv"),
    row.names = FALSE

)


# Table S2 - Class balance in training data -------------------------------
dataset <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-training-data.gpkg")
) |>
    mutate(marsh_type = factor(marsh_type, levels = c("low", "middle", "high")))

data_summary <- dataset |>
    st_drop_geometry() |>
    mutate(train = ifelse(train, "Train", "Test")) |>
    group_by(train, marsh_type, water_system) |>
    summarise(n_water_system = n()) |>
    group_by(train, marsh_type) |>
    mutate(
        perc_water_system = scales::percent(n_water_system / sum(n_water_system)),
        n_marsh_type = sum(n_water_system)
    ) |>
    group_by(train) |>
    mutate(perc_marsh = scales::percent(n_marsh_type / sum(n_water_system))) |>
    ungroup() |>
    select(
        "Data set" = train, "Community" = marsh_type, "N per community" = n_marsh_type,
        "% per community" = perc_marsh, "Study system" = water_system,
        "N per water system" = n_water_system, "% per water system" = perc_water_system
    )

write.csv(
    data_summary,
    here(TALBES_DIR, "S2 - dataset summary.csv"),
    row.names = FALSE
)

# TABLE S3 Final model performance per system  -----------------------------------------

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

table_conf_mat_system <- janitor::tabyl(test_data, .pred_class, marsh_type, water_system) |>
    adorn_totals(c("col", "row")) |>
    adorn_title(
        row_name = "Predicted",
        col_name = "Truth",
        placement = "top"
    )

table_conf_mat_system <- bind_rows(table_conf_mat_system, .id = "water_system")


write.csv(
    table_conf_mat_system,
    here(TABLES_DIR, "S2 - Confusion Matrix per system.csv"),
    row.names = FALSE
)



