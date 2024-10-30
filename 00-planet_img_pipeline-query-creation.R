# Prepare the queries for Planet satellite images
# and export them as a .csv to be used by `planet_img_pipeline`
# Márcio Martins - marciomartinsred@gmail.com
# 2023-03-15

# Setup -------------------------------------------------------------------
if ("PTtidaltools" %in% installed.packages()) {
    library(PTtidaltools)
} else {
    devtools::install_github("https://github.com/MarcioFCMartins/PTtidaltools")
    library(PTtidaltools)
}
library(stringr)
library(tidyr)

# Prepare elements required to combine ------------------------------------
# Tidal height quantiles
tides <- get_tides(19, "2022-01-01", "2022-12-31") # Get tides for 2022
tide_quantiles <- quantile(tides$height, probs = seq(0, 1, 0.1))
tide_ranges <- paste(tide_quantiles, tide_quantiles[-1], sep = ",")[1:10]

# Path for region of interest (relative to the program's working folder)
# First do only alvor - seems like images overlap with arade
rois <- list.files("~/Repos/planet_img_pipeline/inputs/", pattern = "*.geojson")
rois <- paste0("./inputs/", rois)

# Start and end dates
start_date <- c("2022-01-01", "2022-04-01", "2022-07-01", "2022-10-01")
end_date <- c("2022-04-01", "2022-07-01", "2022-10-01", "2022-12-31")

cloud_cover <- "0.1"

asset_type <- "ortho_analytic_8b_sr"

n_layers <- 2


# Create all combinations -------------------------------------------------

combinations <- expand.grid(rois, start_date, cloud_cover, asset_type, tide_ranges, n_layers)
names(combinations) <- c("roi", "start_date", "max_cloud_cover", "asset_type", "min_tide", "n_layers")

combinations$end_date <- end_date[match(combinations$start_date, start_date)]
combinations$port <- ifelse(str_detect(combinations$roi, "ria-formosa"), 19, 18)

combinations <- separate(
    combinations,
    "min_tide",
    into = c("min_tide", "max_tide"),
    sep = ","
)

col_order <- c(
    "roi", "start_date", "end_date", "max_cloud_cover",
    "asset_type", "min_tide", "max_tide", "port", "n_layers"
)
combinations <- combinations[, col_order]

write.csv(
    combinations,
    "~/tmp/image-queries.csv",
    row.names = FALSE
)
