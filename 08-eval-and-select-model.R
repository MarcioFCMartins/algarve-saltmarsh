###################################################
##        Model evaluation and selection         ##
##          Márcio Martins    2024-02-07         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##
###################################################


# Setup -------------------------------------------------------------------

libs <- c(
    "tidymodels", # Wrapper for several machine learning packages
    "sf", # Simple-features for vector spatial data
    "here", # Path normalization
    "tidyverse" # Data wrangling
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T))

source(here("code/helpers-models.R"))

FINAL_MODEL_PATH <- here("outputs/models/final-model.RDS")

# Load models -------------------------------------------------------------

uni_model <- list.files(
    here("outputs/models/universal"),
    pattern = "-crosscv.RDS",
    full.names = T
) |>
    map(readRDS)

names(uni_model) <- str_replace(
    list.files(here("outputs/models/universal"), pattern = "-crosscv.RDS"),
    "-crosscv.RDS",
    ""
)

system_model <- list.files(
    here("outputs/models/per-system"),
    full.names = T,
    pattern = "-crosscv.RDS"
) |>
    map(readRDS)

names(system_model) <- str_replace(
    list.files(here("outputs/models/per-system"), pattern = "-crosscv.RDS"),
    "-crosscv.RDS",
    ""
)

# Combine results from all candidate models  ------------------------------

# Get results from model that used all systems' data
uni_model_results <- map2(
    uni_model,
    names(uni_model),
    \(x, y) get_pred_vs_truth(x, c("system" = "all", "response" = y))
) |>
    bind_rows()

# Get results from individual systems' models
system_model_results <- map2(
    system_model,
    names(system_model),
    function(x, y) {
        system <- str_extract(y, "^[^-]+")
        response <- str_extract(y, "(?<=-)[^-]+$")
        get_pred_vs_truth(x, c("system" = system, "response" = response))
    }
) |>
    bind_rows()

all_results <- bind_rows(uni_model_results, system_model_results)


# Translate models that predict "cluster" to a "marsh_type" response --------

cluster_vegetation_type <- read.csv(
    here("outputs/cluster-community-dict.csv")
) |>
    mutate(
        cluster = as.factor(cluster),
        marsh_type = as.factor(marsh_type)
    )

# Create a new column, where only the truth/predictions will be stored,
# after converting all results to "marsh_type
final_perf_results <- list()
for (i in 1:nrow(all_results)) {
    wf_results <- all_results$results[[i]]
    if (length(wf_results[[1]]) == 1) next # Skip if model training failed

    # If this workflow predicts marsh type, we just need to extract the predictions and truth
    if (all_results$response[i] == "marsh_type") {
        wf_final_results <- map(
            wf_results,
            \(x) x[, c(".pred_class", "marsh_type")]
        )
    } else {
        # If the workflow is prediction cluster, we need to
        # translate both truth and predictions and return a tibble with those translations
        wf_final_results <- map(
            wf_results,
            \(x) {
                # First translate the truth - add the marsh type to the data
                truth <- left_join(
                    x,
                    cluster_vegetation_type,
                    by = c("water_system", "cluster")
                ) |>
                    # And extract that vector
                    pluck("marsh_type")

                # Then translate the predictions
                predictions <- left_join(
                    x,
                    cluster_vegetation_type,
                    by = c("water_system", ".pred_class" = "cluster")
                ) |>
                    pluck("marsh_type")

                return(
                    tibble(
                        ".pred_class" = predictions,
                        "marsh_type" = truth
                    )
                )
            }
        )
    }
    # Add the wf final results to the results of remaining workflows
    final_perf_results[[i]] <- wf_final_results
}

all_results$final_perf_results <- final_perf_results


# Model performance metrics -----------------------------------------------

# Calculate the kappa for every candidate model and all folds
kappa_s <- map(
    all_results$final_perf_results,
    \(x) map_dbl(x, \(y) kap(y, truth = marsh_type, estimate = .pred_class)$.estimate)
)

all_results$kappa <- kappa_s


# Save the results per model for later use
saveRDS(
    all_results,
    here("outputs/models/all_model_results.RDS")
)

all_results |>
    # Now turn every fold's kappa into an individual row for plotting
    unnest(kappa) |>
    mutate(wf_id = fct_reorder(wf_id, kappa, .desc = TRUE)) |>
    ggplot(aes(y = kappa, x = wf_id)) +
    geom_point() +
    geom_line() +
    stat_summary(
        geom = "point",
        shape = "-",
        fun = "median",
        color = "red",
        size = 10
    ) +
    facet_grid(
        cols = vars(response),
        rows = vars(system)
    ) +
    labs(
        title = "Model performance",
        x = "Model",
        y = "Kappa"
    ) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1)
    )



# Select and train the best model -----------------------------------------
full_data_sf <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-training-data.gpkg"),
    quiet = TRUE
) |>
    mutate(cluster = as.factor(cluster))


train_data <- full_data_sf |>
    filter(train == TRUE) |>
    st_drop_geometry()


best_model_parameters <- uni_model$marsh_type |>
    extract_workflow_set_result("raw_rf") |>
    select_best(metric = "kap")

final_wf <- uni_model$marsh_type |>
    extract_workflow("raw_rf") |>
    finalize_workflow(best_model_parameters) |>
    fit(data = train_data)

saveRDS(final_wf, here("outputs/models/final-model.RDS"))


