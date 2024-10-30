###################################################
## title ##
##          Márcio Martins    date           ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##
###################################################


# Plot class balance of a variable in folds -------------------------------

plot_class_balance <- function(folds, variable) {
    folds_df <- map2(
        folds$splits,
        folds$id,
        function(x, y) {
            x <- x$data[x$in_id, ]
            x$id <- y
            return(x)
        }
    )

    folds_df <- do.call(rbind, folds_df)

    return(
        ggplot(folds_df) +
            geom_bar(aes(x = id, fill = cluster)) +
            theme_bw()
    )
}


# Plot performance metrics from a workflow set ----------------------------
plot_metrics <- function(workflow_set) {
    results <- rank_results(model_tuning) |>
        separate_wider_delim(
            cols = wflow_id,
            names = c("rec", "mod"),
            delim = "_",
            too_many = "merge",
            cols_remove = FALSE
        )

    # From each model + recipe, keep the one with highest kappa
    keep_wfs <- results |>
        group_by(rec, model) |>
        filter(mean == max(mean[.metric == "kap"])) |>
        pull(rank)

    results |>
        filter(rank %in% keep_wfs) |>
        mutate(mod = fct_reorder(mod, mean, .desc = TRUE)) |>
        ggplot() +
        geom_pointrange(
            aes(
                x = mod,
                y = mean,
                ymin = mean - std_err,
                ymax = mean + std_err,
                color = rec
            ),
            position = position_dodge(0.5)
        ) +
        facet_grid(
            rows = vars(.metric),
            scale = "free_y"
        )
}


# Extract data.frames with predictions used to evaluate wf splits ---------

# Returns a data.frame with the predictions and truth for each split, with
# additional identifiers, used for different response variables and systems
get_pred_vs_truth <- function(wf, identifiers) {
    # Identifier for the models, which are not included in the workflow
    identifiers_df <- as_tibble(
        sapply(identifiers, \(x) rep(x, nrow(wf)))
    )

    # Names of workflows included in set
    wf_ids <- wf$wflow_id

    # Extract the full results for each workflow
    wf_results <- map(
        wf_ids,
        \(x) extract_workflow_set_result(wf, x)$.predictions
    )

    # Extract the row IDs used to evaluate each fold
    wf_results_rows <- map(
        wf_results,
        \(x) map(x, \(y) y$.row)
    )

    # Extract the full data for each fold
    wf_data <- map(
        wf_ids,
        \(x) extract_workflow_set_result(wf, x)$splits
    )

    # From the fold data, extract the system name for each row used to evaluate the splits
    wf_results_system <- map2(
        wf_data,
        wf_results_rows,
        \(x, y) map2(x, y, \(z, k) z$data$water_system[k])
    )

    # Add the water system identifier to the results data.frames
    wf_results <- map2(
        wf_results,
        wf_results_system,
        \(x, y) map2(
            x,
            y,
            function(z, k) {
                z$water_system <- k
                return(z)
            }
        )
    )

    results <- identifiers_df
    results$wf_id <- wf_ids
    results$results <- wf_results

    return(results)
}
