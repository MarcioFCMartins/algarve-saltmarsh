## 1.1 Plot cluster species composition ----------------------------------------
# Accepts a wide-format presence/absence matrix with clustered quadrats
# and plots species presence per cluster
plot_composition <- function(data, sp_cols, group) {
    # Calculate the relative frequencies of species presence in quadrats per marsh community
    vegetation_composition_community <- data |>
        gather(key = "species", value = "presence", sp_cols) |>
        group_by(water_system, {{ group }}, species) |>
        summarise(prop_presences = sum(presence) / n()) |>
        ungroup()

    # Now for the plot:
    ## Only keep species which are present in at least 10% of obs for a community in a system
    common_sps <- vegetation_composition_community |>
        group_by(species) |>
        summarise(max_presence = max(prop_presences)) |>
        filter(max_presence > 0.10) |>
        pull(species)
    print(names(vegetation_composition_community))
    # And let's calculate a species ordering based on average presence across all systems
    sp_order <- vegetation_composition_community |>
        # Let's give each community a weight to make ordering account for it
        mutate(
            marsh_weight = case_when(
                marsh_type == "low" ~ 6,
                marsh_type == "middle" ~ 3,
                marsh_type == "high" ~ 1
            ),
            weight = prop_presences * marsh_weight
        ) |>
        group_by(species) |>
        summarise(weight = mean(weight)) |>
        arrange(weight) |>
        pull(species)


    vegetation_composition_community |>
        filter(species %in% common_sps) |>
        mutate(species = factor(species, levels = sp_order)) |>
        ggplot(aes(
            x = {{ group }},
            y = species,
            fill = prop_presences
        )) +
        geom_raster() +
        geom_text(
            aes(
                label =
                    ifelse(
                        prop_presences > 0.1,
                        percent(prop_presences, accuracy = 1),
                        NA
                    )
            ),
            color = "white",
            size = 3
        ) +
        facet_grid(
            cols = vars(water_system),
            scales = "free_x",
            space = "free_x"
        ) +
        scale_fill_viridis_c(guide = "none") +
        scale_y_discrete(
            labels = \(x) str_replace(x, "\\.", " "),
            expand = expansion(0)
        ) +
        scale_x_discrete(expand = expansion(0)) +
        labs(
            y = NULL,
            x = "Saltmarsh community"
        ) +
        theme(
            text = element_text(size = 10),
            axis.text.y = element_text(face = "italic")
        )
}


# 1.2 Candidate models performance ----------------------------------------

plot_candidate_mod_performances <- function(data, model_names, max_x = 10) {
    data |>
        # Now turn every fold's kappa into an individual row for plotting
        unnest(kappa) |>
        filter(system != "Arade") |>
        separate(
            wf_id,
            into = c("preprocess", "model"),
            sep = "_",
            extra = "merge",
            remove = FALSE
        ) |>
        mutate(
            wf_id = fct_reorder(wf_id, kappa, .desc = TRUE, .na_rm = TRUE)
        ) |>
        ggplot(
            aes(
                y = kappa,
                x = wf_id,
                color = preprocess
            )
        ) +
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
            rows = vars(system)
        ) +
        scale_x_discrete(
            limits = \(x) x[1:max_x],
            labels = \(x) model_names[str_split_fixed(x, "_", 2)[, 2]]
        ) +
        scale_color_discrete(
            name = "Predictor transformation",
            labels = c("raw" = "Normalized", "pca" = "PCA")
        ) +
        labs(
            x = "Classifier",
            y = "Kappa"
        ) +
        theme(
            text = element_text(size = 10),
            legend.position = "bottom",
            legend.title.position = "top",
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        )
}
