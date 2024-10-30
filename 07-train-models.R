###################################################
##    Train saltmarsh classification models      ##
##          Márcio Martins    2023/08/09         ##
##  email     marciomartinsred@gmail.com         ##
## Github     MarcioFCMartins                    ##
###################################################


# Setup -------------------------------------------------------------------
libs <- c(
    "tidymodels", # Wrapper for several machine learning packages
    "sf", # Simple-features for vector spatial data
    "here", # Path normalization
    "doParallel", # Paralelization methods, used in backend for model training
    "tidyverse" # Data wrangling
)

# Install any missing libraries
not_installed <- !libs %in% installed.packages()
if (any(not_installed)) {
    install.packages(libs[not_installed])
}
# Load all libraries
invisible(lapply(libs, library, character.only = T, warn.conflicts = FALSE))

here::i_am("code/08-model-training.R")

# Setup parallel processing
all_cores <- parallel::detectCores(logical = FALSE)
all_cores <- ifelse(all_cores > 12L, 12L, all_cores)
cl <- makePSOCKcluster(all_cores - 1L)
registerDoParallel(cl)


# Dataset -----------------------------------------------------------------
# Data sets
full_data_sf <- st_read(
    here("data/clean/mapping/all-systems/saltmarsh-training-data.gpkg"),
    quiet = TRUE
) |>
    mutate(cluster = as.factor(cluster))


train_data <- full_data_sf |>
    filter(train == TRUE) |>
    st_drop_geometry()


# Define model variables --------------------------------------------------
# Variables to predict
pred_vars <- c("marsh_type")
# Predictors
train_features <- c(
    "coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge",
    "nir", "ndwi_low", "ndwi_high", "ndvi"
)

# Loop over possible predictors
for (pred_var in pred_vars) {
    # Prepare a model recipe to pre-process data
    recipe_formula <- as.formula(
        paste(
            pred_var,
            "~",
            paste(train_features, collapse = " + ")
        )
    )
    # Model with no pre-processing
    sm_rec <- recipe(
        recipe_formula,
        data = train_data
    ) |>
        step_normalize(all_predictors())

    # Model based on principal components
    sm_rec_pca <-
        recipe(
            recipe_formula,
            data = train_data
        ) |>
        step_normalize(all_predictors()) |>
        step_pca(all_predictors(), threshold = 0.98)


    ## Resampling folds and metrics --------------------------------------------
    # Create folds for cross-validation
    set.seed(99)
    folds <- vfold_cv(train_data, strata = !!pred_var, v = 10)

    # Evaluation metrics that I want to compare across models
    metrics <- metric_set(kap, accuracy)


    ## Models to compare -------------------------------------------------------
    # https://www.tmwr.org/workflow-sets - suggested method to compare models
    lasso_spec <- multinom_reg(
        penalty = tune(),
        mixture = tune()
    ) %>%
        set_engine("glmnet") |>
        set_mode("classification")

    rforest_spec <- rand_forest(
        mtry = tune(),
        min_n = tune(),
        trees = 500
    ) |>
        set_mode("classification") |>
        set_engine("ranger")

    svm_r_spec <- svm_rbf(
        cost = tune(),
        rbf_sigma = tune()
    ) |>
        set_engine("kernlab") |>
        set_mode("classification")

    svm_p_spec <- svm_poly(
        cost = tune(),
        degree = tune()
    ) |>
        set_engine("kernlab") |>
        set_mode("classification")

    knn_spec <- nearest_neighbor(
        neighbors = tune(),
        dist_power = tune(),
        weight_func = tune()
    ) |>
        set_engine("kknn") |>
        set_mode("classification")

    nnet_spec <- mlp(
        hidden_units = tune(),
        penalty = tune(),
        epochs = tune()
    ) |>
        set_mode("classification") |>
        set_engine("nnet")

    xgb_spec <- boost_tree(
        trees = tune(),
        min_n = tune(),
        mtry = tune(),
        learn_rate = 0.01
    ) %>%
        set_engine("xgboost") %>%
        set_mode("classification")


    ## Hyperparameter optimization - models for all systems --------------------
    model_wfs <- workflow_set(
        preproc = list(
            "raw" = sm_rec,
            "pca" = sm_rec_pca
        ),
        models = list(
            "lasso" = lasso_spec,
            "rf" = rforest_spec,
            "svm_rad" = svm_r_spec,
            "svm_pol" = svm_p_spec,
            "knn" = knn_spec,
            "nn" = nnet_spec,
            "xgb" = xgb_spec
        ),
        cross = TRUE
    )

    message("Training universal model for ", pred_var)
    # Tune the combinations of recipes and models
    model_tuning <- model_wfs |>
        workflow_map(
            seed = 99,
            resamples = folds,
            metrics = metrics,
            grid = 25,
            control = tune::control_grid(
                parallel_over = "everything",
                save_pred = TRUE,
                save_workflow = TRUE
            )
        )


    saveRDS(
        model_tuning,
        here(
            paste0("outputs/models/universal/", pred_var, "-crosscv.RDS")
        )
    )
    saveRDS(
        folds,
        here(
            paste0("outputs/models/universal/", pred_var, "-crosscv-data.RDS")
        )
    )


    ## Hyperparameter optimization - 1 model per system ------------------------

    region_names <- c(
        "Alvor",
        "Arade",
        "Ria Formosa",
        "Guadiana"
    )

    for (system in region_names) {
        system_train_data <- full_data_sf |>
            filter(water_system == system, train == TRUE) |>
            st_drop_geometry()


        # Model with no pre-processing
        sm_rec <- recipe(
            recipe_formula,
            data = system_train_data
        ) |>
            step_normalize(all_predictors())

        # Model based on principal components
        sm_rec_pca <- recipe(
            recipe_formula,
            data = system_train_data
        ) |>
            step_normalize(all_predictors()) |>
            step_pca(all_predictors(), threshold = 0.98)

        # Create folds for cross-validation based on spatial clusters
        set.seed(99)
        folds <- vfold_cv(system_train_data, strata = !!pred_var, v = 10)

        model_wfs <- workflow_set(
            preproc = list(
                "raw" = sm_rec,
                "pca" = sm_rec_pca
            ),
            models = list(
                "lasso" = lasso_spec,
                "rf" = rforest_spec,
                "svm_rad" = svm_r_spec,
                "svm_pol" = svm_p_spec,
                "knn" = knn_spec,
                "nn" = nnet_spec,
                "xgb" = xgb_spec
            ),
            cross = TRUE
        )

        message("Training model for ", system, ", response: ", pred_var)
        model_tuning <- model_wfs |>
            workflow_map(
                seed = 99,
                resamples = folds,
                metrics = metrics,
                grid = 25,
                control = tune::control_grid(
                    parallel_over = "everything",
                    save_pred = TRUE,
                    save_workflow = TRUE
                )
            )

        saveRDS(
            model_tuning,
            here(
                "outputs/models/per-system",
                paste0(system, "-", pred_var, "-crosscv.RDS")
            )
        )
        saveRDS(
            folds,
            here(
                "outputs/models/per-system",
                paste0(system, "-", pred_var, "-crosscv-data.RDS")
            )
        )

        rm(model_tuning)
        gc()
    }
}

stopCluster(cl)
