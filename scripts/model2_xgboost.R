#!/usr/bin/env Rscript
# Model 3: Gradient Boosted Decision Trees (XGBoost)
#
# Run from the top level of the repo, after scripts/cv_and_preprocess.R:
#
#     Rscript scripts/model3_xgboost.R
#
# Writes:
#   output/model3_xgb_results.rds          full CV metrics, per-fold tuning
#                                          summaries, OOF predictions, and refits
#   data/processed/model3_xgb_results.rds  copy in data/processed/

suppressPackageStartupMessages({
  library(tidyverse)
  library(xgboost)
  library(pROC)
})

source("scripts/cv_and_preprocess.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)
cv_folds         <- readRDS(file.path("data", "processed", "cv_folds.rds"))
cv_fold_summaries <- read_csv(file.path("data", "processed", "cv_fold_summaries.csv"),
                              show_col_types = FALSE)

# ---- 1. Native-NA Preprocessing & Metric Helpers -----------------------------

calc_log_loss <- function(y_true, p_hat, eps = 1e-15) {
  p_clip <- pmin(pmax(p_hat, eps), 1 - eps)
  -mean(y_true * log(p_clip) + (1 - y_true) * log(1 - p_clip))
}

calc_brier <- function(y_true, p_hat) {
  mean((y_true - p_hat)^2)
}

calc_auroc <- function(y_true, p_hat) {
  as.numeric(pROC::auc(pROC::roc(y_true, p_hat, levels = c(0, 1),
                                 direction = "<", quiet = TRUE)))
}

#' Prepare raw (un-imputed, unscaled) design matrices for XGBoost
#' Preserves NA values so XGBoost learns default branch directions natively.
prep_xgb_fold_data <- function(train_sub, val_sub, window = c("24h", "48h")) {
  window    <- match.arg(window)
  pred_cols <- if (window == "24h") predictors_24h else predictors_48h
  num_cols  <- setdiff(pred_cols, calendar_cols)
  
  # Drop any numeric column that is 100% NA or zero-variance in train_sub
  sds      <- map_dbl(num_cols, ~ sd(train_sub[[.x]], na.rm = TRUE))
  keep_num <- num_cols[!is.na(sds) & sds > 0]
  keep_all <- c(calendar_cols, keep_num)
  
  # Build model.matrix preserving NAs (na.action = "na.pass")
  rhs_formula <- reformulate(keep_all)
  old_na      <- getOption("na.action")
  options(na.action = "na.pass")
  on.exit(options(na.action = old_na), add = TRUE)
  
  X_train <- model.matrix(rhs_formula, data = train_sub)[, -1, drop = FALSE]
  X_val   <- model.matrix(rhs_formula, data = val_sub)[, -1, drop = FALSE]
  
  storage.mode(X_train) <- "double"
  storage.mode(X_val)   <- "double"
  
  list(
    window     = window,
    predictors = keep_all,
    train_df   = train_sub,
    val_df     = val_sub,
    X_train    = X_train,
    y_train    = as.numeric(train_sub$y_exceed),
    X_val      = X_val,
    y_val      = as.numeric(val_sub$y_exceed)
  )
}

# Version-safe helper to extract best iteration & inner CV log-loss from xgb.cv
extract_xgb_cv_best <- function(cv_obj) {
  elog <- cv_obj$evaluation_log
  if (is.null(elog)) elog <- attr(cv_obj, "evaluation_log")
  
  # Identify the test/validation logloss mean column across xgboost v1.7 / v2 / v3
  col_idx <- grep("^(test|valid|eval).*logloss.*mean$", names(elog), value = TRUE)[1]
  losses  <- as.numeric(elog[[col_idx]])
  
  best_iter <- cv_obj$best_iteration
  if (is.null(best_iter) || length(best_iter) == 0L || is.na(best_iter)) {
    best_iter <- attr(cv_obj, "best_iteration")
  }
  if (is.null(best_iter) || length(best_iter) == 0L || is.na(best_iter)) {
    best_iter <- which.min(losses)
  }
  
  list(
    best_nrounds = as.integer(best_iter),
    min_logloss  = losses[best_iter],
    eval_log     = as_tibble(elog)
  )
}

# ---- 2. Inner 5-Fold Temporal CV Tuning & Early Stopping ---------------------
# Hyperparameter search grid & fixed parameters (Checklist 2.6):
#   - max_depth (interaction depth d): {2, 3, 4, 6}
#   - nrounds (number of trees B): up to 1500, chosen via early stopping (50 rounds)
#   - learning_rate (eta): 0.05
#   - min_child_weight (minimum leaf Hessian sum): 5
#   - subsample (row subsampling): 0.80
#   - colsample_bytree (predictor subsampling): 0.80

fit_nested_xgb <- function(X_train, y_train, train_dates,
                           depth_grid = c(2L, 3L, 4L, 6L),
                           eta = 0.05,
                           min_child_weight = 5,
                           subsample = 0.8,
                           colsample_bytree = 0.8,
                           max_nrounds = 1500L,
                           early_stopping_rounds = 50L,
                           n_inner_folds = 5L) {
  
  # Construct 5 contiguous temporal blocks inside the training fold so early
  # stopping and depth selection never see the outer held-out fold or leak adjacent days
  unique_d   <- sort(unique(train_dates))
  d_map      <- tibble(date = unique_d, inner_fold = ntile(seq_along(unique_d), n_inner_folds))
  inner_id   <- tibble(date = train_dates) |>
    left_join(d_map, by = "date") |>
    pull(inner_fold)
  inner_list <- split(seq_along(y_train), inner_id)
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train, missing = NA)
  n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  
  best_logloss <- Inf
  best_depth   <- NA_integer_
  best_nrounds <- NA_integer_
  best_elog    <- NULL
  grid_history <- list()
  
  for (d in depth_grid) {
    params <- list(
      objective        = "binary:logistic",
      eval_metric      = "logloss",
      tree_method      = "hist",
      eta              = eta,
      max_depth        = d,
      min_child_weight = min_child_weight,
      subsample        = subsample,
      colsample_bytree = colsample_bytree,
      nthread          = n_cores
    )
    
    cv_fit <- xgb.cv(
      params                = params,
      data                  = dtrain,
      nrounds               = max_nrounds,
      folds                 = inner_list,
      early_stopping_rounds = early_stopping_rounds,
      verbose               = 0
    )
    
    res <- extract_xgb_cv_best(cv_fit)
    grid_history[[length(grid_history) + 1]] <- tibble(
      max_depth       = d,
      optimal_nrounds = res$best_nrounds,
      inner_deviance  = res$min_logloss
    )
    
    if (res$min_logloss < best_logloss) {
      best_logloss <- res$min_logloss
      best_depth   <- d
      best_nrounds <- res$best_nrounds
      best_elog    <- res$eval_log
    }
  }
  
  # Refit final model on the full outer training fold at (best_depth, best_nrounds)
  final_params <- list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    tree_method      = "hist",
    eta              = eta,
    max_depth        = best_depth,
    min_child_weight = min_child_weight,
    subsample        = subsample,
    colsample_bytree = colsample_bytree,
    nthread          = n_cores
  )
  
  fitted_bst <- xgb.train(
    params  = final_params,
    data    = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  list(
    model          = fitted_bst,
    best_depth     = best_depth,
    best_nrounds   = best_nrounds,
    inner_deviance = best_logloss,
    eval_log       = best_elog,
    grid_results   = bind_rows(grid_history)
  )
}

# ---- 3. Run Outer Cross-Validation Across All 4 Schemes & Both Windows -------

set.seed(707)

windows <- c("24h", "48h")
schemes <- names(cv_folds)

fold_metrics_list <- list()
oof_preds_list    <- list()

message("Running XGBoost across 4 CV schemes and 2 feature windows (24h, 48h)...")

for (w in windows) {
  for (s in schemes) {
    folds <- cv_folds[[s]]
    message(sprintf("  -> Window: %s | Scheme: %s (%d folds)", w, s, length(folds)))
    
    for (f_idx in seq_along(folds)) {
      f_info    <- folds[[f_idx]]
      train_sub <- train_df[f_info$train_idx, ]
      val_sub   <- train_df[f_info$val_idx, ]
      
      # 1. Prepare matrices with native NA handling
      prepped <- prep_xgb_fold_data(train_sub, val_sub, window = w)
      
      # 2. Tune tree depth and nrounds via inner 5-fold temporal CV
      tuned <- fit_nested_xgb(
        X_train     = prepped$X_train,
        y_train     = prepped$y_train,
        train_dates = prepped$train_df$date
      )
      
      # 3. Predict on outer held-out fold
      dval      <- xgb.DMatrix(data = prepped$X_val, missing = NA)
      p_val_xgb <- as.numeric(predict(tuned$model, dval))
      
      # Persistence baseline on the exact same validation fold
      base_obs <- !is.na(prepped$val_df$exceed_today)
      p_base   <- as.numeric(prepped$val_df$exceed_today[base_obs])
      y_base   <- prepped$y_val[base_obs]
      
      fold_metrics_list[[length(fold_metrics_list) + 1]] <- tibble(
        window          = w,
        scheme          = s,
        cv_fold         = f_idx,
        fold_id         = f_info$fold_id,
        n_train         = nrow(prepped$X_train),
        n_val           = nrow(prepped$X_val),
        optimal_depth   = tuned$best_depth,
        optimal_nrounds = tuned$best_nrounds,
        inner_deviance  = tuned$inner_deviance,
        auroc           = calc_auroc(prepped$y_val, p_val_xgb),
        logloss         = calc_log_loss(prepped$y_val, p_val_xgb),
        brier           = calc_brier(prepped$y_val, p_val_xgb),
        auroc_persist   = calc_auroc(y_base, p_base),
        acc_xgb         = mean((p_val_xgb >= 0.5) == prepped$y_val),
        acc_persist     = mean(p_base == y_base)
      )
      
      oof_preds_list[[length(oof_preds_list) + 1]] <- tibble(
        window       = w,
        scheme       = s,
        cv_fold      = f_idx,
        fold_id      = f_info$fold_id,
        train_row_id = f_info$val_idx,
        station      = prepped$val_df$station,
        date         = prepped$val_df$date,
        target_date  = prepped$val_df$target_date,
        y_exceed     = prepped$y_val,
        exceed_today = as.integer(prepped$val_df$exceed_today),
        p_xgb        = p_val_xgb
      )
    }
  }
}

fold_metrics <- bind_rows(fold_metrics_list) |>
  left_join(
    cv_fold_summaries |> select(scheme, cv_fold, descriptor),
    by = c("scheme", "cv_fold")
  )

oof_preds <- bind_rows(oof_preds_list)

scheme_summary <- fold_metrics |>
  group_by(window, scheme) |>
  summarise(
    n_folds        = n(),
    auroc_mean     = mean(auroc),
    auroc_sd       = sd(auroc),
    logloss_mean   = mean(logloss),
    logloss_sd     = sd(logloss),
    brier_mean     = mean(brier),
    brier_sd       = sd(brier),
    persist_auroc  = mean(auroc_persist),
    median_depth   = median(optimal_depth),
    mean_nrounds   = mean(optimal_nrounds),
    .groups        = "drop"
  )

# ---- 4. Refit Final 24h & 48h XGBoost Models on All Training Data ------------

refit_full_xgb <- function(window_choice) {
  prepped_full <- prep_xgb_fold_data(train_df, test_df, window = window_choice)
  tuned_full   <- fit_nested_xgb(
    X_train     = prepped_full$X_train,
    y_train     = prepped_full$y_train,
    train_dates = prepped_full$train_df$date
  )
  imp_table <- as_tibble(xgb.importance(
    feature_names = colnames(prepped_full$X_train),
    model         = tuned_full$model
  ))
  
  list(
    window       = window_choice,
    tuned        = tuned_full,
    importance   = imp_table,
    feature_cols = colnames(prepped_full$X_train)
  )
}

refit_24h <- refit_full_xgb("24h")
refit_48h <- refit_full_xgb("48h")

# ---- 5. Save Results ---------------------------------------------------------

xgb_results <- list(
  fold_metrics   = fold_metrics,
  scheme_summary = scheme_summary,
  oof_preds      = oof_preds,
  refit_24h      = refit_24h,
  refit_48h      = refit_48h
)

saveRDS(xgb_results, file.path("output", "model3_xgb_results.rds"))
saveRDS(xgb_results, file.path("data", "processed", "model3_xgb_results.rds"))

message("Saved XGBoost results to output/model3_xgb_results.rds")
print(scheme_summary |> select(window, scheme, n_folds, auroc_mean, auroc_sd,
                               brier_mean, persist_auroc, median_depth, mean_nrounds))