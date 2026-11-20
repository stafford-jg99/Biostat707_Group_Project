#!/usr/bin/env Rscript
# Model 1: Regularized Logistic Regression (Elastic Net)
#
# Run from the top level of the repo, after scripts/cv_and_preprocess.R:
#
#     pixi run Rscript scripts/model1_elastic_net.R
#
# Writes:
#   output/model1_enet_results.rds          full CV metrics, fold coefficients,
#                                           out-of-fold predictions, and refits
#   data/processed/model1_enet_results.rds  copy in data/processed/ for convenience

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(pROC)
})

source("scripts/cv_and_preprocess.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)
cv_folds <- readRDS(file.path("data", "processed", "cv_folds.rds"))
catalog  <- read_csv(file.path("data", "processed", "predictor_catalog.csv"),
                     show_col_types = FALSE)

# ---- 1. Helper Functions for Metrics & Nested CV Tuning ----------------------

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

#' Fit Elastic Net with nested temporal blocked CV inside a training fold
#'
#' @param X_train Standardized training design matrix
#' @param y_train Binary 0/1 outcome vector
#' @param train_dates Date vector for constructing contiguous inner temporal folds
#' @param alpha_grid Grid of Elastic Net mixing parameters (0 = Ridge, 1 = LASSO)
#' @param n_inner_folds Number of contiguous temporal blocks for inner CV
fit_nested_enet <- function(X_train, y_train, train_dates,
                            alpha_grid = c(0.25, 0.50, 0.75, 1.00),
                            n_inner_folds = 5) {
  # Build contiguous temporal fold IDs inside the training fold so inner
  # hyperparameter tuning does not leak adjacent days across inner splits
  unique_d <- sort(unique(train_dates))
  d_map    <- tibble(date = unique_d, inner_fold = ntile(seq_along(unique_d), n_inner_folds))
  foldid   <- tibble(date = train_dates) |>
    left_join(d_map, by = "date") |>
    pull(inner_fold)
  
  # Tune alpha and lambda using binomial deviance (2 * log loss)
  best_fit  <- NULL
  best_cvm  <- Inf
  best_a    <- NA_real_
  
  for (a in alpha_grid) {
    cv_fit <- cv.glmnet(
      x            = X_train,
      y            = y_train,
      family       = "binomial",
      alpha        = a,
      foldid       = foldid,
      type.measure = "deviance", # Minimizes log loss (Checklist 2.3)
      standardize  = FALSE       # Already standardized inside prep_fold_data()
    )
    min_cvm <- min(cv_fit$cvm, na.rm = TRUE)
    if (min_cvm < best_cvm) {
      best_cvm <- min_cvm
      best_fit <- cv_fit
      best_a   <- a
    }
  }
  
  list(
    cv_fit     = best_fit,
    alpha      = best_a,
    lambda_min = best_fit$lambda.min,
    lambda_1se = best_fit$lambda.1se,
    min_dev    = best_cvm
  )
}

extract_enet_coefs <- function(cv_fit, s = "lambda.min") {
  raw_c <- as.matrix(coef(cv_fit, s = s))
  tibble(
    term = rownames(raw_c),
    beta = as.numeric(raw_c[, 1])
  ) |>
    filter(term != "(Intercept)") |>
    mutate(
      kept    = abs(beta) > 1e-8,
      std_or  = exp(beta)
    )
}

# ---- 2. Run Outer Cross-Validation Across All 4 Schemes & Both Windows -------

set.seed(707)

windows <- c("24h", "48h")
schemes <- names(cv_folds)

fold_metrics_list <- list()
fold_coefs_list   <- list()
oof_preds_list    <- list()

message("Running Elastic Net across 4 CV schemes and 2 feature windows (24h, 48h)...")

for (w in windows) {
  for (s in schemes) {
    folds <- cv_folds[[s]]
    
    for (f_idx in seq_along(folds)) {
      f_info    <- folds[[f_idx]]
      train_sub <- train_df[f_info$train_idx, ]
      val_sub   <- train_df[f_info$val_idx, ]
      
      # 1. Fold-safe imputation and standardization (Checklist 2.3)
      prepped <- prep_fold_data(train_sub, val_sub, window = w)
      
      # 2. Inner CV tuning on training fold only
      tuned <- fit_nested_enet(
        X_train     = prepped$X_train,
        y_train     = prepped$y_train,
        train_dates = prepped$train_df$date
      )
      
      # 3. Predict on held-out validation fold (using lambda.min & lambda.1se)
      p_val_min <- as.numeric(predict(
        tuned$cv_fit, newx = prepped$X_val, s = "lambda.min", type = "response"
      ))
      p_val_1se <- as.numeric(predict(
        tuned$cv_fit, newx = prepped$X_val, s = "lambda.1se", type = "response"
      ))
      
      # Persistence baseline on the exact same validation fold (where observed)
      base_obs <- !is.na(prepped$val_df$exceed_today)
      p_base   <- as.numeric(prepped$val_df$exceed_today[base_obs])
      y_base   <- prepped$y_val[base_obs]
      
      # Extract coefficients for predictor stability tracking (Checklist 2.4)
      coefs_min <- extract_enet_coefs(tuned$cv_fit, s = "lambda.min") |>
        mutate(window = w, scheme = s, fold_id = f_info$fold_id, rule = "lambda.min")
      coefs_1se <- extract_enet_coefs(tuned$cv_fit, s = "lambda.1se") |>
        mutate(window = w, scheme = s, fold_id = f_info$fold_id, rule = "lambda.1se")
      
      fold_coefs_list[[length(fold_coefs_list) + 1]] <- bind_rows(coefs_min, coefs_1se)
      
      # Record fold performance metrics
      fold_metrics_list[[length(fold_metrics_list) + 1]] <- tibble(
        window         = w,
        scheme         = s,
        fold_num       = f_idx,
        fold_id        = f_info$fold_id,
        n_train        = nrow(prepped$X_train),
        n_val          = nrow(prepped$X_val),
        alpha          = tuned$alpha,
        lambda_min     = tuned$lambda_min,
        lambda_1se     = tuned$lambda_1se,
        n_total_terms  = nrow(coefs_min),
        n_kept_min     = sum(coefs_min$kept),
        n_kept_1se     = sum(coefs_1se$kept),
        auroc_min      = calc_auroc(prepped$y_val, p_val_min),
        logloss_min    = calc_log_loss(prepped$y_val, p_val_min),
        brier_min      = calc_brier(prepped$y_val, p_val_min),
        auroc_1se      = calc_auroc(prepped$y_val, p_val_1se),
        logloss_1se    = calc_log_loss(prepped$y_val, p_val_1se),
        brier_1se      = calc_brier(prepped$y_val, p_val_1se),
        auroc_persist  = calc_auroc(y_base, p_base),
        acc_enet_min   = mean((p_val_min >= 0.5) == prepped$y_val),
        acc_persist    = mean(p_base == y_base)
      )
      
      # Store out-of-fold predictions
      oof_preds_list[[length(oof_preds_list) + 1]] <- tibble(
        window       = w,
        scheme       = s,
        fold_id      = f_info$fold_id,
        train_row_id = f_info$val_idx,
        station      = prepped$val_df$station,
        date         = prepped$val_df$date,
        target_date  = prepped$val_df$target_date,
        y_exceed     = prepped$y_val,
        exceed_today = as.integer(prepped$val_df$exceed_today),
        p_enet_min   = p_val_min,
        p_enet_1se   = p_val_1se
      )
    }
  }
}

fold_metrics <- bind_rows(fold_metrics_list)
fold_coefs   <- bind_rows(fold_coefs_list)
oof_preds    <- bind_rows(oof_preds_list)

scheme_summary <- fold_metrics |>
  group_by(window, scheme) |>
  summarise(
    n_folds          = n(),
    auroc_mean       = mean(auroc_min),
    auroc_sd         = sd(auroc_min),
    logloss_mean     = mean(logloss_min),
    logloss_sd       = sd(logloss_min),
    brier_mean       = mean(brier_min),
    brier_sd         = sd(brier_min),
    auroc_1se_mean   = mean(auroc_1se),
    auroc_1se_sd     = sd(auroc_1se),
    persist_auroc    = mean(auroc_persist),
    mean_kept_min    = mean(n_kept_min),
    mean_kept_1se    = mean(n_kept_1se),
    .groups          = "drop"
  )

# ---- 3. Refit Elastic Net on All Training Data (2013-2015) -------------------
# Fits on all of train_df using the same common-scale standardization (Checklist 2.4)

refit_full_train <- function(window_choice) {
  prepped_full <- prep_fold_data(train_df, test_df, window = window_choice)
  tuned_full   <- fit_nested_enet(
    X_train     = prepped_full$X_train,
    y_train     = prepped_full$y_train,
    train_dates = prepped_full$train_df$date
  )
  coefs_min <- extract_enet_coefs(tuned_full$cv_fit, s = "lambda.min")
  coefs_1se <- extract_enet_coefs(tuned_full$cv_fit, s = "lambda.1se")
  
  list(
    window       = window_choice,
    tuned        = tuned_full,
    coefs_min    = coefs_min,
    coefs_1se    = coefs_1se,
    scale_params = prepped_full$scale_params,
    imputation   = prepped_full$imputation
  )
}

refit_24h <- refit_full_train("24h")
refit_48h <- refit_full_train("48h")

# ---- 4. Build Table 2.3 Data (Selection Stability & Standardized Coefs) ------

build_stability_table <- function(window_choice, refit_obj, rule_choice = "lambda.min") {
  refit_c <- if (rule_choice == "lambda.min") refit_obj$coefs_min else refit_obj$coefs_1se
  
  # Summarize selection frequency across the 4 CV schemes
  stab_by_scheme <- fold_coefs |>
    filter(window == window_choice, rule == rule_choice) |>
    group_by(scheme, term) |>
    summarise(
      n_kept  = sum(kept),
      n_folds = n(),
      pct     = 100 * mean(kept),
      label   = sprintf("%d/%d (%.0f%%)", sum(kept), n(), 100 * mean(kept)),
      .groups = "drop"
    ) |>
    select(scheme, term, label) |>
    pivot_wider(names_from = scheme, values_from = label)
  
  # Summarize fold coefficient range in the primary temporal schemes
  beta_range <- fold_coefs |>
    filter(window == window_choice, rule == rule_choice,
           scheme == "blocked_gap_6fold") |>
    group_by(term) |>
    summarise(
      fold_beta_mean = mean(beta),
      fold_beta_min  = min(beta),
      fold_beta_max  = max(beta),
      sign_consistent = (all(beta >= 0) | all(beta <= 0)),
      .groups = "drop"
    )
  
  # Map dummy-encoded factor terms back to base variable names for catalog join
  refit_c |>
    mutate(
      base_var = case_when(
        str_starts(term, "season") ~ "season",
        str_starts(term, "dow")    ~ "dow",
        TRUE                       ~ term
      )
    ) |>
    left_join(catalog |> distinct(variable, group), by = c("base_var" = "variable")) |>
    left_join(stab_by_scheme, by = "term") |>
    left_join(beta_range, by = "term") |>
    arrange(desc(kept), desc(abs(beta)))
}

table_2_3_24h <- build_stability_table("24h", refit_24h, "lambda.min")
table_2_3_48h <- build_stability_table("48h", refit_48h, "lambda.min")

# ---- 5. Save Results ---------------------------------------------------------

enet_results <- list(
  fold_metrics   = fold_metrics,
  scheme_summary = scheme_summary,
  oof_preds      = oof_preds,
  fold_coefs     = fold_coefs,
  refit_24h      = refit_24h,
  refit_48h      = refit_48h,
  table_2_3_24h  = table_2_3_24h,
  table_2_3_48h  = table_2_3_48h
)

saveRDS(enet_results, file.path("output", "model1_enet_results.rds"))
saveRDS(enet_results, file.path("data", "processed", "model1_enet_results.rds"))

message("Saved Elastic Net results to output/model1_enet_results.rds")
print(scheme_summary |> select(window, scheme, n_folds, auroc_mean, auroc_sd,
                               brier_mean, persist_auroc, mean_kept_min))