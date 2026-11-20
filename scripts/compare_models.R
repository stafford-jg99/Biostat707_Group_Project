#!/usr/bin/env Rscript
# Section 2.7: Comparison of Models 1 (Elastic Net), 2 (XGBoost), and Persistence Baseline
#
# Run from the top level of the repo, after scripts/model1_elastic_net.R and scripts/model3_xgboost.R:
#
#     Rscript scripts/compare_models.R
#
# Writes:
#   output/oof_predictions.csv              combined out-of-fold predictions for all models
#   output/model_comparison_summary.rds     per-fold and summary comparison tables + ROC objects

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
})

enet_path <- file.path("output", "model1_enet_results.rds")
xgb_path  <- file.path("output", "model3_xgb_results.rds")
sum_path  <- file.path("data", "processed", "cv_fold_summaries.csv")

stopifnot(
  "Run scripts/model1_elastic_net.R first" = file.exists(enet_path),
  "Run scripts/model3_xgboost.R first"     = file.exists(xgb_path),
  "Run scripts/cv_and_preprocess.R first"  = file.exists(sum_path)
)

enet_res          <- readRDS(enet_path)
xgb_res           <- readRDS(xgb_path)
cv_fold_summaries <- read_csv(sum_path, show_col_types = FALSE)

# ---- 1. Combine Out-of-Fold Predictions Across Models ------------------------

oof_combined <- enet_res$oof_preds |>
  left_join(
    xgb_res$oof_preds |>
      select(window, scheme, fold_id, train_row_id, cv_fold, p_xgb),
    by = c("window", "scheme", "fold_id", "train_row_id")
  ) |>
  relocate(cv_fold, .after = scheme) |>
  arrange(window, scheme, cv_fold, train_row_id)

write_csv(oof_combined, file.path("output", "oof_predictions.csv"))
write_csv(oof_combined, file.path("data", "processed", "oof_predictions.csv"))
message(sprintf("Wrote %d combined out-of-fold prediction rows to output/oof_predictions.csv",
                nrow(oof_combined)))

# ---- 2. Build Per-Fold AUROC Comparison Table Across Schemes -----------------

# Merge per-fold AUROC from Elastic Net (24h & 48h), XGBoost (24h & 48h), and Persistence
enet_folds <- enet_res$fold_metrics |>
  select(window, scheme, cv_fold = fold_num, fold_id,
         auroc_enet = auroc_min, logloss_enet = logloss_min, brier_enet = brier_min,
         auroc_persist)

xgb_folds <- xgb_res$fold_metrics |>
  select(window, scheme, cv_fold, fold_id,
         auroc_xgb = auroc, logloss_xgb = logloss, brier_xgb = brier)

# Count excluded missing benchmark rows per fold (Week 8 slide 7 requirement)
persist_missing <- oof_combined |>
  filter(window == "24h") |>
  group_by(scheme, cv_fold, fold_id) |>
  summarise(
    n_val              = n(),
    n_persist_obs      = sum(!is.na(exceed_today)),
    n_persist_excluded = sum(is.na(exceed_today)),
    .groups            = "drop"
  )

fold_comparison_long <- enet_folds |>
  left_join(xgb_folds, by = c("window", "scheme", "cv_fold", "fold_id"))

# Wide per-fold table comparing 24h and 48h for both models + Persistence on each fold
fold_comparison_wide <- fold_comparison_long |>
  select(window, scheme, cv_fold, fold_id, auroc_enet, auroc_xgb, auroc_persist) |>
  pivot_wider(
    id_cols     = c(scheme, cv_fold, fold_id, auroc_persist),
    names_from  = window,
    values_from = c(auroc_enet, auroc_xgb)
  ) |>
  left_join(
    cv_fold_summaries |> select(scheme, cv_fold, descriptor, val_n, val_exceed_pct),
    by = c("scheme", "cv_fold")
  ) |>
  left_join(
    persist_missing |> select(scheme, cv_fold, n_persist_excluded),
    by = c("scheme", "cv_fold")
  ) |>
  select(
    scheme, cv_fold, fold_id, descriptor, val_n, val_exceed_pct, n_persist_excluded,
    auroc_enet_24h, auroc_xgb_24h, auroc_enet_48h, auroc_xgb_48h, auroc_persist
  ) |>
  arrange(scheme, cv_fold)

# Summary across folds by scheme
scheme_comparison_summary <- fold_comparison_wide |>
  group_by(scheme) |>
  summarise(
    n_folds             = n(),
    total_persist_excl  = sum(n_persist_excluded),
    enet_24h_mean       = mean(auroc_enet_24h),
    enet_24h_sd         = sd(auroc_enet_24h),
    xgb_24h_mean        = mean(auroc_xgb_24h),
    xgb_24h_sd          = sd(auroc_xgb_24h),
    enet_48h_mean       = mean(auroc_enet_48h),
    enet_48h_sd         = sd(auroc_enet_48h),
    xgb_48h_mean        = mean(auroc_xgb_48h),
    xgb_48h_sd          = sd(auroc_xgb_48h),
    persist_mean        = mean(auroc_persist),
    persist_sd          = sd(auroc_persist),
    diff_xgb_vs_enet_24 = mean(auroc_xgb_24h - auroc_enet_24h),
    sd_xgb_vs_enet_24   = sd(auroc_xgb_24h - auroc_enet_24h),
    .groups             = "drop"
  )

# ---- 3. Build Out-of-Fold ROC Curve Coordinates for Each Scheme --------------

build_scheme_roc_df <- function(oof_df, scheme_choice, window_choice = "24h") {
  sub <- oof_df |>
    filter(scheme == scheme_choice, window == window_choice)
  
  roc_enet <- pROC::roc(sub$y_exceed, sub$p_enet_min, levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_xgb  <- pROC::roc(sub$y_exceed, sub$p_xgb,      levels = c(0, 1), direction = "<", quiet = TRUE)
  
  sub_pers <- sub |> filter(!is.na(exceed_today))
  roc_pers <- pROC::roc(sub_pers$y_exceed, sub_pers$exceed_today, levels = c(0, 1), direction = "<", quiet = TRUE)
  
  extract_coords <- function(roc_obj, model_label) {
    tibble(
      scheme      = scheme_choice,
      window      = window_choice,
      model       = model_label,
      pooled_auc  = as.numeric(pROC::auc(roc_obj)),
      fpr         = 1 - roc_obj$specificities,
      tpr         = roc_obj$sensitivities
    ) |>
      arrange(fpr, tpr)
  }
  
  bind_rows(
    extract_coords(roc_enet, "Model 1: Elastic Net"),
    extract_coords(roc_xgb,  "Model 2: XGBoost"),
    extract_coords(roc_pers, "Benchmark: Persistence Baseline")
  )
}

roc_curves_24h <- map_dfr(unique(oof_combined$scheme),
                          ~ build_scheme_roc_df(oof_combined, .x, window_choice = "24h"))

comparison_results <- list(
  fold_comparison_wide      = fold_comparison_wide,
  fold_comparison_long      = fold_comparison_long,
  scheme_comparison_summary = scheme_comparison_summary,
  roc_curves_24h            = roc_curves_24h
)

saveRDS(comparison_results, file.path("output", "model_comparison_summary.rds"))
saveRDS(comparison_results, file.path("data", "processed", "model_comparison_summary.rds"))

message("Saved comparison tables and ROC curves to output/model_comparison_summary.rds")
print(scheme_comparison_summary |>
        select(scheme, enet_24h_mean, enet_24h_sd, xgb_24h_mean, xgb_24h_sd,
               persist_mean, diff_xgb_vs_enet_24))