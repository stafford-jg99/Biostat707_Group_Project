#!/usr/bin/env Rscript
# Cross-validation design and fold-safe preprocessing pipeline.
#
# Run from the top level of the repo, after scripts/build_features.R:
#
#     Rscript scripts/cv_and_preprocess.R
#
# Writes:
#   data/processed/cv_folds.rds           list of CV fold definitions on train
#   data/processed/predictor_catalog.csv  structured table of predictors
#
# Provides helper functions when sourced:
#   prep_fold_data()                      fold-safe median imputation, zero-var
#                                         filtering, scaling, and model matrices

suppressPackageStartupMessages(library(tidyverse))
source("scripts/impute.R")

# ---- 1. Load & Verify Modeling Table -------------------------

in_path  <- file.path("data", "processed", "station_day_features.csv")
out_dir  <- file.path("data", "processed")
stopifnot("Run scripts/build_features.R first" = file.exists(in_path))

station_day <- read_csv(in_path, show_col_types = FALSE) |>
  mutate(
    station = factor(station),
    period  = factor(period, levels = c("train", "test")),
    season  = factor(season, levels = c("spring", "summer", "autumn", "winter")),
    dow     = factor(dow, levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
  ) |>
  arrange(station, date)

# Layout & integrity assertions
stopifnot(
  "Duplicate (station, date) rows detected" =
    anyDuplicated(station_day[, c("station", "date")]) == 0L,
  "Missing outcome values detected" =
    !anyNA(station_day$y_exceed),
  "Target date must be exactly date + 1" =
    all(station_day$target_date == station_day$date + 1),
  "Train/test split must occur strictly at 2016-03-01" =
    all((station_day$target_date < as.Date("2016-03-01")) == (station_day$period == "train"))
)

train_df <- station_day |> filter(period == "train")
test_df  <- station_day |> filter(period == "test")

message(sprintf(
  "Modeling table verified: %d total rows | Train: %d rows (%.1f%% exceed) | Test: %d rows (%.1f%% exceed)",
  nrow(station_day),
  nrow(train_df), 100 * mean(train_df$y_exceed),
  nrow(test_df),  100 * mean(test_df$y_exceed)
))

# ---- 2. Define Predictor Sets & Catalog ----------------------
# Strictly separate IDs, outcomes, benchmark, and 24h vs 48h predictor sets.
# Note: n_hours_24h is always 24 (constant across all rows), so it is excluded
# from the active predictor lists to prevent zero-variance scaling issues.

id_cols        <- c("station", "date", "target_date", "period", "season_year")
outcome_cols   <- c("y_exceed", "y_pm25_next")
benchmark_cols <- c("exceed_today", "PM25_today_mean", "hours_obs_pm25_today")
calendar_cols  <- c("season", "month", "dow")

cols_24h <- names(station_day) |>
  str_subset("_24h$") |>
  setdiff("n_hours_24h")

cols_48h <- names(station_day) |>
  str_subset("_48h$") |>
  setdiff("n_hours_48h")

predictors_24h <- c("PM25_last", calendar_cols, cols_24h)
predictors_48h <- c("PM25_last", calendar_cols, cols_48h)

# Build a structured predictor catalog for the report
predictor_catalog <- tibble(variable = c(predictors_24h, cols_48h)) |>
  mutate(
    window = case_when(
      str_ends(variable, "_24h") ~ "24h window",
      str_ends(variable, "_48h") ~ "48h window",
      TRUE                       ~ "Shared (both windows)"
    ),
    group = case_when(
      variable == "PM25_last"                ~ "Terminal lag reading (23:00)",
      variable %in% calendar_cols            ~ "Calendar & seasonality",
      str_detect(variable, "^(n_obs|n_out)") ~ "Informative missingness / outage",
      str_detect(variable, "^(PM|SO2|NO2|CO|O3).*_mean_") ~ "Pollutant window mean",
      str_detect(variable, "^(PM|SO2|NO2|CO|O3).*_max_")  ~ "Pollutant window peak (max)",
      str_detect(variable, "^(TEMP|PRES|DEWP|WSPM)_mean_") ~ "Meteorological mean",
      str_detect(variable, "^(TEMP|PRES|DEWP|WSPM)_(min|max)_") ~ "Meteorological extremes (min/max)",
      str_detect(variable, "^RAIN_total_")   ~ "Cumulative precipitation",
      str_detect(variable, "^wind_[uv]_")    ~ "Directional wind vector components",
      TRUE                                   ~ "Other"
    ),
    rationale = case_when(
      group == "Terminal lag reading (23:00)"      ~ "Captures immediate pollutant concentration 1 hour prior to target day.",
      group == "Calendar & seasonality"            ~ "Accounts for winter heating emissions and seasonal atmospheric shifts.",
      group == "Informative missingness / outage"  ~ "Flags sensor downtime and reduced observation counts in the window.",
      group == "Pollutant window mean"             ~ "Summarizes baseline aerosol and gas accumulation over the window.",
      group == "Pollutant window peak (max)"       ~ "Captures acute pollution spikes that may persist into the next day.",
      group == "Meteorological mean"               ~ "Reflects stagnation, humidity, and thermal inversion conditions.",
      group == "Meteorological extremes (min/max)" ~ "Captures diurnal weather swings and frontal passages.",
      group == "Cumulative precipitation"          ~ "Measures wet deposition (rain washout) of particulate matter.",
      group == "Directional wind vector components"~ "Separates clean northerly winds from polluted southerly transport."
    )
  )

write_csv(predictor_catalog, file.path(out_dir, "predictor_catalog.csv"))

# ---- 3. Cross-Validation Design on Training Set --------------
# All CV schemes operate exclusively on train_df (target_date < 2016-03-01).
# Row indices below refer to row positions in `train_df` (1 .. nrow(train_df)).

set.seed(707)

# Scheme A: Naive Stratified Random 5-Fold CV (to quantify temporal leakage)
make_random_cv <- function(df, k = 5) {
  # Stratified by y_exceed
  fold_assign <- integer(nrow(df))
  for (val in c(FALSE, TRUE)) {
    idx <- which(df$y_exceed == val)
    fold_assign[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  map(1:k, function(f) {
    val_idx   <- which(fold_assign == f)
    train_idx <- which(fold_assign != f)
    list(
      scheme     = "random_5fold",
      fold_num   = f,
      fold_id    = sprintf("Random_Fold_%d", f),
      descriptor = sprintf(
        "Stratified random 20%% sample across all %d stations (%s to %s)",
        n_distinct(df$station[val_idx]),
        min(df$date[val_idx]),
        max(df$date[val_idx])
      ),
      train_idx  = train_idx,
      val_idx    = val_idx
    )
  })
}

# Scheme B: Blocked Temporal CV with a 7-Day Buffer Gap
# Splits the 3-year training period into 6 contiguous 6-month blocks. When a
# block is held out for validation, a 7-day buffer on either side of the validation
# window is dropped from the training fold so adjacent-day autocorrelation cannot leak.
make_blocked_gap_cv <- function(df, n_blocks = 6, buffer_days = 7) {
  unique_dates <- sort(unique(df$date))
  date_blocks  <- tibble(
    date  = unique_dates,
    block = ntile(seq_along(unique_dates), n_blocks)
  )
  df_blocked <- df |> left_join(date_blocks, by = "date")
  
  map(1:n_blocks, function(b) {
    val_dates <- date_blocks |> filter(block == b) |> pull(date)
    min_val   <- min(val_dates)
    max_val   <- max(val_dates)
    
    val_idx   <- which(df_blocked$block == b)
    train_idx <- which(
      df_blocked$block != b &
        (df_blocked$date < (min_val - buffer_days) | df_blocked$date > (max_val + buffer_days))
    )
    n_buffer  <- nrow(df) - length(val_idx) - length(train_idx)
    
    list(
      scheme     = "blocked_gap_6fold",
      fold_num   = b,
      fold_id    = sprintf("Blocked_Gap_Fold_%d (%s to %s)", b, min_val, max_val),
      descriptor = sprintf(
        "Contiguous 6-mo block: %s to %s (all 12 sites; ±%dd buffer drops %d rows)",
        min_val, max_val, buffer_days, n_buffer
      ),
      train_idx  = train_idx,
      val_idx    = val_idx
    )
  })
}

# Scheme C: Leave-One-Season-Year-Out Temporal CV (3 folds: 2013, 2014, 2015)
# Plus a 7-day buffer at the year boundary to prevent boundary day leakage.
make_year_temporal_cv <- function(df, buffer_days = 7) {
  years <- sort(unique(df$season_year))
  map2(seq_along(years), years, function(f_num, yr) {
    val_idx   <- which(df$season_year == yr)
    min_val   <- min(df$date[val_idx])
    max_val   <- max(df$date[val_idx])
    train_idx <- which(
      df$season_year != yr &
        (df$date < (min_val - buffer_days) | df$date > (max_val + buffer_days))
    )
    train_yrs <- paste(setdiff(years, yr), collapse = " & ")
    n_buffer  <- nrow(df) - length(val_idx) - length(train_idx)
    
    list(
      scheme     = "temporal_year_3fold",
      fold_num   = f_num,
      fold_id    = sprintf("Holdout_Year_%d", yr),
      descriptor = sprintf(
        "Holdout Season-Year %d (%s to %s; trained on %s, %d buffer rows dropped)",
        yr, min_val, max_val, train_yrs, n_buffer
      ),
      train_idx  = train_idx,
      val_idx    = val_idx
    )
  })
}

# Scheme D: Spatio-Temporal Leave-One-Site-Out CV
# For each of the 12 stations: train on the OTHER 11 stations in earlier training
# years (season_year 2013-2014), and validate on the HELD-OUT 12th station in the
# final training year (season_year 2015). Eliminates same-day shared weather leakage.
make_spatiotemporal_loso_cv <- function(df) {
  stations <- levels(df$station)
  map2(seq_along(stations), stations, function(f_num, st) {
    train_idx <- which(df$station != st & df$season_year %in% c(2013, 2014))
    val_idx   <- which(df$station == st & df$season_year == 2015)
    list(
      scheme     = "spatiotemporal_loso_12fold",
      fold_num   = f_num,
      fold_id    = sprintf("LOSO_%s_2015", st),
      descriptor = sprintf(
        "Holdout Site: %s in Season-Year 2015 (Trained on 11 other sites in 2013–2014)",
        st
      ),
      train_idx  = train_idx,
      val_idx    = val_idx
    )
  })
}

cv_folds <- list(
  random_5fold               = make_random_cv(train_df, k = 5),
  blocked_gap_6fold          = make_blocked_gap_cv(train_df, n_blocks = 6, buffer_days = 7),
  temporal_year_3fold        = make_year_temporal_cv(train_df, buffer_days = 7),
  spatiotemporal_loso_12fold = make_spatiotemporal_loso_cv(train_df)
)

# Build a structured summary table of every fold across all 4 CV schemes
summarize_cv_scheme <- function(folds_list, df) {
  map_dfr(folds_list, function(f) {
    val_y   <- df$y_exceed[f$val_idx]
    train_y <- df$y_exceed[f$train_idx]
    tibble(
      scheme            = f$scheme,
      cv_fold           = f$fold_num,
      fold_id           = f$fold_id,
      descriptor        = f$descriptor,
      train_n           = length(f$train_idx),
      train_exceed_pct  = 100 * mean(train_y),
      val_n             = length(f$val_idx),
      val_exceed_n      = sum(val_y),
      val_exceed_pct    = 100 * mean(val_y)
    )
  })
}

cv_fold_summaries <- map_dfr(cv_folds, summarize_cv_scheme, df = train_df)

saveRDS(cv_folds, file.path(out_dir, "cv_folds.rds"))
write_csv(cv_fold_summaries, file.path(out_dir, "cv_fold_summaries.csv"))
message("Saved CV fold definitions and split summaries to data/processed/")

# ---- 4. Fold-Safe Imputation & Preprocessing Helper ----------

# Wrapper around station-median calculation that preserves all 12 factor levels
# of `station` even when a station is completely held out in LOSO CV.
fit_station_medians_safe <- function(train_sub, vars, all_stations = levels(station_day$station)) {
  vars <- unique(vars)
  
  by_station <- train_sub |>
    group_by(station) |>
    summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), .groups = "drop") |>
    pivot_longer(-station, names_to = "variable", values_to = "station_median")
  
  overall <- train_sub |>
    summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE))) |>
    pivot_longer(everything(), names_to = "variable", values_to = "overall_median")
  
  expand_grid(
    station  = factor(all_stations, levels = all_stations),
    variable = vars
  ) |>
    left_join(by_station, by = c("station", "variable")) |>
    left_join(overall, by = "variable") |>
    mutate(fill = coalesce(station_median, overall_median)) |>
    select(station, variable, fill)
}

#' Prepare a single training/validation (or train/test) split without leakage
#'
#' @param train_sub Data frame of training rows for this fold
#' @param val_sub   Data frame of validation (or test) rows for this fold
#' @param window    Either "24h" or "48h"
#' @return A list with:
#'   - train_df, val_df: Imputed data frames on original physical scale (for GAM)
#'   - X_train, X_val:   Standardized dummy-encoded numeric matrices (for glmnet / XGBoost)
#'   - y_train, y_val:   Numeric 0/1 outcome vectors
#'   - scale_params:     Tibble of training means and SDs used for standardization
prep_fold_data <- function(train_sub, val_sub, window = c("24h", "48h")) {
  window <- match.arg(window)
  pred_cols <- if (window == "24h") predictors_24h else predictors_48h
  num_cols  <- setdiff(pred_cols, calendar_cols)
  
  # 1. Fit station-median imputation strictly on train_sub and apply to both
  meds      <- fit_station_medians_safe(train_sub, num_cols)
  train_imp <- impute_station_medians(train_sub, meds)
  val_imp   <- impute_station_medians(val_sub,   meds)
  
  # 2. Drop any numeric predictor with zero variance inside train_sub
  sds <- map_dbl(num_cols, ~ sd(train_imp[[.x]], na.rm = TRUE))
  keep_num <- num_cols[!is.na(sds) & sds > 0]
  keep_all <- c(calendar_cols, keep_num)
  
  # 3. Compute scaling parameters (mean & SD) strictly on train_imp
  scale_params <- tibble(
    variable = keep_num,
    mean     = map_dbl(keep_num, ~ mean(train_imp[[.x]])),
    sd       = map_dbl(keep_num, ~ sd(train_imp[[.x]]))
  )
  
  # 4. Create scaled copies for regularized regression / matrix generation
  train_scaled <- train_imp
  val_scaled   <- val_imp
  for (i in seq_len(nrow(scale_params))) {
    v <- scale_params$variable[i]
    m <- scale_params$mean[i]
    s <- scale_params$sd[i]
    train_scaled[[v]] <- (train_imp[[v]] - m) / s
    val_scaled[[v]]   <- (val_imp[[v]]   - m) / s
  }
  
  # 5. Build model matrices (one-hot/dummy encoding for factors, no intercept column)
  rhs_formula <- reformulate(keep_all)
  X_train <- model.matrix(rhs_formula, data = train_scaled)[, -1, drop = FALSE]
  X_val   <- model.matrix(rhs_formula, data = val_scaled)[, -1, drop = FALSE]
  
  list(
    window       = window,
    predictors   = keep_all,
    num_cols     = keep_num,
    train_df     = train_imp |> select(all_of(c(id_cols, outcome_cols, benchmark_cols, keep_all))),
    val_df       = val_imp   |> select(all_of(c(id_cols, outcome_cols, benchmark_cols, keep_all))),
    X_train      = X_train,
    y_train      = as.integer(train_imp$y_exceed),
    X_val        = X_val,
    y_val        = as.integer(val_imp$y_exceed),
    scale_params = scale_params,
    imputation   = meds
  )
}

# Quick self-test on Fold 1 of the spatio-temporal LOSO scheme when run as script
if (sys.nframe() == 0L) {
  f1 <- cv_folds$spatiotemporal_loso_12fold[[1]]
  prepped_24 <- prep_fold_data(
    train_df[f1$train_idx, ],
    train_df[f1$val_idx, ],
    window = "24h"
  )
  stopifnot(
    "Missing values remain in X_train" = !anyNA(prepped_24$X_train),
    "Missing values remain in X_val"   = !anyNA(prepped_24$X_val)
  )
  message(sprintf(
    "Self-test passed on %s (24h): X_train (%d x %d), X_val (%d x %d), 0 NAs.",
    f1$fold_id,
    nrow(prepped_24$X_train), ncol(prepped_24$X_train),
    nrow(prepped_24$X_val),   ncol(prepped_24$X_val)
  ))
}