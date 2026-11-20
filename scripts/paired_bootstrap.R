#!/usr/bin/env Rscript
# Paired Block-Bootstrap Comparison: Model 2 (XGBoost, 24h) vs. Model 1 (Elastic Net, 24h)
#
# Run from the top level of the repo, after scripts/compare_models.R:
#
#     Rscript scripts/paired_bootstrap.R
#
# Writes:
#   output/paired_bootstrap_results.rds          bootstrap draws & 95% CIs across all 4 CV schemes
#   data/processed/paired_bootstrap_results.rds  copy in data/processed/

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
})

# Prefer data/processed/oof_predictions.csv (guaranteed to hold both models)
oof_path <- if (file.exists(file.path("data", "processed", "oof_predictions.csv"))) {
  file.path("data", "processed", "oof_predictions.csv")
} else {
  file.path("output", "oof_predictions.csv")
}

stopifnot("Run scripts/compare_models.R first" = file.exists(oof_path))
oof_all <- read_csv(oof_path, show_col_types = FALSE)

stopifnot(
  "Both p_enet_min and p_xgb must exist in oof_predictions.csv (re-run scripts/compare_models.R)" =
    all(c("p_enet_min", "p_xgb") %in% names(oof_all))
)

oof_24 <- oof_all |>
  filter(window == "24h") |>
  arrange(scheme, date, station)

# ---- 1. Exact Metric Functions (Verified Against pROC::auc) ------------------

# Exact Hanley-McNeil / Mann-Whitney AUROC (ties counted as half, identical to pROC::auc)
fast_auc <- function(y, p) {
  n1 <- as.numeric(sum(y == 1))
  n0 <- as.numeric(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_brier <- function(y, p) {
  mean((p - y)^2)
}

calc_logloss <- function(y, p, eps = 1e-15) {
  p_clip <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p_clip) + (1 - y) * log(1 - p_clip))
}

# Verify fast_auc() matches pROC::auc() to machine precision
check_sub <- oof_24 |> filter(scheme == unique(oof_24$scheme)[1])
auc_fast  <- fast_auc(check_sub$y_exceed, check_sub$p_enet_min)
auc_proc  <- as.numeric(pROC::auc(check_sub$y_exceed, check_sub$p_enet_min,
                                  levels = c(0, 1), direction = "<", quiet = TRUE))
stopifnot("fast_auc mismatch with pROC::auc" = abs(auc_fast - auc_proc) < 1e-12)

# ---- 2. Paired Temporal Block Bootstrap Engine -------------------------------
# Resamples contiguous multi-day blocks (default: 14 days), keeping all 12
# stations and consecutive days within each block together to respect both
# spatial correlation across sites and multi-day weather autocorrelation.

run_paired_block_bootstrap <- function(df, B = 2000L, block_days = 14L, seed = 707L) {
  set.seed(seed)
  
  df <- df |> arrange(date, station)
  unique_dates <- sort(unique(df$date))
  n_dates      <- length(unique_dates)
  
  # Partition timeline into contiguous 14-day blocks (merge tiny trailing remainder)
  raw_block <- ((seq_len(n_dates) - 1L) %/% block_days) + 1L
  max_b     <- max(raw_block)
  if (max_b > 1L && sum(raw_block == max_b) < (block_days / 2)) {
    raw_block[raw_block == max_b] <- max_b - 1L
  }
  
  date_map <- tibble(date = unique_dates, block_id = raw_block)
  df <- df |> left_join(date_map, by = "date")
  
  block_indices <- split(seq_len(nrow(df)), df$block_id)
  n_blocks      <- length(block_indices)
  
  y   <- as.integer(df$y_exceed)
  p_a <- df$p_enet_min  # Model A: Elastic Net (24h)
  p_b <- df$p_xgb       # Model B: XGBoost (24h)
  
  # Observed pooled out-of-fold metrics
  obs_auc_a   <- fast_auc(y, p_a)
  obs_auc_b   <- fast_auc(y, p_b)
  obs_brier_a <- calc_brier(y, p_a)
  obs_brier_b <- calc_brier(y, p_b)
  obs_ll_a    <- calc_logloss(y, p_a)
  obs_ll_b    <- calc_logloss(y, p_b)
  
  # Pre-allocate bootstrap draw vectors
  auc_a_b   <- numeric(B); auc_b_b   <- numeric(B)
  brier_a_b <- numeric(B); brier_b_b <- numeric(B)
  ll_a_b    <- numeric(B); ll_b_b    <- numeric(B)
  
  for (k in seq_len(B)) {
    drawn <- sample.int(n_blocks, n_blocks, replace = TRUE)
    idx   <- unlist(block_indices[drawn], use.names = FALSE)
    
    y_k  <- y[idx]
    pa_k <- p_a[idx]
    pb_k <- p_b[idx]
    
    # Paired evaluation: both models scored on the exact same resampled blocks
    auc_a_b[k]   <- fast_auc(y_k, pa_k)
    auc_b_b[k]   <- fast_auc(y_k, pb_k)
    brier_a_b[k] <- calc_brier(y_k, pa_k)
    brier_b_b[k] <- calc_brier(y_k, pb_k)
    ll_a_b[k]    <- calc_logloss(y_k, pa_k)
    ll_b_b[k]    <- calc_logloss(y_k, pb_k)
  }
  
  draws <- tibble(
    iter       = seq_len(B),
    auc_enet   = auc_a_b,
    auc_xgb    = auc_b_b,
    d_auc      = auc_b_b - auc_a_b,       # Positive favors XGBoost
    brier_enet = brier_a_b,
    brier_xgb  = brier_b_b,
    d_brier    = brier_b_b - brier_a_b,   # Negative favors XGBoost (lower is better)
    ll_enet    = ll_a_b,
    ll_xgb     = ll_b_b,
    d_logloss  = ll_b_b - ll_a_b          # Negative favors XGBoost (lower is better)
  )
  
  # Helper to pack one metric row
  pack_metric <- function(metric_name, obs_a, obs_b, vec_a, vec_b, diff_vec, higher_is_better) {
    ci_a <- quantile(vec_a,    c(0.025, 0.975), na.rm = TRUE)
    ci_b <- quantile(vec_b,    c(0.025, 0.975), na.rm = TRUE)
    ci_d <- quantile(diff_vec, c(0.025, 0.975), na.rm = TRUE)
    pct_xgb_wins <- if (higher_is_better) 100 * mean(diff_vec > 0, na.rm = TRUE)
    else                  100 * mean(diff_vec < 0, na.rm = TRUE)
    
    tibble(
      metric        = metric_name,
      n_rows        = nrow(df),
      n_dates       = n_dates,
      n_blocks      = n_blocks,
      block_days    = block_days,
      enet_est      = obs_a,
      enet_lo       = ci_a[1],
      enet_hi       = ci_a[2],
      xgb_est       = obs_b,
      xgb_lo        = ci_b[1],
      xgb_hi        = ci_b[2],
      diff_est      = obs_b - obs_a,
      diff_lo       = ci_d[1],
      diff_hi       = ci_d[2],
      excludes_zero = (ci_d[1] > 0) | (ci_d[2] < 0),
      corr_ab       = cor(vec_a, vec_b, use = "complete.obs"),
      pct_xgb_wins  = pct_xgb_wins
    )
  }
  
  summary_tbl <- bind_rows(
    pack_metric("AUROC",       obs_auc_a,   obs_auc_b,   auc_a_b,   auc_b_b,   draws$d_auc,     TRUE),
    pack_metric("Brier Score", obs_brier_a, obs_brier_b, brier_a_b, brier_b_b, draws$d_brier,   FALSE),
    pack_metric("Log Loss",    obs_ll_a,    obs_ll_b,    ll_a_b,    ll_b_b,    draws$d_logloss, FALSE)
  )
  
  list(summary = summary_tbl, draws = draws)
}

# ---- 3. Run Across All 4 Cross-Validation Schemes ----------------------------

schemes <- c(
  "random_5fold",
  "blocked_gap_6fold",
  "temporal_year_3fold",
  "spatiotemporal_loso_12fold"
)

message("Running 2,000 paired 14-day block-bootstrap resamples across all 4 CV schemes...")

boot_list <- map(schemes, function(s) {
  sub <- oof_24 |> filter(scheme == s)
  res <- run_paired_block_bootstrap(sub, B = 2000L, block_days = 14L, seed = 707L)
  list(
    summary = res$summary |> mutate(scheme = s, .before = 1),
    draws   = res$draws   |> mutate(scheme = s, .before = 1)
  )
})

boot_summary <- map_dfr(boot_list, "summary")
boot_draws   <- map_dfr(boot_list, "draws")

boot_results <- list(
  summary    = boot_summary,
  draws      = boot_draws,
  B          = 2000L,
  block_days = 14L
)

saveRDS(boot_results, file.path("output", "paired_bootstrap_results.rds"))
saveRDS(boot_results, file.path("data", "processed", "paired_bootstrap_results.rds"))

message("Saved paired bootstrap results to output/paired_bootstrap_results.rds")
print(
  boot_summary |>
    filter(metric %in% c("AUROC", "Brier Score")) |>
    select(scheme, metric, n_blocks, enet_est, xgb_est, diff_est, diff_lo, diff_hi, corr_ab, pct_xgb_wins)
)