# Biostat707_Group_Project
The Repository for the 707 Project

# Project Goals:
## Project Overview & Predictive Objective
* Predict whether tomorrow’s PM2.5 levels will exceed a predefined health threshold (i.e., 75 µg/m³) by modeling hourly pollutant and weather readings from a fixed 24- or 48-hour observation window.
* Forecast Setup: Predictions are made at 23:00 on day *d* for day *d + 1* using only prior pollutant, meteorological, calendar, and sensor-availability summaries from a fixed **24-hour** or **48-hour** lookback window
* Evaluate and quantify the optimistic performance bias caused by temporal leakage by actively comparing naive random cross-validation results against a strict temporal validation design.

## Dataset Characteristics
* Source: The Beijing Multi-Site Air Quality dataset, sourced from the UCI Machine Learning Repository.
* Structure: Multivariate, time-series data comprising 420,768 total instances. Feature types consist of real and integer values.
* Variables: The data includes 6 main air pollutants and 6 relevant meteorological variables.
* Collection Sites: Hourly data was collected from 12 nationally-controlled air-quality monitoring sites across Beijing. The pollutant readings were supplied by the Beijing Municipal Environmental Monitoring Center and matched with meteorological data from the nearest weather stations.

## Validation & Leakage Prevention Strategy
Because daily PM2.5 is autocorrelated across multi-day weather fronts and strongly correlated across Beijing’s 12 stations on the same day, we evaluate all candidate models on the 2013–2015 training set using a three-tier cross-validation ladder:

1. **Naive Stratified Random 5-Fold CV (Leaking Benchmark):** Randomly assigns individual station-day rows (20% per fold) across all dates and stations, stratified only by the binary exceedance outcome. Intentionally permits adjacent-day and same-day spatial leakage to measure optimistic bias.
2. **Blocked Temporal 6-Fold CV with a ±7-Day Buffer:** Partitions the 3-year training timeline into 6 contiguous ~6-month validation blocks across all 12 stations simultaneously. A ±7-day buffer adjacent to each validation block is dropped from training (82–168 station-days per fold), eliminating same-day spatial leakage and short-term weather autocorrelation.
3. **Spatio-Temporal Leave-One-Site-Out (LOSO) 12-Fold CV:** Loops across all 12 stations, training strictly on the other 11 stations during earlier season-years (2013–2014) and validating on the held-out 12th station in season-year 2015. Combines unseen-site transfer with a forward temporal split to avoid same-day regional weather leakage (measuring within-Beijing site transfer rather than cross-city transfer).

All preprocessing—including station-median imputation and feature standardization—is fit strictly inside each training fold and applied to validation rows without leakage.

## Modeling, Comparison & Final Model Selection

We compared two model families across both observation windows (24h and 48h) against a **Persistence Baseline** (predicting tomorrow’s exceedance from today’s exceedance flag):

* **Model 1 — Elastic Net Logistic Regression:** Regularized linear classifier with mixing parameter ($\alpha$) and penalty strength ($\lambda$) tuned via nested 5-fold temporal inner CV.
* **Model 2 — Gradient Boosted Trees (XGBoost):** Non-linear tree ensemble with tree depth and early-stopping boosting rounds (`nrounds`) tuned via nested 5-fold temporal inner CV.

### Out-of-Fold Performance Summary (24-Hour Window, 2013–2015 Training Set)

| Validation Design | Persistence AUROC | Elastic Net (24h) AUROC | XGBoost (24h) AUROC | Paired ΔAUROC (XGB − ENet, 95% CI)* | Paired ΔBrier (XGB − ENet, 95% CI)* |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **1. Naive Random 5-Fold CV** | 0.686 ± 0.004 | 0.879 ± 0.005 | 0.908 ± 0.007 | +0.0293 (+0.0237, +0.0353) | −0.0187 (−0.0224, −0.0153) |
| **2. Blocked Temporal 6-Fold (±7d Gap)** | 0.684 ± 0.048 | 0.874 ± 0.025 | 0.875 ± 0.030 | +0.0042 (−0.0022, +0.0112) | −0.0046 (−0.0087, −0.0009) |
| **3. Spatio-Temporal LOSO 12-Fold** | 0.722 ± 0.018 | 0.891 ± 0.010 | 0.904 ± 0.010 | +0.0135 (+0.0049, +0.0215) | −0.0133 (−0.0181, −0.0087) |

*\*95% percentile CIs from a paired 14-day block bootstrap ($B = 2,000$ resamples) on pooled out-of-fold predictions, keeping all 12 stations and consecutive days within each 2-week block together. Positive ΔAUROC and negative ΔBrier favor XGBoost.*

### Key Findings & Selected Model

* **24-Hour Window Outperforms 48-Hour Window:** Across every validation design and both model classes, 24-hour features beat 48-hour features on AUROC, log loss, and Brier score (e.g., Blocked Temporal AUROC: `0.874` vs. `0.867` for Elastic Net; `0.875` vs. `0.869` for XGBoost), indicating that averaging over 48 hours dilutes recent atmospheric shifts.
* **Leakage Distorts Both Mean Performance and Uncertainty:** Naive Random CV inflated XGBoost’s AUROC lead to `+0.029` by letting deeper trees (`median depth = 4`) memorize leaked neighbors, and understated temporal AUROC standard deviation four- to fivefold (`±0.005–0.007` vs. `±0.025–0.030` in Blocked Temporal CV). Under stricter validation, XGBoost simplified to shallow `depth = 2` trees across both Blocked Temporal and LOSO designs.
* **Final Model Choice — 24-Hour Elastic Net (`output/final_model.rds`):** Across the full 2013–2015 timeline (Blocked Temporal CV), both 24-hour models tied in discrimination (`0.8744` vs. `0.8749`; paired ΔAUROC 95% CI includes zero) and Elastic Net outperformed XGBoost in 4 of the 6 temporal blocks. Although XGBoost achieved a modest Brier reduction (`−0.0046`) and stronger single-year 2015 LOSO transfer (`+0.0135` AUROC), we selected the **24-hour Elastic Net** for parsimony, stability, and direct interpretability, refitting it once on the complete 2013–2015 training set.
* **Stable Core Drivers:** Out of 40 candidate terms, Elastic Net pruned 22 to zero in the full refit and retained a compact core of **8 predictors in 100% of folds across all three designs**: the 23:00 PM2.5 reading (standardized `OR = 5.36`), 24-hour mean southerly wind component (`OR = 1.75`), minimum barometric pressure (`OR = 1.62`), peak and mean ozone (`OR = 1.33`, `1.22`), peak wind speed (`OR = 0.82`), cumulative rainfall (`OR = 0.97`), and calendar month (`OR = 0.97`).

---

## Data Engineering Challenges
* Missingness: The dataset contains known missing values. The modeling pipeline must address informative missing data resulting from gaps in the monitoring site sensor records.
* Temporal Drift: Manage class imbalance and handle seasonal environmental drift across the longitudinal hourly observations.

# Group Members and Roles (Name, Email)
1. Andreas Tjavaras (andreas.tjavaras@duke.edu):
2. Tim Deng (tim.deng@duke.edu):
3. Joe Stafford (joseph.stafford@duke.edu): 
4. JT Rapp (jt.rapp@duke.edu):

# Tentative Assignment of Roles (subject to change as project advances):
Data Wrangling: Joe, Tim 

Descriptive Analysis: Joe, Tim

Modeling: JT

Modeling Evaluation: Andreas

Project write-up/presentation: Everyone 

# Tasks to be completed: 
Data Cleaning, handling of missing data, descriptive analysis, predictor/Variable Selection, Cross Validation Process, Model fitting, etc.

# Repository layout

```
scripts/download_data.R              downloads the 12 station CSVs to data/raw/ (gitignored)
scripts/build_features.R             cleans hourly data and builds the station-day modeling table
scripts/impute.R                     fold-safe station-median imputation (fit on training rows only)
scripts/cv_and_preprocess.R          generates CV folds and fold-level preprocessing pipelines
scripts/fit_enet.R                   fits nested-CV Elastic Net models (24h & 48h) and full refit*
scripts/fit_xgboost.R                fits nested-CV XGBoost models (24h & 48h)*
scripts/compare_models.R             compiles cross-scheme comparison tables and ROC coordinates
scripts/paired_bootstrap.R           runs B = 2,000 paired 14-day block-bootstrap comparisons
707_finalproj_code.qmd               full analysis report: EDA, cleaning log, CV, modeling & selection
data/raw/                            hourly CSVs from UCI (gitignored, ~32 MB)
data/processed/                      station_day_features.csv, cleaning_log.csv, model RDS objects
output/final_model.rds               locked 24-hour Elastic Net model refit on full 2013–2015 training set
output/oof_predictions.csv           out-of-fold predictions for the selected model across CV schemes
output/model_comparison_summary.rds  compiled comparison tables and ROC curve data
output/paired_bootstrap_results.rds  paired 14-day block-bootstrap draws and 95% CIs
```

To reproduce the analysis table from a fresh clone:

```r
# 1. Download raw hourly data and build the station-day modeling table
Rscript scripts/download_data.R
Rscript scripts/build_features.R

# 2. Generate CV folds, fit models, compare schemes, and run paired block bootstrap
Rscript scripts/cv_and_preprocess.R
Rscript scripts/fit_enet.R
Rscript scripts/fit_xgboost.R
Rscript scripts/compare_models.R
Rscript scripts/paired_bootstrap.R

# 3. Render the full report
quarto render 707_finalproj_code.qmd
```

# Analysis table (`data/processed/station_day_features.csv`)

One row per station and day *d*, 17,035 rows. The prediction is made at the end of day *d* for day *d + 1*.

| Column(s) | Meaning |
|---|---|
| `station`, `date`, `target_date` | station, feature day *d*, outcome day *d + 1* |
| `period`, `season_year`, `season`, `month`, `dow` | train/test flag; March-to-February year of the target day (2013–2016); calendar features |
| `y_pm25_next`, `y_exceed` | mean PM2.5 on day *d + 1*, and whether it is above 75 µg/m³ (the outcome) |
| `PM25_today_mean`, `exceed_today`, `PM25_last` | today's mean and exceedance (persistence baseline), and the 23:00 reading |
| `*_mean_24h`, `*_max_24h`, `*_min_24h`, `RAIN_total_24h`, `wind_u_mean_24h`, `wind_v_mean_24h` | summaries of hours 00–23 of day *d* |
| `*_48h` | the same summaries over days *d − 1* and *d* |
| `n_hours_*`, `n_obs_pm25_*`, `n_outage_*` | hours in the window, hours with PM2.5 observed, hours with all pollutants missing |

Cleaning rules and the number of hourly values each touched are in `data/processed/cleaning_log.csv` and in the qmd.

# Open decisions

* **Prediction window.** The table carries both 24 h and 48 h features. Pick one, or compare them as a secondary analysis; using both at once doubles the feature count with highly correlated columns.
* **Train/test split.** Currently target days before 2016-03-01 are training (three years) and the final year is test. Change `test_start` in `scripts/build_features.R` if the group prefers a different cut. `season_year` allows rolling-origin validation inside the training years.
* **Imputation.** `scripts/impute.R` must be fit on the training rows of each fold and applied to the held-out rows; fitting it on the whole table first is leakage.

#Project Status & Next Steps
##Completed
* Data ingestion, hourly cleaning log, and station-day feature engineering (24h and 48h windows).

* Exploratory data analysis, missingness auditing, and spatio-temporal correlation assessment.

* Leakage-resistant cross-validation framework, nested-CV Elastic Net and XGBoost training, paired 14-day block-bootstrap comparison, and final model lock-in (output/final_model.rds).

## Remaining

* One-time evaluation of output/final_model.rds on the untouched 2016–2017 held-out test set (period == "test").

* Comparison of held-out test discrimination, Brier score, and calibration against the three CV estimates to verify which validation design best tracked true future performance.

* Subgroup evaluation (by station and season), error analysis of confident misclassifications, and final report/presentation synthesis.

# Data source and license

Chen, S. (2017). Beijing Multi-Site Air Quality [Dataset]. UCI Machine Learning Repository. https://doi.org/10.24432/C5RK5G. License: CC BY 4.0.
