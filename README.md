# Biostat707_Group_Project
The Repository for the 707 Project

# Project Goals:
## Project Overview & Predictive Objective
* Predict whether tomorrow’s PM2.5 levels will exceed a predefined health threshold by modeling hourly pollutant and weather readings from a fixed 24- or 48-hour observation window.
* Evaluate and quantify the optimistic performance bias caused by temporal leakage by actively comparing naive random cross-validation results against a strict temporal validation design.

## Dataset Characteristics
* Source: The Beijing Multi-Site Air Quality dataset, sourced from the UCI Machine Learning Repository.
* Structure: Multivariate, time-series data comprising 420,768 total instances. Feature types consist of real and integer values.
* Variables: The data includes 6 main air pollutants and 6 relevant meteorological variables.
* Collection Sites: Hourly data was collected from 12 nationally-controlled air-quality monitoring sites across Beijing. The pollutant readings were supplied by the Beijing Municipal Environmental Monitoring Center and matched with meteorological data from the nearest weather stations.

## Validation & Leakage Prevention Strategy
* Implement strict temporal validation by training models on earlier years (March 2013 to early 2017) and testing on later years to prevent temporal data leakage.
* Execute leave-one-site-out validation across the 12 monitoring stations to evaluate spatial generalizability across different regions of Beijing.

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
scripts/download_data.R        downloads the 12 station CSVs to data/raw/ (not committed)
scripts/build_features.R       cleans the hourly data and writes the station-day table
scripts/impute.R               station-median imputation, to be fit on training rows only
707_finalproj_code.qmd         exploratory analysis, cleaning log, Table 1
data/raw/                      hourly CSVs from UCI (gitignored, ~32 MB)
data/processed/                station_day_features.csv and cleaning_log.csv (gitignored)
```

To reproduce the analysis table from a fresh clone:

```r
Rscript scripts/download_data.R
Rscript scripts/build_features.R
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

# Data source and license

Chen, S. (2017). Beijing Multi-Site Air Quality [Dataset]. UCI Machine Learning Repository. https://doi.org/10.24432/C5RK5G. License: CC BY 4.0.
