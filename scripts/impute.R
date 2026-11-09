# Station-wise median imputation for the station-day feature table.
#
# The medians are estimated on training rows only and then applied to both
# training and test rows, so this can sit inside a cross-validation loop
# without leaking test information:
#
#     source("scripts/impute.R")
#     feature_cols <- grep("_(24h|48h)$|^PM25_last$", names(dat), value = TRUE)
#     imp   <- fit_station_medians(train, feature_cols)
#     train <- impute_station_medians(train, imp)
#     test  <- impute_station_medians(test,  imp)
#
# A station with no training value for some variable falls back to the overall
# training median. Missingness itself is already carried by the n_obs_pm25_*
# and n_outage_* columns, so no extra indicator columns are added here.

fit_station_medians <- function(train, vars) {
  by_station <- train |>
    group_by(station) |>
    summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), .groups = "drop") |>
    pivot_longer(-station, names_to = "variable", values_to = "station_median")

  overall <- train |>
    summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE))) |>
    pivot_longer(everything(), names_to = "variable", values_to = "overall_median")

  by_station |>
    left_join(overall, by = "variable") |>
    mutate(fill = coalesce(station_median, overall_median)) |>
    select(station, variable, fill)
}

impute_station_medians <- function(data, medians) {
  fills <- medians |>
    pivot_wider(names_from = variable, values_from = fill)
  row_of_station <- match(data$station, fills$station)

  for (v in unique(medians$variable)) {
    data[[v]] <- coalesce(data[[v]], fills[[v]][row_of_station])
  }
  data
}
