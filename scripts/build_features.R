#!/usr/bin/env Rscript
# Build the station-day analysis table for the Beijing PM2.5 project.
#
# Run from the top level of the repo, after scripts/download_data.R:
#
#     Rscript scripts/build_features.R
#
# Writes
#   data/processed/station_day_features.csv   one row per station-day: the
#                                             next-day outcome plus 24 h and
#                                             48 h predictor summaries
#   data/processed/cleaning_log.csv           every cleaning rule and the
#                                             number of hourly values it touched
#
# Design (open decisions are listed in the README):
#   * The prediction is made at the end of day d (after the 23:00 reading)
#     for day d + 1. Nothing from day d + 1 enters the features.
#   * 24 h features use hours 00-23 of day d; 48 h features use days d - 1
#     and d. Both are built so the group can decide which window to use.
#   * Outcome: mean PM2.5 on day d + 1 above 75 ug/m3, computed from the raw
#     (not gap-filled) hourly values, and only when >= 20 hours were observed.
#   * Short gaps (<= 3 h) in the hourly predictors are filled by carrying the
#     last value forward. Longer gaps stay NA and show up in the missingness
#     features (n_obs_pm25_*, n_outage_*). NA left in the daily table are for
#     scripts/impute.R, which must be fit on training rows only.
#   * Train / test split: target days before 2016-03-01 are training (three
#     March-to-February years), the last year is the test set. The column
#     season_year supports rolling-origin CV within training.

suppressPackageStartupMessages(library(tidyverse))

# ---- settings ----------------------------------------------------------------

raw_dir <- file.path("data", "raw")
out_dir <- file.path("data", "processed")

pm25_threshold    <- 75                    # ug/m3, China Grade II 24 h standard
min_hours_outcome <- 20                    # hours of PM2.5 needed for a valid daily mean
max_locf_gap      <- 3                     # longest gap (hours) filled by carry-forward
o3_max_plausible  <- 500                   # ug/m3; above this is a sensor fault
test_start        <- as.Date("2016-03-01") # first *target* day of the test period

pollutants <- c("PM2.5", "PM10", "SO2", "NO2", "CO", "O3")
weather    <- c("TEMP", "PRES", "DEWP", "RAIN", "WSPM")
compass    <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- read --------------------------------------------------------------------

files <- list.files(raw_dir, pattern = "^PRSA_Data_.*\\.csv$", full.names = TRUE)
stopifnot("Run scripts/download_data.R first" = length(files) == 12)

hourly <- files |>
  map(read_csv, show_col_types = FALSE) |>
  list_rbind() |>
  select(-No) |>
  mutate(
    datetime = make_datetime(year, month, day, hour, tz = "Asia/Shanghai"),
    date     = as_date(datetime),
    station  = factor(station)
  ) |>
  arrange(station, datetime)

stopifnot(nrow(hourly) == 420768)

# ---- cleaning ----------------------------------------------------------------
# Every rule is logged with how many hourly values it changed or flagged, so
# the write-up can report the numbers directly.

cleaning_log <- tibble(rule = character(), action = character(), n_hours = integer())
log_rule <- function(rule, action, n) {
  cleaning_log <<- add_row(cleaning_log, rule = rule, action = action,
                           n_hours = as.integer(n))
}

# 1. O3: one block of hours stuck at 1071 ug/m3 (Dongsi, 15-17 Aug 2013). The
#    99.9th percentile is 312, so anything above 500 is a sensor fault.
n_o3 <- sum(hourly$O3 > o3_max_plausible, na.rm = TRUE)
hourly <- hourly |>
  mutate(O3 = if_else(O3 > o3_max_plausible, NA, O3))
log_rule(sprintf("O3 > %d ug/m3 (stuck sensor, Dongsi Aug 2013)", o3_max_plausible),
         "set to NA", n_o3)

# 2. PM2.5 / PM10 equal to 999 look like an instrument ceiling, but they fall on
#    real extreme-pollution events (2016-02-08 02:00 is Chinese New Year's Eve
#    fireworks, with neighbouring stations at 941 and 957; 2015-04-15 was a
#    dust storm). They are kept as right-censored values.
log_rule("PM2.5 == 999", "kept (instrument ceiling during a real event)",
         sum(hourly$`PM2.5` == 999, na.rm = TRUE))
log_rule("PM10 == 999", "kept (instrument ceiling during a real event)",
         sum(hourly$PM10 == 999, na.rm = TRUE))

# 3. PM2.5 > PM10 is physically impossible (PM2.5 is a subset of PM10), but it
#    affects about 5% of hours, so it is measurement noise between two
#    instruments rather than a handful of bad rows. Both values are kept.
log_rule("PM2.5 > PM10", "kept (inter-instrument noise, reported)",
         sum(hourly$`PM2.5` > hourly$PM10, na.rm = TRUE))

# 4. Hours where all six pollutants are missing are monitor outages. They are
#    flagged so the window features can count them; the values stay NA.
hourly <- hourly |>
  mutate(outage = if_all(all_of(pollutants), is.na))
log_rule("all six pollutants missing (monitor outage)", "flagged", sum(hourly$outage))

# 5. Carry the last observed value forward across gaps of at most max_locf_gap
#    hours. Only the predictors are filled; the raw PM2.5 is kept for the
#    outcome and for counting observed hours.
fill_short_gaps <- function(x, max_gap) {
  if (!anyNA(x)) return(x)
  runs    <- rle(is.na(x))
  run_len <- rep(runs$lengths, runs$lengths)   # length of the run each hour sits in
  carried <- vctrs::vec_fill_missing(x, direction = "down")
  if_else(is.na(x) & run_len <= max_gap, carried, x)
}

locf_vars <- c(pollutants, weather, "wd")
n_na_before <- sum(is.na(hourly[locf_vars]))

hourly <- hourly |>
  mutate(pm25_raw = `PM2.5`) |>
  group_by(station) |>
  mutate(across(all_of(locf_vars), ~ fill_short_gaps(.x, max_locf_gap))) |>
  ungroup()

log_rule(sprintf("predictor gaps <= %d h", max_locf_gap),
         "last observed value carried forward", n_na_before - sum(is.na(hourly[locf_vars])))

# Wind as a vector (meteorological convention: wd is where the wind comes FROM),
# so that direction can be averaged over a window.
hourly <- hourly |>
  mutate(
    wd_deg = (match(wd, compass) - 1) * 22.5,
    wind_u = -WSPM * sin(wd_deg * pi / 180),
    wind_v = -WSPM * cos(wd_deg * pi / 180)
  )

# ---- daily outcome -----------------------------------------------------------

daily_pm25 <- hourly |>
  group_by(station, date) |>
  summarise(
    hours_obs_pm25 = sum(!is.na(pm25_raw)),
    pm25_mean = if_else(hours_obs_pm25 >= min_hours_outcome,
                        mean(pm25_raw, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  mutate(exceed = pm25_mean > pm25_threshold)

# ---- window features ---------------------------------------------------------

# Summarise all hours whose win_date is the same into one row. Called once with
# win_date = date (24 h) and once with each hour duplicated into the following
# day as well (48 h). min/max of an all-NA window give +/-Inf; those become NA.
summarise_window <- function(h, suffix) {
  h |>
    group_by(station, win_date) |>
    summarise(
      n_hours = n(),
      across(all_of(pollutants),
             list(mean = ~ mean(.x, na.rm = TRUE), max = ~ max(.x, na.rm = TRUE))),
      across(c(TEMP, PRES, DEWP, WSPM),
             list(mean = ~ mean(.x, na.rm = TRUE), min = ~ min(.x, na.rm = TRUE),
                  max = ~ max(.x, na.rm = TRUE))),
      RAIN_total  = sum(RAIN, na.rm = TRUE),
      wind_u_mean = mean(wind_u, na.rm = TRUE),
      wind_v_mean = mean(wind_v, na.rm = TRUE),
      n_obs_pm25  = sum(!is.na(pm25_raw)),
      n_outage    = sum(outage),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA))) |>
    rename_with(~ paste0(.x, "_", suffix), -c(station, win_date))
}

win24 <- hourly |>
  mutate(win_date = date) |>
  summarise_window("24h") |>
  suppressWarnings()

win48 <- bind_rows(hourly |> mutate(win_date = date),
                   hourly |> mutate(win_date = date + 1)) |>
  filter(win_date <= max(hourly$date)) |>
  summarise_window("48h") |>
  suppressWarnings()

# The most recent reading before the prediction is made.
pm25_last <- hourly |>
  filter(hour == 23) |>
  select(station, date, PM25_last = `PM2.5`)

# ---- assemble ----------------------------------------------------------------

station_day <- daily_pm25 |>
  rename(PM25_today_mean = pm25_mean, exceed_today = exceed,
         hours_obs_pm25_today = hours_obs_pm25) |>
  # outcome: the next day's row, shifted back by one day
  left_join(daily_pm25 |>
              transmute(station, date = date - 1,
                        y_pm25_next = pm25_mean, y_exceed = exceed),
            by = c("station", "date")) |>
  left_join(win24,     by = c("station", "date" = "win_date")) |>
  left_join(win48,     by = c("station", "date" = "win_date")) |>
  left_join(pm25_last, by = c("station", "date")) |>
  mutate(
    target_date = date + 1,
    period      = factor(if_else(target_date < test_start, "train", "test"),
                         levels = c("train", "test")),
    # March-to-February years, labelled by the starting calendar year
    season_year = year(target_date) - (month(target_date) <= 2),
    month       = month(target_date),
    season      = case_when(month %in% c(12, 1, 2) ~ "winter",
                            month %in% 3:5         ~ "spring",
                            month %in% 6:8         ~ "summer",
                            TRUE                   ~ "autumn"),
    dow         = wday(target_date, label = TRUE, week_start = 1)
  ) |>
  rename_with(~ str_replace(.x, "PM2\\.5", "PM25")) |>
  select(station, date, target_date, period, season_year, season, month, dow,
         y_pm25_next, y_exceed,
         hours_obs_pm25_today, PM25_today_mean, exceed_today, PM25_last,
         ends_with("_24h"), ends_with("_48h"))

n_all     <- nrow(station_day)
n_no_outc <- sum(is.na(station_day$y_exceed))
station_day <- station_day |> filter(!is.na(y_exceed))

log_rule(sprintf("station-days without a valid next-day outcome (< %d h of PM2.5, or last day)",
                 min_hours_outcome),
         "dropped from the analysis table", n_no_outc)

# ---- write -------------------------------------------------------------------

write_csv(station_day,  file.path(out_dir, "station_day_features.csv"))
write_csv(cleaning_log, file.path(out_dir, "cleaning_log.csv"))

message(sprintf("Wrote %d of %d station-days (%d stations, %s to %s) to %s",
                nrow(station_day), n_all, n_distinct(station_day$station),
                min(station_day$target_date), max(station_day$target_date), out_dir))
print(count(station_day, period, y_exceed) |>
        group_by(period) |> mutate(pct = round(100 * n / sum(n), 1)) |> ungroup())
