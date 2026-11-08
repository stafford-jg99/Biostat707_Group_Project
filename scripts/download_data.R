#!/usr/bin/env Rscript
# Download the Beijing Multi-Site Air Quality data from the UCI repository.
#
# Run from the top level of the repo:
#
#     Rscript scripts/download_data.R
#
# Writes the 12 station CSVs to data/raw/ and checks the file and row counts.
# Uses only base R.
#
# Source:  https://archive.ics.uci.edu/dataset/501
# License: CC BY 4.0
# Cite:    Chen, S. (2017). Beijing Multi-Site Air Quality [Dataset]. UCI Machine
#          Learning Repository. https://doi.org/10.24432/C5RK5G

url <- "https://archive.ics.uci.edu/static/public/501/beijing+multi+site+air+quality+data.zip"
raw_dir <- file.path("data", "raw")
expected_files <- 12L
expected_rows <- 420768L
pattern <- "^PRSA_Data_.*\\.csv$"

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

if (length(list.files(raw_dir, pattern = pattern)) == expected_files) {
  message("Found existing files in ", raw_dir, "; skipping download.")
} else {
  tmp <- tempfile("uci_")
  dir.create(tmp)
  outer_zip <- file.path(tmp, "download.zip")

  options(timeout = max(600, getOption("timeout")))
  download.file(url, outer_zip, mode = "wb")

  # The download holds another zip, so unzip until none are left.
  done <- character()
  repeat {
    zips <- setdiff(
      list.files(tmp, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE),
      done
    )
    if (length(zips) == 0) break
    for (z in zips) unzip(z, exdir = dirname(z))
    done <- c(done, zips)
  }

  csvs <- list.files(tmp, pattern = pattern, recursive = TRUE, full.names = TRUE)
  csvs <- csvs[!grepl("__MACOSX", csvs, fixed = TRUE)]
  file.copy(csvs, raw_dir, overwrite = TRUE)
  unlink(tmp, recursive = TRUE)
}

files <- list.files(raw_dir, pattern = pattern, full.names = TRUE)
rows <- sum(vapply(files, function(f) length(readLines(f, warn = FALSE)) - 1L, integer(1)))

if (length(files) != expected_files || rows != expected_rows) {
  stop(sprintf(
    "Check failed: got %d files and %d rows, expected %d files and %d rows.",
    length(files), rows, expected_files, expected_rows
  ))
}
message(sprintf("OK: %d files, %d rows in %s", length(files), rows, raw_dir))
