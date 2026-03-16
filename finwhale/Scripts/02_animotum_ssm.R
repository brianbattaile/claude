library(dplyr)
library(readr)
library(stringr)
library(aniMotum)

# ─── 1. Load location data ────────────────────────────────────────────────────
# For each whale folder, prefer the *-1-Locations.csv (FastGPS + Argos combined)
# over the plain *-Locations.csv (Argos only). Exclude Mac resource fork files (._).

whale_dirs <- list.dirs("RawTagData", recursive = TRUE, full.names = TRUE) |>
  Filter(f = function(d) grepl("/[0-9]{6}$", d), x = _)

location_files <- lapply(whale_dirs, function(dir) {
  whale_id <- basename(dir)

  # FastGPS-combined files: any *-1-Locations.csv not starting with ._
  fastgps <- list.files(dir, pattern = "^[^.][^_].*-1-Locations\\.csv$",
                        full.names = TRUE)
  # Plain Argos-only file
  plain   <- list.files(dir,
                        pattern   = paste0("^", whale_id, "-Locations\\.csv$"),
                        full.names = TRUE)

  if (length(fastgps) > 0) fastgps[1] else if (length(plain) > 0) plain[1] else NULL
}) |> Filter(f = Negate(is.null), x = _) |> unlist()

cat("Found", length(location_files), "location files\n")

locations <- lapply(location_files, function(f) {
  whale_id <- basename(dirname(f))
  year     <- str_extract(f, "(?<=RawTagData/)[0-9]{4}")
  region   <- str_extract(f, "(?<=/)[A-Za-z_]+(?=/[0-9]{6})")

  df <- read_csv(f, show_col_types = FALSE,
                 col_types = cols(.default = col_character()))

  # Normalize column names: lowercase, spaces/special chars → underscore
  names(df) <- gsub("[^a-z0-9]+", "_", tolower(trimws(names(df))))
  names(df) <- gsub("_+$", "", names(df))   # strip trailing underscores

  df$whale_id <- whale_id
  df$year     <- as.integer(year)
  df$region   <- region
  df
}) |> bind_rows()

cat("Loaded", nrow(locations), "raw records for",
    n_distinct(locations$whale_id), "whales\n")

# ─── 2. Assign location class and parse columns ───────────────────────────────
# User type   → "G"  (GPS deployment fix)
# FastGPS type → "G"  (GPS-quality fix; raw quality "4"/"7" is discarded)
# Argos type  → keep existing quality (0, 1, 2, 3, A, B)

locations <- locations |>
  mutate(
    # %OS handles both "HH:MM:SS" and "HH:MM:SS.ffffff" fractional seconds
    date = as.POSIXct(date, format = "%H:%M:%OS %d-%b-%Y", tz = "UTC"),
    lat  = as.numeric(latitude),
    lon  = as.numeric(longitude),
    smaj = as.numeric(error_semi_major_axis),
    smin = as.numeric(error_semi_minor_axis),
    eor  = as.numeric(error_ellipse_orientation),
    lc   = case_when(
      type == "User"    ~ "G",
      type == "FastGPS" ~ "G",
      type == "Argos"   ~ quality,
      TRUE              ~ quality
    )
  )

# ─── 3. Build aniMotum input ──────────────────────────────────────────────────

ssm_input <- locations |>
  rename(id = whale_id) |>
  filter(!is.na(date), !is.na(lat), !is.na(lon), !is.na(lc), lc != "") |>
  select(id, date, lc, lon, lat, smaj, smin, eor, year, region) |>
  arrange(id, date)

# ─── 4. Trim leading/trailing isolated points (gap > 1 day at each end) ───────
# Removes lone observations separated from the main tracking body at start/end.

trim_end_gaps <- function(df, gap_days = 1) {
  # Trim from start: remove first row while gap to next > 1 day
  while (nrow(df) >= 2 &&
         as.numeric(difftime(df$date[2], df$date[1], units = "days")) > gap_days) {
    df <- df[-1, ]
  }
  # Trim from end: remove last row while gap from previous > 1 day
  while (nrow(df) >= 2) {
    n <- nrow(df)
    if (as.numeric(difftime(df$date[n], df$date[n - 1], units = "days")) > gap_days) {
      df <- df[-n, ]
    } else {
      break
    }
  }
  df
}

ssm_input <- ssm_input |>
  group_by(id) |>
  group_modify(~ trim_end_gaps(.x)) |>
  ungroup()

# ─── 5. Exclude whales with < 1 day of data ───────────────────────────────────

durations <- ssm_input |>
  group_by(id) |>
  summarise(
    n_obs    = n(),
    duration = as.numeric(difftime(max(date), min(date), units = "days")),
    year     = first(year),
    region   = first(region),
    .groups  = "drop"
  ) |>
  arrange(year, region, id)

cat("\nWhale data durations (days):\n")
print(durations, n = Inf)

exclude_ids <- durations |> filter(duration < 1) |> pull(id)

if (length(exclude_ids) > 0) {
  cat("\nExcluding", length(exclude_ids), "whale(s) with < 1 day of data:",
      paste(exclude_ids, collapse = ", "), "\n")
  ssm_input <- ssm_input |> filter(!id %in% exclude_ids)
}

cat("\nRunning SSM on", n_distinct(ssm_input$id), "whale(s),",
    nrow(ssm_input), "total observations\n\n")

# ─── 6. Fit aniMotum state-space model ────────────────────────────────────────
# model   = "crw"  → continuous random walk (directional persistence)
# time.step = 6   → regularised predicted locations every 6 hours
# vmax      = 5   → max plausible travel speed (m/s); fin whales ~3–4 m/s

ssm_data <- ssm_input |> select(id, date, lc, lon, lat, smaj, smin, eor)

fit <- fit_ssm(
  ssm_data,
  model     = "crw",
  time.step = 6,
  vmax      = 5,
  control   = ssm_control(verbose = 1)
)

cat("\n── SSM fit summary ──\n")
print(fit)

# ─── 7. Extract results and save ──────────────────────────────────────────────

predicted <- grab(fit, what = "predicted")   # regularised 6-hour locations
fitted    <- grab(fit, what = "fitted")       # model locations at obs times

# Re-attach year and region metadata
meta <- ssm_input |>
  distinct(id, year, region)

predicted <- left_join(predicted, meta, by = "id")
fitted    <- left_join(fitted,    meta, by = "id")

write_csv(predicted, "Results/ssm_predicted_locations.csv")
write_csv(fitted,    "Results/ssm_fitted_locations.csv")

cat("\nResults saved:\n")
cat("  Results/ssm_predicted_locations.csv  (", nrow(predicted), "rows )\n")
cat("  Results/ssm_fitted_locations.csv     (", nrow(fitted),    "rows )\n")
