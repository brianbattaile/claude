library(dplyr)
library(readr)
library(stringr)
library(aniMotum)

# ─── 1. Load location data ────────────────────────────────────────────────────
# For each whale folder, prefer the FastGPS+Argos combined file over the plain
# Argos-only file. Exclude Mac resource fork files (._).
# Cape Cod whales:    use *-1-Locations.csv
# Bay of Fundy whales: use *-2-Locations.csv

whale_dirs <- list.dirs("RawTagData", recursive = TRUE, full.names = TRUE) |>
  Filter(f = function(d) grepl("/[0-9]{6}$", d), x = _)

location_files <- lapply(whale_dirs, function(dir) {
  whale_id <- basename(dir)
  is_bof   <- grepl("Bay_Of_Fundy", dir)

  # FastGPS-combined files: *-2-Locations.csv for BOF, *-1-Locations.csv for Cape Cod
  suffix  <- if (is_bof) "-2-Locations\\.csv$" else "-1-Locations\\.csv$"
  fastgps <- list.files(dir, pattern = paste0("^[^.][^_].*", suffix),
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

# ─── 4. Split at gaps > 24 hours, keep largest segment, save trimmed sections ──
# For each whale, splits the track at every gap > 24 hours. The largest segment
# (by number of observations) is kept for the SSM. All other segments are stored
# in trimmed_segments as named data frames (name = "whaleID_seg#").

trimmed_segments <- list()

split_on_gaps <- function(df, whale_id, gap_hours = 24) {
  gaps    <- as.numeric(difftime(df$date[-1], df$date[-nrow(df)], units = "hours"))
  breaks  <- which(gaps > gap_hours)

  # Build segment index: start/end row for each segment
  starts <- c(1, breaks + 1)
  ends   <- c(breaks, nrow(df))
  segs   <- mapply(function(s, e) df[s:e, ], starts, ends, SIMPLIFY = FALSE)

  # Keep the largest segment (most observations)
  keep_idx  <- which.max(vapply(segs, nrow, integer(1)))
  kept      <- segs[[keep_idx]]
  discarded <- segs[-keep_idx]

  # Store discarded segments in the global trimmed_segments list
  if (length(discarded) > 0) {
    labels <- seq_along(discarded)
    labels <- labels[labels != keep_idx]   # preserve original segment numbering
    for (i in seq_along(discarded)) {
      key <- paste0(whale_id, "_seg", labels[i])
      trimmed_segments[[key]] <<- discarded[[i]]
    }
  }

  kept
}

ssm_input <- ssm_input |>
  group_by(id) |>
  group_modify(~ split_on_gaps(.x, whale_id = .y$id)) |>
  ungroup()

if (length(trimmed_segments) > 0) {
  cat("\nTrimmed segments saved (", length(trimmed_segments), "total ):\n")
  for (nm in names(trimmed_segments)) {
    cat(" ", nm, "—", nrow(trimmed_segments[[nm]]), "observations\n")
  }
} else {
  cat("\nNo segments trimmed (no gaps > 24 hours found)\n")
}

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
  model     = "mp",
  time.step = 6,
  #vmax      = 5,
  control   = ssm_control(verbose = 1)
)

cat("\n── SSM fit summary ──\n")
print(fit)

# ─── 7. Extract results and save ──────────────────────────────────────────────

predicted <- grab(fit, what = "predicted")   # regularised 6-hour locations
fitted    <- grab(fit, what = "fitted")       # model locations at obs times

# Reroute predicted locations around land
fit_rerouted <- route_path(fit, what = "predicted", map_scale = 10, buffer = 10000)
rerouted     <- grab(fit_rerouted, what = "rerouted")

# prt_trim() inside route_path drops predicted locations that fall on land.
# The code below (currently commented out) snaps those points to the nearest
# water location instead of deleting them.
#
# OPTION A — restore at predicted coordinates (on land):
# dropped <- anti_join(
#   predicted |> select(id, date, lon, lat, x, y, x.se, y.se, logit_g, logit_g.se, g),
#   rerouted  |> select(id, date),
#   by = c("id", "date")
# )
# if (nrow(dropped) > 0) {
#   rerouted <- bind_rows(rerouted, dropped) |> arrange(id, date)
# }
#
# OPTION B — snap dropped points to nearest water (requires sf, rnaturalearth):
# dropped_sf <- dropped |>
#   st_as_sf(coords = c("lon","lat"), crs = 4326) |>
#   st_transform(crs = 3857)
# world_mc <- ne_countries(scale = 10, returnclass = "sf") |>
#   st_transform(crs = 3857) |> st_make_valid()
# land_region <- suppressWarnings(
#   st_buffer(dropped_sf, dist = 10000) |> st_union() |> st_convex_hull() |>
#   st_intersection(world_mc) |> st_collection_extract("POLYGON") |> st_sf()
# )
# land_boundary <- st_cast(st_union(land_region), "MULTILINESTRING")
# nudge <- 200  # metres past the coastline into water
# snapped <- lapply(seq_len(nrow(dropped_sf)), function(i) {
#   pt    <- dropped_sf[i, ]
#   nline <- st_nearest_points(pt, land_boundary) |> st_cast("POINT")
#   coast <- nline[2]
#   dist_m <- as.numeric(st_distance(pt, coast))
#   vec    <- st_coordinates(coast) - st_coordinates(pt)
#   norm   <- sqrt(sum(vec^2))
#   st_sfc(st_point(st_coordinates(pt) + vec / norm * (dist_m + nudge)), crs = 3857)
# })
# snapped_sf <- do.call(c, snapped) |> st_sf() |> st_transform(4326)
# dropped$lon <- st_coordinates(snapped_sf)[,1]
# dropped$lat <- st_coordinates(snapped_sf)[,2]
# rerouted <- bind_rows(rerouted, dropped) |> arrange(id, date)

# Re-attach year and region metadata
meta <- ssm_input |>
  distinct(id, year, region)

predicted <- left_join(predicted, meta, by = "id")
fitted    <- left_join(fitted,    meta, by = "id")
rerouted  <- left_join(rerouted,  meta, by = "id")

write_csv(predicted, "Results/ssm_predicted_locations.csv")
write_csv(fitted,    "Results/ssm_fitted_locations.csv")
write_csv(rerouted,  "Results/ssm_rerouted_locations.csv")

cat("\nResults saved:\n")
cat("  Results/ssm_predicted_locations.csv  (", nrow(predicted), "rows )\n")
cat("  Results/ssm_fitted_locations.csv     (", nrow(fitted),    "rows )\n")
cat("  Results/ssm_rerouted_locations.csv   (", nrow(rerouted),  "rows )\n")
