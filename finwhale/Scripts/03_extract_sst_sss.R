library(ncdf4)
library(readr)
library(dplyr)

# ─── Paths ────────────────────────────────────────────────────────────────────

NC_FILE      <- "E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/SSTS_GRDT_combined.nc"
LOCATIONS_IN <- "Results/ssm_predicted_locations.csv"
OUTPUT_CSV   <- "Results/ssm_predicted_with_env.csv"

# ─── 1. Load whale predicted locations ────────────────────────────────────────

locs <- read_csv(LOCATIONS_IN, show_col_types = FALSE) |>
  mutate(date = as.POSIXct(date, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) |>
  select(id, date, lon, lat, x, y, g, year = year.x, region = region.x)

cat("Loaded", nrow(locs), "whale locations for", n_distinct(locs$id), "whales\n")

# ─── 2. Open netCDF and read coordinate arrays ────────────────────────────────

nc <- nc_open(NC_FILE)

nc_lon  <- ncvar_get(nc, "lon")   # [Y, X]
nc_lat  <- ncvar_get(nc, "lat")   # [Y, X]
nc_time <- ncvar_get(nc, "time")  # seconds since 2006-01-01 00:00:00

# Convert netCDF time to POSIXct
time_origin <- as.POSIXct("2006-01-01 00:00:00", tz = "UTC")
nc_datetimes <- time_origin + nc_time  # POSIXct vector, length = 7098

cat("netCDF time range:", format(min(nc_datetimes)), "to", format(max(nc_datetimes)), "\n")
cat("netCDF grid size: Y =", nrow(nc_lon), ", X =", ncol(nc_lon), "\n\n")

# ─── 3. For each whale location, find nearest time and grid cell ───────────────
# Strategy:
#   - Time: nearest hourly step (snap to closest index in nc_datetimes)
#   - Space: nearest grid cell by minimum great-circle-ish distance (lon/lat degrees)
#     Using Euclidean distance in lon/lat is fine at this scale.

find_nearest_cell <- function(whale_lon, whale_lat) {
  dist2 <- (nc_lon - whale_lon)^2 + (nc_lat - whale_lat)^2
  idx   <- which(dist2 == min(dist2), arr.ind = TRUE)[1, ]
  # R/ncdf4 reverses file dim order: file (Y=106, X=242) -> R array (X=242, Y=106)
  # so idx[1] = X index, idx[2] = Y index
  list(xi = idx[1], yi = idx[2])
}

find_nearest_time <- function(whale_time) {
  which.min(abs(as.numeric(nc_datetimes) - as.numeric(whale_time)))
}

# ─── 4. Extract values ────────────────────────────────────────────────────────

n <- nrow(locs)
out_sst  <- numeric(n)
out_sss  <- numeric(n)
out_sstg <- numeric(n)
out_sssg <- numeric(n)
out_nc_time <- as.POSIXct(rep(NA, n), tz = "UTC")
out_nc_lon  <- numeric(n)
out_nc_lat  <- numeric(n)

cat("Extracting environmental data for", n, "locations...\n")

for (i in seq_len(n)) {
  if (i %% 200 == 0) cat("  ", i, "/", n, "\n")

  ti <- find_nearest_time(locs$date[i])
  cell <- find_nearest_cell(locs$lon[i], locs$lat[i])
  yi <- cell$yi
  xi <- cell$xi

  # ncvar_get: start = c(xi, yi, ti), count = c(1, 1, 1)
  # Note: ncdf4 uses (X, Y, time) order for 3D vars stored as (time, Y, X) in file
  # Check dim order from nc object
  out_sst[i]  <- ncvar_get(nc, "SST",  start = c(xi, yi, ti), count = c(1, 1, 1))
  out_sss[i]  <- ncvar_get(nc, "SSS",  start = c(xi, yi, ti), count = c(1, 1, 1))
  out_sstg[i] <- ncvar_get(nc, "SSTG", start = c(xi, yi, ti), count = c(1, 1, 1))
  out_sssg[i] <- ncvar_get(nc, "SSSG", start = c(xi, yi, ti), count = c(1, 1, 1))

  out_nc_time[i] <- nc_datetimes[ti]
  out_nc_lon[i]  <- nc_lon[xi, yi]
  out_nc_lat[i]  <- nc_lat[xi, yi]
}

nc_close(nc)
cat("Done extracting.\n\n")

# ─── 5. Save results ──────────────────────────────────────────────────────────

results <- locs |>
  mutate(
    nc_time = out_nc_time,
    nc_lon  = out_nc_lon,
    nc_lat  = out_nc_lat,
    SST     = out_sst,
    SSS     = out_sss,
    SSTG    = out_sstg,
    SSSG    = out_sssg
  )

write_csv(results, OUTPUT_CSV)

cat("Results saved to:", OUTPUT_CSV, "\n")
cat("Rows:", nrow(results), "\n")
cat("\nSummary of extracted values:\n")
print(summary(results[, c("SST", "SSS", "SSTG", "SSSG")]))
