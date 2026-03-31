library(ncdf4)
library(readr)
library(dplyr)

NC_FILE   <- "E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/SSTS_GRDT_combined.nc"
ENV_CSV   <- "Results/ssm_predicted_with_env.csv"
LOCS_CSV  <- "Results/ssm_predicted_locations.csv"
WHALE_ID  <- 234798

# ─── 1. First location for whale 234798 from ssm_predicted_with_env.csv ───────

env <- read_csv(ENV_CSV, show_col_types = FALSE) |>
  mutate(date = as.POSIXct(date, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) |>
  filter(id == WHALE_ID) |>
  arrange(date)

first_loc <- env[1, ]
cat("--- First predicted location for whale", WHALE_ID, "---\n")
cat("  date     :", format(first_loc$date), "\n")
cat("  lon/lat  :", first_loc$lon, first_loc$lat, "\n")
cat("  SSS (03) :", first_loc$SSS, "\n\n")

# ─── 2. Re-extract SSS using the animation script method ──────────────────────

nc <- nc_open(NC_FILE)
nc_lon  <- ncvar_get(nc, "lon")
nc_lat  <- ncvar_get(nc, "lat")
nc_time <- ncvar_get(nc, "time")
time_origin  <- as.POSIXct("2006-01-01 00:00:00", tz = "UTC")
nc_datetimes <- time_origin + nc_time

# Nearest time step
ti <- which.min(abs(as.numeric(nc_datetimes) - as.numeric(first_loc$date)))
cat("  Nearest netCDF time:", format(nc_datetimes[ti]), "\n")

# Nearest grid cell
dist2   <- (nc_lon - first_loc$lon)^2 + (nc_lat - first_loc$lat)^2
nearest <- which(dist2 == min(dist2), arr.ind = TRUE)
xi <- nearest[1]; yi <- nearest[2]

sss_anim <- ncvar_get(nc, "SSS", start = c(xi, yi, ti), count = c(1, 1, 1))
cat("  SSS (anim):", sss_anim, "\n\n")

nc_close(nc)

cat("--- Comparison ---\n")
cat("  SSS from 03_extract_sst_sss.R :", first_loc$SSS, "\n")
cat("  SSS from animation script     :", sss_anim, "\n")
cat("  Difference                    :", abs(first_loc$SSS - sss_anim), "\n")
