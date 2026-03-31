library(ncdf4)
library(ggplot2)
library(dplyr)
library(readr)
library(lubridate)
library(viridis)

# ─── Settings ─────────────────────────────────────────────────────────────────

WHALE_ID     <- 234798
NC_FILE      <- "E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/SSTS_GRDT_combined.nc"
LOCATIONS_IN <- "Results/ssm_predicted_locations.csv"
FRAMES_DIR   <- "Results/movie_frames"
OUTPUT_MP4   <- "Results/sss_whale234798.mp4"
PLOT_VAR     <- "SSS"   # change to "SST", "SSTG", or "SSSG" if desired
FFMPEG       <- "ffmpeg"  # assumes ffmpeg is in PATH; use full path if not

# ─── 1. Load whale locations ──────────────────────────────────────────────────

locs <- read_csv(LOCATIONS_IN, show_col_types = FALSE) |>
  mutate(date = as.POSIXct(date, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) |>
  filter(id == WHALE_ID) |>
  select(id, date, lon, lat)

t_start <- min(locs$date)
t_end   <- max(locs$date)
cat("Whale", WHALE_ID, ":", nrow(locs), "locations\n")
cat("Temporal extent:", format(t_start), "to", format(t_end), "\n\n")

# ─── 2. Open netCDF and subset to whale's time window ─────────────────────────

nc <- nc_open(NC_FILE)

nc_lon  <- ncvar_get(nc, "lon")   # R order: (X=242, Y=106)
nc_lat  <- ncvar_get(nc, "lat")
nc_time <- ncvar_get(nc, "time")
time_origin   <- as.POSIXct("2006-01-01 00:00:00", tz = "UTC")
nc_datetimes  <- time_origin + nc_time

# Keep only time steps within the whale's extent, sampled every 6 hours
ti_keep <- which(nc_datetimes >= t_start & nc_datetimes <= t_end)
ti_keep <- ti_keep[seq(1, length(ti_keep), by = 6)]
cat("netCDF time steps in window:", length(ti_keep), "\n\n")

# Flatten full 2D curvilinear lon/lat arrays — as.vector() order matches slab
lon_flat <- as.vector(nc_lon)
lat_flat <- as.vector(nc_lat)

# Uniform tile size from full grid extent, slightly oversized to close gaps
tile_w <- diff(range(lon_flat)) / 241 * 1.1
tile_h <- diff(range(lat_flat)) / 105 * 1.1

# ─── 3. Plot extent ───────────────────────────────────────────────────────────

lon_range <- range(lon_flat) + c(-0.5, 0.5)
lat_range <- range(lat_flat) + c(-0.5, 0.5)

# ─── 4. Create frames directory ───────────────────────────────────────────────

if (!dir.exists(FRAMES_DIR)) dir.create(FRAMES_DIR, recursive = TRUE)

# ─── 5. Render one PNG frame per netCDF time step ────────────────────────────

# Pre-compute colour scale limits from a sample of frames for consistency
cat("Computing colour scale limits...\n")
sample_ti <- ti_keep[seq(1, length(ti_keep), length.out = min(50, length(ti_keep)))]
sample_vals <- sapply(sample_ti, function(ti) {
  # ncvar_get dim order for 3D var in R: (X, Y, time)
  slab <- ncvar_get(nc, PLOT_VAR, start = c(1, 1, ti), count = c(-1, -1, 1))
  quantile(slab, c(0.02, 0.98), na.rm = TRUE)
})
clim <- c(min(sample_vals[1, ]), max(sample_vals[2, ]))
cat("Colour limits:", round(clim, 2), "\n\n")

var_units <- list(SST = "°C", SSS = "psu", SSTG = "°C / 10 km", SSSG = "psu / 10 km")

cat("Rendering", length(ti_keep), "frames...\n")

for (k in seq_along(ti_keep)) {
  ti <- ti_keep[k]

  if (k %% 100 == 0)
    cat("  frame", k, "/", length(ti_keep), "\n")

  # Read one time slice — dim order in R: (X=242, Y=106)
  slab <- ncvar_get(nc, PLOT_VAR, start = c(1, 1, ti), count = c(-1, -1, 1))

  field_df <- data.frame(lon = lon_flat, lat = lat_flat, val = as.vector(slab))

  frame_time <- nc_datetimes[ti]

  # Nearest whale location at this time
  whale_now <- locs[which.min(abs(as.numeric(locs$date) - as.numeric(frame_time))), ]

  # SSS value at nearest grid point to whale
  dist2 <- (nc_lon - whale_now$lon)^2 + (nc_lat - whale_now$lat)^2
  nearest <- which(dist2 == min(dist2), arr.ind = TRUE)
  whale_now$val <- slab[nearest[1], nearest[2]]
# Trail: last 10 whale positions up to this frame
  trail <- locs |>
    filter(date <= frame_time) |>
    tail(10)

  p <- ggplot() +
    # SSS raster
    geom_tile(data = field_df, aes(x = lon, y = lat, fill = val),
              width = tile_w, height = tile_h, na.rm = TRUE) +
    scale_fill_viridis_c(
      name   = paste0(PLOT_VAR, "\n(", var_units[[PLOT_VAR]], ")"),
      limits = clim,
      oob    = scales::squish,
      option = "plasma",
      na.value = NA
    ) +
    # Whale trail
    geom_path(data = trail, aes(x = lon, y = lat),
              colour = "white", linewidth = 0.6, alpha = 0.7) +
    # Current whale position
    geom_point(data = whale_now, aes(x = lon, y = lat, fill = val),
               shape = 21, colour = "white", size = 4, stroke = 1.2) +
    annotate("text", x = Inf, y = -Inf,
             label = format(frame_time, "%Y-%m-%d %H:%M UTC"),
             hjust = 1.05, vjust = -0.5, size = 3.5, colour = "black") +
    # North up, fixed aspect ratio
    scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
    scale_x_continuous(limits = lon_range, expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(
      title = paste0("Whale ", WHALE_ID, "  |  ", PLOT_VAR),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30"),
      panel.grid    = element_blank(),
      axis.text     = element_text(size = 7)
    )

  frame_file <- file.path(FRAMES_DIR, sprintf("frame_%04d.png", k))
  ggsave(frame_file, p, width = 8, height = 6, dpi = 120)
}

nc_close(nc)
cat("All frames rendered.\n\n")

# ─── 6. Combine frames into MP4 with ffmpeg ───────────────────────────────────

cat("Assembling MP4...\n")
ffmpeg_cmd <- sprintf(
  '%s -y -framerate 12 -i "%s/frame_%%04d.png" -vcodec libx264 -pix_fmt yuv420p "%s"',
  FFMPEG, FRAMES_DIR, OUTPUT_MP4
)
ret <- system(ffmpeg_cmd)

if (ret == 0) {
  cat("Movie saved to:", OUTPUT_MP4, "\n")
} else {
  cat("ffmpeg failed (exit code", ret, "). Frames are in:", FRAMES_DIR, "\n")
  cat("Run manually:\n  ", ffmpeg_cmd, "\n")
}
