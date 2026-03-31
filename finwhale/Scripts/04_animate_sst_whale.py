"""
Animate SST (or any variable) from the combined netCDF over a map,
overlaying predicted locations for all whales in a given year.

Frames are rendered in parallel using multiprocessing, then assembled
into an MP4 with ffmpeg.
"""

import os
import sys
import subprocess
import multiprocessing as mp

import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# ─── Settings ─────────────────────────────────────────────────────────────────

YEAR         = int(sys.argv[1]) if len(sys.argv) > 1 else 2022
PLOT_VAR     = sys.argv[2].upper() if len(sys.argv) > 2 else "SST"   # SST, SSS, SSTG, SSSG
NC_FILE      = "E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/SSTS_GRDT_combined.nc"
LOCATIONS_IN = "Results/ssm_predicted_locations.csv"
FRAMES_DIR   = f"Results/movie_frames_py_{YEAR}_{PLOT_VAR}"
OUTPUT_MP4   = f"Results/{PLOT_VAR.lower()}_{YEAR}_whales_py.mp4"
N_WORKERS    = max(1, mp.cpu_count() - 1)
FPS          = 24
TRAIL_STEPS  = 10
COLORMAP     = "RdYlBu_r"  # blue=cold, red=hot
N_TEST_FRAMES = None       # set to an int to limit frames, None for all
ZOOM_FRAMES   = 240        # frames over which to zoom from full extent to whale extent

VAR_LABEL = {
    "SST":  "SST (°C)",
    "SSS":  "SSS (psu)",
    "SSTG": "SST Gradient (°C / 10 km)",
    "SSSG": "SSS Gradient (psu / 10 km)",
}

# ─── Frame rendering function (must be at module level for pickling) ───────────

def render_frame(args):
    k, frame_time, slab, locs, nc_lon, nc_lat, clim, extent, in_zoom = args
    lon_min, lon_max, lat_min, lat_max = extent

    frame_file = os.path.join(FRAMES_DIR, f"frame_{k:04d}.png")
    if os.path.exists(frame_file):
        return  # resume support

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    im = ax.pcolormesh(
        nc_lon, nc_lat, slab,
        cmap=COLORMAP, vmin=clim[0], vmax=clim[1],
        shading="auto", transform=ccrs.PlateCarree()
    )

    ax.add_feature(cfeature.LAND, facecolor="lightgrey")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.STATES, linewidth=0.4, edgecolor="grey")
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)

    # Plot each whale active near this frame time
    for whale_id, whale_locs in locs.groupby("id"):
        dt = np.abs((whale_locs["date"] - frame_time).dt.total_seconds())

        # Skip whale if nearest location is more than 24 hours away
        if dt.min() > 24 * 3600:
            continue

        nearest_idx = dt.idxmin()
        whale_now = whale_locs.loc[nearest_idx]

        # SST at nearest grid point
        dist2 = (nc_lon - whale_now["lon"])**2 + (nc_lat - whale_now["lat"])**2
        nearest_cell = np.unravel_index(dist2.argmin(), dist2.shape)
        whale_sst = float(slab[nearest_cell])

        # Trail
        trail = whale_locs[whale_locs["date"] <= frame_time].tail(TRAIL_STEPS)
        if len(trail) > 1:
            ax.plot(trail["lon"].values, trail["lat"].values,
                    color="black", linewidth=0.8, alpha=0.5,
                    transform=ccrs.PlateCarree(), zorder=4)

        # Whale position — SST fill, black border
        ax.scatter(whale_now["lon"], whale_now["lat"],
                   c=[whale_sst], cmap=COLORMAP, vmin=clim[0], vmax=clim[1],
                   edgecolors="black", linewidths=1.2, s=80, zorder=5,
                   transform=ccrs.PlateCarree())

    # Colorbar
    cb = fig.colorbar(im, ax=ax, orientation="vertical", pad=0.02, shrink=0.85)
    cb.set_label(VAR_LABEL[PLOT_VAR], fontsize=10)

    ax.set_title(f"{YEAR} Fin Whales  |  {PLOT_VAR}", fontsize=13, fontweight="bold")

    # Date/time in lower right
    ax.text(0.99, 0.01, frame_time.strftime("%Y-%m-%d %H:%M UTC"),
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, color="black")

    # Suppress gridline labels during zoom to avoid jitter from repositioning labels
    gl = ax.gridlines(draw_labels=not in_zoom, linewidth=0.3, color="grey", alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    fig.savefig(frame_file, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ─── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    # 1. Load all whale locations for the given year
    print(f"Loading whale locations for {YEAR}...")
    locs = pd.read_csv(LOCATIONS_IN, parse_dates=["date"])
    locs["date"] = pd.to_datetime(locs["date"], utc=True)
    locs = locs[locs["date"].dt.year == YEAR][["id", "date", "lon", "lat"]].reset_index(drop=True)
    whale_ids = sorted(locs["id"].unique())
    t_start, t_end = locs["date"].min(), locs["date"].max()
    print(f"  {len(whale_ids)} whales, {len(locs)} locations  {t_start:%Y-%m-%d} to {t_end:%Y-%m-%d}")
    print(f"  Whale IDs: {whale_ids}")

    # 2. Open netCDF
    print("Loading netCDF metadata...")
    ds = nc.Dataset(NC_FILE)
    nc_lon  = ds.variables["lon"][:]
    nc_lat  = ds.variables["lat"][:]
    nc_time = ds.variables["time"][:]
    time_origin  = pd.Timestamp("2006-01-01", tz="UTC")
    nc_datetimes = pd.to_datetime([time_origin + pd.Timedelta(seconds=float(s)) for s in nc_time])

    # Subset to time window covering all whales, every 6 hours
    mask    = (nc_datetimes >= t_start) & (nc_datetimes <= t_end)
    ti_keep = np.where(mask)[0][::6]
    times_keep = nc_datetimes[ti_keep]
    print(f"netCDF frames in window: {len(ti_keep)}")

    # Apply test limit if set
    n_frames = len(ti_keep) if N_TEST_FRAMES is None else min(N_TEST_FRAMES, len(ti_keep))
    ti_use = ti_keep[:n_frames]

    # 3. Read only the frames we need
    print(f"Reading {PLOT_VAR} data for {n_frames} frame(s)...")
    raw = ds.variables[PLOT_VAR][ti_use, :, :]
    if isinstance(raw, np.ma.MaskedArray):
        data = raw.filled(np.nan).astype(float)
    else:
        data = np.array(raw, dtype=float)
    data[~np.isfinite(data)] = np.nan
    data[(data < -2) | (data > 50)] = np.nan
    ds.close()
    print("Data loaded.\n")

    clim = np.nanpercentile(data, [2, 98])
    print(f"Colour limits: {clim[0]:.2f} – {clim[1]:.2f}")

    # Full extent (the whole netCDF grid)
    full_extent = (float(nc_lon.min()) - 0.5, float(nc_lon.max()) + 0.5,
                   float(nc_lat.min()) - 0.5, float(nc_lat.max()) + 0.5)

    # Target extent: 2x the whale bounding box, centred on the whales
    wlon_c = (locs["lon"].min() + locs["lon"].max()) / 2
    wlat_c = (locs["lat"].min() + locs["lat"].max()) / 2
    wlon_h = (locs["lon"].max() - locs["lon"].min())   # full span → use as half-width for 2x
    wlat_h = (locs["lat"].max() - locs["lat"].min())
    target_extent = (wlon_c - wlon_h, wlon_c + wlon_h,
                     wlat_c - wlat_h, wlat_c + wlat_h)

    # Per-frame extents: smoothstep interpolation over first ZOOM_FRAMES frames
    def smoothstep(t):
        t = max(0.0, min(1.0, t))
        return t * t * t * (t * (t * 6 - 15) + 10)  # smootherstep: zero 1st & 2nd derivatives at endpoints

    def interp_extent(k):
        t = smoothstep(k / max(ZOOM_FRAMES - 1, 1))
        return tuple(f + t * (tgt - f) for f, tgt in zip(full_extent, target_extent))

    frame_extents = [interp_extent(min(k, ZOOM_FRAMES - 1)) for k in range(n_frames)]

    os.makedirs(FRAMES_DIR, exist_ok=True)

    # Pre-download all cartopy shapefiles in the main process before spawning workers
    # (avoids race condition where multiple workers corrupt the same download)
    print("Pre-loading cartopy map data...")
    for scale in ["110m", "50m"]:
        shpreader.natural_earth(resolution=scale, category="physical",  name="land")
        shpreader.natural_earth(resolution=scale, category="physical",  name="coastline")
        shpreader.natural_earth(resolution=scale, category="cultural",  name="admin_1_states_provinces_lakes")
        shpreader.natural_earth(resolution=scale, category="cultural",  name="admin_0_boundary_lines_land")
    print("Map data ready.\n")

    # 4. Render frames
    args_list = [
        (k + 1, times_keep[k], data[k], locs, nc_lon, nc_lat, clim, frame_extents[k], k < ZOOM_FRAMES)
        for k in range(n_frames)
    ]

    print(f"Rendering {n_frames} frames using {min(N_WORKERS, n_frames)} workers...")
    if n_frames == 1:
        render_frame(args_list[0])
    else:
        with mp.Pool(processes=N_WORKERS) as pool:
            for i, _ in enumerate(pool.imap_unordered(render_frame, args_list, chunksize=10)):
                if (i + 1) % 50 == 0:
                    print(f"  {i+1} / {n_frames} frames done")

    print("All frames rendered.\n")

    # 5. Assemble MP4
    print("Assembling MP4...")
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(FPS),
        "-i", os.path.join(FRAMES_DIR, "frame_%04d.png"),
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-vcodec", "libx264",
        "-pix_fmt", "yuv420p",
        OUTPUT_MP4,
    ]
    ret = subprocess.run(cmd, capture_output=True, text=True)
    if ret.returncode == 0:
        print(f"Movie saved to: {OUTPUT_MP4}")
        # Clean up individual frames now that the movie is assembled
        for f in os.listdir(FRAMES_DIR):
            if f.endswith(".png"):
                os.remove(os.path.join(FRAMES_DIR, f))
        os.rmdir(FRAMES_DIR)
        print("Frames deleted.")
    else:
        print(f"ffmpeg failed:\n{ret.stderr}")
        print(f"Frames are in: {FRAMES_DIR}")
