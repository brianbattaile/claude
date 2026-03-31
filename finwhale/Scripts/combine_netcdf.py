"""
Combine hourly netCDF files from SSTS_GRDT_22-24/Hourly into a single .nc file.

Each source file has one time step with variables: SST, SSS, SSTG, SSSG
over a 2D grid (Y=106, X=242). Files are concatenated along the time dimension.

Output: E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/SSTS_GRDT_combined.nc
"""

import glob
import os
import numpy as np
import netCDF4 as nc

INPUT_DIR = "E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/Hourly"
OUTPUT_FILE = "E:/Zhitao/From Zhitao/SSTS_GRDT_22-24/SSTS_GRDT_combined.nc"
DATA_VARS = ["SST", "SSS", "SSTG", "SSSG"]

# Collect and sort input files
files = sorted(glob.glob(os.path.join(INPUT_DIR, "SSTS_GRDT_*.nc")))
n_files = len(files)
print(f"Found {n_files} files to combine.")

if n_files == 0:
    raise FileNotFoundError(f"No .nc files found in {INPUT_DIR}")

# Read spatial grid and metadata from the first file
with nc.Dataset(files[0]) as src:
    X = src.variables["X"][:]
    Y = src.variables["Y"][:]
    lon = src.variables["lon"][:]
    lat = src.variables["lat"][:]
    time_units = src.variables["time"].Units
    var_meta = {
        v: {
            "Units": src.variables[v].Units,
            "Long name": getattr(src.variables[v], "Long name", v),
            "dtype": src.variables[v].dtype,
        }
        for v in DATA_VARS
    }

# Create output file
with nc.Dataset(OUTPUT_FILE, "w", format="NETCDF4") as dst:
    # Dimensions
    dst.createDimension("X", len(X))
    dst.createDimension("Y", len(Y))
    dst.createDimension("time", None)  # unlimited

    # Coordinate variables
    vX = dst.createVariable("X", "f8", ("X",))
    vX[:] = X

    vY = dst.createVariable("Y", "f8", ("Y",))
    vY[:] = Y

    vlon = dst.createVariable("lon", "f8", ("Y", "X"))
    vlon[:] = lon
    vlon.long_name = "longitude"
    vlon.units = "degrees_east"

    vlat = dst.createVariable("lat", "f8", ("Y", "X"))
    vlat[:] = lat
    vlat.long_name = "latitude"
    vlat.units = "degrees_north"

    vtime = dst.createVariable("time", "f8", ("time",))
    vtime.Units = time_units

    # Data variables
    dvars = {}
    for v in DATA_VARS:
        dv = dst.createVariable(v, "f8", ("time", "Y", "X"), zlib=True, complevel=4)
        dv.Units = var_meta[v]["Units"]
        dv.setncattr("Long name", var_meta[v]["Long name"])
        dvars[v] = dv

    # Write each file as one time slice
    for i, fpath in enumerate(files):
        if i % 500 == 0:
            print(f"  Processing file {i+1}/{n_files}: {os.path.basename(fpath)}")
        with nc.Dataset(fpath) as src:
            vtime[i] = src.variables["time"][0]
            for v in DATA_VARS:
                dvars[v][i, :, :] = src.variables[v][0, :, :]

print(f"\nDone. Combined file written to:\n  {OUTPUT_FILE}")
