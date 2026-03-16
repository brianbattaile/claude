library(dplyr)
library(readr)
library(stringr)

# Find all Locations CSV files (exclude Mac resource fork files starting with ._)
location_files <- list.files(
  path       = "RawTagData",
  pattern    = "^[0-9]{6}-Locations\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

# Read and combine all location files, adding columns for whale ID, year, and location
locations <- lapply(location_files, function(f) {
  whale_id <- str_extract(basename(f), "^[0-9]{6}")
  year     <- str_extract(f, "(?<=RawTagData/)[0-9]{4}")
  region   <- str_extract(f, "(?<=/)[A-Za-z_]+(?=/[0-9]{6})")

  df <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))
  df$whale_id <- whale_id
  df$year     <- as.integer(year)
  df$region   <- region
  df
}) |> bind_rows()

# Convert numeric columns and parse the Date column to POSIXct
locations <- locations |>
  mutate(
    datetime  = as.POSIXct(Date, format = "%H:%M:%S %d-%b-%Y", tz = "UTC"),
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

cat("Loaded", nrow(locations), "location records for",
    n_distinct(locations$whale_id), "whales\n")

glimpse(locations)
