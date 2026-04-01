# Finwhale Project

Analysis of New England Fin Whale satellite tag data.

## Repository
- GitHub: https://github.com/brianbattaile/claude.git
- Local root: `C:/Users/Green Sturgeon/claude/`

## Environment
- **R**: `C:/Program Files/R/R-4.4.2/bin/Rscript` — NOT in bash PATH, always use full path
- **Python**: 3.11.5, available in PATH
- **Platform**: Windows 11, shell is bash (use Unix-style paths and syntax)

## Project Structure
```
claude/
└── finwhale/
    ├── RawTagData/          # Satellite tag data organized by year and location
    │   ├── 2022/Cape_Cod/   # Tag IDs: 234798-234805
    │   ├── 2023/Cape_Cod/   # Tag IDs: 243052-243065
    │   └── 2024/            # Cape_Cod and Bay_Of_Fundy deployments
    │       ├── Cape_Cod/    # Tag IDs: 243053-243067
    │       └── Bay_Of_Fundy/ # Tag IDs: 267425-267426
    └── Rasters/
        └── R_Script_7/      # 30arcsecGEBCO.tif (bathymetry raster)
```

## Tag Data Files
Each tag folder contains CSVs exported from Wildlife Computers portal:
- `*-Locations.csv` — Argos/GPS location data
- `*-Behavior.csv` — dive behavior summary
- `*-Series.csv` — time-series depth/temperature data
- `*-Histos.csv` — histogram data
- `*-SST.csv` — sea surface temperature
- `*-All.csv` — combined data
- `*-Argos.csv`, `*-RawArgos.csv` — raw Argos messages
- `*-Status.csv`, `*-Summary.csv`, `*-MinMaxDepth.csv`, `*-RTC.csv`

## Git Workflow
- Always commit from `C:/Users/Green Sturgeon/claude/` (the repo root)
- Push to `origin main` after committing
