{
  "name": "Echogram and CTD Interpolation Project",
  "version": "1.0",
  "description": "A project to read and interpolate echogram and CTD data.",
  "dependencies": {
    "oce": "latest",
    "ggplot2": "latest",
    "dplyr": "latest",
    "tidyverse": "latest",
    "viridis": "latest",
    "akima": "latest",
    "interp": "latest",
    "lubridate": "latest"
  },
  "scripts": {
    "setup": "setup.R",
    "data_reading": "data_reading.R",
    "interpolation": "interpolation.R",
    "visualization": "visualization.R"
  },
  "data": {
    "echogram_data": "path/to/echogram_data.csv",
    "ctd_data": "path/to/ctd_data.csv"
  },
  "output": {
    "interpolated_profiles": "output/interpolated_profiles.csv",
    "plots": "output/plots/"
  }
}