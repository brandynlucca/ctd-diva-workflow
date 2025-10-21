{
  "name": "Echogram and CTD Interpolation Project",
  "version": "1.0",
  "description": "A project to read and interpolate echogram and CTD data.",
  "main": "main.R",
  "scripts": {
    "install": [
      "install.packages(c('oce', 'ggplot2', 'dplyr', 'tidyverse', 'viridis', 'akima', 'interp'))"
    ],
    "run": "Rscript main.R"
  },
  "dependencies": {
    "R": ">= 4.0.0"
  },
  "files": [
    {
      "name": "data",
      "type": "directory",
      "files": [
        {
          "name": "echogram_data.csv",
          "type": "file"
        },
        {
          "name": "ctd_data.csv",
          "type": "file"
        }
      ]
    },
    {
      "name": "src",
      "type": "directory",
      "files": [
        {
          "name": "main.R",
          "type": "file"
        },
        {
          "name": "interpolation.R",
          "type": "file"
        }
      ]
    }
  ]
}