{
  "name": "Echogram and CTD Interpolation Project",
  "version": "1.0",
  "description": "A project to read and interpolate binned acoustic backscatter echogram data and CTD data.",
  "main": "main.R",
  "scripts": {
    "install": [
      "install.packages(c('dplyr', 'ggplot2', 'akima', 'oce', 'lubridate'))"
    ],
    "run": "Rscript main.R"
  },
  "dependencies": {
    "dplyr": "^1.0.7",
    "ggplot2": "^3.3.5",
    "akima": "^0.6-2",
    "oce": "^1.5-0",
    "lubridate": "^1.7.10"
  },
  "files": [
    {
      "name": "data",
      "type": "folder",
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
      "name": "R",
      "type": "folder",
      "files": [
        {
          "name": "main.R",
          "type": "file",
          "content": [
            "library(dplyr)",
            "library(ggplot2)",
            "library(akima)",
            "library(oce)",
            "library(lubridate)",
            "",
            "echogram_data <- read.csv('data/echogram_data.csv')",
            "ctd_data <- read.csv('data/ctd_data.csv')",
            "",
            "interpolate_ctd <- function(ctd_df) {",
            "  interp_result <- interp(ctd_df$Int, ctd_df$Depth, ctd_df$Sigma,",
            "                          linear = FALSE, duplicate = 'mean')",
            "  return(interp_result)",
            "}",
            "",
            "transects <- unique(ctd_data$Transect)",
            "for (transect in transects) {",
            "  ctd_subset <- ctd_data %>% filter(Transect == transect)",
            "  interp_result <- interpolate_ctd(ctd_subset)",
            "  interp_df <- as.data.frame(interp2xyz(interp_result))",
            "  ggplot(interp_df, aes(x = Interval, y = Depth, z = Sigma)) +",
            "    geom_raster(aes(fill = Sigma)) +",
            "    scale_y_reverse() +",
            "    labs(title = paste('CTD Interpolation for Transect', transect)) +",
            "    theme_minimal()",
            "}"
          ]
        }
      ]
    }
  ]
}