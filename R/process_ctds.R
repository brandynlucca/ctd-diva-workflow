library(oce)
library(dplyr)
library(ggplot2)
library(lubridate)
library(fields)      # For Tps (thin plate splines)
library(mgcv)        # For GAM-based interpolation
library(sf)          # For spatial operations
library(purrr)
library(akima)       # For fallback interpolation

#' Read CTD files for a transect and extract station information
#' 
#' @param ctd_dir Directory containing CTD files
#' @param pattern File pattern to match CNV files
#' @return List with CTD data and station information
read_ctd_transect <- function(ctd_dir, pattern = "\\.cnv$") {
  # Find all CNV files
  ctd_files <- list.files(ctd_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  if(length(ctd_files) == 0) {
    stop("No CTD files found matching the pattern in directory: ", ctd_dir)
  }
  
  cat("Found", length(ctd_files), "CTD files in", ctd_dir, "\n")
  
  # Read each CTD file using the updated read_ctd_file function
  ctd_list <- lapply(ctd_files, function(file) {
    tryCatch({
      # Read CTD data
      ctd <- read.ctd(file)
      
      # Extract station ID from filename format like "20210525_C35"
      file_name <- basename(file)
      station_pattern <- ".*C(\\d+)"  # Match "C" followed by digits
      station_match <- regexpr(station_pattern, file_name)
      
      if(station_match != -1) {
        # Extract just the digits after "C"
        station_id <- as.integer(gsub(".*C(\\d+).*", "\\1", file_name))
      } else {
        # Fallback if pattern not found
        station_id <- NA
        warning(paste("Could not extract station number from filename:", file_name))
      }
      
      # Extract cruise info from path
      cruise_pattern <- "NYOS\\d{4}"
      cruise <- regmatches(file_name, regexpr(cruise_pattern, file))
      if(length(cruise) == 0) cruise <- NA_character_
      
      # Extract line info if available
      line_pattern <- "Line\\d{2}"
      line_match <- regexpr(line_pattern, file)
      line <- if(line_match != -1) regmatches(file, line_match) else NA_character_
      
      # Create data frame with CTD data
      data.frame(
        Station = station_id,
        Filename = file_name,
        FilePath = file,
        Depth = swDepth(ctd@data$pressure, latitude = ctd@metadata$latitude),
        Pressure = ctd@data$pressure,
        Temperature = ctd@data$temperature,
        Salinity = ctd@data$salinity,
        Sigma_t = swRho(ctd@data$salinity, ctd@data$temperature, ctd@data$pressure) - 1000,
        Latitude = ctd@metadata$latitude,
        Longitude = ctd@metadata$longitude,
        Datetime = as.POSIXct(ctd@metadata$startTime),
        Cruise = cruise,
        Line = line
      )
    }, error = function(e) {
      warning(paste("Error reading file:", file, "-", e$message))
      return(NULL)
    })
  })
  
  # Remove NULL entries (from errors)
  ctd_list <- ctd_list[!sapply(ctd_list, is.null)]
  
  if(length(ctd_list) == 0) {
    stop("No valid CTD data could be read.")
  }
  
  # Combine all CTD data
  ctd_data <- bind_rows(ctd_list)
  
  # Extract station information
  stations <- ctd_data %>%
    group_by(Station, Filename) %>%
    summarize(
      FilePath = first(FilePath),
      Latitude = first(Latitude),
      Longitude = first(Longitude),
      Datetime = first(Datetime),
      Cruise = first(Cruise),
      Line = first(Line),
      .groups = "drop"
    ) %>%
    arrange(Datetime)  # Sort by time
  
  cat("Processed", nrow(stations), "CTD stations\n")
  
  return(list(
    data = ctd_data,
    stations = stations
  ))
}

#' Generate a synthetic transect line based on CTD station positions
#' 
#' @param stations Data frame with station positions
#' @param resolution Number of points to generate along the transect
#' @return Data frame with transect points and distances
generate_transect_line <- function(stations, distance_resolution_km = 1) {
  # Check if we have enough stations
  if(nrow(stations) < 2) {
    stop("Need at least 2 stations to generate a transect")
  }
  
  # Ensure Station column exists
  if(!"Station" %in% names(stations)) {
    warning("No 'Station' column found. Using row numbers.")
    stations$Station <- 1:nrow(stations)
  }
  
  # Sort stations (by time if available, otherwise by station)
  if("Datetime" %in% names(stations)) {
    stations <- stations %>% arrange(Datetime)
  } else {
    stations <- stations %>% arrange(Station)
  }
  
  cat("Generating transect line through", nrow(stations), "stations...\n")
  
  # Calculate total distance
  total_dist_km <- 0
  for(i in 2:nrow(stations)) {
    lat1 <- stations$Latitude[i-1] * pi/180
    lat2 <- stations$Latitude[i] * pi/180
    dlon <- (stations$Longitude[i] - stations$Longitude[i-1]) * pi/180
    dlat <- (stations$Latitude[i] - stations$Latitude[i-1]) * pi/180
    
    a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    dist_km <- 6371 * c  # Earth radius in km
    
    total_dist_km <- total_dist_km + dist_km
  }
  
  # Calculate number of points based on resolution
  # No minimum enforced - use exactly what's requested
  num_points <- ceiling(total_dist_km / distance_resolution_km) + 1
  
  cat("Transect length:", round(total_dist_km, 2), "km\n")
  cat("Generating", num_points, "points with", 
      round(total_dist_km/(num_points-1), 3), "km spacing\n")
  
  # Generate evenly spaced points along the full transect
  # This is the key change - create points along the entire path, not segment by segment
  
  # First, create a cumulative distance for each station
  station_dists <- numeric(nrow(stations))
  for(i in 2:nrow(stations)) {
    lat1 <- stations$Latitude[i-1] * pi/180
    lat2 <- stations$Latitude[i] * pi/180
    dlon <- (stations$Longitude[i] - stations$Longitude[i-1]) * pi/180
    dlat <- (stations$Latitude[i] - stations$Latitude[i-1]) * pi/180
    
    a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    dist_km <- 6371 * c
    
    station_dists[i] <- station_dists[i-1] + dist_km
  }
  
  # Create evenly spaced distances
  even_distances <- seq(0, total_dist_km, length.out = num_points)
  
  # Initialize transect dataframe
  transect <- data.frame(
    TransectPoint = 1:num_points,
    Distance = even_distances,
    Longitude = NA,
    Latitude = NA
  )
  
  # Interpolate lat/lon for each point
  for(i in 1:num_points) {
    # Find which segment this point belongs to
    current_dist <- even_distances[i]
    seg_idx <- max(which(station_dists <= current_dist), 1)
    
    # If it's the last station, special case
    if(seg_idx == nrow(stations)) {
      transect$Longitude[i] <- stations$Longitude[seg_idx]
      transect$Latitude[i] <- stations$Latitude[seg_idx]
      next
    }
    
    # Calculate proportion along segment
    seg_start_dist <- station_dists[seg_idx]
    seg_end_dist <- station_dists[seg_idx + 1]
    seg_prop <- (current_dist - seg_start_dist) / (seg_end_dist - seg_start_dist)
    
    # If proportion is very close to 0 or 1, snap to station
    if(seg_prop < 0.001) {
      transect$Longitude[i] <- stations$Longitude[seg_idx]
      transect$Latitude[i] <- stations$Latitude[seg_idx]
      next
    } else if(seg_prop > 0.999) {
      transect$Longitude[i] <- stations$Longitude[seg_idx + 1]
      transect$Latitude[i] <- stations$Latitude[seg_idx + 1]
      next
    }
    
    # Linear interpolation along great circle path
    # Using the spherical interpolation formula
    lat1 <- stations$Latitude[seg_idx] * pi/180
    lon1 <- stations$Longitude[seg_idx] * pi/180
    lat2 <- stations$Latitude[seg_idx + 1] * pi/180
    lon2 <- stations$Longitude[seg_idx + 1] * pi/180
    
    d <- 2 * asin(sqrt(sin((lat2-lat1)/2)^2 + 
                         cos(lat1) * cos(lat2) * sin((lon2-lon1)/2)^2))
    
    A <- sin((1-seg_prop) * d) / sin(d)
    B <- sin(seg_prop * d) / sin(d)
    
    x <- A * cos(lat1) * cos(lon1) + B * cos(lat2) * cos(lon2)
    y <- A * cos(lat1) * sin(lon1) + B * cos(lat2) * sin(lon2)
    z <- A * sin(lat1) + B * sin(lat2)
    
    lat_interp <- atan2(z, sqrt(x^2 + y^2))
    lon_interp <- atan2(y, x)
    
    transect$Latitude[i] <- lat_interp * 180/pi
    transect$Longitude[i] <- lon_interp * 180/pi
  }
  
  # Calculate station distances along the transect
  stations$Distance <- NA
  stations$OffsetDistance_km <- NA
  
  for(i in 1:nrow(stations)) {
    # Find closest point on the transect
    dists <- geosphere::distGeo(
      cbind(stations$Longitude[i], stations$Latitude[i]),
      cbind(transect$Longitude, transect$Latitude)
    ) / 1000  # Convert to km
    
    closest_idx <- which.min(dists)
    stations$Distance[i] <- transect$Distance[closest_idx]
    stations$OffsetDistance_km[i] <- dists[closest_idx]
  }
  
  # Print station distances for verification
  cat("Station distances along transect (km):\n")
  stations_summary <- stations %>% 
    dplyr::select(Station, Distance, OffsetDistance_km) %>%
    arrange(Distance)
  
  for(i in 1:nrow(stations_summary)) {
    cat("  Station", stations_summary$Station[i], "at", 
        round(stations_summary$Distance[i], 2), 
        "km (offset:", round(stations_summary$OffsetDistance_km[i], 3), "km)\n")
  }
  
  return(list(
    transect = transect,
    stations = stations
  ))
}

#' Find bathymetry along the transect
#' 
#' @param transect Transect line data frame
#' @param bathy Bathymetry data frame with Longitude, Latitude, and Depth columns
#' @param search_radius Radius to search for bathymetry points (km)
#' @return Transect with interpolated bathymetry
find_bathymetry_along_transect <- function(transect, bathy, search_radius = 1) {
  # Convert search radius to degrees (approximate)
  search_radius_deg <- search_radius / 111  # ~111 km per degree
  
  # Find bathymetry for each transect point
  transect$BottomDepth <- NA
  
  for(i in 1:nrow(transect)) {
    # Find nearby bathymetry points
    nearby <- bathy %>%
      filter(
        Longitude >= transect$Longitude[i] - search_radius_deg,
        Longitude <= transect$Longitude[i] + search_radius_deg,
        Latitude >= transect$Latitude[i] - search_radius_deg,
        Latitude <= transect$Latitude[i] + search_radius_deg
      )
    
    if(nrow(nearby) > 0) {
      # Calculate distances
      dists <- geosphere::distGeo(
        c(transect$Longitude[i], transect$Latitude[i]),
        cbind(nearby$Longitude, nearby$Latitude)
      ) / 1000  # Convert to km
      
      # Use inverse distance weighting
      weights <- 1 / (dists^2)
      transect$BottomDepth[i] <- sum(nearby$Depth * weights) / sum(weights)
    }
  }
  
  # Fill any remaining NA values using linear interpolation
  if(any(is.na(transect$BottomDepth))) {
    transect$BottomDepth <- approx(
      x = transect$Distance[!is.na(transect$BottomDepth)],
      y = transect$BottomDepth[!is.na(transect$BottomDepth)],
      xout = transect$Distance,
      rule = 2
    )$y
  }
  
  # Smooth the bathymetry profile
  if(length(unique(transect$BottomDepth)) > 3) {
    bath_smooth <- try({
      smooth_model <- mgcv::gam(BottomDepth ~ s(Distance, bs = "cs"), data = transect)
      predict(smooth_model, transect)
    }, silent = TRUE)
    
    if(!inherits(bath_smooth, "try-error")) {
      transect$BottomDepth <- bath_smooth
    }
  }
  
  return(transect)
}

#' Perform DIVA-like interpolation that respects bathymetry
#' 
#' @param ctd_data CTD data frame
#' @param transect Transect with bathymetry
#' @param variable Variable to interpolate
#' @param nx Number of grid points in x dimension (distance)
#' @param ny Number of grid points in y dimension (depth)
#' @param correlation_x Correlation length in x dimension (km)
#' @param correlation_y Correlation length in y dimension (m)
#' @return List with interpolated data
# Function to interpolate CTD data along transect with bathymetry constraints
# Improved interpolation function to respect CTD depth limits
interpolate_ctd_along_transect <- function(
    ctd_data, 
    transect_with_bathy,
    variables = c("Temperature", "Salinity", "Sigma_t"),
    nx = 100,  # Number of grid points along distance
    ny = 100,  # Number of grid points along depth
    depth_buffer_m = 10  # How far beyond max CTD depth to interpolate
) {
  
  # Extract CTD station data with distances along transect
  ctd_stations <- ctd_data$stations
  
  # Add full CTD profile data with distances
  ctd_profiles <- ctd_data$data %>%
    left_join(dplyr::select(ctd_stations, Station, Distance), by = "Station")
  
  cat("Interpolating CTD data for", paste(variables, collapse = ", "), "\n")
  
  # Find maximum CTD measurement depth for limiting interpolation
  max_ctd_depth <- max(ctd_profiles$Depth, na.rm = TRUE)
  cat("Maximum CTD measurement depth:", round(max_ctd_depth, 1), "m\n")
  
  # Limit interpolation to slightly deeper than deepest measurement
  depth_max <- min(max_ctd_depth + depth_buffer_m, 
                   max(transect_with_bathy$BottomDepth, na.rm = TRUE))
  
  # Setup grid boundaries
  x_range <- range(transect_with_bathy$Distance, na.rm = TRUE)
  y_range <- c(0, depth_max)
  
  x_grid <- seq(x_range[1], x_range[2], length.out = nx)
  y_grid <- seq(y_range[1], y_range[2], length.out = ny)
  
  # Create bathymetry function for masking
  bathy_func <- approxfun(
    x = transect_with_bathy$Distance, 
    y = transect_with_bathy$BottomDepth,
    rule = 2
  )
  
  # Create prediction grid
  pred_grid <- expand.grid(Distance = x_grid, Depth = y_grid)
  
  # Check each grid point against bathymetry AND max CTD depth
  for(i in 1:nrow(pred_grid)) {
    bottom_depth <- bathy_func(pred_grid$Distance[i])
    pred_grid$below_bottom[i] <- pred_grid$Depth[i] > bottom_depth
    pred_grid$beyond_data[i] <- pred_grid$Depth[i] > max_ctd_depth + depth_buffer_m
  }
  
  # Results list for each variable
  results <- list()
  
  # Process each variable
  for(var in variables) {
    if(!var %in% names(ctd_profiles)) {
      warning(paste("Variable", var, "not found in CTD data. Skipping."))
      next
    }
    
    cat("  Interpolating", var, "\n")
    
    # Remove NAs for this variable
    var_data <- ctd_profiles[!is.na(ctd_profiles[[var]]), ]
    
    if(nrow(var_data) < 4) {
      warning(paste("Not enough valid data points for", var, "interpolation"))
      next
    }
    
    # Try GAM interpolation with tensor product smooth
    tryCatch({
      # Fit GAM model - tensor product allows anisotropy
      gam_formula <- as.formula(paste0(var, " ~ te(Distance, Depth, k=c(5, 10))"))
      
      gam_model <- mgcv::gam(
        gam_formula,
        data = var_data,
        method = "REML"
      )
      
      # Predict only for valid grid points (above seafloor AND within CTD depth range)
      valid_grid <- pred_grid[!pred_grid$below_bottom & !pred_grid$beyond_data, ]
      valid_grid$Value <- predict(gam_model, valid_grid)
      
      # Create full results matrix (with NAs below bottom and beyond data)
      result_matrix <- matrix(NA, nrow = length(x_grid), ncol = length(y_grid))
      
      # Fill in predicted values
      for(i in 1:nrow(valid_grid)) {
        x_idx <- which.min(abs(x_grid - valid_grid$Distance[i]))
        y_idx <- which.min(abs(y_grid - valid_grid$Depth[i]))
        result_matrix[x_idx, y_idx] <- valid_grid$Value[i]
      }
      
      # Store results
      results[[var]] <- list(
        x = x_grid,
        y = y_grid,
        z = result_matrix,
        method = "gam",
        min_value = min(valid_grid$Value, na.rm = TRUE),
        max_value = max(valid_grid$Value, na.rm = TRUE)
      )
      
      cat("    Range:", round(min(valid_grid$Value, na.rm = TRUE), 3), "to", 
          round(max(valid_grid$Value, na.rm = TRUE), 3), "\n")
      
    }, error = function(e) {
      warning(paste("GAM interpolation failed for", var, ":", e$message))
      
      # Fallback to Thin Plate Spline interpolation
      tryCatch({
        tps_model <- fields::Tps(
          cbind(var_data$Distance, var_data$Depth),
          var_data[[var]],
          lambda = NULL  # Auto-determine smoothing parameter
        )
        
        # Predict on grid
        valid_grid <- pred_grid[!pred_grid$below_bottom & !pred_grid$beyond_data, ]
        valid_grid$Value <- predict(tps_model, as.matrix(valid_grid[, c("Distance", "Depth")]))
        
        # Create full results matrix
        result_matrix <- matrix(NA, nrow = length(x_grid), ncol = length(y_grid))
        
        # Fill in predicted values
        for(i in 1:nrow(valid_grid)) {
          x_idx <- which.min(abs(x_grid - valid_grid$Distance[i]))
          y_idx <- which.min(abs(y_grid - valid_grid$Depth[i]))
          result_matrix[x_idx, y_idx] <- valid_grid$Value[i]
        }
        
        # Store results
        results[[var]] <- list(
          x = x_grid, 
          y = y_grid,
          z = result_matrix,
          method = "tps",
          min_value = min(valid_grid$Value, na.rm = TRUE),
          max_value = max(valid_grid$Value, na.rm = TRUE)
        )
        
        cat("    Range:", round(min(valid_grid$Value, na.rm = TRUE), 3), "to", 
            round(max(valid_grid$Value, na.rm = TRUE), 3), "\n")
        
      }, error = function(e2) {
        warning(paste("All interpolation methods failed for", var))
        results[[var]] <- NULL
      })
    })
  }
  
  # Return results with bathymetry and grid information
  return(list(
    interpolated = results,
    bathymetry = transect_with_bathy,
    stations = ctd_stations,
    grid = list(x = x_grid, y = y_grid),
    max_ctd_depth = max_ctd_depth
  ))
}

interpolate_ctd_with_diva <- function(
    ctd_data, 
    transect_with_bathy,
    variables = c("Temperature", "Salinity", "Sigma_t"),
    nx = 100,
    ny = 100,
    len_x = 40,   # Correlation length in x direction (km)
    len_z = 20,   # Correlation length in z direction (m)
    epsilon = 1   # Signal-to-noise ratio
) {
  # Extract CTD data
  ctd_stations <- ctd_data$stations
  ctd_profiles <- ctd_data$data %>%
    left_join(dplyr::select(ctd_stations, Station, Distance), by = "Station")
  
  # Setup grid
  x_range <- range(transect_with_bathy$Distance, na.rm = TRUE)
  y_range <- c(0, max(transect_with_bathy$BottomDepth, na.rm = TRUE))
  
  x_grid <- seq(x_range[1], x_range[2], length.out = nx)
  y_grid <- seq(y_range[1], y_range[2], length.out = ny)
  
  # Create mask for bathymetry
  mask <- matrix(TRUE, nrow = nx, ncol = ny)
  for(i in 1:nx) {
    bottom_depth <- approx(
      x = transect_with_bathy$Distance, 
      y = transect_with_bathy$BottomDepth,
      xout = x_grid[i],
      rule = 2
    )$y
    
    mask[i, y_grid > bottom_depth] <- FALSE
  }
  
  # Results list
  results <- list()
  
  # Process each variable
  for(var in variables) {
    if(!var %in% names(ctd_profiles)) {
      warning(paste("Variable", var, "not found in CTD data. Skipping."))
      next
    }
    
    cat("  Interpolating", var, "with DIVA\n")
    
    # Prepare data
    var_data <- ctd_profiles[!is.na(ctd_profiles[[var]]), ]
    
    if(nrow(var_data) < 4) {
      warning(paste("Not enough valid data points for", var, "interpolation"))
      next
    }
    
    # Data points for DIVA
    obs_x <- var_data$Distance
    obs_y <- var_data$Depth
    obs_f <- var_data[[var]]
    
    # Implement DIVA-based interpolation
    # This uses optimal interpolation with anisotropic covariance
    result_matrix <- matrix(NA, nrow = nx, ncol = ny)
    
    # Calculate weights for each grid point
    for(i in 1:nx) {
      for(j in 1:ny) {
        if(mask[i, j]) {  # Only process points above seafloor
          # Grid point coordinates
          x0 <- x_grid[i]
          y0 <- y_grid[j]
          
          # Skip if no nearby observations
          if(length(obs_x) == 0) next
          
          # Distance matrix between observations
          n_obs <- length(obs_x)
          dist_xx <- matrix(0, n_obs, n_obs)
          
          for(ii in 1:n_obs) {
            for(jj in 1:n_obs) {
              # Anisotropic distance metric
              dx <- (obs_x[ii] - obs_x[jj]) / len_x
              dy <- (obs_y[ii] - obs_y[jj]) / len_z
              dist_xx[ii, jj] <- sqrt(dx^2 + dy^2)
            }
          }
          
          # Covariance matrix between observations
          cov_xx <- exp(-dist_xx) + diag(1/epsilon, n_obs)
          
          # Distance vector between grid point and observations
          dist_x0 <- numeric(n_obs)
          for(ii in 1:n_obs) {
            dx <- (x0 - obs_x[ii]) / len_x
            dy <- (y0 - obs_y[ii]) / len_z
            dist_x0[ii] <- sqrt(dx^2 + dy^2)
          }
          
          # Covariance vector
          cov_x0 <- exp(-dist_x0)
          
          # Calculate optimal weights
          weights <- try(solve(cov_xx, cov_x0), silent = TRUE)
          
          # Calculate interpolated value
          if(!inherits(weights, "try-error")) {
            result_matrix[i, j] <- sum(weights * obs_f)
          }
        }
      }
    }
    
    # Store results
    results[[var]] <- list(
      x = x_grid,
      y = y_grid,
      z = result_matrix,
      method = "diva"
    )
    
    # Output value range to verify reasonability
    cat("    Range:", round(min(result_matrix, na.rm = TRUE), 3), "to", 
        round(max(result_matrix, na.rm = TRUE), 3), "\n")
  }
  
  # Return results
  return(list(
    interpolated = results,
    bathymetry = transect_with_bathy,
    stations = ctd_stations,
    grid = list(x = x_grid, y = y_grid),
    mask = mask
  ))
}

interpolate_ctd_with_diva_fast <- function(
    ctd_data, 
    transect_with_bathy,
    variables = c("Temperature", "Salinity", "Sigma_t"),
    nx = 100,
    ny = 100,
    len_x = 40,   # Correlation length in x direction (km)
    len_z = 20,   # Correlation length in z direction (m)
    epsilon = 1,  # Signal-to-noise ratio
    local_radius = 3.0,  # Use only observations within this many correlation lengths
    use_parallel = TRUE  # Use parallel processing if available
) {
  # Set up parallel processing if requested
  if(use_parallel && requireNamespace("parallel", quietly = TRUE)) {
    num_cores <- max(1, parallel::detectCores() - 1)
    cat("Using parallel processing with", num_cores, "cores\n")
    cl <- parallel::makeCluster(num_cores)
    
    # Fix: Properly export the variables to the cluster
    parallel::clusterExport(cl, varlist = c("len_x", "len_z", "epsilon", "local_radius"), 
                            envir = environment())
    
    on.exit(parallel::stopCluster(cl))
  } else {
    use_parallel <- FALSE
  }
  
  # Extract CTD data
  ctd_stations <- ctd_data$stations
  ctd_profiles <- ctd_data$data %>%
    left_join(dplyr::select(ctd_stations, Station, Distance), by = "Station")
  
  # Setup grid
  x_range <- range(transect_with_bathy$Distance, na.rm = TRUE)
  y_range <- c(0, max(transect_with_bathy$BottomDepth, na.rm = TRUE))
  
  x_grid <- seq(x_range[1], x_range[2], length.out = nx)
  y_grid <- seq(y_range[1], y_range[2], length.out = ny)
  
  # Create bathymetry mask matrix - much faster as a vectorized operation
  grid_x_mat <- matrix(rep(x_grid, each = ny), nrow = nx, ncol = ny)
  grid_y_mat <- matrix(rep(y_grid, nx), nrow = nx, ncol = ny, byrow = TRUE)
  
  # Get bottom depth at each x position in a vectorized way
  bottom_depths <- approx(
    x = transect_with_bathy$Distance, 
    y = transect_with_bathy$BottomDepth,
    xout = x_grid,
    rule = 2
  )$y
  
  # Create mask matrix
  mask <- matrix(FALSE, nrow = nx, ncol = ny)
  for(i in 1:nx) {
    mask[i, y_grid <= bottom_depths[i]] <- TRUE
  }
  
  # Results list
  results <- list()
  
  # Process each variable
  for(var in variables) {
    if(!var %in% names(ctd_profiles)) {
      warning(paste("Variable", var, "not found in CTD data. Skipping."))
      next
    }
    
    cat("  Interpolating", var, "with optimized DIVA\n")
    
    # Prepare data
    var_data <- ctd_profiles[!is.na(ctd_profiles[[var]]), ]
    
    if(nrow(var_data) < 4) {
      warning(paste("Not enough valid data points for", var, "interpolation"))
      next
    }
    
    # Data points 
    obs_x <- var_data$Distance
    obs_y <- var_data$Depth
    obs_f <- var_data[[var]]
    n_obs <- length(obs_x)
    
    # Pre-compute the covariance matrix between observations - this is a major speedup
    dist_xx <- matrix(0, n_obs, n_obs)
    for(i in 1:n_obs) {
      for(j in i:n_obs) {
        dx <- (obs_x[i] - obs_x[j]) / len_x
        dy <- (obs_y[i] - obs_y[j]) / len_z
        dist_xx[i, j] <- dist_xx[j, i] <- sqrt(dx^2 + dy^2)
      }
    }
    cov_xx <- exp(-dist_xx) + diag(1/epsilon, n_obs)
    
    # Pre-compute the inverse of the covariance matrix - huge speedup
    cov_xx_inv <- try(solve(cov_xx), silent = TRUE)
    if(inherits(cov_xx_inv, "try-error")) {
      warning("Covariance matrix inversion failed. Using regularization.")
      # Add regularization
      cov_xx_inv <- solve(cov_xx + diag(1e-8, n_obs))
    }
    
    # Function to process a single grid point - used for both serial and parallel
    process_grid_point <- function(idx) {
      i <- (idx - 1) %/% ny + 1
      j <- (idx - 1) %% ny + 1
      
      if(!mask[i,j]) return(NA)
      
      x0 <- x_grid[i]
      y0 <- y_grid[j]
      
      # Use local observations only - major speedup for large datasets
      dx <- (obs_x - x0) / len_x
      dy <- (obs_y - y0) / len_z
      dist_x0 <- sqrt(dx^2 + dy^2)
      
      # Only use observations within local_radius correlation lengths
      local_idx <- which(dist_x0 <= local_radius)
      
      if(length(local_idx) < 3) {
        # Not enough nearby observations
        return(NA)
      }
      
      # Use only local observations - huge speedup
      local_dist <- dist_x0[local_idx]
      local_cov <- exp(-local_dist)
      
      # Extract relevant parts of the inverse covariance matrix
      local_inv_cov <- cov_xx_inv[local_idx, local_idx]
      
      # Calculate optimal weights
      weights <- local_inv_cov %*% local_cov
      
      # Calculate interpolated value
      return(sum(weights * obs_f[local_idx]))
    }
    
    # Create result matrix
    result_matrix <- matrix(NA, nrow = nx, ncol = ny)
    
    # Process all grid points - either in parallel or serial
    valid_idx <- which(mask)
    if(use_parallel) {
      # Parallel processing
      results_list <- parallel::parLapply(cl, valid_idx, process_grid_point)
      for(k in seq_along(valid_idx)) {
        idx <- valid_idx[k]
        i <- (idx - 1) %/% ny + 1
        j <- (idx - 1) %% ny + 1
        result_matrix[i,j] <- results_list[[k]]
      }
    } else {
      # Serial processing with progress indicator
      pb <- txtProgressBar(min = 0, max = length(valid_idx), style = 3)
      for(k in seq_along(valid_idx)) {
        idx <- valid_idx[k]
        i <- (idx - 1) %/% ny + 1
        j <- (idx - 1) %% ny + 1
        result_matrix[i,j] <- process_grid_point(idx)
        setTxtProgressBar(pb, k)
      }
      close(pb)
    }
    
    # Store results
    results[[var]] <- list(
      x = x_grid,
      y = y_grid,
      z = result_matrix,
      method = "diva_fast"
    )
    
    # Output value range to verify reasonability
    cat("    Range:", round(min(result_matrix, na.rm = TRUE), 3), "to", 
        round(max(result_matrix, na.rm = TRUE), 3), "\n")
  }
  
  # Return results
  return(list(
    interpolated = results,
    bathymetry = transect_with_bathy,
    stations = ctd_stations,
    grid = list(x = x_grid, y = y_grid),
    mask = mask
  ))
}

# Akima interpolation for oceanographic data
interpolate_ctd_with_akima <- function(
    ctd_data, 
    transect_with_bathy,
    variables = c("Temperature", "Salinity", "Sigma_t"),
    nx = 100,
    ny = 100
) {
  # Extract CTD data
  ctd_stations <- ctd_data$stations
  ctd_profiles <- ctd_data$data %>%
    left_join(dplyr::select(ctd_stations, Station, Distance), by = "Station")
  
  # Setup grid
  x_range <- range(transect_with_bathy$Distance, na.rm = TRUE)
  y_range <- c(0, max(transect_with_bathy$BottomDepth, na.rm = TRUE))
  
  x_grid <- seq(x_range[1], x_range[2], length.out = nx)
  y_grid <- seq(y_range[1], y_range[2], length.out = ny)
  
  # Create mask for bathymetry
  mask <- matrix(TRUE, nrow = nx, ncol = ny)
  for(i in 1:nx) {
    bottom_depth <- approx(
      x = transect_with_bathy$Distance, 
      y = transect_with_bathy$BottomDepth,
      xout = x_grid[i],
      rule = 2
    )$y
    
    mask[i, y_grid > bottom_depth] <- FALSE
  }
  
  # Results list
  results <- list()
  
  # Process each variable
  for(var in variables) {
    if(!var %in% names(ctd_profiles)) {
      warning(paste("Variable", var, "not found in CTD data. Skipping."))
      next
    }
    
    cat("  Interpolating", var, "with Akima\n")
    
    # Prepare data
    var_data <- ctd_profiles[!is.na(ctd_profiles[[var]]), ]
    
    if(nrow(var_data) < 4) {
      warning(paste("Not enough valid data points for", var, "interpolation"))
      next
    }
    
    # Use Akima interpolation
    if(!requireNamespace("akima", quietly = TRUE)) {
      stop("Package 'akima' is required for interpolation. Install with install.packages('akima')")
    }
    
    # Interpolate with Akima
    akima_result <- akima::interp(
      x = var_data$Distance,
      y = var_data$Depth,
      z = var_data[[var]],
      xo = x_grid,
      yo = y_grid,
      linear = FALSE  # Use bivariate spline
    )
    
    # Apply bathymetry mask
    akima_result$z[!mask] <- NA
    
    # Store results
    results[[var]] <- list(
      x = x_grid,
      y = y_grid,
      z = akima_result$z,
      method = "akima"
    )
    
    # Output value range to verify reasonability
    cat("    Range:", round(min(akima_result$z, na.rm = TRUE), 3), "to", 
        round(max(akima_result$z, na.rm = TRUE), 3), "\n")
  }
  
  # Return results
  return(list(
    interpolated = results,
    bathymetry = transect_with_bathy,
    stations = ctd_stations,
    grid = list(x = x_grid, y = y_grid),
    mask = mask
  ))
}

#' Plot interpolated section with bathymetry
#' 
#' @param interp_result Interpolation result from diva_like_interpolation
#' @param transect Transect with bathymetry
#' @param stations Station information
#' @param variable Name of the variable (for labeling)
#' @param custom_colors Optional color palette to use
#' @return ggplot object
# Function to plot an interpolated variable
# Improved plotting function with variable-specific color scales
# Clean plotting function without seafloor extension
plot_interpolated_section <- function(
    interp_result, 
    variable = "Temperature",
    plot_title = NULL,
    label_contours = TRUE
) {
  # Check if variable exists
  if(!variable %in% names(interp_result$interpolated)) {
    stop(paste("No interpolated data for variable:", variable))
  }
  
  var_result <- interp_result$interpolated[[variable]]
  bathymetry <- interp_result$bathymetry
  stations <- interp_result$stations
  
  # Convert matrix to data frame for ggplot
  grid_data <- expand.grid(
    Distance = var_result$x,
    Depth = var_result$y
  )
  grid_data$Value <- as.vector(var_result$z)
  
  # Set up basic title and legend
  if(is.null(plot_title)) {
    line <- unique(stations$Line)
    plot_title <- paste(variable, "Section", line)
  }
  
  # Choose palette based on variable
  if(variable == "Temperature") {
    palette_option <- "plasma"
    legend_title <- "Temperature (°C)"
  } else if(variable == "Salinity") {
    palette_option <- "viridis" 
    legend_title <- "Salinity (PSU)"
  } else if(variable == "Sigma_t") {
    palette_option <- "magma"
    legend_title <- expression(sigma[t]~"(kg/m"^3*")")
  } else {
    palette_option <- "turbo"
    legend_title <- variable
  }
  
  # Create contour breaks
  value_range <- range(grid_data$Value, na.rm = TRUE)
  contour_breaks <- pretty(value_range, n = 8)
  
  # Create plot
  p <- ggplot() +
    # Data fill
    geom_raster(
      data = grid_data,
      aes(x = Distance, y = Depth, fill = Value),
      na.rm = TRUE
    ) +
    # Contour lines
    geom_contour(
      data = grid_data,
      aes(x = Distance, y = Depth, z = Value),
      breaks = contour_breaks,
      color = "black", alpha = 0.7,
      na.rm = TRUE
    )
  
  # Add contour labels if requested and metR is available
  if(label_contours && requireNamespace("metR", quietly = TRUE)) {
    p <- p + metR::geom_text_contour(
      data = grid_data,
      aes(x = Distance, y = Depth, z = Value, label = ..level..),
      breaks = contour_breaks,
      stroke = 0.2,
      size = 3,
      skip = 1,
      check_overlap = TRUE
    )
  }
  
  # Complete the plot
  p <- p +
    # Add bathymetry line
    geom_line(
      data = bathymetry,
      aes(x = Distance, y = BottomDepth),
      color = "black", linewidth = 1.2
    ) +
    # Gray polygon for seafloor
    geom_polygon(
      data = data.frame(
        x = c(min(bathymetry$Distance), bathymetry$Distance, max(bathymetry$Distance)),
        y = c(max(bathymetry$BottomDepth)*1.5, bathymetry$BottomDepth, max(bathymetry$BottomDepth)*1.5)
      ),
      aes(x = x, y = y),
      fill = "gray20"
    ) +
    # Add stations
    geom_point(
      data = stations,
      aes(x = Distance, y = 0),
      shape = 25, fill = "yellow", color = "black", size = 3
    ) +
    geom_text(
      data = stations,
      aes(x = Distance, y = 0, label = Station),
      vjust = 2, color = "white", fontface = "bold"
    ) +
    # Set scales and theme
    scale_fill_viridis_c(option = palette_option) +
    scale_y_reverse(expand = c(0, 0)) +
    coord_cartesian(ylim = c(max(var_result$y) * 1.05, 0)) +
    theme_bw() +
    labs(
      title = plot_title,
      x = "Distance (km)",
      y = "Depth (m)",
      fill = legend_title
    )
  
  return(p)
}

# Comprehensive wrapper function
run_ctd_transect_interpolation <- function(
    ctd_dir,                   # Directory with CTD files
    bathy_data,                # Bathymetry data (tibble with Longitude, Latitude, Depth)
    distance_resolution_km = 1,# Transect resolution
    buffer_km = 5,             # Bathymetry search buffer
    variables = c("Temperature", "Salinity", "Sigma_t"),  # Variables to interpolate
    nx = 100,                  # Interpolation grid x resolution
    ny = 100,                  # Interpolation grid y resolution
    depth_buffer_m = 10,       # Buffer below max CTD depth
    save_plots = FALSE,        # Whether to save plot files
    save_data = TRUE,          # Whether to save intermediate data
    output_dir = NULL          # Directory to save outputs (default: current dir)
) {
  
  # Set output directory if needed
  if((save_plots || save_data) && is.null(output_dir)) {
    output_dir <- getwd()
    cat("Using current directory for outputs:", output_dir, "\n")
  }
  
  # Create output directory if it doesn't exist
  if((save_plots || save_data) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Extract line name from directory path
  line_name <- basename(ctd_dir)
  if(!grepl("Line", line_name)) {
    # If directory doesn't contain "Line", try to extract from path
    line_parts <- regmatches(ctd_dir, regexpr("Line\\d+", ctd_dir))
    line_name <- if(length(line_parts) > 0) line_parts[1] else "Transect"
  }
  
  # 1. Load CTD data
  cat("Reading CTD data from", ctd_dir, "\n")
  ctd_result <- read_ctd_transect(ctd_dir)
  
  if(save_data) {
    ctd_file <- file.path(output_dir, paste0(line_name, "_ctd_data.rds"))
    cat("Saving CTD data to", ctd_file, "\n")
    saveRDS(ctd_result, ctd_file)
  }
  
  # 2. Generate transect line
  cat("Generating transect line with", distance_resolution_km, "km resolution\n")
  transect_result <- generate_transect_line(
    ctd_result$stations, 
    distance_resolution_km = distance_resolution_km
  )
  
  if(save_data) {
    transect_file <- file.path(output_dir, paste0(line_name, "_transect.rds"))
    cat("Saving transect data to", transect_file, "\n")
    saveRDS(transect_result, transect_file)
  }
  
  # 3. Filter bathymetry to area around transect
  cat("Filtering bathymetry data within", buffer_km, "km of transect\n")
  filtered_bathy <- filter_bathy_for_transect(
    bathy_data = bathy_data, 
    transect_df = transect_result$transect,
    buffer_km = buffer_km
  )
  
  if(save_data) {
    bathy_file <- file.path(output_dir, paste0(line_name, "_filtered_bathy.rds"))
    cat("Saving filtered bathymetry to", bathy_file, "\n")
    saveRDS(filtered_bathy, bathy_file)
  }
  
  # 4. Extract bathymetry along transect
  cat("Extracting bathymetry along transect line\n")
  transect_with_bathy <- extract_bathy_for_transect(
    transect_df = transect_result$transect,
    filtered_bathy = filtered_bathy
  )
  
  if(save_data) {
    transect_bathy_file <- file.path(output_dir, paste0(line_name, "_transect_with_bathy.rds"))
    cat("Saving transect with bathymetry to", transect_bathy_file, "\n")
    saveRDS(transect_with_bathy, transect_bathy_file)
  }
  
  # 5. Add distance information to CTD stations
  cat("Connecting CTD stations with transect distances\n")
  ctd_result$stations <- ctd_result$stations %>%
    left_join(
      dplyr::select(transect_result$stations, Station, Distance, OffsetDistance_km),
      by = "Station"
    )
  
  if(save_data) {
    ctd_with_dist_file <- file.path(output_dir, paste0(line_name, "_ctd_with_distances.rds"))
    cat("Saving CTD data with distances to", ctd_with_dist_file, "\n")
    saveRDS(ctd_result, ctd_with_dist_file)
  }
  
  # 6. Run interpolation
  cat("Interpolating variables:", paste(variables, collapse=", "), "\n")
  if (interp_method=="diva"){
    interp_result <- interpolate_ctd_with_diva(
      ctd_data = ctd_result,
      transect_with_bathy = transect_with_bathy,
      variables = variables,
      nx = nx,
      ny = ny
    )
 
  }else if (interp_method=="akima"){
    interp_result <- interpolate_ctd_with_akima(
      ctd_data = ctd_result,
      transect_with_bathy = transect_with_bathy,
      variables = variables,
      nx = nx,
      ny = ny
    )
    
  }else{
    interp_result <- interpolate_ctd_along_transect(
      ctd_data = ctd_result,
      transect_with_bathy = transect_with_bathy,
      variables = variables,
      nx = nx,
      ny = ny,
      depth_buffer_m = depth_buffer_m
    )       
  }

  
  if(save_data) {
    interp_file <- file.path(output_dir, paste0(line_name, "_interpolation.rds"))
    cat("Saving interpolation results to", interp_file, "\n")
    saveRDS(interp_result, interp_file)
  }
  
  # 7. Create plots
  plots <- list()
  for(var in variables) {
    if(!var %in% names(interp_result$interpolated)) {
      cat("Warning: Variable", var, "not found in interpolation results. Skipping plot.\n")
      next
    }
    
    cat("Creating plot for", var, "\n")
    plots[[var]] <- plot_interpolated_section(interp_result, var)
    
    # Save plot if requested
    if(save_plots) {
      filename <- file.path(output_dir, paste0(line_name, "_", var, "_section.png"))
      cat("Saving plot to", filename, "\n")
      ggsave(filename, plots[[var]], width = 10, height = 6, dpi = 300)
    }
  }
  
  # Save combined plot data if requested
  if(save_data) {
    plots_file <- file.path(output_dir, paste0(line_name, "_plots.rds"))
    cat("Saving plot objects to", plots_file, "\n")
    saveRDS(plots, plots_file)
  }
  
  # Save summary report
  if(save_data) {
    # Create a summary data frame
    summary_data <- data.frame(
      Line = line_name,
      NumStations = nrow(ctd_result$stations),
      TransectLength_km = max(transect_result$transect$Distance),
      MaxDepth_m = max(transect_with_bathy$BottomDepth, na.rm=TRUE),
      Variables = paste(intersect(variables, names(interp_result$interpolated)), collapse=", "),
      DateTime = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    
    summary_file <- file.path(output_dir, paste0(line_name, "_summary.csv"))
    cat("Saving summary data to", summary_file, "\n")
    write.csv(summary_data, summary_file, row.names = FALSE)
  }
  
  # 8. Return results
  return(list(
    plots = plots,
    interpolation = interp_result,
    transect = transect_with_bathy,
    ctd = ctd_result,
    filtered_bathymetry = filtered_bathy,
    output_dir = if(save_data || save_plots) output_dir else NULL
  ))
}

# Filter bathymetry to just the area around the transect
filter_bathy_for_transect <- function(bathy_data, transect_df, buffer_km = 5) {
  # Convert buffer from km to approximate degrees
  buffer_deg <- buffer_km / 111  # ~111 km per degree
  
  # Get transect bounds with buffer
  lon_min <- min(transect_df$Longitude) - buffer_deg
  lon_max <- max(transect_df$Longitude) + buffer_deg
  lat_min <- min(transect_df$Latitude) - buffer_deg
  lat_max <- max(transect_df$Latitude) + buffer_deg
  
  # Filter the data
  filtered_bathy <- bathy_data %>%
    filter(
      Longitude >= lon_min,
      Longitude <= lon_max,
      Latitude >= lat_min,
      Latitude <= lat_max
    )
  
  cat("Filtered from", nrow(bathy_data), "to", nrow(filtered_bathy), "bathymetry points\n")
  return(filtered_bathy)
}

# Extract bathymetry for transect points using the filtered data
extract_bathy_for_transect <- function(transect_df, filtered_bathy, search_radius_km = 0.5) {
  transect_df$BottomDepth <- NA
  
  # Convert search radius from km to degrees (approximate)
  search_radius_deg <- search_radius_km / 111
  
  # For each transect point, find nearest bathymetry
  for(i in 1:nrow(transect_df)) {
    # Find nearby points
    nearby <- filtered_bathy %>%
      filter(
        Longitude >= transect_df$Longitude[i] - search_radius_deg,
        Longitude <= transect_df$Longitude[i] + search_radius_deg,
        Latitude >= transect_df$Latitude[i] - search_radius_deg,
        Latitude <= transect_df$Latitude[i] + search_radius_deg
      )
    
    if(nrow(nearby) > 0) {
      # Calculate distances
      dists <- sqrt(
        ((nearby$Longitude - transect_df$Longitude[i]) * 111 * cos(transect_df$Latitude[i] * pi/180))^2 +
          ((nearby$Latitude - transect_df$Latitude[i]) * 111)^2
      )
      
      # Use nearest points within radius
      valid <- dists <= search_radius_km
      if(sum(valid) > 0) {
        # Inverse distance weighting
        weights <- 1 / (dists[valid]^2)
        transect_df$BottomDepth[i] <- sum(nearby$Depth[valid] * weights) / sum(weights)
      }
    }
  }
  
  # Fill any gaps
  if(any(is.na(transect_df$BottomDepth))) {
    transect_df$BottomDepth <- approx(
      x = transect_df$Distance[!is.na(transect_df$BottomDepth)],
      y = transect_df$BottomDepth[!is.na(transect_df$BottomDepth)],
      xout = transect_df$Distance,
      rule = 2
    )$y
  }
  
  return(transect_df)
}

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Source the script
# source("c:/Users/Brandyn/OneDrive - UW/echogram-ctd-interpolation/diva_like_interpolation.R")

# Path to your CTD files for a specific transect/line
ctd_dir <- "C:/Users/Brandyn/OneDrive - UW/NYOS_all_CTD/NYOS2105/Line06"

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Source the script
# source("c:/Users/Brandyn/OneDrive - UW/echogram-ctd-interpolation/diva_like_interpolation.R")

# Path to your CTD files for a specific transect/line
ctd_dir <- "C:/Users/Brandyn Lucca/OneDrive - UW/NYOS_all_CTD/NYOS2105/Line06"

# Load only the needed data - first extract CTD positions
ctd_result <- read_ctd_transect(ctd_dir)
# Use a distance resolution of 0.1 km (100 meters)
transect_result <- generate_transect_line(
  ctd_result$stations, 
  distance_resolution_km = 0.1 # 100 meter spacing between points
)

# Read bathymetry data (example format - adjust as needed for your actual data)
# This should be a tibble/dataframe with Longitude, Latitude, and Depth columns
# Replace this with your actual bathymetry data loading
# bathy_data <- read.csv("path/to/your/bathymetry/data.csv")
bathy_data <- NE_bathy_tib %>% 
  rename("Longitude"="long", "Latitude"="lat", "Depth"="depth") %>% 
  mutate(Depth=abs(Depth))

# Assuming your bathymetry data is in the variable 'bathymetry_data'
# and has columns Longitude, Latitude, and Depth

# 1. Filter bathymetry to just the area around your transect
filtered_bathy <- filter_bathy_for_transect(
  bathy_data = bathy_data, 
  transect_df = transect_result$transect,
  buffer_km = 5
)

# 2. Extract bathymetry for your transect points
transect_with_bathy <- extract_bathy_for_transect(
  transect_df = transect_result$transect,
  filtered_bathy = filtered_bathy
)

transect_result$stations %>% 
  left_join(ctd_result$data) %>% 
  ggplot() + 
  geom_point(
    mapping = aes(
      x=Distance,
      y=Depth,
      color=Sigma_t
    )
  ) +
  geom_path(
    mapping=aes(x=Distance, y=BottomDepth),
    data=transect_with_bathy
  ) +
  scale_y_reverse() +
  coord_cartesian(ylim=c(325, 0), expand=0,
                  xlim=c(-5, 205)) +
  scico::scale_color_scico(
    palette="lapaz",
    end=0.8,
    limits=c(22, 28),
    oob=scales::squish,
    breaks=seq(22, 28, 2)
  ) +
  labs(
    x="Distance (km)",
    y="Depth (m)",
    color=expression(sigma[t]~"(kg/m"^3*")")
  ) +
  theme_bw() +
  theme(
    text=element_text(color="black", size=15.5),
    axis.text=element_text(color="black", size=15.0),
    panel.grid=element_blank()
  ) 


bathy_data %>% 
  dplyr::filter(
    between(Latitude, 38.7, 40.6) & 
      between(Longitude, -73.5, -72.7)
  ) -> bathy_subset
# 1153 x 351
transect_result$stations %>% 
  ggplot() +
  geom_contour(
    mapping=aes(x=Longitude, y=Latitude, z=Depth, color="A"),
    data=bathy_subset,
    breaks=seq(0,1400,50)
  ) +
  scale_color_manual(
    name=NULL,
    labels=c("A"="Isobath (50m res.)"),
    values=c("A"="gray80")
  ) +
  ggnewscale::new_scale_color() +
  geom_path(
    mapping=aes(x=Longitude, y=Latitude),
    data=transect_result$transect
  ) +
  geom_point(
    mapping=aes(x=Longitude, 
                y=Latitude,
                color=msigma),
    data=ctd_result$data %>% 
      group_by(Station, Longitude, Latitude) %>% 
      reframe(msigma=mean(Sigma_t)),
    size=5,
    show.legend=F
  ) +
  coord_cartesian(ylim=c(38.7, 40.6),
                  xlim=c(-73.5, -72.7),
                  expand=0) +
  scico::scale_color_scico(
    palette="lapaz",
    end=0.8,
    limits=c(22, 28),
    oob=scales::squish,
    breaks=seq(22, 28, 2)
  ) +
  labs(
    x=expression(Longitude~(degree*E)),
    y=expression(Latitude~(degree*N)),
    color=NA
  ) +
  theme_bw() +
  theme(
    text=element_text(color="black", size=15.5),
    axis.text=element_text(color="black", size=15.0),
    panel.grid=element_blank(),
    legend.position="top"
  )


readr::read_csv(
  "C:/Users/Brandyn Lucca/OneDrive - UW/echogram-ctd-interpolation/output/interpolated_sigma_t.csv"
) -> intermediate_dat

intermediate_dat %>% 
  ggplot() +
  geom_raster(
    mapping=aes(x=Distance*1e-3, y=Depth, fill=Value)
  ) +
  geom_path(
    mapping=aes(x=Distance, y=BottomDepth),
    data=transect_with_bathy
  ) +
  geom_contour(
    aes(x = Distance*1e-3, y = Depth, z = Value),
    breaks = c(22, 24, 26, 28),
    color = "black", alpha = 0.7,
    na.rm = TRUE
  ) +
  metR::geom_text_contour(
    mapping=aes(x=Distance*1e-3, y=Depth, z=Value, label=after_stat(level)),
    breaks=c(22, 24, 26, 28),
    stroke=0.2,
    size=3,
    skip=0,
    check_overlap=T
  ) +
  scale_y_reverse() +
  coord_cartesian(ylim=c(325, 0), expand=0,
                  xlim=c(-5, 205)) +
  scico::scale_fill_scico(
    palette="lapaz",
    end=0.8,
    limits=c(22, 28),
    oob=scales::squish,
    breaks=seq(22, 28, 2)
  ) +
  labs(
    x="Distance (km)",
    y="Depth (m)",
    fill=expression(sigma[t]~"(kg/m"^3*")")
  ) +
  theme_bw() +
  theme(
    text=element_text(color="black", size=15.5),
    axis.text=element_text(color="black", size=15.0),
    panel.grid=element_blank()
  ) 



readr::write_csv(
  x = transect_with_bathy,
  file = "C:/Users/Brandyn Lucca/OneDrive - UW/echogram-ctd-interpolation/data/test_transect_bathy.csv"
)

ctd_result$stations <- ctd_result$stations %>%
  left_join(
    dplyr::select(transect_result$stations, Station, Distance, OffsetDistance_km),
    by = "Station"
  )

ctd_result$data

readr::write_csv(
  x = ctd_result$data,
  file = "C:/Users/Brandyn Lucca/OneDrive - UW/echogram-ctd-interpolation/data/test_hydrography.csv"
)

# 3. Run the interpolation with the transect that has bathymetry
interp_result <- interpolate_ctd_along_transect(
  ctd_data = transect_result,
  transect_with_bathy = transect_with_bathy,
  variables = c("Temperature", "Salinity", "Sigma_t"),
  nx = 100,  # Points along distance axis
  ny = 100   # Points along depth axis
)

# 4. Create plots for each variable
temp_plot <- plot_interpolated_section(interp_result, "Temperature")
salt_plot <- plot_interpolated_section(interp_result, "Salinity")
density_plot <- plot_interpolated_section(interp_result, "Sigma_t")

# Display the plots
print(temp_plot)
print(salt_plot)
print(density_plot)

# If you have echogram bottom data, you can use that instead:
# Extract bottom from echogram data
bathy_from_echo <- echo6 %>%
  mutate(Datetime = ymd_hms(paste0(Date_M, " ", Time_M)),
         Depth = Layer * 5) %>%
  group_by(Int) %>%
  summarize(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    Depth = max(Depth, na.rm = TRUE),  # Assuming max depth = bottom
    .groups = "drop"
  )

# Run the interpolation
result <- run_ctd_transect_interpolation(
  ctd_dir = ctd_dir,
  bathy_data = bathy_from_echo,  # Or use your bathymetry tibble
  variables = c("Temperature", "Salinity", "Sigma_t"),
  resolution = c(100, 100),
  output_dir = "C:/Users/Brandyn/OneDrive - UW/echogram-ctd-interpolation/output"
)

# Using optimized DIVA interpolation
diva_result <- interpolate_ctd_with_diva_fast(
  ctd_data = ctd_result,
  transect_with_bathy = transect_with_bathy,
  variables = c("Temperature", "Salinity", "Sigma_t"),
  use_parallel = TRUE
)

max(diva_result$interpolated$Temperature$z)

akima_result <- interpolate_ctd_with_akima(
  ctd_data = ctd_result,
  transect_with_bathy = transect_with_bathy,
  variables = variables,
  nx = 100,
  ny = 100
)

# Plot the results
temp_plot_diva <- plot_interpolated_section(diva_result, "Temperature")

results <- run_ctd_transect_interpolation(
  ctd_dir = "C:/Users/Brandyn/OneDrive - UW/NYOS_all_CTD/NYOS2105/Line06",
  bathy_data = NE_bathy_tib %>% 
    rename("Longitude"="long", "Latitude"="lat", "Depth"="depth") %>% 
    mutate(Depth=abs(Depth)),
  distance_resolution_km = 0.5,  # Higher resolution transect
  buffer_km = 10,                # Wider bathymetry search
  variables = c("Temperature", "Salinity", "Sigma_t", "Oxygen"),  # Add more variables
  nx = 200,                      # Higher resolution interpolation
  ny = 150,
  save_plots = TRUE,
  output_dir = "C:/Users/Brandyn/OneDrive - UW/echogram-ctd-interpolation/output"
)

# View the temperature plot
result$plots$Temperature

# Access the interpolated data
temp_interp <- result$interpolation$Temperature

# Read bathymetry data (example format - adjust as needed for your actual data)
# This should be a tibble/dataframe with Longitude, Latitude, and Depth columns
# Replace this with your actual bathymetry data loading
# bathy_data <- read.csv("path/to/your/bathymetry/data.csv")
bathy_data <- NE_bathy_tib %>% 
  rename("Latitude"="lat", "Longitude"="long", "Depth"="depth") %>% 
  mutate(Depth=abs(Depth))

# If you have echogram bottom data, you can use that instead:
# Extract bottom from echogram data
bathy_from_echo <- echo6 %>%
  mutate(Datetime = ymd_hms(paste0(Date_M, " ", Time_M)),
         Depth = Layer * 5) %>%
  group_by(Int) %>%
  summarize(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    Depth = max(Depth, na.rm = TRUE),  # Assuming max depth = bottom
    .groups = "drop"
  )

# Run the interpolation
results <- run_ctd_transect_interpolation(
  ctd_dir = "C:/Users/Brandyn/OneDrive - UW/NYOS_all_CTD/NYOS2207/Line06",
  bathy_data = NE_bathy_tib %>% 
    rename("Longitude"="long", "Latitude"="lat", "Depth"="depth") %>% 
    mutate(Depth=abs(Depth)),
  distance_resolution_km = 1,
  buffer_km = 5,
  variables = c("Temperature", "Salinity", "Sigma_t"),
  save_plots = TRUE,
  save_data = TRUE,
  output_dir = "C:/Users/Brandyn/OneDrive - UW/echogram-ctd-interpolation/output/NYOS2207/Line06"
)

results$ctd$data %>% 
  as_tibble() %>% reframe(U=max(Temperature))

# View the temperature plot
results$plots$Temperature
results$interpolation$interpolated$Temperature

# Access the interpolated data
temp_interp <- result$interpolation$Temperature

