module DIVAInterp

using CSV, DataFrames, DIVAnd, Measures, Statistics, Interpolations, Plots, Random
export interpolate_field, optimize_diva_parameters, prepare_observations, prepare_diva_grid, visualize_field, process_diva_variable

"""
    prepare_observations(hydro_df, variable="Temperature", km_to_m=true, verbose=true)

Extract and process observation data from hydrography DataFrame.

# Arguments
- `hydro_file::String`: Hydrography data with columns 'Distance', 'Depth', and variable
- `variable::String="Temperature"`: Name of the variable to interpolate
- `km_to_m::Bool`: Convert distance units from km to m 
- `verbose::Bool`: Whether to print progress messages

# Returns
- `NamedTuple`: Contains processed observation data:
  - `x::Vector{Float64}`: X-coordinates of observations in meters
  - `y::Vector{Float64}`: Y-coordinates of observations in meters
  - `x_obs::Tuple{Vector,Vector}`: Coordinates tuple for DIVAnd
  - `values::Vector{Float64}`: Original observation values
  - `anomaly::Vector{Float64}`: Observation values with mean subtracted
  - `background::Float64`: Mean value of observations
  - `typical_station_spacing::Float64`: Median horizontal distance between stations

# Notes
- Automatically removes rows with missing values in the specified variable
- Calculates background field as the mean of all observations
- Computes anomalies by subtracting background from observations
- Determines typical station spacing for correlation length estimation
"""
function prepare_observations(hydro_file, variable="Temperature"; km_to_m=true, verbose=true)
    
    # Load CSV file
    hydro_df = CSV.read(hydro_file, DataFrame)

    # Convert distances, if needed
    if km_to_m
        hydro_df.Distance = hydro_df.Distance * 1e3
    end

    # Remove missing values
    valid_data = dropmissing(hydro_df, Symbol(variable))
    
    # Extract coordinates and values
    obs_x = valid_data.Distance
    obs_y = valid_data.Depth
    obs_values = valid_data[!, Symbol(variable)]
    
    # Calculate typical station spacing
    station_distances = [];
    unique_stations = unique(obs_x);
    sort!(unique_stations);
    for i in 2:lastindex(unique_stations)
        push!(station_distances, unique_stations[i] - unique_stations[i-1])
    end
    typical_station_spacing = median(station_distances)
    
    # Calculate background field
    background = mean(obs_values)
    anomaly = obs_values .- background

    # Calculate maximum station depths
    deepest_stations = combine(groupby(hydro_df, :Station), sdf -> sdf[argmax(sdf.Depth), :])

    # Verbosity
    if verbose
        println("$variable data from hydrographic file ($hydro_file) was successfully ingested and pre-processed!")
    end
    
    return (x=obs_x, y=obs_y, x_obs=(obs_x, obs_y), 
            values=obs_values, anomaly=anomaly, 
            background=background, typical_station_spacing=typical_station_spacing, station_depths=deepest_stations)
end

"""
    prepare_diva_grid(bathy_file, min_x, max_x, horizontal_res, vertical_res, km_to_m=true, verbose=true)

Prepare a bathymetry-aware grid for DIVAnd interpolation.

# Arguments
- `bathy_file::String`: Bathymetry data with columns 'Distance' and 'BottomDepth'
- `min_x::Float64`: Minimum x-coordinate value in meters
- `max_x::Float64`: Maximum x-coordinate value in meters
- `horizontal_res::Float64`: Horizontal resolution of the grid in meters
- `vertical_res::Float64`: Vertical resolution of the grid in meters
- `km_to_m::Bool`: Convert km to m 
- `verbose::Bool`: Whether to print progress messages

# Returns
- `NamedTuple`: Contains all grid elements:
  - `x_grid::StepRangeLen`: X-coordinates of grid points
  - `y_grid::StepRangeLen`: Y-coordinates of grid points
  - `X_grid::Matrix{Float64}`: 2D matrix of X-coordinates
  - `Y_grid::Matrix{Float64}`: 2D matrix of Y-coordinates
  - `xi::Tuple{Matrix,Matrix}`: Tuple of coordinate matrices for DIVAnd
  - `pmn::Tuple{Matrix,Matrix}`: Scale factors for DIVAnd
  - `mask::BitMatrix`: Boolean mask for valid points (true above seafloor)
  - `n_horizontal::Int`: Number of horizontal grid points
  - `n_vertical::Int`: Number of vertical grid points
  - `depth_extent::Float64`: Maximum depth in meters
  - `distance_extent::Float64`: Horizontal extent in meters

# Notes
- Creates a mask where points below the seafloor are marked as false
- Grid starts at depth=0 (surface) and extends to maximum bottom depth
- Uses linear interpolation of bathymetry to determine seafloor location
"""
function prepare_diva_grid(bathy_file, min_x, max_x, horizontal_res, vertical_res; km_to_m=true, verbose=true)
    # Load CSV
    bathy_df = CSV.read(bathy_file, DataFrame)

    # Convert distances, if needed
    if km_to_m
        bathy_df.Distance = bathy_df.Distance * 1e3
    end

    # Calculate grid dimensions
    distance_extent = max_x - min_x
    depth_extent = maximum(bathy_df.BottomDepth)
    
    n_horizontal = Int(ceil(distance_extent / horizontal_res))
    n_vertical = Int(ceil(depth_extent / vertical_res))
    
    # Create grid
    x_regular = range(min_x, max_x, length=n_horizontal)
    y_regular = range(0, depth_extent, length=n_vertical)
    
    # Coordinate grids
    X_grid = [x for x in x_regular, y in y_regular]
    Y_grid = [y for x in x_regular, y in y_regular]
    xi = (X_grid, Y_grid)
    
    # Scale factors (inverse of grid resolution)
    pm = ones(length(x_regular), length(y_regular)) / horizontal_res
    pn = ones(length(x_regular), length(y_regular)) / vertical_res
    pmn = (pm, pn)
    
    # Create mask based on bathymetry
    mask = trues(length(x_regular), length(y_regular));
    for (i, x) in enumerate(x_regular)
        bathy_depth = LinearInterpolation(bathy_df.Distance, 
                                        bathy_df.BottomDepth, 
                                        extrapolation_bc=Flat())(x);
        for (j, y) in enumerate(y_regular)
            if y > bathy_depth
                mask[i, j] = false;
            end
        end
    end

    # Verbosity
    if verbose
        println("The bathymetric file ($bathy_file) was successfully ingested and pre-processed!")
    end
    
    return (x_grid=x_regular, y_grid=y_regular, 
            X_grid=X_grid, Y_grid=Y_grid, 
            xi=xi, pmn=pmn, mask=mask, 
            n_horizontal=n_horizontal, n_vertical=n_vertical,
            depth_extent=depth_extent, 
            distance_extent=distance_extent,
            bottom_depth=bathy_df)
end

"""
    optimize_diva_parameters(mask, pmn, xi, x_obs, obs_values, typical_station_spacing, depth_extent)

Optimize DIVAnd parameters (correlation lengths and noise-to-signal ratio) using cross-validation.

# Arguments
- `mask::BitMatrix`: Boolean mask of valid grid points (true = valid, false = invalid)
- `pmn::Tuple{Matrix, Matrix}`: Scale factors for X and Y dimensions
- `xi::Tuple{Matrix, Matrix}`: Grid coordinates in X and Y dimensions
- `x_obs::Tuple{Vector, Vector}`: Observation coordinates in X and Y dimensions
- `obs_values::Vector`: Observation values (should be anomalies from mean)
- `typical_station_spacing::Float64`: Median distance between observation stations in meters
- `depth_extent::Float64`: Maximum depth of the domain in meters

# Returns
- `NamedTuple`: Contains optimized parameters:
  - `len_x`: Optimal X correlation length in meters
  - `len_y`: Optimal Y correlation length in meters
  - `epsilon2`: Optimal noise-to-signal ratio
  - `cvval`: Cross-validation error metric

# Notes
- Uses two-phase approach with coarse grid for efficiency
- First phase uses 5×5 parameter test points
- Second phase uses 7×7 parameter test points with refined starting point
- Memory usage is limited to 4GB with MEMTOFIT parameter
"""

"""
    optimize_diva_parameters_fast(mask, pmn, xi, x_obs, obs_values, typical_station_spacing, depth_extent; 
                                 downsample=10, nl_coarse=3, ne_coarse=3, nl_fine=5, ne_fine=5, 
                                 memtofit=4, verbose=true, timeout=60)

Faster version of parameter optimization with more aggressive downsampling, fewer test points,
timeout protection, and error handling.

# Arguments
- `mask::BitMatrix`: Boolean mask of valid grid points (true = valid, false = invalid)
- `pmn::Tuple{Matrix, Matrix}`: Scale factors for X and Y dimensions
- `xi::Tuple{Matrix, Matrix}`: Grid coordinates in X and Y dimensions
- `x_obs::Tuple{Vector, Vector}`: Observation coordinates in X and Y dimensions
- `obs_values::Vector`: Observation values (should be anomalies from mean)
- `typical_station_spacing::Float64`: Median distance between observation stations in meters
- `depth_extent::Float64`: Maximum depth of the domain in meters

# Keyword Arguments
- `downsample::Int=10`: Factor by which to downsample the grid for optimization
- `nl_coarse::Int=3`: Number of correlation length factors to test in phase 1
- `ne_coarse::Int=3`: Number of epsilon2 factors to test in phase 1
- `nl_fine::Int=5`: Number of correlation length factors to test in phase 2
- `ne_fine::Int=5`: Number of epsilon2 factors to test in phase 2
- `memtofit::Int=4`: Memory limit in GB for DIVAnd_cv
- `verbose::Bool=true`: Whether to print progress messages
- `timeout::Int=60`: Maximum seconds to allow for optimization before falling back

# Returns
- `NamedTuple`: Contains optimized parameters and status:
  - `len_x`: Optimal X correlation length in meters
  - `len_y`: Optimal Y correlation length in meters
  - `epsilon2`: Optimal noise-to-signal ratio
  - `cvval`: Cross-validation error metric or NaN if optimization failed
  - `status`: String indicating result - "success", "timeout_phase1", or "error"

# Notes
- Safe to use in production with fallback mechanisms
- Aggressively reduces grid size to speed up computation
- Will return initial parameter estimates if optimization fails
"""
function optimize_diva_parameters(obs, grid; 
                                  epsilon2_initial=0.05,
                                  downsample=10, # More aggressive downsampling
                                  ne=7, # Minimal test points in first phase
                                  memtofit=8, # Memory limit
                                  verbose=true)

    # Assign variables
    mask = grid.mask;
    pmn = grid.pmn
    x_obs = obs.x_obs
    xi = grid.xi
    obs_values = obs.anomaly
    typical_station_spacing = obs.typical_station_spacing           
    distance_extent = grid.distance_extent           
    station_depths = obs.station_depths                 
    
    if verbose
        println("Starting parameter optimization...")
        println("Grid size: $(size(mask))")
        println("Downsampling grid by factor: $downsample")
        flush(stdout)
    end

    # Get the minimum and maximum depths with appropriate offsets
    min_depth = minimum(station_depths.Depth) * 0.75
    max_depth = mean(station_depths.Depth) * 1.25

    # Get better initial length scale estimates using DIVAnd's built-in functions
    # For 2D transect data, we need a workaround since fithorz/vert assume 3D data
    # Create a pseudo-3D dataset with constant latitude
    fake_lat = zeros(length(obs.x)) # Constant latitude (we're on a transect)
    x_obs_3d = (obs.x, fake_lat, obs.y) # (lon, lat, depth) format expected by fithorzlen

    # Get depth levels that will be used
    depth_levels = sort(unique(obs.y))

    # Downsample and delimit
    depth_levels_ds = depth_levels[1:(downsample*1):end]
    test_depths = depth_levels_ds[(depth_levels_ds .>= min_depth) .& (depth_levels_ds .<= max_depth)]

    # Add error handling to the limitfun
    limitfun = function(z, len)
        if isnan(len)
            return 100_000  # Return a default value if length is NaN
        else
            return max(min(len, distance_extent), typical_station_spacing)
        end
    end

    if verbose
        println("Phase 1: Coarse parameter search...")
        flush(stdout)
    end

    # Fit the horizontal length correlation
    lenxy, dbinfo = fithorzlen(
            x_obs_3d, 
            obs.values,
            test_depths, 
            limitfun = limitfun,
            maxnsamp=50,
            limitlen=true,
            smoothz=50,
            min_rqual=0.05
        )

    # Fit the vertical length correlation
    lenz, dbinfo_vert = fitvertlen(
        x_obs_3d,
        obs.values,
        test_depths,
        smoothz=50,
        maxnsamp=50,
        min_rqual=0.05
    )

    # Get the median entries
    best_len_x = median(lenxy)
    best_len_y = median(lenz)
    len = (best_len_x, best_len_y)

    # Create a coarser grid
    mask_coarse = mask[1:downsample:end, 1:downsample:end]
    pm_coarse = pmn[1][1:downsample:end, 1:downsample:end]
    pn_coarse = pmn[2][1:downsample:end, 1:downsample:end]
    pmn_coarse = (pm_coarse, pn_coarse)
    
    # Downsample coordinate grid
    X_coarse = xi[1][1:downsample:end, 1:downsample:end]
    Y_coarse = xi[2][1:downsample:end, 1:downsample:end]
    xi_coarse = (X_coarse, Y_coarse)
    
    if verbose
        println("Downsampled grid size: $(size(mask_coarse))")
        flush(stdout)
    end
        
    if verbose
        println("Phase 2: Fit epsilon2...")
        println("Initial correlation lengths (len_xy, len_z): $len")
        println("Initial epsilon2: $epsilon2_initial")
        flush(stdout)
    end
    
    # First phase with very few test points
    try
        bestfactorl_coarse, bestfactore_coarse, cvval = DIVAnd_cv(
            mask_coarse, 
            pmn_coarse, 
            xi_coarse, 
            x_obs, 
            obs_values, 
            len, 
            epsilon2_initial, 
            0, 
            ne, 
            0;
            MEMTOFIT=memtofit
        )
        
        # Calculate refined starting points for phase 2
        epsilon2_refined = epsilon2_initial * bestfactore_coarse
        
        if verbose
            println("Optimization complete!")
            println("Refined parameters: len_x = $best_len_x, len_y = $best_len_y, epsilon2 = $epsilon2_refined")
            flush(stdout)
        end
                
        return (len_x=best_len_x, len_y=best_len_y, 
               epsilon2=epsilon2_refined, cvval=cvval,
               status="success")
               
    catch e
        # If optimization fails, return the initial values
        if verbose
            println("Optimization failed: $e")
            println("Using initial values as fallback.")
        end
        
        return (len_x=typical_station_spacing * 1.5, len_y=mean(station_depths.Depth), 
               epsilon2=epsilon2_initial, cvval=NaN,
               status="error")
    end
end

"""
    interpolate_field(grid, observations, len, epsilon2; velocity=nothing, verbose=true)

Perform DIVAnd interpolation using provided grid, observations, and parameters.

# Arguments
- `grid::NamedTuple`: Grid information from prepare_diva_grid()
- `observations::NamedTuple`: Observation data from prepare_observations()
- `len::Tuple{Float64,Float64}`: Correlation lengths (len_x, len_y) in meters
- `epsilon2::Float64`: Noise-to-signal ratio for DIVAnd

# Keyword Arguments
- `velocity::Union{Nothing,Tuple{Matrix,Matrix}}=nothing`: Optional advection 
  constraint as tuple of velocity components; if nothing, zero velocity is used
- `verbose::Bool=true`: Whether to print progress messages

# Returns
- `NamedTuple`: Contains interpolation results:
  - `field::Matrix{Float64}`: Interpolated field values (with background added)
  - `field_anomaly::Matrix{Float64}`: Interpolated anomalies
  - `s::DIVAnd.DIVAnd_constrained`: DIVAnd state information

# Notes
- Performs interpolation on anomalies, then adds background back
- Uses DIVAndrun with velocity constraint if provided
- Masked grid points will contain NaN values
"""
function interpolate_field(grid, observations, len, epsilon2; velocity=nothing, verbose=true)
    # Extract required parameters
    mask = grid.mask
    pmn = grid.pmn
    xi = grid.xi
    x_obs = observations.x_obs
    obs_anomaly = observations.anomaly
    background = observations.background
    
    # Prepare velocity if not provided
    if isnothing(velocity)
        velocity = (zeros(size(grid.X_grid)), zeros(size(grid.Y_grid)))
    end
    
    # Run interpolation on anomalies
    field_anomaly, s = DIVAndrun(mask, pmn, xi, x_obs, obs_anomaly, len, epsilon2, velocity=velocity)
    
    # Add background back
    field = field_anomaly .+ background

    # Verbosity
    if verbose
        println("Interpolation complete!")
    end
    
    return (field=field, field_anomaly=field_anomaly, s=s)
end

"""
    process_variable(bathy_file, hydro_file, variable; 
                     horizontal_res=100, vertical_res=1, 
                     max_depth=nothing, save_plot=true, 
                     output_file=nothing, km_to_m=true, 
                     epsilon2_initial=0.05, downsample=10, ne=7,
                     memtofit=8, contour_levels=[6, 7, 8, 9, 10],
                     midpoint=10, plot_size=(1200, 800), verbose=true)

Process a single oceanographic variable using DIVAnd interpolation.

# Arguments
- `bathy_file::String`: Path to bathymetry CSV file
- `hydro_file::String`: Path to hydrography CSV file
- `variable::String`: Variable name to process (e.g., "Temperature", "Salinity")

# Keyword Arguments
- `horizontal_res::Float64=100`: Horizontal resolution in meters
- `vertical_res::Float64=1`: Vertical resolution in meters
- `max_depth::Union{Nothing,Float64}=nothing`: Maximum depth for visualization
- `save_plot::Bool=true`: Whether to save plot to file
- `output_file::Union{Nothing,String}=nothing`: Base output filename (will append variable name)
- `km_to_m::Bool=true`: Convert distance units from km to m
- `epsilon2_initial::Float64=0.05`: Initial noise-to-signal ratio
- `downsample::Int=10`: Factor by which to downsample the grid for optimization
- `ne::Int=7`: Number of epsilon2 factors to test
- `memtofit::Int=8`: Memory limit in GB for DIVAnd_cv
- `contour_levels::Vector{Float64}=[6,7,8,9,10]`: Temperature contour levels to display
- `midpoint::Float64=10`: Reference point for color scale
- `plot_size::Tuple{Int,Int}=(1200,800)`: Size of output plot
- `verbose::Bool=true`: Whether to print progress messages

# Returns
- `NamedTuple`: Contains processed results and parameters
"""
function process_diva_variable(bathy_file, hydro_file, variable; 
                             horizontal_res=100, vertical_res=1, 
                             output_file=nothing, km_to_m=true,
                             epsilon2_initial=0.05, downsample=10, ne=7,
                             memtofit=8, verbose=true)
    
    # Set default output directory and base filename if not provided
    if isnothing(output_file)
        output_dir = joinpath(dirname(bathy_file), "..", "output")
        mkpath(output_dir)
        output_base = joinpath(output_dir, "interpolated")
    else
        output_dir = dirname(output_file)
        mkpath(output_dir)
        output_base = output_file
    end
    
    # Process the variable
    if verbose
        println("Processing variable: $variable")
    end
    
    # Prepare observations for this variable
    obs = prepare_observations(hydro_file, variable, km_to_m=km_to_m, verbose=verbose)
    
    # Prepare grid
    grid = prepare_diva_grid(bathy_file, minimum(obs.x), maximum(obs.x), horizontal_res, vertical_res, km_to_m=km_to_m, verbose=verbose)
    
    # Optimize parameters
    optimized_params = optimize_diva_parameters(
        obs, grid, 
        epsilon2_initial=epsilon2_initial,
        downsample=downsample,
        ne=ne,
        memtofit=memtofit,
        verbose=verbose
    )
    
    # Perform interpolation
    result = interpolate_field(
        grid, 
        obs, 
        (optimized_params.len_x, optimized_params.len_y),
        optimized_params.epsilon2,
        verbose=verbose
    )
    
    # Export interpolated field to CSV
    export_interpolated_field(grid, result, variable, "$(output_base)_$(lowercasefirst(variable)).csv", verbose=verbose)
        
    return (grid=grid, result=result, params=optimized_params)
end

function batch_process_diva(bathy_file, hydro_file, variables::Vector{String}; 
                              horizontal_res=100, vertical_res=1, 
                              output_file=nothing, km_to_m=true,
                              epsilon2_initial=0.05, downsample=10, ne=7,
                              memtofit=8, verbose=true)

    # Process each variable and collect results
    results = Dict{String, NamedTuple}()
    start_time = time()

    for (i, variable) in enumerate(variables)
        if verbose
            println("\n" * "="^50)
            println("[$i/$(length(variables))] Processing variable: $variable")
            println("="^50)
        end
        
        # Process this variable
        results[variable] = process_diva_variable(
            bathy_file, hydro_file, variable;
            horizontal_res=horizontal_res,
            vertical_res=vertical_res,
            output_file=output_file, 
            km_to_m=km_to_m,
            epsilon2_initial=epsilon2_initial,
            downsample=downsample,
            ne=ne,
            memtofit=memtofit,
            verbose=verbose
        )
    end

    # Calculate total processing time
    elapsed = time() - start_time
    minutes = floor(Int, elapsed / 60)
    seconds = round(elapsed % 60, digits=1)
    
    if verbose
        println("\n" * "="^50)
        println("Processing complete for $(length(variables)) variables")
        println("Total time: $minutes min $seconds sec")
        println("="^50)
    end
    
    return results
end

"""
    export_interpolated_field(grid, result, variable, output_file; include_masked=false, verbose=true)

Export interpolated field to CSV file with proper coordinates.

# Arguments
- `grid`: Grid information from prepare_diva_grid()
- `result`: Result from interpolate_field()
- `variable::String`: Variable name
- `output_file::String`: Path to save the CSV file

# Keyword Arguments
- `include_masked::Bool=false`: Whether to include masked (NaN) values in output
- `verbose::Bool=true`: Whether to print progress messages
"""
function export_interpolated_field(grid, result, variable, output_file; 
                                   include_masked=false, verbose=true)
    # Create dataframe to store results
    df = DataFrame(
        Distance = Float64[],
        Depth = Float64[],
        Variable = String[],
        Value = Float64[]
    )
    
    # Fill with all grid points and interpolated values
    for i in 1:length(grid.x_grid)
        for j in 1:length(grid.y_grid)
            # Only include valid points (not masked) unless include_masked=true
            if include_masked || grid.mask[i, j]
                push!(df, (
                    grid.x_grid[i],         # Distance in meters
                    grid.y_grid[j],         # Depth in meters
                    variable,               # Variable name
                    result.field[i, j]      # Interpolated value
                ))
            end
        end
    end
    
    # Write to CSV
    CSV.write(output_file, df)
    
    if verbose
        println("Exported interpolated $variable field to $output_file")
    end
    
    return df
end

"""
    visualize_field(grid, observations, result; max_depth=nothing, midpoint=10)

Create a visualization of the interpolated field.
"""
function visualize_field(grid, observations, bathy_df, result; max_depth=nothing, midpoint=10)
    if isnothing(max_depth)
        max_depth = maximum(observations.y)
    end
    
    # Find which y indices correspond to 0-{max_depth}
    depth_mask = grid.y_grid .<= max_depth
    y_subset = grid.y_grid[depth_mask]
    temp_subset = result.field'[depth_mask, :]
    
    # Update color palette for divergence
    temp_min, temp_max = extrema(filter(!isnan, result.field))
    temp_range = max(abs(temp_min - midpoint), abs(temp_max - midpoint))
    clim_min = midpoint - temp_range
    clim_max = midpoint + temp_range
    
    # Plot interpolated field
    p = heatmap(
        collect(grid.x_grid), 
        collect(y_subset),
        color=:vik,
        temp_subset,        
        xlabel="Distance (km)", ylabel="Depth (m)", 
        title="Interpolated Field", 
        colorbar_title="Temperature (°C)",
        xlims=(minimum(grid.x_grid), maximum(grid.x_grid)),
        clims=(clim_min, clim_max),
        size=(1200, 800),
        margin=5Measures.mm
    )

    # Add bathymetry
    plot!(p, bathy_df.Distance, bathy_df.BottomDepth, color=:black, linewidth=2, 
      label="Bottom Depth", yflip=true,
      ylim=(0, max_depth),
      legend=:bottomright)

    # Add contour lines at specific isotherms inside the function
    contour!(p, collect(grid.x_grid), collect(grid.y_grid), 
             result.field', 
             levels=[6, 7, 8, 9, 10], 
             color=:gray50, 
             colorbar_entry=false,
             linewidth=2)

    yflip!(p)
    
    return p
end

end