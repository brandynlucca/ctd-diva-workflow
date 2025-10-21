include("diva_processing.jl");
using .DIVAInterp;
using DIVAnd, Statistics
# using CSV, DataFrames, DIVAnd, Statistics;
# using Measures

# 0. Define target files.

println("Testing DIVAInterp Module");
println("========================\n");

# 1. Define paths with @__DIR__ to ensure correct path resolution
bathy_file = joinpath(@__DIR__, "..", "data", "test_transect_bathy.csv");
hydro_file = joinpath(@__DIR__, "..", "data", "test_hydrography.csv");
output_file=joinpath(@__DIR__, "..", "data", "test_output_new");
variable = "Temperature";

using Dates


# Header information
println("=" * repeat("=", 50))
println("| DIVAnd Oceanographic Interpolation")
println("| $(now())")
println("=" * repeat("=", 50))

# File metadata
println("\nProcessing Files:")
println("- Bathymetry: $(basename(bathy_file))")
println("- Hydrography: $(basename(hydro_file))")
println("- Variable: $variable")
if !isnothing(output_file)
    println("- Output: $(basename(output_file))")
end

proc = process_diva_variable(bathy_file, hydro_file, variable; output_file=output_file);

res = batch_process_diva(bathy_file, hydro_file, ["Temperature", "Salinity", "Sigma_t"])

# visualize_field(res["Salinity"].grid, res["Salinity"].obs 

# function batch_process_diva(bathy_file, hydro_file, variables::Vector{String}; 
#                               horizontal_res=100, vertical_res=1, 
#                               output_file=nothing, km_to_m=true,
#                               epsilon2_initial=0.05, downsample=10, ne=7,
#                               memtofit=8, verbose=true)

#     # Process each variable and collect results
#     results = Dict{String, NamedTuple}()
#     start_time = time()

#     for (i, variable) in enumerate(variables)
#         if verbose
#             println("\n" * "="^50)
#             println("[$i/$(length(variables))] Processing variable: $variable")
#             println("="^50)
#         end
        
#         # Process this variable
#         results[variable] = process_diva_variable(
#             bathy_file, hydro_file, variable;
#             horizontal_res=horizontal_res,
#             vertical_res=vertical_res,
#             output_file=output_file, 
#             km_to_m=km_to_m,
#             epsilon2_initial=epsilon2_initial,
#             downsample=downsample,
#             ne=ne,
#             memtofit=memtofit,
#             verbose=verbose
#         )
#     end

#     # Calculate total processing time
#     elapsed = time() - start_time
#     minutes = floor(Int, elapsed / 60)
#     seconds = round(elapsed % 60, digits=1)
    
#     if verbose
#         println("\n" * "="^50)
#         println("Processing complete for $(length(variables)) variables")
#         println("Total time: $minutes min $seconds sec")
#         println("="^50)
#     end
    
#     return results
# end

# Visualize result
# p = visualize_field(grid, obs, bathy_df, result, max_depth=150);
# display(p)
