__precompile__()

module SeasonsAnalysis

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using DelimitedFiles
using LinearAlgebra
using Unicode
using StatsBase
using PlotUtils
using GeneralizedGenerated
using DataFrames
using XLSX
using Dates

photochemistry_source_dir = "$(@__DIR__)/../Photochemistry/src/"
println("loading Photochemistry.jl into SeasonsAnalysis from $photochemistry_source_dir")
push!(LOAD_PATH, photochemistry_source_dir)
using Photochemistry  # custom module

export 
make_equilibrium_plots_temp, make_seasonal_cycle_plots_temp, make_seasonal_cycle_plots_inclusive,
make_equilibrium_plots_water, make_seasonal_cycle_plots_water, make_insolation_plots, plot_limiting_flux, 
make_3panel_figure, DH_6panel, f_vs_time, flux_vs_time, seasonal_cycling_figure_original, 
seasonal_cycling_figure_skeleton, make_temperature_panel, DH_alt_profile_singlepanels, make_density_panel, 
make_water_panel, make_insolation_panel, esc_flux_vs_time, DH_vs_time, density_vs_time, generate_indvar_vs_time_array,
parent_folders_from_full_path, draw_DH_lines, collect_atmospheres


# ~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~        

code_loc = "$(@__DIR__)/../../"

# Basic files
include(code_loc*"CUSTOMIZATIONS.jl") 
include(code_loc*"CONSTANTS.jl")

# Submodules
include("MakePlots.jl")
include("Support.jl")

end
