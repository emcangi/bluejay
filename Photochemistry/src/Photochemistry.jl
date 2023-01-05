__precompile__()

module Photochemistry

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils
using GeneralizedGenerated
using DataFrames
using XLSX
using Random
using Printf

export 
# Atmosphere.jl
atm_dict_to_matrix, 
atm_matrix_to_dict, 
column_density, 
column_density_above, 
compile_ncur_all,
electron_density, 
find_exobase, 
flatten_atm, 
get_deuterated, 
meanmass, 
ncur_with_boundary_layers, 
n_tot, 
scaleH,
unflatten_atm, 

# BasicUtilities.jl
## Standard, miscellaneous functions
    deletefirst, 
    find_nonfinites, 
    fluxsymbol, 
    generate_code, 
    getpos, 
    input, 
    get_paramfile, 
    logrange, 
    nans_present, 
    next_in_loop, 
    searchsortednearest, 
    subtract_difflength,

## String Manipulation
    charge_type, 
    decompose_chemistry_string, 
    format_chemistry_string, 
    format_scin, 
    format_sec_or_min, 
    string_to_latexstr,

# Chemistry functions
calculate_stiffness, 
check_jacobian_eigenvalues, 
chemical_jacobian, 
eval_rate_coef, 
get_column_rates, 
get_volume_rates,
getrate,
loss_equations, 
loss_rate, 
make_chemjac_key, 
make_net_change_expr, 
production_equations, 
production_rate, 
reactant_density_product, 
rxns_where_species_is_observer, 
volume_rate_wrapper,

# Crosssections.jl
padtosolar, 
populate_xsect_dict, 

# FileIO.jl
create_folder,
get_ncurrent, 
searchdir, 
search_subfolders, 
write_atmosphere,
write_final_state,

## Logging functions
    get_param, 
    load_bcdict_from_paramdf, 
    load_from_paramlog, 
    write_to_log,

# PhotochemicalEquilibrium.jl
choose_solutions, 
construct_quadratic, 
group_terms, 
loss_coef, 
linear_in_species_density,

# Plotting functions
get_colors,
get_grad_colors, 
plot_atm, 
plot_bg, 
plot_extinction, 
plot_Jrates, 
plot_production_and_loss,
plot_rxns,
plot_temp_prof, 
plot_water_profile, 
top_mechanisms,    

# ReactionNetwork.jl
filter_network, 
format_Jrates, 
load_reaction_network, 
log_reactions, 

# Temperature.jl
T,

# TransportAndBoundaryConditions.jl
## Boundary conditions and flux
    boundaryconditions, 
    effusion_velocity, 
    escape_probability, 
    escaping_hot_atom_production, 
    nonthermal_escape_flux,
    flux_pos_and_neg,
    get_flux, 
    limiting_flux, 
    limiting_flux_molef, 

## Vertical transport
    binary_dcoeff_inCO2, 
    Dcoef_neutrals, 
    Dcoef!, 
    diffparams, 
    fluxcoefs, 
    flux_param_arrays,  
    get_transport_PandL_rate, 
    Keddy, 
    update_diffusion_and_scaleH, 
    update_transport_coefficients,  


# UnitConversions.jl
GEL_to_molecule, molec_to_GEL, 

# Water profile functions  
precip_microns, 
Psat, 
Psat_HDO, 
setup_water_profile!, 
water_tanh_prof

# ~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~        

code_loc = "$(@__DIR__)/../../"
# println("Photochemistry.jl code_loc = $code_loc")
include(code_loc*"CUSTOMIZATIONS.jl") 
include(code_loc*"CONSTANTS.jl")
                                                                                
include("Atmosphere.jl")
include("BasicUtilities.jl")
include("Chemistry.jl")
include("Crosssections.jl")
include("FileIO.jl")
include("PhotochemicalEquilibrium.jl")
include("Plotting.jl")
include("ReactionNetwork.jl")
include("Temperature.jl")
include("TransportAndBoundaryConditions.jl")
include("UnitConversions.jl")
include("Water.jl")

end
