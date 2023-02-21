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

# Core.jl
    ## Basic mechanical functions
    atm_dict_to_matrix, 
    atm_matrix_to_dict, 
    compile_ncur_all,
    flatten_atm, 
    solve_sparse,
    subtract_difflength,
    unflatten_atm, 

    ## Chemistry functions
    calculate_stiffness, 
    check_jacobian_eigenvalues, 
    chemical_jacobian, 
    eval_rate_coef, 
    getrate,
    loss_equations, 
    loss_rate, 
    production_equations, 
    production_rate, 
    update_Jrates!,

    ## Escape
    effusion_velocity, 
    escape_probability, 
    escaping_hot_atom_production, 
    nonthermal_escape_flux,

    ## Photochemical equilibrium
    choose_solutions, 
    construct_quadratic, 
    group_terms, 
    loss_coef, 
    linear_in_species_density,
    setup_photochemical_equilibrium,

    ## Transport and boundary conditions
    binary_dcoeff_inCO2, 
    boundaryconditions, 
    Dcoef_neutrals, 
    Dcoef!, 
    diffparams, 
    fluxcoefs, 
    Keddy, 
    update_diffusion_and_scaleH, 
    update_transport_coefficients,  

# Atmosphere.jl
column_density, 
column_density_above, 
column_density_species,
electron_density, 
find_exobase, 
meanmass, 
ncur_with_boundary_layers, 
n_tot, 
optical_depth,
scaleH,

# BasicUtilities.jl
## Standard, miscellaneous functions
deletefirst, 
df_lookup,
find_nonfinites, 
fluxsymbol, 
generate_code, 
get_deuterated, 
get_paramfile, 
getpos, 
input, 
logrange, 
nans_present, 
next_in_loop, 
searchsortednearest, 

## String Manipulation
charge_type, 
decompose_chemistry_string, 
format_chemistry_string, 
format_scin, 
format_sec_or_min, 
string_to_latexstr,

# AnalyzeChemAndTransport.jl
## Chemistry
chemical_lifetime,
get_column_rates, 
get_volume_rates,
make_chemjac_key, 
make_net_change_expr, 
reactant_density_product, 
volume_rate_wrapper,

## Transport and escape
diffusion_timescale,
final_escape,
fractionation_factor,
get_transport_PandL_rate, 
flux_pos_and_neg,
#get_flux, 
limiting_flux, 
limiting_flux_molef, 

# Crosssections.jl
padtosolar, 
populate_xsect_dict, 

# FileIO.jl
create_folder,
get_elapsed_time, 
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
load_network_and_make_functions,
load_reaction_network, 
log_reactions, 
modify_rxn_spreadsheet,
rxns_where_species_is_observer, 

# Temperature.jl
T,

# UnitConversions.jl
GEL_to_molecule, 
molec_to_GEL, 
DH_conversion, 

# Water profile functions  
precip_microns, 
Psat, 
Psat_HDO, 
set_h2oinitfrac_bySVP,
setup_water_profile!, 
water_tanh_prof

# ~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~        

code_loc = "$(@__DIR__)/../../"
# println("Photochemistry.jl code_loc = $code_loc")
include(code_loc*"CUSTOMIZATIONS.jl") 
include(code_loc*"CONSTANTS.jl")

include("AnalyzeChemAndTransport.jl")                                                                                
include("Atmosphere.jl")
include("BasicUtilities.jl")
include("Core.jl")
include("Crosssections.jl")
include("FileIO.jl")
include("Plotting.jl")
include("ReactionNetwork.jl")
include("Temperature.jl")
include("UnitConversions.jl")
include("Water.jl")

end
