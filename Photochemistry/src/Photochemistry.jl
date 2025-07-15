__precompile__()

module Photochemistry

using DataFrames
using DelimitedFiles
using GeneralizedGenerated
using HDF5, JLD
using IncompleteLU 
using Interpolations
using LaTeXStrings
using LinearAlgebra
using PlotUtils
using Printf
using PyPlot # Do not change the order of PyPlot and PyCall, it causes problems.
using PyCall
using Random
using SparseArrays
using XLSX
using Latexify

export 

# Note: There are some functions in the module which don't appear in this list.
# That's because they haven't ever needed to be called outside the file in which they
# are created.
# Core.jl
    ## Atmosphere attributes
    column_density, 
    column_density_above, 
    column_density_species,
    electron_density, 
    find_exobase, 
    meanmass, 
    n_tot, 
    optical_depth,
    reduced_mass,
    scaleH,
    scaleH_lowerboundary,

    ## Atmospheric object manipulation
    atm_dict_to_matrix, 
    atm_matrix_to_dict, 
    compile_ncur_all,
    flatten_atm, 
    ncur_with_boundary_layers, 
    unflatten_atm, 
    TooManyIterationsException,

    ## Chemistry functions
    calculate_stiffness, 
    charge_type, 
    check_jacobian_eigenvalues, 
    chemical_jacobian, 
    eval_rate_coef, 
    getrate,
    loss_equations, 
    loss_rate, 
    make_net_change_expr, 
    production_equations, 
    production_rate, 
    solve_sparse,
    subtract_difflength,
    update_Jrates!,

    ## Escape
    effusion_velocity, 
    escape_probability, 
    escaping_hot_atom_production, 
    nonthermal_escape_flux,

    ## Temperature
    T_Mars,
    T_Venus,

    ## Transport and boundary conditions
    binary_dcoeff_inCO2, 
    boundaryconditions, 
    boundaryconditions_horiz,
    Dcoef_neutrals, 
    Dcoef!, 
    diffparams, 
    fluxcoefs, 
    fluxcoefs_horiz,
    Keddy, 
    update_diffusion_and_scaleH, 
    update_transport_coefficients,
    update_horiz_transport_coefficients,

    ## Water
    precip_microns, 
    colabund_from_prum,
    Psat, 
    Psat_HDO, 
    set_h2oinitfrac_bySVP,
    setup_water_profile!, 
    water_tanh_prof,

    ## Photochemical equilibrium
    choose_solutions, 
    construct_quadratic, 
    group_terms, 
    loss_coef, 
    linear_in_species_density,
    setup_photochemical_equilibrium,

    ## Helper
    generate_code, 
    get_deuterated,
    get_paramfile,

# BasicUtilities.jl
    ## Standard, miscellaneous functions
    check_requirements,
    deletefirst, 
    df_lookup,
    find_nonfinites, 
    fluxsymbol, 
    generate_code, 
    get_counts,
    getpos, 
    input, 
    logrange,
    nans_present,
    arrays_equal_with_nan,
    next_in_loop,
    searchsortednearest,

    ## String Manipulation
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
    reactant_density_product, 
    volume_rate_wrapper,

    ## Transport and escape
    diffusion_timescale,
    final_escape,
    fractionation_factor,
    get_transport_PandL_rate, 
    get_directional_fluxes,
    flux_pos_and_neg,
    limiting_flux, 
    limiting_flux_molef, 
    limiting_flow_velocity,

# Crosssections.jl
populate_xsect_dict, 
padtosolar, 

# FileIO.jl
add_column,
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
plot_directional_flux,
plot_extinction, 
plot_Jrates, 
plot_net_volume_change,
plot_production_and_loss,
plot_reaction_on_demand,
plot_rxns,
plot_species_on_demand,
plot_temp_prof, 
plot_tophot_lineandbar,
plot_water_profile, 
set_rc_params,
top_mechanisms,    
turn_off_borders,

# ReactionNetwork.jl
    ## Load and manipulate
    calculate_and_write_column_rates,
    enthalpy_of_reaction,
    filter_network, 
    find_duplicates,
    load_network_and_make_functions,
    load_reaction_network, 
    log_reactions, 
    make_k_expr,
    rxns_where_species_is_observer, 

    ## Formatting the reaction network object
    format_Jrates, 
    get_Jrate_symb,

    ## Enthalpy calculations
    modify_rxn_spreadsheet,

# UnitConversions.jl
DH_conversion, 
GEL_to_molecule, 
molec_to_GEL, 
total_escape_to_GEL,
total_escape_to_area_escape,
prum_to_ppm,
rayleigh_fractionation


# ~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~        

code_loc = "$(@__DIR__)/../../"
# println("Photochemistry.jl code_loc = $code_loc")
include(code_loc*"CONSTANTS.jl")
include(code_loc*"PLOT_STYLES.jl")

include("AnalyzeChemAndTransport.jl")                                                                                
include("BasicUtilities.jl")
include("Core.jl")
include("Crosssections.jl")
include("FileIO.jl")
include("Plotting.jl")
include("ReactionNetwork.jl")
include("UnitConversions.jl")

end
