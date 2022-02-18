################################################################################
# plot_rate_balances.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Plots the chemical reaction rates, transport rates, and sum of both
# to check that the atmosphere is equilibrated.
# 
# Eryn Cangi
# Created March 2021
# Last edited: 23 March 2021
# Currently tested for Julia: 1.5.3
################################################################################

using Photochemistry: input, search_subfolders, T_updated, create_folder, get_ncurrent, Psat, Psat_HDO, plot_rxns, effusion_velocity, charge_type, scaleH, meanmass
using PyPlot
using PyCall

include("CONSTANTS.jl")
include("CUSTOMIZATIONS.jl")

user_input_paramfile = input("Enter a parameter file or press enter to use default (PARAMETERS.jl): ")
paramfile = user_input_paramfile == "" ? "PARAMETERS.jl" : user_input_paramfile*".jl"
include(paramfile)

# Load the new standard reaction network from file if it the parameter file doesn''t have its own network.
if !@isdefined reactionnet 
    include("$(@__DIR__)/reaction_network.jl")
end

println("Found the folder name: $(sim_folder_name)")
simfolder = results_dir*sim_folder_name*"/"

converged_file = final_atm_file == "" ? "final_atmosphere.h5" : final_atm_file
if !isfile(converged_file)
    converged_file = input("No file found, please enter file to use including .h5: ")
end

println("Using final atmosphere file: $(converged_file)")

# Do the actual plotting =====================================================================================
create_folder("chemeq_plots", simfolder)

filelist = search_subfolders(simfolder, "$(converged_file)", type="files")
println(filelist)
for f in filelist
    ncur = get_ncurrent(simfolder*f)

    # Now plot the reactions ------------------------------------------------------------
    # looks for a pattern in the filename that looks like X{.XeX}, where curly braces are optional. 
    # i.e. it will match "100", "100.0", "0.001", "1e2" or "1.0e2".
    dtmatch = match(r"\d+\.\d*e*\d*", f) 
    if typeof(dtmatch) != Nothing
        dtval = dtmatch.match
    else
        dtval = input("Please enter dt: ")
    end

    println("Working on dt=$(dtval)")

    # Plotting for the old version of the code
    if @isdefined fullspecieslist
        println("ALERT: You have to go back to the parameter file, $(paramfile), and change all the regular T's in the reaction network to Tn.")
        println("Surely there is a better way to do this but for right now this is what I got.")
        pritnln("If you've left them in as Tn, go back and change them to T later. For posterity and such.")
        for sp in fullspecieslist
            plot_rxns(sp, ncur, Tn_arr, Ti_arr, Te_arr, speciesbclist, reactionnet, fullspecieslist, [:CO2pl], transportspecies, speciesmolmasslist, alt, n_alt_index, dz, nothing, 
                      num_layers, plot_grid, results_dir, Tprof_for_Hs, Tprof_for_diffusion, subfolder=sim_folder_name, plotsfolder="chemeq_plots", num=dtval, extra_title="dt=$(dtval)")
        end
    else  # The normal case, the new code, where we use all_species and not fullspecieslist. 
        for sp in all_species
            println("Working on $(sp)")

                (sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, Tp, # params to pass
                   bcdict::Dict{Symbol, Matrix{Any}}, rxnnet, 
                   allsp, ionsp, transsp, chemsp, mmass::Dict{Symbol, Int64}, alts, n_alt_index::Dict, dz, 
                   polar::Dict{Symbol, Float64}, numlyrs::Int64, plotalts, 
                   results_dir::String, T_for_Hs::Dict{String, Vector{Any}}, T_for_diff::Dict{String, Vector{Any}}; 
                   shown_rxns=nothing, subfolder="", plotsfolder="", dt=nothing, num="", extra_title="", 
                   plot_timescales=false, plot_total_rate_coefs=false, showonly=false)
    #=

            plot_rxns(sp, ncur, Tn_arr, Ti_arr, Te_arr, Tplasma_arr, speciesbclist, reactionnet, 
                      all_species, ion_species, transport_species, chem_species, 
                      molmass, alt, n_alt_index, dz, polarizability, 
                      num_layers, plot_grid, results_dir, Tprof_for_Hs, Tprof_for_diffusion, 
                      subfolder=sim_folder_name, plotsfolder="chemeq_plots", num=dtval, extra_title="dt=$(dtval)")
        end
    end

end
