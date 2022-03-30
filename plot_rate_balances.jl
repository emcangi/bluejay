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

# Modules and critical files ================================================================================
using Revise
photochemistry_source_dir = "$(@__DIR__)/Photochemistry/src/"
println("loading Photochemistry.jl from $photochemistry_source_dir")
push!(LOAD_PATH, photochemistry_source_dir)
using Photochemistry: create_folder, get_ncurrent, get_paramfile, input, search_subfolders, T_updated, Psat, Psat_HDO, plot_rxns, 
                      effusion_velocity, charge_type, scaleH, meanmass, load_reaction_network
using PyPlot
using PyCall

# **************************************************************************** #
#                                                                              #
#           LOAD PARAMETERS, REACTION NETWORK, ATMOSPHERIC STATE               #
#                                                                              #
# **************************************************************************** #

include("CONSTANTS.jl")
include("CUSTOMIZATIONS.jl")
paramfile = get_paramfile(code_dir)
include(paramfile)

# Load the new standard reaction network from spreadsheet
reaction_network = load_reaction_network(reaction_network_spreadsheet; Jratelist, absorber, photolysis_products, all_species, ions_on=ions_included)

println("Found the folder name: $(sim_folder_name)")
simfolder = results_dir*sim_folder_name*"/"

converged_file = "final_atmosphere.h5"
if !isfile(simfolder*converged_file)
    converged_file = input("No file found, please enter file to use: ")
    if !occursin(".h5", converged_file)
        converged_file = converged_file * ".h5"
    end
end

println("Using final atmosphere file: $(converged_file)")

# **************************************************************************** #
#                                                                              #
#                       DO THE ACTUAL PLOTTING                                 #
#                                                                              #
# **************************************************************************** #
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
            plot_rxns(sp, ncur, Tn_arr, Ti_arr, Te_arr, speciesbclist, reaction_network, fullspecieslist, [:CO2pl], transportspecies, speciesmolmasslist, alt, n_alt_index, dz, nothing, 
                      num_layers, plot_grid, results_dir, Tprof_for_Hs, Tprof_for_diffusion, subfolder=sim_folder_name, plotsfolder="chemeq_plots", num=dtval, extra_title="dt=$(dtval)")
        end
    else  # The normal case, the new code, where we use all_species and not fullspecieslist. 
        for sp in all_species
            println("Working on $(sp)")

            # Uncomment this section in case you want to show individual reactions on the chemistry panel of the plot.
            # theserxns = []
            # for rxn in reaction_network
            #     if occursin(string(sp), string(rxn[1])) | occursin(string(sp), string(rxn[2]))
            #         push!(theserxns, rxn)
            #     end
            # end

            plot_rxns(sp, ncur, results_dir; subfolder=sim_folder_name, plotsfolder="chemeq_plots", num=dtval, extra_title="dt=$(dtval)", 
                      Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, bcdict=speciesbclist, rxnnet=reaction_network, 
                      all_species, neutral_species, ion_species, transport_species, chem_species, 
                      molmass, polarizability, Tprof_for_Hs, Tprof_for_diffusion, Hs_dict,
                      alt, n_alt_index, dz, num_layers, n_all_layers, plot_grid, upper_lower_bdy_i, upper_lower_bdy, q)#, shown_rxns=theserxns) # Uncomment to show individual reactions.
        end
    end
end
