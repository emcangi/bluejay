################################################################################
# plot_rate_balances.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Plots the chemical reaction rates, transport rates, and sum of both
# to check that the atmosphere is equilibrated. Stand-alone. These plots are also
# now part of the main simulation routine as of 7 June 2022, so this file is 
# likely to be used more rarely.
# 
# Eryn Cangi
# Created March 2021
# Last edited: 7 June 2022
# Currently tested for Julia: 1.7.1
################################################################################

# Modules and critical files ================================================================================
using Revise
photochemistry_source_dir = "$(@__DIR__)/Photochemistry/src/"
# println("loading Photochemistry.jl from $photochemistry_source_dir")
push!(LOAD_PATH, photochemistry_source_dir)
using Photochemistry: create_folder, generate_code, get_ncurrent, get_paramfile, input, search_subfolders, T, Psat, Psat_HDO, plot_rxns, 
                      effusion_velocity, charge_type, scaleH, meanmass, load_reaction_network, format_Jrates, load_from_paramlog
using PyPlot
using PyCall
using GeneralizedGenerated

include("CONSTANTS.jl")
include("CUSTOMIZATIONS.jl")

# **************************************************************************** #
#                                                                              #
#                  SELECT FOLDER AND LOAD THE ATMOSPHERIC STATE                #
#                                                                              #
# **************************************************************************** #

simfolder = input("Enter folder name: ")

file_to_use = input("Enter filename to use or press enter to look for: ")
full_file = results_dir*simfolder*"/"*file_to_use
println("Searching for $(full_file)")
while !isfile(full_file)
    if !occursin(".h5", file_to_use)
        global file_to_use = file_to_use * ".h5"
    end
    global file_to_use = input("$(full_file) not found, please enter file name again: ")
end

# **************************************************************************** #
#                                                                              #
#           LOAD PARAMETERS, REACTION NETWORK, ATMOSPHERIC STATE               #
#                                                                              #
# **************************************************************************** #

vardict = load_from_paramlog(results_dir*simfolder*"/")

ions_included = vardict["ions_included"]
hrshortcode = vardict["hrshortcode"]
rshortcode = vardict["rshortcode"]
neutral_species = vardict["neutral_species"]
ion_species = vardict["ion_species"]
all_species = vardict["all_species"]
transport_species = vardict["transport_species"]
chem_species = vardict["chem_species"]
Tn_arr = vardict["Tn_arr"]
Ti_arr = vardict["Ti_arr"]
Te_arr = vardict["Te_arr"]
Tplasma_arr = vardict["Tplasma_arr"]
Tprof_for_Hs = vardict["Tprof_for_Hs"]
Tprof_for_diffusion = vardict["Tprof_for_diffusion"]
Hs_dict = vardict["Hs_dict"]
speciesbclist = vardict["speciesbclist"]
reaction_network_spreadsheet = vardict["rxn_spreadsheet"]
water_bdy = vardict["water_bdy"]

# Load the new standard reaction network from spreadsheet
if ions_included==true
    ions_included=true
    reaction_network, hot_H_network, hot_D_network, hot_H2_network, hot_HD_network = load_reaction_network(reaction_network_spreadsheet; ions_on=ions_included, 
                                                                           get_hot_rxns=true, all_species)
    const hot_H_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_H_network]);
    const hot_D_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_D_network]);
    const hot_H2_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_H2_network]);
    const hot_HD_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_HD_network]);
    const functions_for_rxns = [hot_H_rc_funcs, hot_D_rc_funcs]
else 
    reaction_network, hot_H_network, hot_D_network = load_reaction_network(reaction_network_spreadsheet; ions_on=ions_included, all_species)
    const functions_for_rxns = []
end 

println("Using final atmosphere file: $(file_to_use)")

# extra_file = input("Enter any additional files you would like to make plots for or press enter to skip: ")

# **************************************************************************** #
#                                                                              #
#                       DO THE ACTUAL PLOTTING                                 #
#                                                                              #
# **************************************************************************** #
create_folder("chemeq_plots", results_dir*simfolder*"/")

filelist = search_subfolders(results_dir*simfolder, "$(file_to_use)", type="files")
# if extra_file != ""
    # push!(filelist, extra_file)
# end

for f in filelist
    ncur = get_ncurrent(results_dir*simfolder*"/"*f)

    Jratedict = Dict([j=>ncur[j] for j in keys(ncur) if occursin("J", string(j))])

    # Plot  ---------------------------------------------------------------------------------
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

            plot_rxns(sp, ncur, results_dir; nonthermal=ions_included, subfolder=simfolder, plotsfolder="chemeq_plots", 
                      num="$(f[1:end-3])", all_species, alt, chem_species, collision_xsect, 
                      dz, hot_D_rc_funcs, hot_H_rc_funcs, hot_H2_rc_funcs, hot_HD_rc_funcs, Hs_dict, hot_H_network, hot_D_network, hot_H2_network, hot_HD_network, hrshortcode, ion_species, Jratedict,
                      molmass, neutral_species, non_bdy_layers, num_layers, n_all_layers, n_alt_index, polarizability, 
                      plot_grid, q, rshortcode, reaction_network, speciesbclist, Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, 
                      Tprof_for_Hs, Tprof_for_diffusion, transport_species, upper_lower_bdy_i=n_alt_index[water_bdy], upper_lower_bdy=water_bdy, zmax)#, shown_rxns=theserxns) # Uncomment to show individual reactions.
        end
    end
end
