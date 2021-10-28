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

using Photochemistry
using PyPlot
using PyCall

user_input_paramfile = input("Enter a parameter file or press enter to use default (PARAMETERS.jl): ")
paramfile = user_input_paramfile == "" ? "PARAMETERS.jl" : user_input_paramfile*".jl"
include(paramfile)

parentfolder = input("Which folder? (no trailing slash please): ")

fnstr = input("Enter the exact name of the file you want to plot for (without .h5): ")


# Set up the temperature stuff that will be used for boundary conditions and plotting.
# Get the control temperatures from the simulation logfile
simfile = search_subfolders(foldr, "simulation_params"; type="files")[1]
logfile = open(foldr*simfile) do file
             read(file, String)
          end
controltemps = [parse(Float64, x) for x in [match(r"(?<=T_s\=)\d+\.\d", logfile).match, match(r"(?<=T_tropo\=)\d+\.\d", logfile).match, match(r"(?<=T_exo\=)\d+\.\d", logfile).match]]

Temp_n(z::Float64) = T_all(z, controltemps[1], controltemps[2], controltemps[3], "neutral")
Temp_keepSVP(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral")


# Do the actual plotting =====================================================================================
create_folder("chemeq_plots", results_dir*parentfolder*"/")

filelist = search_subfolders(results_dir*parentfolder, "$(fnstr).h5", type="files")#r"ncurrent_\d+\.\d+e*14\.h5", type="files")
println(filelist)
for f in filelist

    # Reload the boundary conditions for this file -------------------------------------
    ncur = get_ncurrent(results_dir*parentfolder*"/"*f)

    fix_SVP = true
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt))
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))

    global const speciesbclist=Dict(
                    :CO2=>["n" 2.1e17; "f" 0.],
                    :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                    :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                    :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                    :HDO=>["n" HDOsat[1]; "f" 0.],
                    :O=>["f" 0.; "f" 1.2e8],
                    :H2=>["f" 0.; "v" effusion_velocity(Temp_n(zmax), 2.0, zmax)],  # velocities are in cm/s
                    :HD=>["f" 0.; "v" effusion_velocity(Temp_n(zmax), 3.0, zmax)],
                    :H=>["f" 0.; "v" effusion_velocity(Temp_n(zmax), 1.0, zmax)],
                    :D=>["f" 0.; "v" effusion_velocity(Temp_n(zmax), 2.0, zmax)],
                   );

    # Now plot the reactions ------------------------------------------------------------
    # looks for a pattern in the filename that looks like X{.XeX}, where curly braces are optional. 
    # i.e. it will match "100", "100.0", "0.001", "1e2" or "1.0e2".
    dtmatch = match(r"\d+\.*\d*e*\d*", f) 
    if typeof(dtmatch) != Nothing
        dtval = dtmatch.match
    else
        dtval = "unknown"
    end

    println("Working on dt=$(dtval)")
    for sp in all_species
        plot_rxns(sp, ncur, controltemps, speciesbclist; subfolder=parentfolder, plotsfolder="chemeq_plots", num=dtval, extra_title="dt=$(dtval)")
    end
end