###############################################################################
# converge_atmo_instantiate_ions.jl
# TYPE: (1) Model files - required
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere; this file adds in initial abundances for ions and converges
# the atmosphere so that they will take on their equilibrium profiles
#
# Eryn Cangi
# Created 2018
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
###############################################################################

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Photochemistry  # custom module for this project


include("PARAMETERS.jl")

################################################################################
#                                   FUNCTIONS                                  #
################################################################################

#=
These functions are required to be in this file for one of two reasons:
1) Because they call a function that is couched in an @eval statement, which
   cannot be relocated,
2) They are the main routine functions and will not be shared by any other scripts.

For additional functions, see the Photochemistry module.
=#

# chemistry functions ==========================================================

# TODO: This function isn't currently used. It was used with get_rates_and_fluxes in 
# mike's original code, which was responsible for printing reaction rates into the h5 file.
# It might be a good idea to write some new code to store this information in the converged
# atmosphere files just to have it.
# function reactionrates(n_current)
#     #=
#     Creates an array of size length(intaltgrid) x (number of reactions).
#     Populated with chemical reaction rates for each reaction based on species
#     populations.
#     =#
    
#     theserates = fill(convert(Float64, NaN), (length(intaltgrid), length(reactionnet)))
#     for ialt in 1:length(intaltgrid)
#         theserates[ialt,:] = reactionrates_local([[n_current[sp][ialt] for sp in specieslist];  # non-vectorized function call
#                                                 [n_current[J][ialt] for J in Jratelist];
#                                                 Temp_n(alt[ialt+1]); Temp_i(alt[ialt+1]); Temp_e(alt[ialt+1])]...)
#                                                 # n_tot(n_current, alt[ialt+1])]...)   # M is no longer pssed in.

                                                
#     end
#     return theserates
# end

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    # println("ratefn has been called") # DEBUG: remove
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive,(length(inactivespecies),length(intaltgrid)))
    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] = ratefn_local([nthismat[:,1]; nthismat[:,2];
                                    fill(1.0, length(activespecies));
                                    inactivemat[:,1]; Jrates[:,1]; Tn[1]; Ti[1]; Te[1];
                                    tup[:,1]; tlower[:,1]; tdown[:,2];
                                    tlower[:,2]]...)


    # iterate through other altitudes except the last level, filling the info in
    for ialt in 2:(length(intaltgrid)-1)
        returnrates[:,ialt] = ratefn_local([nthismat[:,ialt];
                                          nthismat[:,ialt+1];
                                          nthismat[:,ialt-1];
                                          inactivemat[:,ialt];
                                          Jrates[:,ialt];
                                          Tn[ialt]; Ti[ialt]; Te[ialt];
                                          tup[:,ialt]; tdown[:,ialt];
                                          tdown[:,ialt+1]; tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude
    returnrates[:,end] = ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       Tn[end]; Ti[end]; Te[end];
                                       tupper[:,1]; tdown[:,end];
                                       tupper[:,2]; tup[:,end-1]]...)
    return [returnrates...;]
end

function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper, dt)
    #=
    TODO: docstring
    =#

    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive, (length(inactivespecies), length(intaltgrid)))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    # tc___ are the coordinate tuples containing (I, J, V) to be use to fill a sparse matrix.
    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1]; nthismat[:,2]; fill(1.0, length(activespecies));
                                                  inactivemat[:,1]; Jrates[:,1];
                                                  Tn[1]; Ti[1]; Te[1]; 
                                                  tup[:,1]; tlower[:,1];
                                                  tdown[:,2]; tlower[:,2]; dt]...) 

    # add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])

    # and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(length(intaltgrid)-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt]; Tn[ialt]; Ti[ialt]; Te[ialt];
                                                      tup[:,ialt];
                                                      tdown[:,ialt];
                                                      tdown[:,ialt+1];
                                                      tup[:,ialt-1]; dt]...)
        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclocal[2].+(ialt-1)*length(activespecies))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tcupper[2].+(ialt  )*length(activespecies))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclower[2].+(ialt-2)*length(activespecies))
        append!(chemJval, tclower[3])
    end

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,end];
                                              fill(1.0, length(activespecies));
                                              nthismat[:,end-1];
                                              inactivemat[:,end];
                                              Jrates[:,end];
                                              Tn[end]; Ti[end]; Te[end]; 
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1]; dt]...)

    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclower[2].+(length(intaltgrid)-2)*length(activespecies))
    append!(chemJval, tclower[3])

    # make sure to add 1's along the diagonal
    append!(chemJi,[1:length(nthis);])
    append!(chemJj,[1:length(nthis);])
    append!(chemJval, fill(1.0, length(nthis)))

    return sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +);
end

# main routine functions =======================================================
function update_Jrates!(n_current::Dict{Symbol, Array{Float64, 1}}; plot_opt_depth=false, dt=nothing, iter=nothing, jrate_ext_to_plot=nothing)
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#

    # Initialize an array, length=num of altitude levels - 2.
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, num_layers)
    for i in range(1, length=num_layers)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_current[species][ialt]*dz
            # if jspecies==:JO2toOpO
            #     println(string("At alt = ",alt[ialt+1],
            #                    ", n_",species," = ",n_current[species][ialt],
            #                    ", jcolumn = ",jcolumn))
            #     println("and solarabs[ialt] is $(solarabs[ialt]) before we do axpy")
            #     readline(STDIN)
            # end

            # add the total extinction to solarabs:
            # multiplies air column density (N, #/cm^2) at all wavelengths by crosssection (σ)
            # to get optical depth (τ). This is an override of axpy! to use the
            # full arguments. For the equation Y' = alpha*X + Y:
            # ARG 1: n (length of arrays in ARGS 3, 5)
            # ARG 2: alpha, a scalar.
            # ARG 3: X, an array of length n.
            # ARG 4: the increment of the index values of X, maybe?
            # ARG 5: Y, an array of length n
            # ARG 6: increment of index values of Y, maybe?
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1,
                       solarabs[ialt],1)
        end
    end

    # solarabs now records the total optical depth of the atmosphere at
    # each wavelength and altitude

    # Some plotting functions. 
    if plot_opt_depth
        if dt==nothing || iter==nothing || jrate_ext_to_plot==nothing
            throw("Error! A plot of optical depth has been requested, but one of its necessary arguments is not provided. dt, iter, and jrate_ext_to_plot must have values. 
                   Currently: dt=$(dt), iter=$(iter), jrate_ext_to_plot=$(jrate_ext_to_plot). First two should be integers, last should be a list of species symbols.")
        end
        plot_extinction(solarabs, results_dir, dt, iter, tauonly=true) # plots tau (optical depth)
        # This one plots some Jrates and their absorption band features: 
        for jrate in jrate_ext_to_plot
            plot_extinction(solarabs, results_dir, dt, iter, xsect_info=Any[crosssection[jrate][1:124], string(absorber[jrate])], solflux=solarflux[:,2])
        end
    end

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2].*exp.(-solarabs[ialt])
    end

    # each species absorbs according to its cross section at each
    # altitude times the actinic flux.
    # BLAS.dot includes an integration (sum) across wavelengths, i.e:
    # (a·b) = aa + ab + ab + bb etc that kind of thing
    for j in Jratelist
        for ialt in [1:nalt;]
            n_current[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1,
                                          crosssection[j][ialt+1], 1)
        end
    end
end

function timeupdate(mytime, imgpath)
    # plot rates to inspect them! These are only used while slowly converging an atmosphere.
    # for s in new_ions
    #     plot_rxns(s, n_current, [T_surf, T_tropo, T_exo], speciesbclist; plot_indiv_rxns=false, subfolder="temp_216_130_205/", dt=mytime, num=PN)
    # end
    # for s in new_neutrals
    #     plot_rxns(s, n_current, [T_surf, T_tropo, T_exo], speciesbclist; plot_indiv_rxns=false, subfolder="temp_216_130_205/", dt=mytime, num=PN)
    # end
    # global PN += 1 # increment the number to append to plot filenames 

    write_ncurrent(n_current, results_dir*FNext*"/ncurrent_$(mytime).h5")  # used while slowly converging an atmosphere

    numiters = 15
    for i = 1:numiters
        update!(n_current, mytime)
        plot_atm(n_current, [neutrallist, ionlist], imgpath, t=mytime, iter=i)  # only plot at the end of the iterations, speeds up code given the way I changed plot_atm.
    end
    
end

function next_timestep(nstart::Array{Float64, 1}, nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1}, activespecies, inactivespecies,
                       Jrates::Array{Float64, 2},
                       Tn::Array{Float64, 1}, Ti::Array{Float64, 1}, Te::Array{Float64, 1},
                       tup::Array{Float64, 2}, tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2}, tupper::Array{Float64, 2},
                       dt::Float64)
    #=
    moves to the next timestep using Newton's method on the linearized
    coupled transport and chemical reaction network.
    =#

    eps = 1.0 # ensure at least one iteration
    iter = 0
    while eps>1e-8 
        nold = deepcopy(nthis)

        # stuff concentrations into update function and jacobian.
        ratefn_output = dt*ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te,  
                                          tup, tdown, tlower, tupper)  #Eulerian (P*L) * dt
        # println("Eulerian update (dt*ratefn / nthis): $(ratefn_output ./ nthis)")
        fval = nthis - nstart - ratefn_output
        # println("ratefn_output: $(ratefn_output)")
        
        updatemat = chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, 
                             tup, tdown, tlower, tupper, dt)
        # println("type updatemat: $(typeof(updatemat))")

        # update
        # updatemat_divided_by_fval = updatemat \ fval  # backslash - matrix solution operator

        # for i in range(1, length=length(nthis))
        #     if ratefn_output[i]/nthis[i] > 1e-6  # Eulerian update  
        #         nthis[i] = nthis[i] - updatemat_divided_by_fval[i]
        #     end
        # end
        nthis = nthis - (updatemat \ fval)#updatemat_divided_by_fval #
        # check relative size of update
        eps = maximum(abs.(nthis-nold)./nold)
        iter += 1
        if iter>2e3
            throw("too many iterations in next_timestep! eps: $(eps)")
        end
        # throw("Breakpoint within next_timestep to examine values")
    end

    return nthis
end

function update!(n_current::Dict{Symbol, Array{Float64, 1}}, dt; plotJratesflag=false, iter=nothing) 
    #=
    update n_current using the coupled reaction network, moving to
    the next timestep

    plotJratesflag: if set to true, it will plot the Jrates.
    =#

    # set auxiliary (not solved for in chemistry) species values, photolysis rates
    inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(intaltgrid)]...])
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(intaltgrid)])

    # extract concentrations
    nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(intaltgrid)]...])

    # plot J rates if required
    if plotJratesflag 
        plot_Jrates(n_current, [T_surf, T_tropo, T_exo], speciesbclist, filenameext="time=$(dt)")
    end

    # set temperatures
    Tn = Float64[Temp_n(a) for a in non_bdy_layers]
    Ti = Float64[Temp_i(a) for a in non_bdy_layers]
    Te = Float64[Temp_e(a) for a in non_bdy_layers]

    # println("The length of the temperature arrays is $(length(Tn))")

    # take initial guess
    nthis = deepcopy(nstart)

    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    tup = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current, [T_surf, T_tropo, T_exo])[2] for sp in specieslist, a in non_bdy_layers]
    tdown = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current, [T_surf, T_tropo, T_exo])[1] for sp in specieslist, a in non_bdy_layers]

    # put the lower layer and upper layer boundary conditions in separate arrays; but they are not the
    # right shape!
    tlower_temporary = [boundaryconditions(sp, dz, n_current, [T_surf, T_tropo, T_exo], speciesbclist)[1,:] for sp in specieslist]
    tupper_temporary = [boundaryconditions(sp, dz, n_current, [T_surf, T_tropo, T_exo], speciesbclist)[2,:] for sp in specieslist]

    # reshape tlower and tupper into 2x2 arrays
    tlower = zeros(Float64, length(tlower_temporary), 2)
    tupper = zeros(Float64, length(tupper_temporary), 2)

    # tlower_temporary & tupper_temporary have same length; OK to use lower for the range
    for r in range(1, length=length(tlower_temporary))
        tlower[r, :] = tlower_temporary[r]
        tupper[r, :] = tupper_temporary[r]
    end

    nthis = next_timestep(nstart, nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te,
                          tup, tdown, tlower, tupper, dt)
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))

    # write found values out to n_current
    for s in 1:length(activespecies)
        # set the updated altitude range
        altrange = get(update_rules, activespecies[s], 1:num_layers)
        for ia in altrange
            tn = nthismat[s, ia]
            n_current[activespecies[s]][ia] = tn > 0. ? tn : 0.  # prevents negative concentrations
        end
    end

    update_Jrates!(n_current)
end

################################################################################
#                                  MAIN SETUP                                  #
################################################################################

# Set up simulation files and experiment type ==================================

# Note: directory paths are in PARAMETERS.jl
# get command line arguments for experiment type. format:
# <experiment type> <parameters> <solar cycle type>
# examples:
# temp Tsurf Ttropo Texo mean; water <mixing ratio> mean; dh <multiplier> mean;
# Oflux <cm^-2s^-1> mean
# examples: temp 190 110 200 mean; water 1e-3 mean; dh 8 mean; Oflux 1.2e8 mean.
# last argument is solar cycle: min, mean, or max.
args = Any["temp", "216", "130", "205", "mean"]
# args = Any[ARGS[i] for i in 1:1:length(ARGS)]  #TODO: revert

# Establish a pattern for filenames. FNext = filename extension
# TODO: We probably won't be running all these experiments this time, so this section can probably go away.
# we may not even have to run from command line at all...
if args[1]=="temp"
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    FNext = "water_$(args[2])"
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
else
    throw("Error! Bad experiment type")
end

# Set up the folder if it doesn't exist
create_folder(FNext, results_dir)

# Case where this file will be used to converge an atmosphere of a different extent
# make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
make_new_alt_grid = "n" # TODO: revert
if make_new_alt_grid=="y"
    readfile = research_dir*"converged_200km_atmosphere.h5"
    const alt = convert(Array, (0:2e5:200e5))
    n_current = get_ncurrent(readfile)

    new_zmax = parse(Int64, input("Enter the new top of the atmosphere in km: "))
    extra_entries = Int64((new_zmax - 200)/(dz/1e5))

    # Extend the grid
    for (k,v) in zip(keys(n_current), values(n_current))
       append!(v, fill(v[end], extra_entries))  # repeats the last value in the array for the upper atmo as an initial value.
    end

    const alt = convert(Array, (0:2e5:new_zmax*10^5))

    if new_zmax != 250
        println("Warning: You entered $(new_zmax) for the new max altitude but
                 250 is hard-coded in the PARAMETERS.jl file. I haven't made this
                 general yet. So the code is probably about to break")
    end
elseif make_new_alt_grid=="n"
    # Set up the converged file to read from and load the simulation state at init.
    file_to_use = input("Enter the name of a file containing a converged, 250 km atmosphere to use (press enter to use default): ")   # TODO: revert
    readfile = file_to_use == "" ? "converged_250km_atmosphere.h5" : file_to_use
    n_current = get_ncurrent(readfile)

    converge_which = input("Converging ions or neutrals? (Used for setting timesteps): ")
    while converge_which != "neutrals" && converge_which != "ions" && converge_which != "both"
        println("Bad entry! Please enter ions or neturals or both.")
        global converge_which = input("Converging ions or neutrals? (Used for setting timesteps): ")
    end
    println()

    if converge_which == "neutrals" #|| converge_which == "both"
        mindt = -3
        maxdt = 14
        println("ALERT: Did you make sure that the ions are added to nochemspecies and notransportspecies?")
    elseif converge_which == "ions"
        mindt = -3
        maxdt = 5
        println("ALERT: Did you make sure that the neutrals are added to nochemspecies and notransportspecies?")
    elseif converge_which == "both"
        mindt = -5
        maxdt = 5
        println("ALERT: Did you make sure nochemspecies and notransportspecies are only :Ar and :N2?")
    else
        throw("Uncaught exception")
    end
else
    error("Didn't understand response")
end

# Whether to initialize new species as zeros or not
use_nonzero_initial_profiles = true

# Set solar cycle file and alert user
cycle = args[end]
solar_data_file = Dict("max"=>"marssolarphotonflux_solarmax.dat",
                       "mean"=>"marssolarphotonflux_solarmean.dat",
                       "min"=>"marssolarphotonflux_solarmin.dat")
solarfile = solar_data_file[cycle]

# Convert the arguments to numbers so we can use them to do maths
for i in 2:1:length(args)-1
    args[i] = parse(Float64, args[i])
end

# Let the user know what is being done
println("ALERT: running sim for $(FNext)")
println("ALERT: Using converged atmosphere file: ", readfile)
if cycle != "mean"
    println("ALERT: Solar $(cycle) data being used")
end

# Plot styles ==================================================================
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 18
rcParams["axes.labelsize"]= 20
rcParams["xtick.labelsize"] = 18
rcParams["ytick.labelsize"] = 18

################################################################################
#                             SET UP SPECIES PROFILES                          #
################################################################################

if args[1] == "dh"
    DH = parse(Float64, args[2]) * 1.6e-4
end

# new ions and neutrals from Roger's model
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
# establish zero profiles
for ion in new_ions
    n_current[ion] = zeros(num_layers)
end

# for nn in new_neutrals
#     n_current[nn] = zeros(num_layers)
# end

# if we want to use the output profiles from Roger's model as an initial guess, we can do that by just overwriting the zeros
if use_nonzero_initial_profiles
    println("Initializing profiles as non-zero")
    for ion in new_ions
        n_current[ion] = reshape(readdlm("../Resources/initial_profiles/$(string(ion))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
    end

    for nn in new_neutrals
        n_current[nn] = reshape(readdlm("../Resources/initial_profiles/$(string(nn))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
    end
end

# NEW: photoionization + photodissociation reactions from Roger
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
for Jrate in newJrates
    n_current[Jrate] = zeros(num_layers)
end

# new D neutrals: TODO: Find?

# new D ions:    # TODO: turn this on after testing that Roger's network works.
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
# n_current[:ArDpl] = n_current[:ArHpl] * DH
# n_current[:CDpl] = n_current[:CHpl] * DH
# n_current[:HDpl] = n_current[:H2pl] * DH
# n_current[:H2Dpl] = n_current[:H3pl] * DH
# n_current[:HD2pl] = n_current[:H3pl] * DH
# n_current[:DCOpl] = n_current[:HCOpl] * DH
# n_current[:Dpl] = n_current[:Hpl] * DH
# n_current[:N_2Dpl] = n_current[:N2Hpl] * DH
# n_current[:NDpl] = n_current[:NHpl] * DH
# n_current[:ODpl] = n_current[:OHpl] * DH
# n_current[:NHDpl] = zeros(length(alt))
# n_current[:NHD2pl] = zeros(length(alt))
# n_current[:ND3pl] = zeros(length(alt))
# n_current[:DOCOpl] = zeros(length(alt))

################################################################################
#                       TEMPERATURE/PRESSURE PROFILES                          #
################################################################################

#If changes to the temperature are needed, they should be made here
if args[1]=="temp"
    global T_surf = args[2]
    global T_tropo = args[3]
    global T_exo = args[4]
    Temp_n(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "neutral")
    Temp_i(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "ion")
    Temp_e(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "electron")
    Temp_keepSVP(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral") # for testing temp without changing SVP. #TODO: adjust if needed for ions?
else
    global T_surf = meanTs
    global T_tropo = meanTt
    global T_exo = meanTe
    Temp_n(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral")
    Temp_i(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "ion")
    Temp_e(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "electron")
end

# plot all 3 profiles on top of each other
plot_nie_temp_profs([Temp_n(a) for a in alt], [Temp_i(a) for a in alt], [Temp_e(a) for a in alt], results_dir*FNext, alt)

################################################################################
#                               WATER PROFILES                                 #
################################################################################

# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
else
    fix_SVP = false
    H2Osat = map(x->Psat(x), map(Temp_n, alt)) # array in #/cm^3 by altitude
    HDOsat = map(x->Psat_HDO(x), map(Temp_n, alt))
end

# H2O Water Profile ============================================================
surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo
# set H2O SVP fraction to minimum for all alts above first time min is reached
H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
               fill(minimum(H2Osatfrac), num_layers-length(H2Oinitfrac))]

# make profile constant in the lower atmosphere (well-mixed).
# when doing water experiments, the temperature profile is such that the minimum
# in the SVP curve occurs at about 55 km alt. This was found by manual tweaking.
# thus we need to only set the lower atmo mixing ratio below that point, or
# there will be a little spike in the water profile.
if args[1] == "water"
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= args[2]
    MR = args[2] # mixing ratio
else  # i believe 30 km is supposed to be approximately the hygropause.
    MR = MR_mean_water
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR # 10 pr μm
end

for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
end

# HDO water profile ============================================================
HDOsatfrac = HDOsat./map(z->n_tot(n_current, z), alt)
# use D/H ratio to set population of HDO
HDOinitfrac = H2Oinitfrac * DH  # initial profile for HDO

# Compute total water column in pr μm starting with mixing ratio array =========
# H2O #/cm^3 = sum(MR * n_tot) for each alt
H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))

# pr μm = (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23,
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

################################################################################
#                             BOUNDARY CONDITIONS                              #
################################################################################

H_veff = effusion_velocity(Temp_n(zmax), 1.0, zmax)
H2_veff = effusion_velocity(Temp_n(zmax), 2.0, zmax)
D_veff = effusion_velocity(Temp_n(zmax), 2.0, zmax)
HD_veff = effusion_velocity(Temp_n(zmax), 3.0, zmax)

#=
    boundary conditions for each species (mostly from Nair 1994, Yung 1988). For
    most species, default boundary condition is zero (f)lux at top and bottom.
    Atomic/molecular hydrogen and deuterated analogues have a nonzero effusion
    (v)elocity at the upper layer of the atmosphere. Some species have a (n)umber
    density boundary condition.
=#

if args[1] == "Oflux"
    global const Of = args[2]
else
    global const Of = 1.2e8
end

# This has to be defined here, because it uses as a boundary condition the H2O
# and HDO saturation at the surface.
global const speciesbclist=Dict(
                :CO2=>["n" 2.1e17; "f" 0.],
                :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                :HDO=>["n" HDOsat[1]; "f" 0.],
                :O=>["f" 0.; "f" Of],
                :H2=>["f" 0.; "v" H2_veff],  # velocities are in cm/s
                :HD=>["f" 0.; "v" HD_veff],
                :H=>["f" 0.; "v" H_veff],
                :D=>["f" 0.; "v" D_veff],
                # TODO: Ion boundary conditions?
               );

################################################################################
#                       COMBINED CHEMISTRY AND TRANSPORT                       #
################################################################################

#=
    We now have objects that return the list of indices and coefficients
    for transport, assuming no other species in the atmosphere
    (transportmat), and for chemistry, assuming no other altitudes
    (chemical_jacobian). We need to perform a kind of outer product on
    these operators, to determine a fully coupled set of equations for
    all species at all altitudes.

    the rates at each altitude can be computed using the reaction network
    already in place, plus additional equations describing the transport
    to and from the cells above and below:
=#
upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
          Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
          for s in specieslist]

downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
            Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
            for s in specieslist]

local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in specieslist]
                          [Symbol("t"*string(s)*"_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_above_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_below_up") for s in specieslist]]...;]

transportnet = [[upeqns...;]; [downeqns...;]]

# define names for all the species active in the coupled rates:
active_above = [Symbol(string(s)*"_above") for s in activespecies]
active_below = [Symbol(string(s)*"_below") for s in activespecies]

# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x), activespecies)...);

# TODO: Turn these notes into a test that looks for "+()" in rates_local instead of debugging lines. 
# -----------------------------------------------------------------------------------------------
# These lines are useful for troubleshooting if you're getting weird errors with ratefn.
# "no method matching +()" means the chemical system is unbalanced and you've got a
# species that has production but no consumption, or vice versa.
# other errors may occur. Uncomment these 3 lines to inspect what goes into rates_local.
# println("The contents of rates_local: ")
# println(rates_local)
# println("That was rates_local")

const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

arglist_local = [activespecies; active_above; active_below; inactivespecies;
                 Jratelist; :Tn; :Ti; :Te; local_transport_rates; :dt]

arglist_local_typed = [:($s::Float64) for s in arglist_local]

Mexpr = Expr(:call, :+, fullspecieslist...)
Eexpr = Expr(:call, :+, ionlist...)

# NOTE: These functions within @eval cannot be moved. Do not move them.
@eval begin
    function ratefn_local($(arglist_local_typed[1:end-1]...))  # mike's suggestion to vectorize: (arglist_vector_typed::Array{Float64, 1})
        #=
        Produces the symbolic expression for production and loss of every species
        via both transport and chemistry
        =#

        # stuff mike suggested for debugging compiler errors:
        # use this for the argument: (arglist_vector_typed::Array{Float64, 1})
        # unpack vector 
        #$(arglist_local_typed[1:end-1]...) = arglist_vector_typed
        # $(rates_local) 
        # end mike's suggestions

        # M and E are calculated here to ensure that the right number of ions/electrons
        # is used. It is for only the altitude at which this function was called 
        # (i.e. all the arguments to the function, when it's called, are for only one altitude)
        M = $Mexpr
        E = $Eexpr
        $rates_local # evaluates the rates_local expression
    end
end

@eval begin
    function chemJmat_local($(arglist_local_typed...))

        M = $Mexpr
        E = $Eexpr

        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])

        # section to populate localchemJval element by element,
        # which is having problems with compiler errors when ----->
        # localchemJval = similar($chemJ_local[3], Float64)
        #     # unroll the loop and evaluate each expression 1 by 1 thank you 张实唯 of SE
        #     $((quote
        #         localchemJval[$i] = $(chemJ_local[3][i])
        #     end for i in 1:length(chemJ_local[3]))...)
        # localchemJval = -dt*localchemJval  # the old code that throws compiler error:*$(chemJ_local[3])  
        # <----- end section
        localchemJval = -dt*$(Expr(:vcat, chemJ_local[3]...)) 

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = -dt*$(Expr(:vcat, chemJ_above[3]...))#(chemJ_above[3])#

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(Expr(:vcat, chemJ_below[3]...))#(chemJ_below[3])#

        # return the actual values of I, J, and V (indices and numerical value):
        return ((localchemJi, localchemJj, localchemJval),
                (abovechemJi, abovechemJj, abovechemJval),
                (belowchemJi, belowchemJj, belowchemJval))
    end
end

# TODO: this function is not currently used. It's used with reactionrates to 
# write reaction rates to the .h5 file, originally done in Mike's version of the 
# model. It would be useful to include this functionality in this version of the code.
# @eval begin
#     function reactionrates_local($(specieslist...), $(Jratelist...), Tn, Ti, Te)#, M, E)  # vectorized: (spclist_vector, Jratelist_vector, Tn, Ti, Te)
         
#           a function to return chemical reaction rates by multiplying the reaction rate coefficient
#           (represented by x[3]) by the concentrations of all reactants (represented by x[1]...).
#           the concentrations and Jrates are contained in specieslist and Jratelist, which are flattened over altitude
#           using the splats (I think).
           

#         # unpack the vectors
#         #$(specieslist...) = spclist_vector
#         #$(Jratelist...) = Jratelist_vector

#         M = eval(Expr(:call, :+, $(fullspecieslist...)))
#         E = eval(Expr(:call, :+, $(ionlist...)))
#         $(Expr(:vcat, map(x->Expr(:call,:*,x[1]..., x[3]), reactionnet)...))
#     end
# end

###############################################################################  
#                        PHOTOCHEMICAL CROSS SECTIONS                         #
###############################################################################

xsecfolder = research_dir * "uvxsect/";

# Crosssection Files ===========================================================
co2file = "CO2.dat"
# co2exfile = "binnedCO2e.csv" # added to shield short λ of sunlight in upper atmo
h2ofile = "h2oavgtbl.dat"
hdofile = "HDO.dat"#"HDO_250K.dat"#
h2o2file = "H2O2.dat"
hdo2file = "H2O2.dat" #TODO: do HDO2 xsects exist?
o3file = "O3.dat"
o3chapfile = "O3Chap.dat"
o2file = "O2.dat"
o2_130_190 = "130-190.cf4"
o2_190_280 = "190-280.cf4"
o2_280_500 = "280-500.cf4"
h2file = "binnedH2.csv"
hdfile = "binnedH2.csv" # TODO: change this to HD file if xsects ever exist
ohfile = "binnedOH.csv"
oho1dfile = "binnedOHo1D.csv"
odfile = "OD.csv"

# NEW: files for crosssections of specific photodissociation/photoionization
# Source: Roger Yelle
# NOTE: Order has changed compared to the order in the Jrates list and the absorber dict,
# for ease of reading.
H2O_ionize_file = "JH2OtoH2Opl.csv"
H2O_ionOdiss_file = "JH2OtoOplpH2.csv"
H2O_ionHdiss_file = "JH2OtoHplpOH.csv"
H2O_ionOHdiss_file = "JH2OtoOHplpH.csv"

CO_diss_file = "JCOtoCpO.csv"
CO_ionize_file = "JCOtoCOpl.csv"
CO_ionOdiss_file = "JCOtoCpOpl.csv"
CO_ionCdiss_file = "JCOtoOpCpl.csv"

N2_diss_file = "JN2OtoN2pO1D.csv"
N2_ionize_file = "JN2toN2pl.csv"
N2_iondiss_file = "JN2toNplpN.csv"
NO2_diss_file = "JNO2toNOpO.csv"
NO2_ionize_file = "JNO2toNO2pl.csv"
NO_diss_file = "JNOtoNpO.csv"
NO_ionize_file = "JNOtoNOpl.csv"
N2O_ionize_file = "JN2OtoN2Opl.csv"

H_ionize_file = "JHtoHpl.csv"
H2_ion_file = "JH2toH2pl.csv"
H2_iondiss_file = "JH2toHplpH.csv"

CO2_totaldiss_file = "JCO2toCpOpO.csv"
CO2_diss_file = "JCO2toCpO2.csv"
CO2_ionize_file = "JCO2toCO2pl.csv"
CO2_doubleion_file = "JCO2toCO2plpl.csv"
CO2_ionC2diss_file = "JCO2toCplplpO2.csv"
CO2_ionCdiss_file = "JCO2toCplpO2.csv"
CO2_ionCOandOdiss_file = "JCO2toCOplpOpl.csv"
CO2_ionCOdiss_file = "JCO2toCOplpO.csv"
CO2_ionOdiss_file = "JCO2toOplpCO.csv"
CO2_ionCandOdiss_file = "JCO2toOplpCplpO.csv"

H2O2_ionize_file = "JH2O2toH2O2pl.csv"

O_ionize_file = "JOtoOpl.csv"
O2_ionize_file = "JO2toO2pl.csv"#"O2-photoionize-swri.csv"#
O3_ionize_file = "JO3toO3pl.csv"

# Loading Data =================================================================
# CO2 photodissociation --------------------------------------------------------
# temperature-dependent between 195-295K
co2xdata = readdlm(xsecfolder*co2file,'\t', Float64, comments=true, comment_char='#')

# CO2 photoionization (used to screen high energy sunlight)
# co2exdata = readdlm(xsecfolder*co2exfile,',',Float64, comments=true, comment_char='#')

# H2O & HDO --------------------------------------------------------------------
h2oxdata = readdlm(xsecfolder*h2ofile,'\t', Float64, comments=true, comment_char='#')

# These crosssections for HDO are for 298K.
hdoxdata = readdlm(xsecfolder*hdofile,'\t', Float64, comments=true, comment_char='#')

# H2O2 + HDO2 ------------------------------------------------------------------
# the data in the following table cover the range 190-260nm
h2o2xdata = readdlm(xsecfolder*h2o2file,'\t', Float64, comments=true, comment_char='#')
hdo2xdata = readdlm(xsecfolder*hdo2file,'\t', Float64, comments=true, comment_char='#')

# O3 ---------------------------------------------------------------------------
# including IR bands which must be resampled from wavenumber
o3xdata = readdlm(xsecfolder*o3file,'\t', Float64, comments=true, comment_char='#')
o3ls = o3xdata[:,1]
o3xs = o3xdata[:,2]
o3chapxdata = readdlm(xsecfolder*o3chapfile,'\t', Float64, comments=true, comment_char='#')
o3chapxdata[:,1] = map(p->1e7/p, o3chapxdata[:,1])
for i in [round(Int, floor(minimum(o3chapxdata[:,1]))):round(Int, ceil(maximum(o3chapxdata))-1);]
    posss = getpos(o3chapxdata, x->i<x<i+1)
    dl = diff([map(x->o3chapxdata[x[1],1], posss); i])
    x = map(x->o3chapxdata[x[1],2],posss)
    ax = reduce(+,map(*,x, dl))/reduce(+,dl)
    global o3ls = [o3ls; i+0.5]
    global o3xs = [o3xs; ax]
end
o3xdata = reshape([o3ls; o3xs],length(o3ls),2)

# O2 ---------------------------------------------------------------------------
o2xdata = readdlm(xsecfolder*o2file,'\t', Float64, comments=true, comment_char='#')
o2schr130K = readdlm(xsecfolder*o2_130_190,'\t', Float64, comments=true, comment_char='#')
o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
o2schr130K = binupO2(o2schr130K)
o2schr190K = readdlm(xsecfolder*o2_190_280,'\t', Float64, comments=true, comment_char='#')
o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
o2schr190K = binupO2(o2schr190K)
o2schr280K = readdlm(xsecfolder*o2_280_500,'\t', Float64, comments=true, comment_char='#')
o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
o2schr280K = binupO2(o2schr280K)

# HO2 & DO2 --------------------------------------------------------------------
ho2xsect = [190.5:249.5;]
ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)
do2xsect = deepcopy(ho2xsect)

# H2 & HD ----------------------------------------------------------------------
h2xdata = readdlm(xsecfolder*h2file,',',Float64, comments=true, comment_char='#')
hdxdata = readdlm(xsecfolder*hdfile,',',Float64, comments=true, comment_char='#')

# OH & OD ----------------------------------------------------------------------
ohxdata = readdlm(xsecfolder*ohfile,',',Float64, comments=true, comment_char='#')
ohO1Dxdata = readdlm(xsecfolder*oho1dfile,',',Float64, comments=true, comment_char='#')
odxdata = readdlm(xsecfolder*odfile,',',Float64, comments=true, comment_char='#')

# PHOTODISSOCIATION ============================================================

const solarflux=readdlm(research_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]  # cut off at 2000 nm
solarflux[:,2] = solarflux[:,2]/2  # To roughly put everything at an SZA=60° (from a Kras comment)

absorber = Dict(#:JCO2ion =>:CO2,
                :JCO2toCOpO=>:CO2,
                :JCO2toCOpO1D=>:CO2,
                :JO2toOpO=>:O2,
                :JO2toOpO1D=>:O2,
                :JO3toO2pO=>:O3,
                :JO3toO2pO1D=>:O3,
                :JO3toOpOpO=>:O3,
                :JH2toHpH=>:H2,
                :JHDtoHpD=>:HD,
                :JOHtoOpH=>:OH,
                :JOHtoO1DpH=>:OH,
                :JODtoOpD=>:OD,
                :JODtoO1DpD=>:OD,
                :JHO2toOHpO=>:HO2,
                :JDO2toODpO=>:DO2,
                :JH2OtoHpOH=>:H2O,
                :JH2OtoH2pO1D=>:H2O,
                :JH2OtoHpHpO=>:H2O,
                :JH2O2to2OH=>:H2O2,
                :JH2O2toHO2pH=>:H2O2,
                :JH2O2toH2OpO1D=>:H2O2,
                :JHDO2toHDOpO1D=>:HDO2,
                :JHDOtoHpOD=>:HDO,
                :JHDO2toOHpOD=>:HDO2,
                :JHDO2toDO2pH=>:HDO2,
                :JHDO2toHO2pD=>:HDO2,
                :JHDOtoDpOH=>:HDO,
                :JHDOtoHpDpO=>:HDO,
                :JHDOtoHDpO1D=>:HDO,
                # NEW: reactions from Roger's model. 
                :JH2OtoH2Opl=>:H2O,
                :JH2OtoOplpH2=>:H2O,
                :JCOtoCpO=>:CO,
                :JCOtoCOpl=>:CO,
                :JN2OtoN2pO1D =>:N2O,
                :JH2toH2pl=>:H2,
                :JCOtoCpOpl=>:CO,
                :JNO2toNOpO=>:NO2,
                :JCO2toCpO2=>:CO2,
                :JCO2toCplplpO2=>:CO2,
                :JNOtoNOpl=>:NO,
                :JH2toHplpH=>:H2,
                :JH2OtoHplpOH=>:H2O,
                :JH2O2toH2O2pl=>:H2O2,
                :JN2toN2pl=>:N2,
                :JCO2toCOplpOpl=>:CO2,
                :JCOtoOpCpl=>:CO,
                :JCO2toOplpCplpO=>:CO2,
                :JNOtoNpO=>:NO,
                :JCO2toCplpO2=>:CO2,
                :JCO2toCO2pl=>:CO2,
                :JOtoOpl=>:O,
                :JH2OtoOHplpH=>:H2O,
                :JNO2toNO2pl=>:NO2,
                :JCO2toCOplpO=>:CO2,
                :JN2toNplpN=>:N2,
                :JCO2toCpOpO=>:CO2,
                :JCO2toCO2plpl=>:CO2,
                :JCO2toOplpCO=>:CO2,
                :JO2toO2pl=>:O2,
                :JHtoHpl=>:H,
                :JN2OtoN2Opl=>:N2O,
                :JO3toO3pl=>:O3
                );

# Enter crosssections into master dictionary ==================================
#=
    this is a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, crosssection has already folded in quantum
    efficiency considerations.
=#
crosssection = Dict{Symbol, Array{Array{Float64}}}()

# CO2 photodissociation --------------------------------------------------------
# setindex!(crosssection, fill(co2exdata, length(alt)), :JCO2ion)  # photoionization

#CO2+hv->CO+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))),
          map(t->co2xsect(co2xdata, t),map(Temp_n, alt))), :JCO2toCOpO)
#CO2+hv->CO+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))),
          map(t->co2xsect(co2xdata, t),map(Temp_n, alt))), :JCO2toCOpO1D)

# O2 photodissociation ---------------------------------------------------------
#O2+hv->O+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), map(Temp_n, alt))),
          :JO2toOpO)
#O2+hv->O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), map(Temp_n, alt))),
          :JO2toOpO1D)

# O3 photodissociation ---------------------------------------------------------
# O3+hv->O2+O
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1-(1.37e-2*193-2.16)),
                               (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
                               (l->225<=l<306, 0.1),
                               (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
                               (l->328<=l<340, 0.92),
                               (l->340<=l, 1.0)
                              )), map(Temp_n, alt)), :JO3toO2pO)
# O3+hv->O2+O1D
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1.37e-2*193-2.16),
                               (l->193<=l<225, l->(1.37e-2*l-2.16)),
                               (l->225<=l<306, 0.9),
                               (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                               (l->328<=l<340, 0.08),
                               (l->340<=l, 0.0)
                              )), map(Temp_n, alt)), :JO3toO2pO1D)
# O3+hv->O+O+O
setindex!(crosssection,
          fill(quantumyield(o3xdata,((x->true, 0.),)),length(alt)),
          :JO3toOpOpO)

# H2 and HD photodissociation --------------------------------------------------
# H2+hv->H+H
setindex!(crosssection, fill(h2xdata, length(alt)), :JH2toHpH)
# HD+hν -> H+D
setindex!(crosssection, fill(hdxdata, length(alt)), :JHDtoHpD)

# OH and OD photodissociation --------------------------------------------------
# OH+hv->O+H
setindex!(crosssection, fill(ohxdata, length(alt)), :JOHtoOpH)
# OH+hv->O1D+H
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JOHtoO1DpH)
# OD + hv -> O+D
setindex!(crosssection, fill(odxdata, length(alt)), :JODtoOpD)
# OD + hν -> O(¹D) + D
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JODtoO1DpD)

# HO2 and DO2 photodissociation ------------------------------------------------
# HO2 + hν -> OH + O
setindex!(crosssection, fill(ho2xsect, length(alt)), :JHO2toOHpO)
# DO2 + hν -> OD + O
setindex!(crosssection, fill(do2xsect, length(alt)), :JDO2toODpO)

# H2O and HDO photodissociation ------------------------------------------------
# H2O+hv->H+OH
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),length(alt)),
          :JH2OtoHpOH)

# H2O+hv->H2+O1D
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JH2OtoH2pO1D)

# H2O+hv->H+H+O
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->true, 0),)),length(alt)),
          :JH2OtoHpHpO)

# HDO + hν -> H + OD
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoHpOD)

# HDO + hν -> D + OH
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoDpOH)

# HDO + hν -> HD + O1D
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JHDOtoHDpO1D)

# HDO + hν -> H + D + O
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->true, 0),)),length(alt)),
          :JHDOtoHpDpO)

# H2O2 and HDO2 photodissociation ----------------------------------------------
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->h2o2xsect(h2o2xdata, t), map(Temp_n, alt))), :JH2O2to2OH)

# H2O2+hv->HO2+H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))),
          map(t->h2o2xsect(h2o2xdata, t), map(Temp_n, alt))), :JH2O2toHO2pH)

# H2O2+hv->H2O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(h2o2xdata, t),
          map(Temp_n, alt))), :JH2O2toH2OpO1D)

# HDO2 + hν -> OH + OD
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alt))), :JHDO2toOHpOD)

# HDO2 + hν-> DO2 + H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alt))), :JHDO2toDO2pH)

# HDO2 + hν-> HO2 + D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alt))), :JHDO2toHO2pD)

# HDO2 + hν -> HDO + O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(hdo2xdata, t),
          map(Temp_n, alt))), :JHDO2toHDOpO1D)

# NEW: CO2 photodissociation ---------------------------------------------------------
# Source: Roger Yelle
# CO₂ + hν -> C + O + O
CO2_totaldiss_data = readdlm(xsecfolder*CO2_totaldiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_totaldiss_data, length(alt)), :JCO2toCpOpO)

# CO2 + hν -> C + O₂
CO2_diss_data = readdlm(xsecfolder*CO2_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_diss_data, length(alt)), :JCO2toCpO2)

# NEW: CO photodissociation ---------------------------------------------------------
# Source: Roger Yelle

# CO + hν -> C + O
CO_diss_data = readdlm(xsecfolder*CO_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_diss_data, length(alt)), :JCOtoCpO)


# NEW: Nitrogen species photodissociation --------------------------------------------
# Source: Roger Yelle

# N₂ + hν -> N₂ + O(¹D)
N2_diss_data = readdlm(xsecfolder*N2_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2_diss_data, length(alt)), :JN2OtoN2pO1D)

# NO₂ + hν -> NO + O
NO2_diss_data = readdlm(xsecfolder*NO2_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO2_diss_data, length(alt)), :JNO2toNOpO)

# NO + hν -> N + O
NO_diss_data = readdlm(xsecfolder*NO_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO_diss_data, length(alt)), :JNOtoNpO)

# Photoionization or ionizing dissociation reactions ============================================

# NEW: CO₂ ionization -----------------------------------------------------------------
# Source: Roger Yelle

# CO₂ + hν -> CO₂⁺
CO2_ionize_data = readdlm(xsecfolder*CO2_ionize_file, ',', Float64, comments=true, comment_char='#')  # NOTE: replaced with Mike's file 19-Jan-2021.
setindex!(crosssection, fill(CO2_ionize_data, length(alt)), :JCO2toCO2pl)

# CO₂ + hν -> CO₂²⁺  (even though we don't track doubly ionized CO₂)
CO2_doubleion_data = readdlm(xsecfolder*CO2_doubleion_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_doubleion_data, length(alt)), :JCO2toCO2plpl)

# CO₂ + hν -> C²⁺ + O₂
CO2_ionC2diss_data = readdlm(xsecfolder*CO2_ionC2diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionC2diss_data, length(alt)), :JCO2toCplplpO2)

# CO₂ + hν -> C⁺ + O₂
CO2_ionCdiss_data = readdlm(xsecfolder*CO2_ionCdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCdiss_data, length(alt)), :JCO2toCplpO2)

# CO₂ + hν -> CO⁺ + O⁺
CO2_ionCOandOdiss_data = readdlm(xsecfolder*CO2_ionCOandOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCOandOdiss_data, length(alt)), :JCO2toCOplpOpl)

# CO₂ + hν -> CO⁺ + O
CO2_ionCOdiss_data = readdlm(xsecfolder*CO2_ionCOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCOdiss_data, length(alt)), :JCO2toCOplpO)

# CO₂ + hν -> CO + O⁺
CO2_ionOdiss_data = readdlm(xsecfolder*CO2_ionOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionOdiss_data, length(alt)), :JCO2toOplpCO)

# CO₂ + hν -> C⁺ + O⁺ + O
CO2_ionCandOdiss_data = readdlm(xsecfolder*CO2_ionCandOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCandOdiss_data, length(alt)), :JCO2toOplpCplpO)

# NEW: H2O ionization --------------------------------------------------------------
# Source: Roger Yelle

# H2O + hν -> H2O⁺
h2o_ionize_data = readdlm(xsecfolder*H2O_ionize_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionize_data, length(alt)), :JH2OtoH2Opl)

# H2O + hν -> O⁺ + H2
h2o_ionOdiss_data = readdlm(xsecfolder*H2O_ionOdiss_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionOdiss_data, length(alt)), :JH2OtoOplpH2)

# # H2O + hν -> H⁺ + OH
h2o_ionHdiss_data = readdlm(xsecfolder*H2O_ionHdiss_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionHdiss_data, length(alt)), :JH2OtoHplpOH)

# # H2O + hν -> OH⁺ + H
h2o_ionOHdiss_data = readdlm(xsecfolder*H2O_ionOHdiss_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionOHdiss_data, length(alt)), :JH2OtoOHplpH)

# # NEW: CO ionization ----------------------------------------------------------------
# # Source: Roger Yelle

# # CO + hν -> CO⁺
CO_ionize_data = readdlm(xsecfolder*CO_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_ionize_data, length(alt)), :JCOtoCOpl)

# # CO + hν -> C + O⁺
CO_ionOdiss_data = readdlm(xsecfolder*CO_ionOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_ionOdiss_data, length(alt)), :JCOtoCpOpl)

# # CO + hν -> C⁺ + O
CO_ionCdiss_data = readdlm(xsecfolder*CO_ionCdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_ionCdiss_data, length(alt)), :JCOtoOpCpl)

# NEW: Nitrogen species ionization --------------------------------------------------
# Source: Roger Yelle

# # N₂ + hν -> N₂⁺
N2_ionize_data = readdlm(xsecfolder*N2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2_ionize_data, length(alt)), :JN2toN2pl)

# # N₂ + hν -> N⁺ + N
N2_iondiss_data = readdlm(xsecfolder*N2_iondiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2_iondiss_data, length(alt)), :JN2toNplpN)

# # NO₂ + hν -> NO₂⁺
NO2_ionize_data = readdlm(xsecfolder*NO2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO2_ionize_data, length(alt)), :JNO2toNO2pl)

# # NO + hν -> NO⁺
NO_ionize_data = readdlm(xsecfolder*NO_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO_ionize_data, length(alt)), :JNOtoNOpl)

# # N₂O + hν -> N₂O⁺
N2O_ionize_data = readdlm(xsecfolder*N2O_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2O_ionize_data, length(alt)), :JN2OtoN2Opl)

# # NEW: Molecular and atomic hydrogen ionization -------------------------------------
# # Source: Roger Yelle

# # H + hν -> H⁺
H_ionize_data = readdlm(xsecfolder*H_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H_ionize_data, length(alt)), :JHtoHpl)

# # H₂ + hν -> H₂⁺
H2_ion_data = readdlm(xsecfolder*H2_ion_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H2_ion_data, length(alt)), :JH2toH2pl)

# # H₂ + hν -> H⁺ + H
H2_iondiss_data = readdlm(xsecfolder*H2_iondiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H2_iondiss_data, length(alt)), :JH2toHplpH)

# # H₂O₂ + hν -> H₂O₂⁺
H2O2_ionize_data = readdlm(xsecfolder*H2O2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H2O2_ionize_data, length(alt)), :JH2O2toH2O2pl)

# # NEW: Oxygen and ozone ionization --------------------------------------------------
# # Source: Roger Yelle

# O + hν -> O⁺
O_iondiss_data = readdlm(xsecfolder*O_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(O_iondiss_data, length(alt)), :JOtoOpl)

# O₂ + hν -> O₂⁺
O2_ionize_data = readdlm(xsecfolder*O2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(O2_ionize_data, length(alt)), :JO2toO2pl)

# # O₃ + hν -> O₃⁺
O3_ionize_data = readdlm(xsecfolder*O3_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(O3_ionize_data, length(alt)), :JO3toO3pl)


# Solar Input ==================================================================
lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

if !(setdiff(solarflux[:,1],lambdas)==[])
    throw("Need a broader range of solar flux values!")
end

# pad all cross-sections to solar
for j in Jratelist, ialt in 1:length(alt)
    crosssection[j][ialt] = padtosolar(solarflux, crosssection[j][ialt])
end

# this is the unitialized array for storing values
solarabs = fill(fill(0.,size(solarflux, 1)), num_layers);

################################################################################
#                                 LOGGING                                      #
################################################################################

# crosssection dict for logging purposes =======================================
xsect_dict = Dict("CO2"=>[co2file],
                  "H2O, HDO"=>[h2ofile, hdofile],
                  "H2O2, HDO2"=>[h2o2file, hdo2file],
                  "O3"=>[o3file, o3chapfile],
                  "O2"=>[o2file, o2_130_190, o2_190_280, o2_280_500],
                  "H2, HD"=>[h2file, hdfile],
                  "OH, OD"=>[ohfile, oho1dfile, odfile],
                  "CO2pl"=>[CO2_diss_file, CO2_ionOdiss_file],
                  "O2pl"=>[O2_ionize_file],
                  "Opl"=>[O_ionize_file])

# Log temperature and water parameters =========================================
if args[1]=="temp"
    input_string = "T_s=$(args[2]), T_tropo=$(args[3]), T_exo=$(args[4])" *
                   "\nwater init=$(MR)\nDH=5.5 \nOflux=1.2e8"
elseif args[1]=="water"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\n" *
                   "water init=$(args[2])\nDH=5.5\nOflux=1.2e8\n"
elseif args[1]=="dh"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)\n" *
                   "DH=$(args[2]) \nOflux=1.2e8\n"
elseif args[1]=="Oflux"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)" *
                   "\nDH=5.5\nOflux=$(Of)\n"
end

# Write the log ================================================================
f = open(results_dir*FNext*"/simulation_params_"*FNext*".txt", "w")
write(f, "$(args[1]) experiment: \n")
write(f, input_string)
write(f, "\n\n")

# Mean temperatures
write(f, "Mean temperatures used:\n")
write(f, "Surface: $(meanTs) K, Tropopause: $(meanTt) K, Exobase: $(meanTe) K\n\n")

# which species are turned on
write(f, "All species: $(join([string(i) for i in fullspecieslist], ", "))\n")
write(f, "No-chem species: $(join([string(i) for i in nochemspecies], ", "))\n")
write(f, "No-transport species: $(join([string(i) for i in notransportspecies], ", "))\n\n")

# Which ion reactions are turned on 
write(f, "New J rates that are turned on: \n")
for rxn in [string(i) for i in newJrates]
    write(f, rxn*"\n")
end
write(f, "\n")

# gets just chemistry reactions with ions as reactants or products
ion_chem_rxns = filter(x->(occursin("pl", string(x[1])) || occursin("pl", string(x[2]))), filter(x->!occursin("J", string(x[3])), reactionnet)) 
write(f, "Ion chemistry reactions: \n")
for rxn in [string(i) for i in ion_chem_rxns]
    write(f, rxn*"\n")
end
write(f, "\n")

# Whether nonzero initial profiles were used
write(f, "Non-zero initial profiles used: $(use_nonzero_initial_profiles)\n\n")

# Whether we fixed things below a certain altitude
write(f, "Altitudes able to have ion densities updated: $(get(update_rules, "ion", "all"))\n\n")

# Solar cycle and SVP stuff
write(f, "\nSVP fixed: $(fix_SVP)\n\n")
write(f, "Solar cycle status: solar $(cycle)\n\n")

# Water profile information
write(f, "Water profile information: \n")
write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
write(f, "Total H2O (pr μm): $(H2Oprum)\n")
write(f, "Total HDO (pr μm): $(HDOprum)\n")
write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)\n\n")


# cross sections
write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
write(f, "\n")

# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="),
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end
write(f, "\n")

# timestep etc
write(f, "Initial atmosphere state: $(readfile)\n")

close(f)

################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# set the water profiles =======================================================
n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers)
n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers)

# Plot initial water profile ===================================================
plot_water_profile(H2Oinitfrac, HDOinitfrac, n_current[:H2O], n_current[:HDO], results_dir*FNext, watersat=[H2Osatfrac, HDOsatfrac])

# do the convergence ===========================================================

println("Plotting the initial condition")
plot_atm(n_current, [neutrallist, ionlist], results_dir*FNext*"/initial_atmosphere.png")

img_savepath = results_dir*FNext*"/converged_atm_"*FNext*".png"
println("Beginning Convergence")

@showprogress 0.1 "Converging..." [timeupdate(t, results_dir*FNext*"/converged_atm_"*FNext*".png") for t in [10.0^(1.0*i) for i in mindt:maxdt]]

if converge_which == "neutrals"
    @showprogress 0.1 "Last convergence steps..." for i in 1:100 
        plot_atm(n_current, [neutrallist, ionlist], img_savepath, t="1e14", iter=i)
        # println("dt: 1e14 iter $(i)")
        update!(n_current, 1e14)
    end
end

for s in fullspecieslist
    plot_rxns(s, n_current, [T_surf, T_tropo, T_exo], speciesbclist; plot_indiv_rxns=false, subfolder=FNext, extra_title="converged")
end

# write out the new converged file to matching folder.
towrite = results_dir*FNext*"/converged_"*FNext*".h5"
write_ncurrent(n_current, towrite)
println("Wrote $(towrite)")

# save the figure
# savefig(img_savepath, bbox_inches="tight")
println("Saved figure to same folder")

println("ALERT: Finished")
println()
