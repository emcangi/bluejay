###############################################################################
# converge_new_file.jl
# TYPE: (1) Model files - required
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere.
# 
# Eryn Cangi
# Created 2018
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
###############################################################################

t1 = time()

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

t2 = time()

println("Time to load modules: $(round(t2-t1, digits=1)) seconds")

user_input_paramfile = input("Enter a parameter file or press enter to use default (PARAMETERS.jl): ")
paramfile = user_input_paramfile == "" ? "PARAMETERS.jl" : user_input_paramfile
t3 = time()
include(paramfile)
t4 = time()
println("Time to load PARAMETERS: $(round(t4-t3, digits=1)) seconds")

t5 = time()

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

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] = ratefn_local([nthismat[:,1]; nthismat[:,2];
                                    fill(1.0, length(activespecies));
                                    inactivemat[:,1]; Jrates[:,1]; Tn[1]; Ti[1]; Te[1];
                                    tup[:,1]; tlower[:,1]; tdown[:,2];
                                    tlower[:,2]]...)


    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(num_layers-1)
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

    # NEW: Overwrite the entries for water in the lower atmosphere with 0s so that it will behave as fixed.
    # Only runs when water is in the activespecies list. If neutrals are set to inactive, it will be taken care of already.
    if in(:H2O, activespecies) && in(:HDO, activespecies)
        returnrates[H2Oi, 1:upper_lower_bdy_i] .= 0
        returnrates[HDOi, 1:upper_lower_bdy_i] .= 0
    end

    return [returnrates...;]
end

function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper, dt)
    #=
    TODO: docstring
    =#

    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    # tc___ are the coordinate tuples containing (I, J, V) to be used to fill a sparse matrix.
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

    for ialt in 2:(num_layers-1)
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
    append!(chemJi, tclocal[1].+(num_layers-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(num_layers-1)*length(activespecies))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(num_layers-1)*length(activespecies))
    append!(chemJj, tclower[2].+(num_layers-2)*length(activespecies))
    append!(chemJval, tclower[3])

    # NEW: fix water below whatever we set as upper/lower atmosphere boundary.
    # This only runs if water is designated as an active species; if it's in inactivespecies, this doesn't need to run,
    # and in fact, CAN'T run. When it is active, this finds all the H2O and HDO indices for the lower atmosphere. 
    # It's like above where we add (ialt-1)*length(activespecies), but this way it's outside the loop.
    if in(:H2O, activespecies) && in(:HDO, activespecies)
        H2Opositions = H2Oi .+ length(activespecies)*collect(0:upper_lower_bdy_i-1) # if H2Oi is 5, this returns 5, 69, 133...
        HDOpositions = HDOi .+ length(activespecies)*collect(0:upper_lower_bdy_i-1)
        water_positions = sort(union(H2Opositions, HDOpositions))

        i_remove = findall(x->in(x, water_positions), chemJi)
        j_remove = findall(x->in(x, water_positions), chemJj)
        remove_these = sort(union(i_remove, j_remove)) # This makes a set, since it describes the locations where the H2O and HDO indices are.
                                                       # Kinda confusing since we're talking about indices of indices.
        deleteat!(chemJi, remove_these)
        deleteat!(chemJj, remove_these)
        deleteat!(chemJval, remove_these)
    end

    # make sure to add 1's along the diagonal
    append!(chemJi, [1:length(nthis);])
    append!(chemJj, [1:length(nthis);])
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

function timeupdate(mytime, thefolder)
    #=
    Key function that runs the update! function and steps the simulation forward.

    mytime: a timestep size, logarithmic, in seconds
    thefolder: location in which to save the atmospheric state files and plots.
    =#

    for i = 1:numiters
        update!(n_current, mytime)   
        # plot_atm(n_current, [neutrallist, ionlist], thefolder*"/converged_atm_$(FNext).png", t=mytime, iter=i)
        if mytime == mindt
            println("iter $(i)")
        end
    end
    plot_atm(n_current, [neutrallist, ionlist], thefolder*"/converged_atm_$(FNext).png", t=mytime)  # only plot at the end of the iterations, speeds up code given the way I changed plot_atm.

    # Write out the current atmospheric state for making plots later 
    write_ncurrent(n_current, thefolder*"/ncurrent_$(mytime).h5") 
    # Uncomment if trying to examine detailed changes timestep to timestep.
    # plot_atm(n_current, [neutrallist, ionlist], thefolder*"/atm_plot_$(PN).png", t=mytime)

    # global PN += 1 # increment the number to append to plot filenames 
end

function next_timestep(nstart::Array{Float64, 1}, nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1}, activespecies, inactivespecies,
                       Jrates::Array{Float64, 2},
                       Tn::Array{Float64, 1}, Ti::Array{Float64, 1}, Te::Array{Float64, 1},
                       tup::Array{Float64, 2}, tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2}, tupper::Array{Float64, 2},
                       dt::Float64, abs_tol::Array{Float64, 1})
    #=
    moves to the next timestep using Newton's method on the linearized
    coupled transport and chemical reaction network.
    =#

    eps = 1.0 # ensure at least one iteration
    iter = 0

    abs_criterion = abs_tol # This criterion is 1 ppt, adjusted for density at different altitudes.
    rel_criterion = 1e-8

    # Set up the array to track whether all elements have met criterion
    met_criterion = BitArray(undef, length(nthis))

    while !all(met_criterion) # will quit when all elements of met_criterion are true.     # eps>1e-8 # old condition 
        nold = deepcopy(nthis)

        # stuff concentrations into update function and jacobian. This is the Eulerian (P*L) * dt
        fval = nthis - nstart - dt*ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te,  
                                          tup, tdown, tlower, tupper) 
        
        updatemat = chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, 
                             tup, tdown, tlower, tupper, dt)

        nthis = nthis - (updatemat \ fval)

        # Check whether error tolerance is met         
        abs_eps = abs.(nthis-nold) ./ abs_criterion  # calculated for every species.
        rel_eps = (abs.(nthis-nold)./nold) ./ rel_criterion
        # eps = maximum(abs.(nthis-nold)./nold) # old epsilon

        # Check each element to see if < 1. 
        abs_bool = abs_eps .< 1
        rel_bool = rel_eps .< 1

        # Assign values to the array that governs when the loop runs.
        met_criterion = abs_bool .| rel_bool

        # if all(met_criterion)
        #     println("All elements of nthis have met the tolerance after $(iter) iterations")
        # end

        iter += 1
        if iter>1e3
            throw("too many iterations in next_timestep! number of elements that met the criteria: $(count(met_criterion))/$(length(met_criterion))")
        end
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
    inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(non_bdy_layers)]...])
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(non_bdy_layers)])

    # extract concentrations
    nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(non_bdy_layers)]...])
    # Next line calculates 1 ppt of the total density at each altitude - used for absolute error tolerance.
    # This gets it into the same shape as nstart.
    ppt_vals = 1e-12 .* [[n_tot(n_current, a) for sp in activespecies, a in non_bdy_layers]...]

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
                          tup, tdown, tlower, tupper, dt, ppt_vals)
    nthismat = reshape(nthis, (length(activespecies), length(non_bdy_layers)))

    # write found values out to n_current
    for s in 1:length(activespecies)
        for ia in 1:num_layers
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
args = Any[ARGS[i] for i in 1:1:length(ARGS)]

# Establish a pattern for filenames. FNext = filename extension
# TODO: We probably won't be running all these experiments this time, so this section can probably go away.
# we may not even have to run from command line at all...
if args[1]=="temp" # TODO: this needs to handle the ion and electron temps.
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    FNext = "water_$(args[2])"
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
    DH = parse(Float64, args[2]) * 1.6e-4
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
else
    throw("Error! Bad experiment type")
end

user_input_folder_name = input("Enter a name for the results folder or press enter to use default: ")
sim_folder_name = user_input_folder_name == "" ? FNext : user_input_folder_name
# Set up the folder if it doesn't exist
create_folder(sim_folder_name, results_dir)

# Case where this file will be used to converge an atmosphere of a different extent
make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
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
    defaultatm = "4wdDae_notquiteconverged.h5"
    file_to_use = input("Enter the name of a file containing a converged, 250 km atmosphere to use (press enter to use default: $(defaultatm)): ")   # TODO: revert
    readfile = file_to_use == "" ?  defaultatm : file_to_use
    n_current = get_ncurrent(readfile)
else
    error("Didn't understand response")
end

# Whether to initialize new species as zeros or not
use_nonzero_initial_profiles = true

# Set whether neutrals or ions are being converged. The other group must be entered in nochemspecies, notransportspecies.
converge_which = input("Converging ions, neutrals or both?: ")
while converge_which != "neutrals" && converge_which != "ions" && converge_which != "both"
    println("Bad entry! Please enter ions or neturals or both.")
    global converge_which = input("Converging ions, neutrals or both?: ")
end
println()

# Set up initial profiles and time steps 
if converge_which == "neutrals"
    println("ALERT: Did you make sure that the ions are added to nochemspecies and notransportspecies?")

    for nn in new_neutrals
        n_current[nn] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        if length(new_neutrals) != 0
            println("Initializing non-zero profiles for $(new_neutrals)")
        end
        for nn in new_neutrals
            n_current[nn] = reshape(readdlm("../Resources/initial_profiles/$(string(nn))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
        end
    end
elseif converge_which == "ions"
    println("ALERT: Did you make sure that the neutrals are added to nochemspecies and notransportspecies?")

    for ni in new_ions
        n_current[ni] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        if length(new_ions) != 0
            println("Initializing non-zero profiles for $(new_ions)")
        end
        for ni in new_ions
            try
                n_current[ni] = reshape(readdlm("../Resources/initial_profiles/$(string(ni))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
            catch
                println("Skipping ion $(ni) for now because we need to set it after its analogue has been set...")
            end  
        end
        for Dion in keys(D_H_analogues)
            if in(Dion, new_ions)
                n_current[Dion] = DH .* n_current[D_H_analogues[Dion]]
                println("Set a profile for $(Dion): $(n_current[Dion][length(n_current[Dion])-10:end])...etc")
            end
        end
    end
elseif converge_which == "both"  # Currently this doesn't really work.
    println("ALERT: Did you make sure nochemspecies and notransportspecies are only :Ar and :N2?")
    for nn in new_neutrals
        n_current[nn] = zeros(num_layers)
    end
    for ni in new_ions
        n_current[ni] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        if length(new_ions) != 0 || length(new_neutrals) != 0
            println("Initializing non-zero profiles for $(new_neutrals) and $(new_ions)")
        end
        for nn in new_neutrals
            n_current[nn] = reshape(readdlm("../Resources/initial_profiles/$(string(nn))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
        end

        for ni in new_ions
            try
                n_current[ni] = reshape(readdlm("../Resources/initial_profiles/$(string(ni))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
            catch
                println("Skipping ion $(ni) for now because we need to set it after its analogue has been set...")
            end  
        end
        for Dion in keys(D_H_analogues)
            if in(Dion, new_ions)
                n_current[Dion] = DH .* n_current[D_H_analogues[Dion]]
                println("Set a profile for $(Dion): $(n_current[Dion][length(n_current[Dion])-10:end])...etc")
            end
        end
    end
else
    throw("Uncaught exception")
end

# Zero out the new Jrates
for nj in newJrates
    n_current[nj] = zeros(num_layers)
end

# Set up the timesteps
mindt = dt_min_and_max[converge_which][1]
maxdt = dt_min_and_max[converge_which][2]
total_timesteps = abs(mindt) + maxdt + 1 # the 1 is to account for 10^0. Used in some plotting codes. TODO: Check if not used anymore.

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
if cycle != "mean"
    println("ALERT: Solar $(cycle) data being used")
end

# Plot styles ==================================================================
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22

################################################################################
#                       TEMPERATURE/PRESSURE PROFILES                          #
################################################################################
println("Setting up temperature functions...")
# If changes to the temperature are needed, they should be made here 
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
plot_temp_prof([Temp_n(a) for a in alt], savepath=results_dir*sim_folder_name, i_temps=[Temp_i(a) for a in alt], e_temps=[Temp_e(a) for a in alt])


################################################################################
#                               WATER PROFILES                                 #
################################################################################
println("Setting up the water profile...")
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
else
    MR = MR_mean_water
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR # 10 pr μm
end

for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
end

# set the water profiles =======================================================
n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers)
n_current[:HDO] = 2 * DH * n_current[:H2O] # This should be the correct way to set the HDO profile.
# n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers) # OLD WAY that is wrong.


# We still have to calculate the HDO initial fraction in order to calculate the pr um 
# and make water plots.
HDOinitfrac = n_current[:HDO] ./ map(z->n_tot(n_current, z), non_bdy_layers)  

# ADD EXCESS WATER AS FOR DUST STORMS. TODO: revert
# H2Oppm = 1e-6*map(x->250 .* exp(-((x-42)/12.5)^2), non_bdy_layers/1e5) + H2Oinitfrac  # 250 ppm at 42 km (peak)
# HDOppm = 1e-6*map(x->0.350 .* exp(-((x-38)/12.5)^2), non_bdy_layers/1e5) + HDOinitfrac  # 350 ppb at 38 km (peak)
# n_current[:H2O] = H2Oppm .* map(z->n_tot(n_current, z), non_bdy_layers)
# n_current[:HDO] = HDOppm .* map(z->n_tot(n_current, z), non_bdy_layers)

# Compute total water column for logging and checking that we did things right =
# H2O #/cm^3 (whole atmosphere) = sum(MR * n_tot) for each alt
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
println("Setting boundary conditions and creating metaprogramming...")
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
                :H2=>["f" 0.; "v" H2_veff],
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
# DO NOT CHANGE!!
const active_above = [Symbol(string(s)*"_above") for s in activespecies]
const active_below = [Symbol(string(s)*"_below") for s in activespecies]

# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x), activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

# TODO: Turn these notes into a test that looks for "+()" in rates_local instead of debugging lines. 
# -----------------------------------------------------------------------------------------------
# These lines are useful for troubleshooting if you're getting weird errors with ratefn.
# "no method matching +()" means the chemical system is unbalanced and you've got a
# species that has production but no consumption, or vice versa.
# other errors may occur. Uncomment these 3 lines to inspect what goes into rates_local.
# println("The contents of rates_local: ")
# println(rates_local)
# println("That was rates_local")

arglist_local = [activespecies; active_above; active_below; inactivespecies;
                 Jratelist; :Tn; :Ti; :Te; local_transport_rates; :dt]

arglist_local_typed = [:($s::Float64) for s in arglist_local]

# These expressions which are evaluated below enable a more accurate assessment of M and E values
# by calculating at the time of being called rather than only at each timestep.
Mexpr = Expr(:call, :+, fullspecieslist...)
Eexpr = Expr(:call, :+, ionlist...)

# NOTE: These functions within @eval cannot be moved. Do not move them.
@eval begin
    function ratefn_local($(arglist_local_typed[1:end-1]...))

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
        abovechemJval = -dt*$(Expr(:vcat, chemJ_above[3]...))

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(Expr(:vcat, chemJ_below[3]...))

        # return the actual values of I, J, and V (indices and numerical value):
        return ((localchemJi, localchemJj, localchemJval),
                (abovechemJi, abovechemJj, abovechemJval),
                (belowchemJi, belowchemJj, belowchemJval))
    end
end


################################################################################
#                         PHOTOCHEMICAL CROSS SECTIONS                         #
################################################################################
println("Populating cross section dictionary...")

crosssection = populate_xsect_dict([T_surf, T_tropo, T_exo])#, o3xdata)

# Solar Input ==================================================================

const solarflux=readdlm(research_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]
solarflux[:,2] = solarflux[:,2]/2  # To roughly put everything at an SZA=60° (from a Kras comment)

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
println("Creating the simulation log file...")
# crosssection dict for logging purposes =======================================
xsect_dict = Dict("CO2"=>[co2file],
                  "H2O, HDO"=>[h2ofile, hdofile],
                  "H2O2, HDO2"=>[h2o2file, hdo2file],
                  "O3"=>[o3file, o3chapfile],
                  "O2"=>[o2file, o2_130_190, o2_190_280, o2_280_500],
                  "H2, HD"=>[h2file, hdfile],
                  "OH, OD"=>[ohfile, oho1dfile, odfile])

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
f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "w")
write(f, "$(args[1]) experiment: \n")
write(f, input_string)
write(f, "\n")

# Mean temperatures
write(f, "Mean temperatures used:\n")
write(f, "Surface: $(meanTs) K, Tropopause: $(meanTt) K, Exobase: $(meanTe) K\n\n")

# which species are turned on
write(f, "All species: $(join([string(i) for i in fullspecieslist], ", "))\n")
write(f, "No-chem species: $(join([string(i) for i in nochemspecies], ", "))\n")
write(f, "No-transport species: $(join([string(i) for i in notransportspecies], ", "))\n")
write(f, "Active species: $(join(sort([string(i) for i in activespecies]), ", "))\n\n")

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
write(f, "Initial atmosphere state: $(readfile)\n\n")

close(f)

################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# Uncomment this line if you'd like to add an extra parcel to some species. You must specify the species.
# n_current[:D] = map(x->1e5*exp(-((x-184)/20)^2), non_bdy_layers/1e5) + n_current[:D]

# write initial atmospheric state ==============================================
write_ncurrent(n_current, results_dir*sim_folder_name*"/initial_state.h5")

# Plot initial water profile ===================================================
plot_water_profile(H2Oinitfrac, HDOinitfrac, n_current[:H2O], n_current[:HDO], results_dir*sim_folder_name, watersat=H2Osatfrac)

# do the convergence ===========================================================
println("Plotting the initial condition")
plot_atm(n_current, [neutrallist, ionlist], results_dir*sim_folder_name*"/initial_atmosphere.png")

img_savepath = results_dir*sim_folder_name*"/converged_atm_"*FNext*".png"
println("Beginning Convergence")

@showprogress 0.1 "Converging..." [timeupdate(t, results_dir*sim_folder_name) for t in [10.0^(1.0*i) for i in mindt:maxdt]]

if converge_which == "neutrals" || converge_which == "both"
    @showprogress 0.1 "Last convergence steps..." for i in 1:100 
        plot_atm(n_current, [neutrallist, ionlist], img_savepath, t="1e14", iter=i)
        # println("dt: 1e14 iter $(i)")
        update!(n_current, 1e14)
    end
end

# write out the new converged file to matching folder.
towrite = results_dir*sim_folder_name*"/converged_"*FNext*".h5"
write_ncurrent(n_current, towrite)
println("Wrote $(towrite)")
println("Saved figure to same folder")

println("ALERT: Finished")
println()

f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "a")
t6 = time()
write(f, "Total runtime $((t6-t5) / 60) minutes\n")
write(f, "Module load time $(t2-t1) seconds\n")
write(f, "Parameter file load time $(t4-t3) seconds\n")
close(f)