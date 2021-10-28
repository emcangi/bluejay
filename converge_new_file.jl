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

#logging witchcraft
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Revise
using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Photochemistry  # custom module
using DifferentialEquations

t2 = time()

println("Time to load modules: $(format_sec_or_min(t2-t1))")

user_input_paramfile = input("Enter a parameter file or press enter to use default (PARAMETERS.jl): ")
paramfile = user_input_paramfile == "" ? "PARAMETERS.jl" : user_input_paramfile*".jl"
t3 = time()
include(paramfile)

# Some user-specified files have their own networks due to needing certain reactions to be on or off.
# For all others, load the default network.
load_rxn_net = ["PARAMETERS.jl", "PARAMETERS-conv3.jl", "PARAMETERS-convall.jl", "PARAMETERS-conv2a.jl", "PARAMETERS-conv4.jl"]
if any(x->occursin(x, paramfile), load_rxn_net)
    include("/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/reaction_network.jl")
end

t4 = time()
println("Time to load PARAMETERS and reaction network: $(format_sec_or_min(t4-t3))")

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

# Chemical Jacobian functions ==================================================
function chemJmat(n_active_longlived, n_active_shortlived, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper; dt=1., check_eigen=false)
    #=
    Collects coordinate tuples of (I, J, V) [row index, column index, value] for a sparse matrix
    representing the chemical jacobian of the atmospheric system. 

    n_active_longlived: The atmospheric densities array, but flattened, in the form [n_CO(z=0), n_CO2(z=0)...n_N2Dpl(z=0), n_CO(z=2)...n_N2Dpl(z=250)], for active 
                        and chemically long-lived species.
    n_active_shortlived: active shortlived species densities, necessary to to calculations for longlived species.
    n_inactive: A flattened array of the atmospheric densities of any inactive species, same format as nthis. Functionally constant.
    activellsp, activeslsp, inactivesp: Lists of active and longlived, active shortlived, and ianctive species (const)
    Jrates: Flattened array of Jrates, same format as nthis.
    Tn, Ti, Te: Temperature-vs-altitude arrays for neutrals, ions, electrons. (const)
    tup, tdown: Transport coefficients
    tlower, tupper: Transport coefficients

    check_eigen: Will check the eigenvalues of the jacobian for non-real or real/positive values if true. Code currently commented out.

    dt: Can be specified to manually control the timestep by which the derivatives are multiplied.
        If not supplied, then dt=1 so that the external solver can manage the multiplication by time.
    =#              

    nmat_llsp = reshape(n_active_longlived, (length(activellsp), num_layers))
    nmat_slsp = reshape(n_active_shortlived, (length(activeslsp), num_layers))
    nmat_inactive = reshape(n_inactive, (length(inactivesp), num_layers))
    
    # For storing the jacobian indices and values
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    # tc___ are the coordinate tuples containing (I, J, V) to be used to fill a sparse matrix.
    (tclocal, tcupper, tclower) = chemJmat_local([nmat_llsp[:, 1];  # active_longlived; 
                                                  nmat_llsp[:, 2]; # active_longlived_above; 
                                                  fill(1.0, length(activellsp)); # active_longlived_below;
                                                  nmat_slsp[:, 1]; # active_shortlived; 
                                                  nmat_inactive[:,1];  # inactive_species; 
                                                  Jrates[:,1]; # Jratelist; 
                                                  Tn[1]; Ti[1]; Te[1];  #:Tn; :Ti; :Te; 
                                                  tup[:,1]; tlower[:,1]; # local_transport_rates
                                                  tdown[:,2]; tlower[:,2]]...) 


    # add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])

    # and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activellsp))
    append!(chemJval, tcupper[3])

    for ialt in 2:(num_layers-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nmat_llsp[:, ialt];
                                                      nmat_llsp[:, ialt+1];
                                                      nmat_llsp[:, ialt-1];
                                                      nmat_slsp[:, ialt];
                                                      nmat_inactive[:, ialt];
                                                      Jrates[:, ialt]; 
                                                      Tn[ialt]; Ti[ialt]; Te[ialt];
                                                      tup[:, ialt];
                                                      tdown[:, ialt];
                                                      tdown[:, ialt+1];
                                                      tup[:, ialt-1]]...)


        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(activellsp))
        append!(chemJj, tclocal[2].+(ialt-1)*length(activellsp))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(activellsp))
        append!(chemJj, tcupper[2].+(ialt  )*length(activellsp))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(activellsp))
        append!(chemJj, tclower[2].+(ialt-2)*length(activellsp))
        append!(chemJval, tclower[3])
    end

    (tclocal, tcupper, tclower) = chemJmat_local([nmat_llsp[:,end];
                                                 fill(1.0, length(activellsp));
                                                 nmat_llsp[:,end-1];
                                                 nmat_slsp[:, end];
                                                 nmat_inactive[:,end];
                                                 Jrates[:,end];
                                                 Tn[end]; Ti[end]; Te[end]; 
                                                 tupper[:,1]; tdown[:,end];
                                                 tupper[:,2]; tup[:,end-1]]...)

    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(num_layers-1)*length(activellsp))
    append!(chemJj, tclocal[2].+(num_layers-1)*length(activellsp))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(num_layers-1)*length(activellsp))
    append!(chemJj, tclower[2].+(num_layers-2)*length(activellsp))
    append!(chemJval, tclower[3])

    # fix water below whatever we set as upper/lower atmosphere boundary.
    # This only runs if water is designated as an active species; if it's in inactive_species, this won't run,
    # When it is active, this finds all the H2O and HDO indices for the lower atmosphere. 
    # It's like above where we add (ialt-1)*length(active_species), but this way it's outside the loop.
    if in(:H2O, activellsp) && in(:HDO, activellsp)
        H2Opositions = H2Oi .+ length(activellsp)*collect(0:upper_lower_bdy_i-1)
        HDOpositions = HDOi .+ length(activellsp)*collect(0:upper_lower_bdy_i-1)
        water_positions = sort(union(H2Opositions, HDOpositions))

        i_remove = findall(x->in(x, water_positions), chemJi)
        j_remove = findall(x->in(x, water_positions), chemJj)
        remove_these = sort(union(i_remove, j_remove)) # This makes a set, since it describes the locations where the H2O and HDO indices are.
                                                       # Kinda confusing since we're talking about indices of indices.
        chemJval[remove_these] .= 0 
    end

    # J = sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +) 
    # println("checking eigenvalues")
    # if check_eigen==true
    #     check_jacobian_eigenvalues(J, results_dir*sim_folder_name)
    #     append!(stiffness, calculate_stiffness(J))
    # end

    return sparse(chemJi, chemJj, chemJval, length(n_active_longlived), length(n_active_longlived), +)
end

function ratefn(n_active_longlived, n_active_shortlived, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    nmat_llsp = reshape(n_active_longlived, (length(activellsp), num_layers))
    nmat_slsp = reshape(n_active_shortlived, (length(activeslsp), num_layers))
    nmat_inactive = reshape(n_inactive, (length(inactivesp), num_layers))
    returnrates = zeros(size(nmat_llsp))
    
    # fill the first altitude entry with information for all species
    returnrates[:,1] .= ratefn_local([nmat_llsp[:,1]; # densities for active_longlived;
                                      nmat_llsp[:,2]; # active_longlived_above; 
                                      fill(1.0, length(activellsp)); #active_longlived_below;
                                      nmat_slsp[:, 1]; # active_shortlived;
                                      nmat_inactive[:,1];  # inactive_species; 
                                      Jrates[:,1];  # Jratelist; 
                                      Tn[1]; Ti[1]; Te[1];  # :Tn; :Ti; :Te;
                                      tup[:,1]; tlower[:,1]; tdown[:,2]; tlower[:,2]]...) # local_transport_rates

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(num_layers-1)
        returnrates[:,ialt] .= ratefn_local([nmat_llsp[:, ialt]; # active_longlived;
                                            nmat_llsp[:, ialt+1];
                                            nmat_llsp[:, ialt-1];
                                            nmat_slsp[:, ialt]; # active_shortlived;
                                            nmat_inactive[:,ialt];
                                            Jrates[:,ialt];
                                            Tn[ialt]; Ti[ialt]; Te[ialt];
                                            tup[:, ialt]; 
                                            tdown[:, ialt];
                                            tdown[:, ialt+1]; 
                                            tup[:, ialt-1]]...)
    end

    # fill in the last level of altitude
    returnrates[:,end] .= ratefn_local([nmat_llsp[:, end];
                                       fill(1.0, length(activellsp));
                                       nmat_llsp[:, end-1];
                                       nmat_slsp[:, end]; # active_shortlived;
                                       nmat_inactive[:,end];
                                       Jrates[:,end];
                                       Tn[end]; Ti[end]; Te[end];
                                       tupper[:,1]; 
                                       tdown[:,end];
                                       tupper[:,2]; 
                                       tup[:,end-1]]...)


    # NEW: Overwrite the entries for water in the lower atmosphere with 0s so that it will behave as fixed.
    # Only runs when water is in the active_species list. If neutrals are set to inactive, it will be taken care of already.
    if in(:H2O, activellsp) && in(:HDO, activellsp)
        returnrates[H2Oi, 1:upper_lower_bdy_i] .= 0
        returnrates[HDOi, 1:upper_lower_bdy_i] .= 0
    end

    return [returnrates...;]
end

function get_transport_and_J_rates(n, n_short, n_inactive, activellsp, activeslsp, inactivesp, neutral_temps::Vector{Float64}, plasma_temps::Vector{Float64}, 
                                   D_arr::Vector{Float64}; compute_Jrates=true)
    #=
    This takes the current densities of all active species (n), transforms  back into an atmospheric state dictionary, calculates the current Jrates, 
    and then returns the transport coefficients, boundary layer transport coefficients, and Jrates needed to run ratefn and chemical_jacobian.

    n: densities by altitude for active, long-lived species 
    n_short: the same but for active short-lived species assumed to be in photochemical equilibrium.
    n_inactive: the same but for inactive species
    inactivesp, activellsp, activeslsp: species lists 
    neutral_temps, plasma_temps: temperature by altitude for neutrals and plasmas
    D_arr: An empty array to store diffusion coefficients, makes initializing the Dcoef dictionary faster.

    =#
    
    # transform n vectors back into dictionary so we can operate on it -------------------------------------------
    # Time for this block is ~0.0002 seconds.
    n_cur_active_long = unflatten_atm(n, activellsp)
    n_cur_active_short = unflatten_atm(n_short, activeslsp)
    n_cur_inactive = unflatten_atm(n_inactive, inactivesp)
    n_cur_all = Dict(vcat([k=>n_cur_active_long[k] for k in keys(n_cur_active_long)],
                          [k=>n_cur_active_short[k] for k in keys(n_cur_active_short)],
                          [k=>n_cur_inactive[k] for k in keys(n_cur_inactive)]))
    
    # Retrieve diffusion coefficients and mean scale height profiles for current atmospheric state ----------------
    # Keddy: 1D by altitude (array) independent of species
    # Molecular diffusion: species=>[D by altitude] (dictionary)
    # Mean atmospheric scale height (dictionary of 1D arrays)
    # Neither meal atmospheric scale height nor eddy diffusion are functions of species, but both depend on the total species density of the whole atmosphere.
    Keddy_arr = zeros(size(alt))
    Keddy_arr .= map(z->Keddy(z, n_tot(n_cur_all, z)), alt) # Eddy diffusion: K coefficient by altitude (array)
    H0_dict = Dict{String, Vector{Float64}}("neutral"=>map((z,t)->scaleH(z, t, n_cur_all), alt, neutral_temps),
                                            "ion"=>map((z,t)->scaleH(z, t, n_cur_all), alt, plasma_temps))

    # Molecular diffusion is only needed for transport species, though.
    Dcoef_dict = Dict{Symbol, Vector{Float64}}([s=>deepcopy(Dcoef!(D_arr, whichtemps[charge_type(s)], s, n_cur_all, speciesbclist)) for s in transport_species])

    # check for NaNs in these important inputs to fluxcoefs_all
    # report_NaNs(Keddy_arr, name="eddy diffusion coefficients")
    # report_NaNs(Dcoef_dict, name="diffusion coefficients")
    # report_NaNs(H0_dict, name="mean atmospheric scale height")
    # report_NaNs(Hs_dict, name="species-specific scale heights")

    # Calculate things which will be passed into ratefn() -------------------------------
    if compute_Jrates==true
        update_Jrates!(n_cur_all)
        # copy all the Jrates into an external dictionary for storage
        for jr in Jratelist                # time for this is ~0.000005 s
            global external_storage[jr] = n_cur_all[jr]
        end
    end

    # time to flatten Jrates is ~0.00014 s
    Jrates = deepcopy(Float64[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:length(non_bdy_layers)])
    
    # transport coefficients for bulk layers:
    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    # Temperature arrays and Hs dict are global scope. Rest are created within this function bc they depend on densities.
    fluxcoefs_all = fluxcoefs(neutral_temps, plasma_temps, Keddy_arr, Dcoef_dict, H0_dict, Hs_dict)
    tup = fill(-999., length(activellsp), num_layers)
    tdown = fill(-999., length(activellsp), num_layers)
    for (i, s) in enumerate(activellsp)
        tup[i, :] .= fluxcoefs_all[s][2:end-1, 2]
        tdown[i, :] .= fluxcoefs_all[s][2:end-1, 1]
    end
    
    # transport coefficients for boundary layers
    bc_dict = boundaryconditions(fluxcoefs_all, speciesbclist)
    tlower = permutedims(reduce(hcat, [bc_dict[sp][1,:] for sp in activellsp]))
    tupper = permutedims(reduce(hcat, [bc_dict[sp][2,:] for sp in activellsp]))
    
    return Jrates, tup, tdown, tlower, tupper
end

function make_jacobian(n, p, t)
    #=
    Constructs the chemical jacobian in the normal way, including stuff to calculate parameters for chemJmat.
    =#

    # Unpack the parameters ---------------------------------------------------------------
    n_inactive, inactivesp, activesp, activellsp, activeslsp, Tn::Vector{Float64}, Ti::Vector{Float64}, Te::Vector{Float64}, Tp::Vector{Float64}, D_arr::Vector{Float64} = p 

    # get the concentrations of species assumed to be in photochemical equilibrium. 
    n_shortlived = flatten_atm(external_storage, activeslsp)  # retrieve the shortlived species from their storage and flatten them

    # update Jrates and transport coefficients
    Jrates, tup, tdown, tlower, tupper = get_transport_and_J_rates(n, n_shortlived, inactive, activellsp, activeslsp, inactivesp, Tn, Tp, D_arr, compute_Jrates=true)

    # and update the shortlived species with the new Jrates 
    # n_short_updated = set_concentrations!(external_storage, n, n_shortlived, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)

    return chemJmat(n, n_short_updated, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper; dt=t)
end

function jacobian_wrapper(J, n, p, t)
    #=
    A wrapper since the Jacobian is sometimes generated outside of the normal simulation run
    =#

    J .= make_jacobian(n, p, t)
end
# Photochemical equilibrium functions ==========================================

function set_concentrations!(external_storage, n_active_long, n_active_short, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)
    #=
    at each altitude, sets the concentrations for short-lived species assumed to be in photochemical equilibrium
    and sends them back into the storage dictionary, external_storage

    external_storage: dictionary storing densities for short-lived and inactive species, as well as Jrates.
    n_active_short, n_active_long, n_inactive: density of short-lived, long-lived, and inactive species
    activeslsp, activellsp, inactivesp:: list of short- and long-lived species names
    Jrates: Jrates for each species, for a particular altitude
    Tn, Ti, Te: temperature arrays
    =#
    
    # rows = species, columns = altitudes. 
    nmat_shortlived = reshape(n_active_short, (length(activeslsp), num_layers))
    
    # auxiliary information that is needed. 
    nmat_longlived = reshape(n_active_long, (length(activellsp), num_layers))
    nmat_inactive = reshape(n_inactive, (length(inactivesp), num_layers))
    
    # storage for the updated concentrations
    new_densities = zeros(size(nmat_shortlived))
    # dist_zero = zeros(length(nmat_shortlived)) # TODO: Figure out how to make useful

    # fill the first altitude entry with information for all species   
    new_densities[:,1] .= set_concentrations_local([nmat_shortlived[:,1]; nmat_longlived[:, 1]; nmat_inactive[:, 1]; Jrates[:, 1]; Tn[1]; Ti[1]; Te[1]]...)
    # dist_zero[1] = check_zero_distance([nmat_shortlived[:,1]; nmat_longlived[:, 1]; nmat_inactive[:, 1]; Jrates[:, 1]; Tn[1]; Ti[1]; Te[1]]...)

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(num_layers-1)
        new_densities[:,ialt] .= set_concentrations_local([nmat_shortlived[:,ialt]; nmat_longlived[:, ialt]; nmat_inactive[:, ialt]; Jrates[:, ialt]; Tn[ialt]; Ti[ialt]; Te[ialt]]...)
        # dist_zero[ialt] = check_zero_distance([nmat_shortlived[:,ialt]; nmat_longlived[:, ialt]; nmat_inactive[:, ialt]; Jrates[:, ialt]; Tn[ialt]; Ti[ialt]; Te[ialt]]...)
    end

    # fill in the last level of altitude
    new_densities[:,end] .= set_concentrations_local([nmat_shortlived[:, end]; nmat_longlived[:, end]; nmat_inactive[:, end]; Jrates[:, end]; Tn[end]; Ti[end]; Te[end]]...)
    # dist_zero[end] = check_zero_distance([nmat_shortlived[:, end]; nmat_longlived[:, end]; nmat_inactive[:, end]; Jrates[:, end]; Tn[end]; Ti[end]; Te[end]]...)

    # Keep track of how many elements have experienced a change of > 1e-6 (1 ppm)
    change = abs.(vec(new_densities) .- n_active_short) ./ (n_active_short)
    ppm = 1e-6*ones(size(n_active_short)) # parts per thousand
    global num_beyond += count(x->x==true, change .> ppm)
    global num_times_called += 1
    
    # write out the new densities for shortlived species to the external storage
    for (s, ssp) in enumerate(activeslsp)
        external_storage[ssp] .= new_densities[s, :]
    end

    # if any(x->x>100, dist_zero)
    #     println("elements >100 from zero: $(dist_zero[findall(x->x>100, dist_zero)])")
    # end

    return vec(new_densities)
end

# Production and loss equation functions =======================================
function PnL_eqn(dndt, n, p, t)
    #=
    This is the primary function that is solved by the solver. For our purposes that means it calls ratefn().

    dndt is written as such because the output of this function is silently multiplied by the timestep within the solver's code.
    n is a vector of active species densities by altitude.
    =#

    # Unpack the parameters needed to call ratefn 
    n_inactive, inactivesp, activesp, activellsp, activeslsp, Tn::Vector{Float64}, Ti::Vector{Float64}, Te::Vector{Float64}, Tp::Vector{Float64}, D_arr::Vector{Float64} = p 

    if t / timestorage >= 10
        f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "a")
        progress_alert = "reached timestep $(t) in $(format_sec_or_min(time() - ti))\n"
        write(f, progress_alert)
        close(f)
        println(progress_alert) # prints present timestep, but only if it's a factor of 10 greater than the last stored timestep (so we don't get a trillion outputs)
        global timestorage = t # updates the last stored timestep to repeat process
        
        # Write out atmospheric state to a file to keep track as things run.
        ncur = merge(external_storage, unflatten_atm(n, activellsp))
        plot_atm(ncur, [neutral_species, ion_species], results_dir*sim_folder_name*"/atm_peek_dt=$(t).png", t="$(round(t, digits=4))")
        write_ncurrent(ncur, results_dir*sim_folder_name*"/atm_dt=$(t).h5")
    end

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, activeslsp)  

    # Use the current state to retrieve the appropriate J rates and transport coefficients for the long-lived species.
    Jrates, tup, tdown, tlower, tupper = get_transport_and_J_rates(n, n_short, n_inactive, activellsp, activeslsp, inactivesp, Tn, Tp, D_arr, compute_Jrates=false)
        
    # set the concentrations of species assumed to be in photochemical equilibrium. Runs multiple times to ensure it gets closer to the true solution.
    n_short_updated = set_concentrations!(external_storage, n, n_short, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)
    for i in 1:9
        n_short_updated = set_concentrations!(external_storage, n, n_short_updated, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)
    end 

    dndt .= ratefn(n, n_short_updated, inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
end

function update_Jrates!(n_cur_densities::Dict{Symbol, Array{Float64, 1}})
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#

    # Initialize an array, length=number of active layers
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
            jcolumn += n_cur_densities[species][ialt]*dz

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
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1, solarabs[ialt], 1)
        end
    end

    # solarabs now records the total optical depth of the atmosphere at
    # each wavelength and altitude

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2] .* exp.(-solarabs[ialt])
    end

    # each species absorbs according to its cross section at each
    # altitude times the actinic flux.
    # BLAS.dot includes an integration (sum) across wavelengths, i.e:
    # (a·b) = aa + ab + ab + bb etc that kind of thing
    for j in Jratelist
        n_cur_densities[j] = zeros(size(n_cur_densities[:CO2]))
        for ialt in [1:nalt;]
            n_cur_densities[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1, crosssection[j][ialt+1], 1)
        end
    end
end

# The equation that actually calls the ODE solver ==============================
function evolve_atmosphere(atm_init::Dict{Symbol, Array{Float64, 1}}, log_t_start, log_t_end; t_to_save=[])
    #=
    Sets up the initial conditions for the simulation and calls the ODE solver. 

    atm_init: dictionary of species densities by altitude for the current sttae.

    log_t_start: Starting time will be 10^log_t_start.
    log_t_end: Similar, but end time.
    t_to_save: timesteps at which to save a snapshot of the atmosphere.
    =#

    # Set up the initial state and check for any problems 
    nstart = flatten_atm(atm_init, active_longlived)
    find_nonfinites(nstart, collec_name="nstart")
    # Next line calculates 1 ppt of the total density at each altitude for absolute error tolerance,
    # vector in same shape as nstart.
    abs_tol_vec = 1e-12 .* [[n_tot(atm_init, a) for sp in active_longlived, a in non_bdy_layers]...] 


    # Set up simulation solver parameters and log them
    odesolver = "QNDF"
    saveall = false
    abst = "1 ppt for active long-lived species"
    relt = rel_tol
    sys_size = length(nstart)
    startdt = 10.0^log_t_start
    tspan = (10.0^log_t_start, 10.0^log_t_end)
    
    Dcoef_arr_template = zeros(size(Tn_arr))  # For making diffusion coefficient calculation go faster 
    params = [inactive, inactive_species, active_species, active_longlived, active_shortlived, Tn_arr::Vector{Float64}, Ti_arr::Vector{Float64}, Te_arr::Vector{Float64}, 
              Tplasma_arr::Vector{Float64}, Dcoef_arr_template::Vector{Float64}]
    params_exjac = [inactive, inactive_species, active_species, active_longlived, active_shortlived, Tn_arr::Vector{Float64}, Ti_arr::Vector{Float64}, Te_arr::Vector{Float64}, 
                    Tplasma_arr::Vector{Float64}, Dcoef_arr_template::Vector{Float64}]

    # Set up an example Jacobian and its preconditioner
    example_jacobian = make_jacobian(nstart, params_exjac, tspan[1])
    find_nonfinites(example_jacobian, collec_name="example_jacobian")
    println("Finished the example jacobian")

    # Log solver options
    f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "a")
    write(f, "SOLVER OPTIONS: \nsystem size: $(sys_size)\njacobian total elements: $(sys_size^2)\n")
    write(f, "timespan=$(tspan)\nsolver=$(odesolver)\nsaveat=$(t_to_save)\nsave_everystep=$(saveall)\nabstol=$(abst)\nreltol=$(relt)\nstarting dt=$(startdt)\n\n")
    
    sparsity = round(length(example_jacobian.nzval)*100/(size(example_jacobian)[1] * size(example_jacobian)[2]), digits=2)
    if sparsity > 1
        println("Warning! Sparsity of the jacobian is rather high: $(sparsity)%")
    end
    write(f, "Sparsity: $(sparsity)\n\n")
    close(f)

    # Now define the problem and solve it
    f = ODEFunction(PnL_eqn, jac=jacobian_wrapper, jac_prototype=example_jacobian)
    prob = ODEProblem(f, nstart, tspan, params)
    
    println("Starting the solver...")
    sol = solve(prob, QNDF(), saveat=t_to_save, progress=true, save_everystep=saveall, dt=startdt, #progress_steps=1, 
                abstol=abs_tol_vec, reltol=relt, isoutofdomain=(u,p,t)->any(x->x<0,u))
end

################################################################################
#                                  MAIN SETUP                                  #
################################################################################

# Set up simulation files and experiment type ==================================

# get command line arguments for experiment type. format:
# <experiment type> <parameters> 
# examples: 
# temp Tsurf Ttropo Texo; water <mixing ratio>; dh <multiplier>; 
# Oflux <cm^-2s^-1>
# examples: temp 190 110 200; water 1e-3; dh 8; Oflux 1.2e8.
args = Any[ARGS[i] for i in 1:1:length(ARGS)]

# Establish a pattern for filenames. FNext = filename extension
# TODO: We probably won't be running all these experiments this time, so this section can probably go away.
# we may not even have to run from command line at all...
if args[1]=="temp"
    const FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    const FNext = "water_$(args[2])"
elseif args[1]=="dh"
    const FNext = "dh_$(args[2])"
    const DH = parse(Float64, args[2]) * 1.6e-4
elseif args[1]=="Oflux"
    const FNext = "Oflux_$(args[2])"
else
    throw("Error! Bad experiment type")
end

# Create the results folder
sim_folder_name = sim_folder_name == "" ? FNext : sim_folder_name  # allows to not specify a name and use the default
create_folder(sim_folder_name, results_dir)

# Establish the final write-out file
final_atm_file = final_atm_file == "" ? "final_atmosphere_$(FNext).h5" : final_atm_file 

# Convert the arguments to numbers so we can use them to do maths
for i in 2:1:length(args)
    args[i] = parse(Float64, args[i])
end

# Code used to change the vertical extent (altitudes) ================================================
if make_new_alt_grid==true # TODO: rework this because it is deprecated.
    throw("The code for extending the altitude grid needs to be redone.")
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
end

# Load the starting atmosphere ======================================================================== 
const n_current = get_ncurrent(initial_atm_file)

# Setup of new species profiles and density storage ===================================================

if converge_which == "neutrals"
    println("Converging neutrals only. The following readout should contain the ions and N-bearing neutrals: $(inactive_species)")

    for nn in new_neutrals
        n_current[nn] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        println("Initializing non-zero profiles for $(new_neutrals)")
        for nn in new_neutrals
            n_current[nn] = reshape(readdlm("../Resources/initial_profiles/$(string(nn))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
        end
    end
elseif converge_which == "ions"
    println("Converging ions. The following readout should contain the non-N-bearing neutrals: $(inactive_species)")

    for ni in new_ions
        n_current[ni] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        println("Initializing non-zero profiles for $(new_ions)")
        # first fill in the H-bearing ions from data-inspired profiles
        for ni in setdiff(new_ions, keys(D_H_analogues))
            n_current[ni] = reshape(readdlm("../Resources/initial_profiles/$(string(ni))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
        end
        # Then create profiles for the D-bearing analogues based on the H-bearing species profiles
        for ni in intersect(new_ions, keys(D_H_analogues))
            n_current[ni] = DH .* n_current[ni]
        end
    end
elseif converge_which == "both" 
    if occursin("PARAMETERS-conv3", paramfile)
        println("Converging N-bearing neutrals and ions together. This list readout of inactive_species should contain non-N-bearing neutrals: $(inactive_species)")
    else
        println("Converging neutrals and ions together. This list readout of inactive_species should be empty: $(inactive_species)")
    end

    for nn in new_neutrals
        n_current[nn] = zeros(num_layers)
    end
    for ni in new_ions
        n_current[ni] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        println("Initializing non-zero profiles for $(new_neutrals) and $(new_ions)")
        for nn in new_neutrals
            try
                n_current[nn] = reshape(readdlm("../Resources/initial_profiles/$(string(nn))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
            catch 
                println("No initial guess found for $(nn). Initial profile will be zero everywhere.")
            end
        end

        for ni in setdiff(new_ions, keys(D_H_analogues))
            try
                n_current[ni] = reshape(readdlm("../Resources/initial_profiles/$(string(ni))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
            catch 
                println("No initial guess found for $(ni). Initial profile will be zero everywhere.")
            end
        end
        
        for ni in intersect(new_ions, keys(D_H_analogues))
            n_current[ni] = DH .* n_current[ni]
        end
    end
else
    throw("Uncaught exception")
end

# Zero out the new Jrates
for nj in newJrates
    n_current[nj] = zeros(num_layers)
end

# Create a separate dictionary for: short lived species, inactive species, and Jrates which will be updated separately from
# the main solving routine because these values cannot be/need not be part of the solver algorithm. 
const external_storage = Dict([j=>n_current[j] for j in union(short_lived_species, inactive_species, Jratelist)])

# Densities of inactive species, which by definition do not change for entire simulation
const inactive = flatten_atm(n_current, inactive_species)

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
    global const T_surf = args[2]
    global const T_tropo = args[3]
    global const T_exo = args[4]
    Temp_n(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "neutral")
    Temp_i(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "ion")
    Temp_e(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "electron")
    Temp_keepSVP(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral") # for testing temp without changing SVP. #TODO: adjust if needed for ions?
else
    global const T_surf = meanTs
    global const T_tropo = meanTt
    global const T_exo = meanTe
    Temp_n(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral")
    Temp_i(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "ion")
    Temp_e(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "electron")
end

const controltemps = [T_surf, T_tropo, T_exo]

# set temperature arrays to be used everywhere
const Tn_arr = [Temp_n(a) for a in alt];
const Ti_arr = [Temp_i(a) for a in alt];
const Te_arr = [Temp_e(a) for a in alt];
const Tplasma_arr = 0.5 * (Ti_arr .+ Te_arr);
const whichtemps = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)

# Species-specific scale heights - has to be done here once the control temps are set
Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>map(z->scaleH(z, sp, controltemps), alt) for sp in all_species])

plot_temp_prof([Temp_n(a) for a in alt], savepath=results_dir*sim_folder_name, i_temps=[Temp_i(a) for a in alt], e_temps=[Temp_e(a) for a in alt])

################################################################################
#                               WATER PROFILES                                 #
################################################################################

println("Setting up the water profile...")

# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    const H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    const HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
else
    fix_SVP = false
    const H2Osat = map(x->Psat(x), map(Temp_n, alt)) # array in #/cm^3 by altitude
    const HDOsat = map(x->Psat_HDO(x), map(Temp_n, alt))
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

# ADD EXCESS WATER AS FOR DUST STORMS.
dust_storm_on = "no"#input("Add excess water to simulate a dust storm? (y/n)")
if dust_storm_on == "y" || dust_storm_on == "yes"
    H2Oppm = 1e-6*map(x->250 .* exp(-((x-42)/12.5)^2), non_bdy_layers/1e5) + H2Oinitfrac  # 250 ppm at 42 km (peak)
    HDOppm = 1e-6*map(x->0.350 .* exp(-((x-38)/12.5)^2), non_bdy_layers/1e5) + HDOinitfrac  # 350 ppb at 38 km (peak)
    n_current[:H2O] = H2Oppm .* map(z->n_tot(n_current, z), non_bdy_layers)
    n_current[:HDO] = HDOppm .* map(z->n_tot(n_current, z), non_bdy_layers)
end

# We still have to calculate the HDO initial fraction in order to calculate the pr um 
# and make water plots.
HDOinitfrac = n_current[:HDO] ./ map(z->n_tot(n_current, z), non_bdy_layers)  

# Compute total water column for logging and checking that we did things right.
# H2O #/cm^3 (whole atmosphere) = sum(MR * n_tot) for each alt. Then convert to 
# pr μm: (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23, 
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

################################################################################
#                             BOUNDARY CONDITIONS                              #
################################################################################
println("Setting boundary conditions and creating metaprogramming...")
const H_veff = effusion_velocity(Temp_n(zmax), 1.0, zmax)
const H2_veff = effusion_velocity(Temp_n(zmax), 2.0, zmax)
const D_veff = effusion_velocity(Temp_n(zmax), 2.0, zmax)
const HD_veff = effusion_velocity(Temp_n(zmax), 3.0, zmax)

#=
    boundary conditions for each species (mostly from Nair 1994, Yung 1988). For 
    most species, default boundary condition is zero (f)lux at top and bottom. 
    Atomic/molecular hydrogen and deuterated analogues have a nonzero effusion 
    (v)elocity at the upper layer of the atmosphere. Some species have a (n)umber 
    density boundary condition.
=#

if args[1] == "Oflux"
    const Of = args[2]
else
    const Of = 1.2e8
end

# This has to be defined here, because it uses as a boundary condition the H2O 
# and HDO saturation at the surface.
const speciesbclist=Dict(
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
    to and from the cells above and below: (CONTENT MOVED TO PARAMETERS FILE)
=#

# Create symbolic expressions for the chemical jacobian at a local layer with influence from that same layer, 
# the one above, and the one below
const chemJ_local = chemical_jacobian(reactionnet, transportnet, active_longlived, active_longlived);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, active_longlived, active_longlived_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, active_longlived, active_longlived_below);

# Active long-lived and short-lived species expression arrays =====================

# Long-lived species expression array ---------------------------------------------
# An array to store the rate equations for active, long-lived species, which are 
# solved for in the production and loss equation.
# each row is for each species; each column is for chemical production, chemical loss, 
# transport production, transport loss, in that order.
active_longlived_species_rates = Array{Array{Expr}}(undef, length(active_longlived), 4)
for (i, sp) in enumerate(active_longlived)
    active_longlived_species_rates[i, :] .= getrate(reactionnet, transportnet, sp, sepvecs=true)
end

# Short-lived species expression array --------------------------------------------
# Similarly, this array stores expressions for the concentrations of
# active, short-lived species, which are assumed to be in photochemical equilibrium. 
# Instead of rates, the expressions calculate the new density, n_s.
#
# The stored value is the equation P - L = 0 for photochemical equilibrium, where each 
# row corresponds to the equation for a different species. 
# If the loss term is linear in n_s, the solution is n_s = P/L.
# If it's quadratic in n_s, the solution is a quadratic equation for which two solutions
# are possible. 
# Column 1: Solution 1 for n_s
# Column 2: Either 0 if loss term is linear in species density, or Solution 2 if quadratic.
short_lived_density_eqn = Array{Expr}(undef, length(short_lived_species), 2)

# This array just stores the same P - L = 0, but in the format P - nL = 0 for linear 
# and an^2 + bn + c = 0 for quadratic. This allows us to check how good a job the 
# densities solved for do in actually getting the system in equilibrium.
# The reason it's possible for P - nL =/= 0 or an^2 + bn + c =/= 0 is that 
# we solve for the density of each species one-by-one rather than solving it
# as a vector system. This may eventually need to be changed. It may also not be
# possible to solve it as a vector system and satisfy the constraints... which would 
# just be the entire original problem that we were trying to avoid by assuming photochemical
# equilibrium. 
equilibrium_eqn_terms = Array{Expr}(undef, length(short_lived_species), 1)

for (i, sp) in enumerate(active_shortlived)
    # Removes from consideration any reaction where a species appears on both sides of the equation (as an observer)
    ret = rxns_where_species_is_observer(sp, reactionnet)
    if ret == nothing
        chemnet = reactionnet
    else
        chemnet = filter(x->!in(x, ret), reactionnet)        
    end
    
    # Get the species production rate and loss rate by chemistry. These are obtained as vectors of reaction vectors
    # in the form [[:R1, :R2], [:P1, :P2], :(rate)]
    chem_prod_rate = production_rate(chemnet, sp, return_peqn_unmapped=true)
    chem_loss_rate = loss_rate(chemnet, sp, return_leqn_unmapped=true)
    
    if linear_in_species_density(sp, chem_loss_rate)
        # factors out the n_s so we can do n_s = P/L and converts to a big expression
        Lcoef_val = make_net_change_expr(loss_coef(chem_loss_rate, sp)) 
        P_val = make_net_change_expr(chem_prod_rate) # convert production to a big expression
        short_lived_density_eqn[i, 1] = :($P_val / $Lcoef_val)
        short_lived_density_eqn[i, 2] = :(0+0) # no second solution for linear.

        equilibrium_eqn_terms[i, 1] = :($P_val - $sp*($Lcoef_val))
    else # if it's quadratic in the species in question, i.e. the species appears twice on the LHS
        # Get the quadratic coefficients A, B, C for P - L = A(n_s^2) + B(n_s) + C = 0
        println("Note: $(sp) is not linear in density")
        qc = construct_quadratic(sp, chem_prod_rate, chem_loss_rate)

        # store the two solutions in the array. Warning: HORRIFYING!
        short_lived_density_eqn[i, 1] = :((-$(qc["B"]) + sqrt($(qc["B"])^2 - 4*$(qc["A"])*$(qc["C"])))/(2*$(qc["A"])) )
        short_lived_density_eqn[i, 2] = :((-$(qc["B"]) - sqrt($(qc["B"])^2 - 4*$(qc["A"])*$(qc["C"])))/(2*$(qc["A"])) )
        
        # Populate the array that lets us check if the densities give us 0
        equilibrium_eqn_terms[i, 1] = :($(qc["A"])*(($sp)^2) + $(qc["B"])*($sp) + $(qc["C"]))
    end
end

# Arguments for ratefn_local
const ratefn_arglist = [active_longlived; active_longlived_above; active_longlived_below; active_shortlived; inactive_species; Jratelist; 
                        :Tn; :Ti; :Te; local_transport_rates];
const ratefn_arglist_typed = [:($s::Float64) for s in ratefn_arglist];

# Arguments for set_concentrations (photochemical equilibrium)
const set_concentration_arglist = [active_shortlived; active_longlived; inactive_species; Jratelist; :Tn; :Ti; :Te];
const set_concentration_arglist_typed = [:($s::Float64) for s in set_concentration_arglist];

# Expressions for more accurate values of M (total third-bodies) and E (total electrons)
const Mexpr = Expr(:call, :+, all_species...)
const Eexpr = Expr(:call, :+, ion_species...)

# NOTE: These functions within @eval cannot be moved. Do not move them.
@eval begin
    function ratefn_local($(ratefn_arglist_typed...))
        #=
        How this works is that every time ratefn_local is called, this block gets evaluated
        and the contents of active_longlived_arglist_typed are evaluated numerically. All this does 
        is provide numerical values for all the symbols such as :CO2, :CO2pl, :H_above, and 
        etc. Then the symbolic expressions for production and loss of each species, which are 
        contained in active_longlived_species_rates, are evaluated.
        
        Note that for a given species, this function only has information for a given atmospheric
        layer, the layer above, and the layer below. So this function is called as part of a loop 
        which is written out in ratefn(). 
        
        The return value is the net change (dndt) to the species concentrations. The order of the array
        is the same as the order of arguments, so in that sense, "which altitude we are calculating for"
        is encoded in the order of terms in the return array. 
        =#

        # M and E are calculated here to ensure that the right number of ions/electrons
        # is used. It is for only the altitude at which this function was called 
        # (i.e. all the arguments to the function, when it's called, are for only one altitude)
        M = $Mexpr
        E = $Eexpr

        # stack overflow - answer (Cite)
        # create a result array for evaluating the production and loss expressions
        result = map(_ -> Float64[], $active_longlived_species_rates) 

        $( (quote
            # i = row, j = vector, k = specific expression
                push!(result[$i], $(expr))
            end 
            for (i, expr_list) in enumerate(active_longlived_species_rates)                                                                                                                 
                for expr in expr_list
            )...  
         )

        # these track the changes for each species
        net_chem_change = zeros(size(result, 1))
        net_trans_change = zeros(size(result, 1))
        
        for r in 1:size(result,1)
            net_chem_change[r] = subtract_difflength(sort(result[r, :][1], rev=true), sort(result[r, :][2], rev=true))
            net_trans_change[r] = subtract_difflength(sort(result[r, :][3], rev=true), sort(result[r, :][4], rev=true))
        end

        return net_chem_change .+ net_trans_change

    end
end

@eval begin
    function check_zero_distance($(set_concentration_arglist_typed...))

        M = $Mexpr
        E = $Eexpr
        
        # Make a result array in which to store evaluated expressions 
        result = map(_ -> 0., $equilibrium_eqn_terms)

        # then we evaluate the expressions in the unrolled loop as below 
        $( (quote
                result[$i] = $(expr)
            end 
            for (i, expr) in enumerate(equilibrium_eqn_terms)                                                                                                                 
            )...  
         )
        
        # Now return the distance from zero
        return sqrt(sum(result .^ 2))
    end
end

@eval begin
    function set_concentrations_local($(set_concentration_arglist_typed...))
        #=
        Calculates the best possible solution for the new density of each short-lived species.
        =#

        M = $Mexpr
        E = $Eexpr
        
        # Make a result array in which to store evaluated expressions 
        result = map(_ -> 0., $short_lived_density_eqn)  # will hold evaluated expressions.

        # then we evaluate the expressions in the unrolled loop as below and sum up each so we have total chem prod and total chem loss.
        # (and also transport)
        $( (quote
                result[$i] = $(expr)
            end 
            for (i, expr) in enumerate(short_lived_density_eqn)                                                                                                                 
            )...  
         )
        
        # get the previous densities before we evaluated the expressions. They already exist as the first set of arguments to this
        # function, but they're hard to pull out specifically. This makes it easier. Same logic as above.
        prev_densities = zeros(size(active_shortlived))

        $( (quote
                prev_densities[$j] = $(ex)
            end 
            for (j, ex) in enumerate(active_shortlived)                                                                                                                 
            )...  
         )

        # Handle negatives and quadratic cases and pick out one best density solution for each species
        best_solutions = choose_solutions(result, prev_densities)

        return best_solutions
    end
end

@eval begin
    function chemJmat_local($(ratefn_arglist_typed...))
        #=
        Generates a matrix of I, J, and V values for a sparse matrix, for the 
        local layer, for influences from the layer above, and influences 
        from the layer below. 
        =#

        M = $Mexpr
        E = $Eexpr

        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = $(Expr(:vcat, chemJ_local[3]...)) 

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = $(Expr(:vcat, chemJ_above[3]...))

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = $(Expr(:vcat, chemJ_below[3]...))

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

const crosssection = populate_xsect_dict(controltemps)

# Solar Input ==================================================================

solarflux=readdlm(research_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]
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
    input_string = "T_surf=$(args[2]), T_tropo=$(args[3]), T_exo=$(args[4])" *
                   "\nwater init=$(MR)\nDH=5.5 \nOflux=1.2e8"
elseif args[1]=="water"
    input_string = "T_surf=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\n" *
                   "water init=$(args[2])\nDH=5.5\nOflux=1.2e8\n"
elseif args[1]=="dh"
    input_string = "T_surf=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)\n" *
                   "DH=$(args[2]) \nOflux=1.2e8\n"
elseif args[1]=="Oflux"
    input_string = "T_surf=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)" *
                   "\nDH=5.5\nOflux=$(Of)\n"
end

# Write the log ================================================================
f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "w")

# Basic parameters written out 
write(f, "$(args[1]) experiment: \n")
write(f, input_string*"\n\n")
write(f, "Initial atmosphere state: $(initial_atm_file)\n\n")
write(f, "Final atmosphere state: $(final_atm_file)\n\n")
write(f, "Mean temperatures used:\n")
write(f, "Surface: $(meanTs) K, Tropopause: $(meanTt) K, Exobase: $(meanTe) K\n\n")
write(f, "SVP fixed: $(fix_SVP)\n\n")
write(f, "Solar cycle status: $(solarfile)\n\n")
write(f, "Non-zero initial profiles used: $(use_nonzero_initial_profiles)\n\n")

# which species are turned on
write(f, "All species: $(join([string(i) for i in all_species], ", "))\n\n")
write(f, "No-chem species: $(join([string(i) for i in no_chem_species], ", "))\n\n")
write(f, "Chem species: $(join([string(i) for i in chem_species], ", "))\n\n")
write(f, "No-transport species: $(join([string(i) for i in no_transport_species], ", "))\n\n")
write(f, "Transport species: $(join([string(i) for i in transport_species], ", "))\n\n")
write(f, "Active long-lived species: $(join(sort([string(i) for i in active_longlived]), ", "))\n\n")
write(f, "Active short-lived species: $(join(sort([string(i) for i in active_shortlived]), ", "))\n\n")

# Water profile information
write(f, "Water profile information: \n")
write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
write(f, "Total H2O (pr μm): $(H2Oprum)\n")
write(f, "Total HDO (pr μm): $(HDOprum)\n")
write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)\n\n")

# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="),
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end
write(f, "\n")

# Which photodissociation reactions are turned on 
# write(f, "New J rates that are turned on: \n")
# for rxn in [string(i) for i in newJrates]
#     write(f, rxn*"\n")
# end
# write(f, "\n")

# gets just chemistry reactions with ions as reactants or products
# TODO: Write out the name of the reaction network spreadsheet, once network is stored in spreadsheets.
# ion_chem_rxns = filter(x->(occursin("pl", string(x[1])) || occursin("pl", string(x[2]))), filter(x->!occursin("J", string(x[3])), reactionnet)) 
# write(f, "Ion chemistry reactions: \n")
# for rxn in [string(i) for i in ion_chem_rxns]
#     write(f, rxn*"\n")
# end
# write(f, "\n")

# cross sections
write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
write(f, "\n")

################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# Uncomment this line if you'd like to add an extra parcel to some species. You must specify the species.
# n_current[:D] = map(x->1e5*exp(-((x-184)/20)^2), non_bdy_layers/1e5) + n_current[:D]

# write initial atmospheric state ==============================================
write_ncurrent(n_current, results_dir*sim_folder_name*"/initial_state.h5")

# Plot initial water profile ===================================================
plot_water_profile(H2Oinitfrac, HDOinitfrac, n_current[:H2O], n_current[:HDO], results_dir*sim_folder_name, watersat=H2Osatfrac)

# Plot initial atmosphere condition  ===========================================
println("Plotting the initial condition")
plot_atm(n_current, [neutral_species, ion_species], results_dir*sim_folder_name*"/initial_atmosphere.png")

# Create a list to keep track of stiffness ratio ===============================
# const stiffness = []

# Record setup time
t6 = time()

write(f, "Setup time $(format_sec_or_min(t6-t5))\n")
close(f)

# do the convergence ===========================================================
println("Beginning Convergence")

# Timesteps
const mindt = dt_min_and_max[converge_which][1]
const maxdt = dt_min_and_max[converge_which][2]
times_to_save = [10.0^t for t in mindt:maxdt]
timestorage = times_to_save[1] / 10  # This variable lets us periodically print the timestep being worked on for tracking purposes

# Keep track of number of elements greater than 1 ppm change when calculating photochemical equilibrium
num_beyond = 0
num_times_called = 0

# Pack the global variables and send them through for faster code!
ti = time()
atm_soln = evolve_atmosphere(n_current, mindt, maxdt, t_to_save=times_to_save) # dz, transport_species, alt, all_species, num_layers, Jratelist, non_bdy_layers, 
tf = time() 

f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "a")
write(f, "Simulation (active convergence) runtime $(format_sec_or_min(tf-ti))\n")

println("Finished convergence in $((tf-ti)/60) minutes")


# Save the results =============================================================

L = length(atm_soln.u)
i = 1
# Write all states to individual files.
for (timestep, atm_state) in zip(atm_soln.t, atm_soln.u)
    # This is the contents of unflatten_atm.

    nc_all = merge(external_storage, unflatten_atm(atm_state, active_longlived))
    
    if i == L
        for jsp in Jratelist  # TODO: This needs to be worked so Jrates are tracked for each timestep. right now, only the last one is kept. :/
            nc_all[jsp] = external_storage[jsp]
        end
        plot_atm(nc_all, [neutral_species, ion_species], results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state")
        write_ncurrent(nc_all, results_dir*sim_folder_name*"/"*final_atm_file)   # write out the final state to a specially named file
    end
    global i += 1

    filepath = results_dir*sim_folder_name*"/atm_state_t_$(timestep).h5"
    write_ncurrent(nc_all, filepath)   
end

t6 = time()

println("Wrote all states to files")

t9 = time()

# write(f, "Stiffness ratio evolution:\n")
# write(f, "$(stiffness)")
# write(f, "\n\n")
# write(f, "Total runtime $(format_sec_or_min(t9-t1))\n")

close(f)



println("Total runtime $(format_sec_or_min(t9-t1))")