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
# using LSODA
using ProgressMeter
# using IncompleteLU
using Photochemistry  # custom module


using DifferentialEquations
using Sundials

t2 = time()

println("Time to load modules: $(format_sec_or_min(t2-t1))")

user_input_paramfile = input("Enter a parameter file or press enter to use default (PARAMETERS.jl): ")
paramfile = user_input_paramfile == "" ? "PARAMETERS.jl" : user_input_paramfile*".jl"
t3 = time()
include(paramfile)
t4 = time()
println("Time to load PARAMETERS: $(format_sec_or_min(t4-t3))")

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
function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper, dt=1.; time=nothing, check_eigen=true)
    #=
    Collects coordinate tuples of (I, J, V) [row index, column index, value] for a sparse matrix
    representing the chemical jacobian of the atmospheric system. 

    nthis: The atmospheric densities array, but flattened, in the form [n_CO(z=0), n_CO2(z=0)...n_N2Dpl(z=0), n_CO(z=2)...n_N2Dpl(z=250)]
    inactive: A flattened array of the atmospheric densities of any inactive species, same format as nthis. Functionally constant.
    activespecies: List of active species (const)
    inactivespecies: List of inactivespecies (const)
    Jrates: Flattened array of Jrates, same format as nthis.
    Tn, Ti, Te: Temperature-vs-altitude arrays for neutrals, ions, electrons. (const)
    tup, tdown: Transport coefficients
    tlower, tupper: Transport coefficients

    dt: Can be specified to manually control the timestep by which the derivatives are multiplied.
        If not supplied, then dt=1 so that the external solver can manage the multiplication by time.
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
                                                  tdown[:,2]; tlower[:,2]]...) 


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
                                                      tup[:,ialt-1]]...)


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
                                              tupper[:,2]; tup[:,end-1]]...)

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
                                                       # Kinda confusing since we're talking about indices of indices.s
        chemJval[remove_these] .= 0 
    end

    # J = sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +) 
    # println("checking eigenvalues")
    # if check_eigen==true
    #     check_jacobian_eigenvalues(J, results_dir*sim_folder_name)
    #     append!(stiffness, calculate_stiffness(J))
    # end

    return sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +) #J
end

function get_transport_and_J_rates(n, inactive, activesp, inactivesp, neutral_temps::Vector{Float64}, plasma_temps::Vector{Float64}, 
                                   D_arr::Vector{Float64}; compute_Jrates=true)
    #=
    This takes the current densities of all active species (n), transforms
    back into an atmospheric state dictionary, calculates the current Jrates, 
    and then returns the transport coefficients, boundary layer transport
    coefficients, and Jrates needed to run ratefn and chemical_jacobian.
    =#
    
    # transform n vector back into n_current so we can operate on it ---------------------------------------------
    # Time for this block is ~0.0002 seconds.
    n_cur_active = unflatten_atm(n, activesp)
    n_cur_inactive = unflatten_atm(inactive, inactivesp)
    n_cur_all = Dict(vcat([k=>n_cur_active[k] for k in keys(n_cur_active)],
                           [k=>n_cur_inactive[k] for k in keys(n_cur_inactive)]))
    # Retrieve diffusion coefficients and mean scale height profiles for current atmospheric state ----------------
    # Keddy: 1D by altitude (array) independent of species
    # Molecular diffusion: species=>[D by altitude] (dictionary)
    # Mean atmospheric scale height (dictionary of 1D arrays)
    Keddy_arr = zeros(size(alt))
    Keddy_arr .= map(z->Keddy(z, n_tot(n_cur_all, z)), alt) # Eddy diffusion: K coefficient by altitude (array)
    Dcoef_dict = Dict{Symbol, Vector{Float64}}([s=>deepcopy(Dcoef!(D_arr, whichtemps[charge_type(s)], s, n_cur_all, speciesbclist)) for s in fullspecieslist])
    H0_dict = Dict{String, Vector{Float64}}("neutral"=>map((z,t)->scaleH(z, t, n_cur_all), alt, neutral_temps),
                                            "ion"=>map((z,t)->scaleH(z, t, n_cur_all), alt, plasma_temps))

    # check for NaNs in these important inputs to fluxcoefs_all
    report_NaNs(Keddy_arr, name="eddy diffusion coefficients")
    report_NaNs(Dcoef_dict, name="diffusion coefficients")
    report_NaNs(H0_dict, name="mean atmospheric scale height")
    report_NaNs(Hs_dict, name="species-specific scale heights")

    # Calculate things which will be passed into ratefn() -------------------------------
    if compute_Jrates==true
        update_Jrates!(n_cur_all)
        # copy all the Jrates into an external dictionary for storage
        for jr in Jratelist                # time for this is ~0.000005 s
            global Jrate_storage[jr] = n_cur_all[jr]
        end
    end

    # time to flatten Jrates is ~0.00014 s
    Jrates = deepcopy(Float64[Jrate_storage[jr][ialt] for jr in Jratelist, ialt in 1:length(non_bdy_layers)])
    # Jrates = deepcopy(Float64[n_cur_all[sp][ialt] for sp in Jratelist, ialt in 1:length(non_bdy_layers)]) # can't use flatten_atm because it splats, this can't splat
    
    # transport coefficients for bulk layers:
    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    # Temperature arrays and Hs dict are global scope. Rest are created within this function bc they depend on densities.
    # println("construct the fluxcoef dict")
    # @time 
    fluxcoefs_all = fluxcoefs(neutral_temps, plasma_temps, Keddy_arr, Dcoef_dict, H0_dict, Hs_dict)
    # println("timed fluxcoef dict")
    tup = fill(-999., length(fullspecieslist), num_layers)
    tdown = fill(-999., length(fullspecieslist), num_layers)
    for (s, i) in zip(fullspecieslist, collect(1:length(fullspecieslist)))
        tup[i, :] .= fluxcoefs_all[s][2:end-1, 2]
        tdown[i, :] .= fluxcoefs_all[s][2:end-1, 1]
    end
    
    # time for these 3 lines is ~0.0003 seconds.
    # transport coefficients for boundary layers
    bc_dict = boundaryconditions(fluxcoefs_all, speciesbclist)
    tlower = permutedims(reduce(hcat, [bc_dict[sp][1,:] for sp in fullspecieslist]))
    tupper = permutedims(reduce(hcat, [bc_dict[sp][2,:] for sp in fullspecieslist]))
    
    return Jrates, tup, tdown, tlower, tupper
end

function make_jacobian(n, p, t)
    #=
    Constructs the chemical jacobian in the normal way, including stuff to calculate parameters for chemJmat.
    =#

    # n[n .< 0] .= 1e-50
    # println("enter make_jacobian")
    # Unpack the parameters ---------------------------------------------------------------
    inactive, inactivesp, activesp, Tn::Vector{Float64}, Ti::Vector{Float64}, Te::Vector{Float64}, Tp::Vector{Float64}, D_arr::Vector{Float64}, check_eigen = p

    Jrates, tup, tdown, tlower, tupper = get_transport_and_J_rates(n, inactive, activesp, inactivesp, Tn, Tp, D_arr, compute_Jrates=true)
    return chemJmat(n, inactive, activesp, inactivesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper; check_eigen=check_eigen, time=t)
end

function jacobian_wrapper(J, n, p, t)#; write_jac_out=true)
    #=
    dear god
    =#

    J .= make_jacobian(n, p, t)
end

# Production and loss equation functions =======================================
function PnL_eqn(dn, n, p, t)
    #=
    Only exists to call ratefn(), but has the proper arguments requested by the ODE solver. 
    So, this is the function the ODE solver will solve.
    =#

    # n[n .< 0] .= 1e-50

    # Unpack the parameters needed to call ratefn 
    inactive, inactivesp, activesp, Tn::Vector{Float64}, Ti::Vector{Float64}, Te::Vector{Float64}, Tp::Vector{Float64}, D_arr::Vector{Float64} = p 
    
    Jrates, tup, tdown, tlower, tupper = get_transport_and_J_rates(n, inactive, activesp, inactivesp, Tn, Tp, D_arr, compute_Jrates=false)

    println(t) # is this the time we are at? probably.
 
    dn .= ratefn(n, inactive, inactivesp, activesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)  # todo: update dn to dndt
end

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    # println("enter ratefn")
    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] .= ratefn_local([nthismat[:,1]; nthismat[:,2];
                                     fill(1.0, length(activespecies));
                                     inactivemat[:,1]; Jrates[:,1]; 
                                     Tn[1]; Ti[1]; Te[1];
                                     tup[:,1]; tlower[:,1]; tdown[:,2]; tlower[:,2]]...)

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(num_layers-1)
        returnrates[:,ialt] .= ratefn_local([nthismat[:,ialt];
                                            nthismat[:,ialt+1];
                                            nthismat[:,ialt-1];
                                            inactivemat[:,ialt];
                                            Jrates[:,ialt];
                                            Tn[ialt]; Ti[ialt]; Te[ialt];
                                            tup[:,ialt]; 
                                            tdown[:,ialt];
                                            tdown[:,ialt+1]; 
                                            tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude
    returnrates[:,end] .= ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       Tn[end]; Ti[end]; Te[end];
                                       tupper[:,1]; 
                                       tdown[:,end];
                                       tupper[:,2]; 
                                       tup[:,end-1]]...)


    # NEW: Overwrite the entries for water in the lower atmosphere with 0s so that it will behave as fixed.
    # Only runs when water is in the activespecies list. If neutrals are set to inactive, it will be taken care of already.
    if in(:H2O, activespecies) && in(:HDO, activespecies)
        returnrates[H2Oi, 1:upper_lower_bdy_i] .= 0
        returnrates[HDOi, 1:upper_lower_bdy_i] .= 0
    end
    # println("exit ratefn")
    return [returnrates...;]
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
    t_to_save: timesteps at which to save a snapshot of the atmosphere.
    
    plotJratesflag: if set to true, it will plot the Jrates.
    =#

    # Set up the initial state and check for any problems 
    nstart = flatten_atm(atm_init, activespecies)
    find_nonfinites(nstart, collec_name="nstart")
    # Next line calculates 1 ppt of the total density at each altitude for absolute error tolerance,
    # vector in same shape as nstart.
    # absfloor = 1e-3
    abs_tol_vec = 1e-12 .* [[n_tot(atm_init, a) for sp in activespecies, a in non_bdy_layers]...] #0.001#
    # spvec = [[sp for sp in activespecies, a in non_bdy_layers]...]
    # abs_tol_vec[spvec .== "neutral"] .= 1e-2

    # Set up simulation solver parameters and log them
    odesolver = "QNDF"
    # linsolv = :KLU
    saveall = false
    abst = "1 ppt"#, 1e-4 for neutrals"#abs_tol_vec
    relt = rel_tol
    sys_size = length(nstart)
    startdt = 10.0^log_t_start
    tspan = (10.0^log_t_start, 10.0^log_t_end)
    
    Dcoef_arr_template = zeros(size(Tn_arr))  # For making diffusion coefficient calculation go faster 
    params = [inactive, inactivespecies, activespecies, Tn_arr::Vector{Float64}, Ti_arr::Vector{Float64}, Te_arr::Vector{Float64}, 
              Tplasma_arr::Vector{Float64}, Dcoef_arr_template::Vector{Float64}, false]
    params_exjac = [inactive, inactivespecies, activespecies, Tn_arr::Vector{Float64}, Ti_arr::Vector{Float64}, Te_arr::Vector{Float64}, 
                    Tplasma_arr::Vector{Float64}, Dcoef_arr_template::Vector{Float64}, false]

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
                abstol=abs_tol_vec, reltol=relt, 
                isoutofdomain=(u,p,t)->any(x->x<0,u))
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
args = Any[ARGS[i] for i in 1:1:length(ARGS)]#Any["temp", "216", "130", "205", "mean"]#

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

user_input_folder_name = input("Enter a name for the results folder or press enter to use default: ")
const sim_folder_name = user_input_folder_name == "" ? FNext : user_input_folder_name
# Set up the folder if it doesn't exist
create_folder(sim_folder_name, results_dir)

# Converging an altitude grid of a new extent ==================================
make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
while make_new_alt_grid != "y" && make_new_alt_grid != "n"
    println("Bad entry!")
    global make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
end
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
    defaultatm = "converged_neutrals_sundials.h5"
    file_to_use = input("Enter the name of a file containing a converged, 250 km atmosphere to use (press enter to use default: $(defaultatm)): ")   # TODO: revert
    readfile = file_to_use == "" ?  defaultatm : file_to_use * ".h5"
    const n_current = get_ncurrent(readfile)
end

# Setup of new species profiles =======================================================================
converge_which = input("Converging ions, neutrals or both?: ")
while converge_which != "neutrals" && converge_which != "ions" && converge_which != "both"
    println("Bad entry! Please enter ions or neturals or both.")
    global converge_which = input("Converging ions, neutrals or both?: ")
end
# Whether to initialize new species as zeros or not
use_nonzero_initial_profiles = true#false 

# Set up initial profiles and time steps 
if converge_which == "neutrals"
    println("ALERT: Did you make sure that the ions are added to nochemspecies and notransportspecies?")

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
    println("ALERT: Did you make sure that the neutrals are added to nochemspecies and notransportspecies?")

    for ni in new_ions
        n_current[ni] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        println("Initializing non-zero profiles for $(new_ions)")
        for ni in new_ions
            try
                n_current[ni] = reshape(readdlm("../Resources/initial_profiles/$(string(ni))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
            catch
                nothing
            end  
        end
        for Dion in keys(D_H_analogues)
            if in(Dion, new_ions)
                n_current[Dion] = DH .* n_current[D_H_analogues[Dion]]
            end
        end
    end
elseif converge_which == "both" 
    println("ALERT: Did you make sure nochemspecies and notransportspecies are only :Ar and :N2?")
    for nn in new_neutrals
        n_current[nn] = zeros(num_layers)
    end
    for ni in new_ions
        n_current[ni] = zeros(num_layers)
    end

    if use_nonzero_initial_profiles
        println("Initializing non-zero profiles for $(new_neutrals) and $(new_ions)")
        for nn in new_neutrals
            n_current[nn] = reshape(readdlm("../Resources/initial_profiles/$(string(nn))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
        end

        for ni in new_ions
            try
                n_current[ni] = reshape(readdlm("../Resources/initial_profiles/$(string(ni))_initial_profile.txt", '\r', comments=true, comment_char='#'), (num_layers,))
            catch
                nothing
            end  
        end
        for Dion in keys(D_H_analogues)
            if in(Dion, new_ions)
                n_current[Dion] = DH .* n_current[D_H_analogues[Dion]]
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

# Create a separate n_current for inactive species and Jrates which will be updated separately from
# the main solving routine, then merged back in.
const Jrate_storage = Dict([j=>n_current[j] for j in keys(n_current) if !in(j, activespecies)])

# Here you can choose whether to turn chemistry or transport on and off for testing ==================
#do_chem = input("Allow chemistry? (on/off): ")
#while do_chem != "on" && do_chem != "off" 
#    println("Bad entry! Please enter on or off")
#    global do_chem = input("Allow chemistry? (on/off): ")
#end
const do_chem = true#do_chem=="on" ? true : false
#println()

#do_trans = input("Allow transport? (on/off): ")
#while do_trans != "on" && do_trans != "off" 
#    println("Bad entry! Please enter on or off")
#    global do_chem = input("Allow transport? (on/off): ")
#end
const do_trans = true#do_trans=="on" ? true : false
#println()

# PROBLEMS? TODO: Make sure this is set to what you want!
# Use this to zero out all neutrals. Used to attempt to get things to converge all together.
# for n in setdiff(fullspecieslist, [:CO2])
#     n_current[n] = zeros(num_layers)
# end

# Various other needful things ===============================================================================

# Densities of inactive species, which won't change by definition
const inactive = flatten_atm(n_current, inactivespecies)

# Set up the timesteps
const mindt = dt_min_and_max[converge_which][1]
const maxdt = dt_min_and_max[converge_which][2]

# Set solar cycle file
cycle = args[end]
solar_data_file = Dict("max"=>"marssolarphotonflux_solarmax.dat", 
                       "mean"=>"marssolarphotonflux_solarmean.dat", 
                       "min"=>"marssolarphotonflux_solarmin.dat")
const solarfile = solar_data_file[cycle]
# Let the user know what is being done 
if cycle != "mean"
    println("ALERT: Solar $(cycle) data being used")
end

# Convert the arguments to numbers so we can use them to do maths
for i in 2:1:length(args)-1
    args[i] = parse(Float64, args[i])
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

global const controltemps = [T_surf, T_tropo, T_exo]

# set temperature arrays to be used everywhere
const Tn_arr = [Temp_n(a) for a in alt];
const Ti_arr = [Temp_i(a) for a in alt];
const Te_arr = [Temp_e(a) for a in alt];
const Tplasma_arr = 0.5 * (Ti_arr .+ Te_arr);
const whichtemps = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)

# Species-specific scale heights - has to be done here once the control temps are set
Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>map(z->scaleH(z, sp, controltemps), alt) for sp in fullspecieslist])

# plot all 3 profiles on top of each other
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

# obtain the rates and jacobian for each altitude
# TODO: rates_local must also be mapped over a new list which designates whether the species is in photochemical equilibrium or not, 
# and that argument must be passed to getrate.
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x, chem_on=do_chem, trans_on=do_trans), activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

# DEBUG 
# make_chemjac_key("chemjac_local_key.txt", results_dir*sim_folder_name, activespecies, activespecies)
# make_chemjac_key("chemjac_above_key.txt", results_dir*sim_folder_name, activespecies, active_above)
# make_chemjac_key("chemjac_below_key.txt", results_dir*sim_folder_name, activespecies, active_below)

# define a replacement for rates_local
# each row is for each species; each column is for chemical production, chemical loss, 
# transport production, transport loss, in that order.

# pchem_eq_by_species = Dict(:CO2pl=>true);   # TODO: this is to start. change to ionlist
prod_loss_array = Array{Array{Expr}}(undef, length(activespecies), 4)
for (i, asp) in enumerate(activespecies)
    prod_loss_array[i, :] .= getrate(reactionnet, transportnet, asp, sepvecs=true)#pchem_eq=get(pchem_eq_by_species, asp, false), sepvecs=false)
end

# TODO: Turn these notes into a test that looks for "+()" in rates_local instead of debugging lines. 
# -----------------------------------------------------------------------------------------------
# These lines are useful for troubleshooting if you're getting weird errors with ratefn.
# "no method matching +()" means the chemical system is unbalanced and you've got a
# species that has production but no consumption, or vice versa.
# other errors may occur. Uncomment these 3 lines to inspect what goes into rates_local.
# println("The contents of rates_local: ")
# println(rates_local)
# println("That was rates_local")

const arglist_local = [activespecies; active_above; active_below; inactivespecies;
                      Jratelist; :Tn; :Ti; :Te; local_transport_rates]

const arglist_local_typed = [:($s::Float64) for s in arglist_local]

# These expressions which are evaluated below enable a more accurate assessment of M and E values
# by calculating at the time of being called rather than only at each timestep.
const Mexpr = Expr(:call, :+, fullspecieslist...)
const Eexpr = Expr(:call, :+, ionlist...)

# NOTE: These functions within @eval cannot be moved. Do not move them.
# @eval begin
#     function ratefn_local_eq_ok($(arglist_local_typed...))

#         # M and E are calculated here to ensure that the right number of ions/electrons
#         # is used. It is for only the altitude at which this function was called 
#         # (i.e. all the arguments to the function, when it's called, are for only one altitude)
#         M = $Mexpr
#         E = $Eexpr

#         # Here: calculate new concentrations for species in photochemical equilibrium.
#         # We can pass an argument to this function telling us whether to do photochemical equilibrium ornot.
#         # We can also identify the position in arglist_local(_typed) of whatever species we calculate for.

#         # first here, we need to get seaprate vectors of chemical production and loss expressions just like in the function below.
#         # the loss expression should come out ready to go for photochemical equilibrium for designated species--i.e. with n_s removed. 
#         result = map(_ -> Float64[], $prod_loss_array)   # will hold evaluated expressions.

#         # then we evaluate the expressions in the unrolled loop as below and sum up each so we have total chem prod and total chem loss.
#         # (and also transport)
#         $( (quote
#             # i = row, j = vector, k = specific expression
#                 push!(result[$i], $(expr))
#             end 
#             for (i, expr_list) in enumerate(prod_loss_array)                                                                                                                 
#                 for expr in expr_list
#             )...  
#          )

#         # then we take n = P/L. This is the value it should be.
#         # ----working----
#         # --- do special P+L addition

#         # but rates_local is basically dn/dt, so we need to subtract off the current density of the species from n here.
#         # the result is our dn/dt for that species. So we assign it to net_chem_change. 
#         net_chem_change = zeros(size(result, 1))
#         net_trans_change = zeros(size(result, 1))
        
#         for r in 1:size(result,1)
#             net_chem_change[r] = subtract_difflength(sort(result[r, :][1], rev=true), sort(result[r, :][2], rev=true))
#             net_trans_change[r] = subtract_difflength(sort(result[r, :][3], rev=true), sort(result[r, :][4], rev=true))
#         end

#         # lastly, we return net_chem_change + net_trans_change.

#         $rates_local # evaluates the rates_local expression
#     end
# end


@eval begin
    function ratefn_local_original($(arglist_local_typed...))

        # M and E are calculated here to ensure that the right number of ions/electrons
        # is used. It is for only the altitude at which this function was called 
        # (i.e. all the arguments to the function, when it's called, are for only one altitude)
        M = $Mexpr
        E = $Eexpr

        $rates_local # evaluates the rates_local expression
    end
end


@eval begin
    function ratefn_local($(arglist_local_typed...))

        # M and E are calculated here to ensure that the right number of ions/electrons
        # is used. It is for only the altitude at which this function was called 
        # (i.e. all the arguments to the function, when it's called, are for only one altitude)
        M = $Mexpr
        E = $Eexpr

        # stack overflow - answer (Cite)
        # create a result array for evaluating the production and loss expressions
        result = map(_ -> Float64[], $prod_loss_array) 

        $( (quote
            # i = row, j = vector, k = specific expression
                push!(result[$i], $(expr))
            end 
            for (i, expr_list) in enumerate(prod_loss_array)                                                                                                                 
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
    function chemJmat_local($(arglist_local_typed...))

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

# PROBLEMS? TODO: revert 
# if !(setdiff(solarflux[:,1],lambdas)==[])
#     throw("Need a broader range of solar flux values!")
# end

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
plot_atm(n_current, [neutrallist, ionlist], results_dir*sim_folder_name*"/initial_atmosphere.png")

# Create a list to keep track of stiffness ratio ===============================
const stiffness = []

# Record setup time
t6 = time()

write(f, "Setup time $(format_sec_or_min(t6-t5))\n")
close(f)

# do the convergence ===========================================================
println("Beginning Convergence")
ti = time()
times_to_save = [10.0^t for t in mindt:maxdt]
# Pack the global variables and send them through for faster code!
atm_soln = evolve_atmosphere(n_current, mindt, maxdt, t_to_save=times_to_save) # dz, transportspecies, alt, fullspecieslist, num_layers, Jratelist, non_bdy_layers, 
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

    nc = unflatten_atm(atm_state, activespecies)
    nc_inactive = unflatten_atm(inactive, inactivespecies)
    
    # add back in the inactive species and Jrates
    for isp in inactivespecies
        nc[isp] = nc_inactive[isp]
    end

    if i == L
        for jsp in Jratelist  # TODO: This needs to be worked so Jrates are tracked for each timestep. right now, only the last one is kept. :/
            nc[jsp] = Jrate_storage[jsp]
        end
        plot_atm(nc, [neutrallist, ionlist], results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state")
    end
    global i += 1

    filepath = results_dir*sim_folder_name*"/atm_state_t_$(timestep).h5"
    write_ncurrent(nc, filepath)   
end
t6 = time()

println("Wrote all states to files")

t9 = time()

write(f, "Stiffness ratio evolution:\n")
write(f, "$(stiffness)")
write(f, "\n\n")
write(f, "Total runtime $(format_sec_or_min(t9-t1))\n")

close(f)



println("Total runtime $(format_sec_or_min(t9-t1))")