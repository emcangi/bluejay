###############################################################################
# converge_new_file.jl
# TYPE: (1) Model files - required
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere.
# 
# Eryn Cangi
# Created 2018
# Last edited:December 2021 
# Currently tested for Julia: 1.6.1
###############################################################################

# **************************************************************************** #
#                                                                              #
#                              MODULE LOADING                                  #
#                                                                              #
# **************************************************************************** #

using Dates

t1 = time()
println("$(Dates.format(now(), "(HH:MM:SS)")) Start")

#logging witchcraft
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Revise
using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using DelimitedFiles
using SparseArrays
using LinearAlgebra
photochemistry_source_dir = "$(@__DIR__)/Photochemistry/src/"
println("loading Photochemistry.jl from $photochemistry_source_dir")
push!(LOAD_PATH, photochemistry_source_dir)
using Photochemistry  # custom module
using DifferentialEquations
using IterativeSolvers
using IncompleteLU
using Dates

t2 = time()
println("$(Dates.format(now(), "(HH:MM:SS)")) Loaded modules in $(format_sec_or_min(t2-t1))")

# **************************************************************************** #
#                                                                              #
#                           LOAD PARAMETERS                                    #
#                                                                              #
# **************************************************************************** #

include("CONSTANTS.jl")
include("CUSTOMIZATIONS.jl")

paramfile = get_paramfile(code_dir)
include(paramfile)

t3 = time()
println("Running a $(simtype) experiment with parameters $(controltemps) in the form of a $(problem_type == "SS" ? "steady state" : problem_type) problem")
println("$(Dates.format(now(), "(HH:MM:SS)")) Loaded parameters in $(format_sec_or_min(t3-t2))")

t4 = time()

# **************************************************************************** #
#                                                                              #
#                              PLOT STYLING                                    #
#                                                                              #
# **************************************************************************** #
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22

# **************************************************************************** #
#                                                                              #
#                                FUNCTIONS                                     #
#                                                                              #
# **************************************************************************** #

#=
These functions are required to be in this file for one of two reasons:
1) Because they call a function that is couched in an @eval statement, which 
   cannot be relocated,
2) They are the main routine functions and will not be shared by any other scripts.
=#

# Basic functions ============================================================

function evolve_atmosphere(atm_init::Dict{Symbol, Array{ftype_ncur, 1}}, log_t_start, log_t_end; t_to_save=[], abstol=1e-12, reltol=1e-2, globvars...)
    #=
    Sets up the initial conditions for the simulation and calls the ODE solver. 

    Input:
        atm_init: dictionary of species densities by altitude for the current state.
        log_t_start: Starting time will be 10^log_t_start.
        log_t_end: Similar, but end time.
        t_to_save: timesteps at which to save a snapshot of the atmosphere.

    Output:
        sol: Solver object, the solution of the atmospheric system at the end of the simulation time.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :neutral_species, :ion_species, :transport_species,
                                   :inactivesp, :active_species, :activellsp, :activeslsp, 
                                   :n_inactive, 
                                   :Tn, :Ti, :Te, :Tp, 
                                   :zmax, :alt, :num_layers, :non_bdy_layers, :dz, 
                                   :plot_grid, :speciescolor, :speciesstyle,
                                   :absorber, :Jratelist, :crosssection, 
                                   :electron_val, :gearmode, 
                                   :Hs_dict, :bcdict, :H2Oi, :HDOi, :upper_lower_bdy_i, 
                                   :molmass, :n_alt_index, :polarizability, :n_all_layers, :Tprof_for_diffusion,
                                   :Dcoef_arr_template, :q])

    println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up initial state")

    # Get the start and stop time information 
    startdt = 10.0^log_t_start
    tspan = (10.0^log_t_start, 10.0^log_t_end)
    
    # Set up the initial state and check for any problems 
    M = sum([n_current[sp] for sp in GV.all_species])
    if GV.electron_val=="constant" # E FIX ATTEMPT 
        E = [1e5 for i in GV.non_bdy_layers]
    elseif electron_val=="quasineutral"
        E = sum([n_current[sp] for sp in GV.ion_species])
    else
        throw("Unhandled electron profile specification: $(elecval)")
    end

    nstart = flatten_atm(atm_init, GV.activellsp; GV.num_layers)
    find_nonfinites(nstart, collec_name="nstart")

    # Set up parameters
    println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up parameters")
    # Dcoef_arr_template = zeros(size(GV.Tn))  # For making diffusion coefficient calculation go faster 
    params = [ # THis one used for ODE solver. 
              #GV.inactive, 
              #GV.inactive_species, GV.active_species, GV.active_longlived, GV.active_shortlived, 
              #GV.Tn, GV.Ti, GV.Te, GV.Tp, 
              GV.Dcoef_arr_template, 
              M, E] # E FIX ATTEMPT
    params_exjac = deepcopy(params)  # I think this is so the Dcoef doesn't get filled in with the wrong info?

    # Set up an example Jacobian 
    println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up example jacobian")
    dndt, example_jacobian = get_rates_and_jacobian(nstart, params_exjac, 0.0; globvars...)
    find_nonfinites(example_jacobian, collec_name="example_jacobian")

    sparsity = round(length(example_jacobian.nzval)*100/(size(example_jacobian)[1] * size(example_jacobian)[2]), digits=2)
    if sparsity > 1
        println("Warning! Sparsity of the jacobian is rather high: $(sparsity)%")
    end

    # Now define the problem function to be solved
    println("$(Dates.format(now(), "(HH:MM:SS)")) Defining ODE function")
    odefunc = ODEFunction(PnL_eqn, jac=jacobian_wrapper, jac_prototype=example_jacobian)

    # Log solver options and run the solver
    more_to_write = ["\nJACOBIAN INFORMATION:", 
                     "Jacobian sparsity: $(sparsity)", 
                     "jacobian total elements: $(length(nstart)^2)", 
                     "\nSOLVER OPTIONS:",
                     "system size: $(length(nstart))",
                     "timespan=$(tspan)",
                     "starting dt=$(startdt)",
                     "saveat=$(t_to_save)", 
                     "abstol=$(abs_tol)",
                     "reltol=$(rel_tol)"]
    write_to_log(logfile, more_to_write, mode="a")

    solver = QNDF

    write_time_stuff = ["\nTIMING",
                        "Start time: $(Dates.format(now(), "yyyy-mm-dd at HH:MM:SS"))"]
    
    if problem_type == "SS"
        write_to_log(logfile, "\nODE solver: DynamicSS/$(solver)", mode="a")
        write_to_log(logfile, write_time_stuff, mode="a")
        println("$(Dates.format(now(), "(HH:MM:SS)")) Defining SS problem")
        probSS = SteadyStateProblem{isinplace}(odefunc, nstart, params)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Calling the solver")
        sol = solve(probSS, DynamicSS(QNDF(), abstol=abs_tol, reltol=rel_tol), saveat=t_to_save, progress=true, save_everystep=false, dt=startdt,
                    abstol=abstol, reltol=reltol, isoutofdomain=(u,p,t)->any(x->x<0,u))

    elseif problem_type == "ODE"
        write_to_log(logfile, "\nODE solver: $solver", mode="a")
        write_to_log(logfile, write_time_stuff, mode="a")
        println("$(Dates.format(now(), "(HH:MM:SS)")) Defining ODE problem")
        probODE = ODEProblem(odefunc, nstart, tspan, params)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Calling the solver")
        sol = solve(probODE, QNDF(), saveat=t_to_save, progress=true, save_everystep=false, dt=startdt,
                    abstol=abstol, reltol=reltol, isoutofdomain=(u,p,t)->any(x->x<0,u))
    elseif problem_type == "Gear"
        write_to_log(logfile, ["\nODE solver: Custom First Order Gear", "Gear timestep type: $(GV.gearmode)", "n_steps for static Gear timesteps: $(GV.n_steps)"], mode="a")
        write_to_log(logfile, write_time_stuff, mode="a")
        println("$(Dates.format(now(), "(HH:MM:SS)")) Calling the solver")
        sol = converge(atm_init, log_t_start, log_t_end; abstol=abstol, reltol=reltol, Dcoef_arr_template, globvars...)
    else
        throw("Invalid problem_type")
    end 

    return sol 
end

function record_atmospheric_state(t, n, actively_solved; globvars...)
    #=
    Input:
        t: current simulation time
        n: current densities of active longlived species 
        actively_solved: species which are actively solved for. Required to compile the full atmospheric state
    Output:
        plot of the atmospheric densities at the present time
        .h5 file containing the present atmospheric state 
        Also prints a progress alert to the terminal and writes out the time to a logfile. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:neutral_species, :ion_species, :plot_grid, :speciescolor, :speciesstyle, :zmax, :alt, :num_layers])

    # This is just to change how many decimal places to include depending if t >= 1.
    rounding_digits = t <= 1 ? Int64(ceil(abs(log10(t)))) : 0 

    progress_alert = "$(Dates.format(now(), "(HH:MM:SS)")) reached timestep $(t), total elapsed time=$(format_sec_or_min(time() - ti))"
    write_to_log(logfile, progress_alert)
    println(progress_alert) # prints present timestep, but only if it's a factor of 10 greater than the last stored timestep (so we don't get a trillion outputs)
    global timestorage = t # updates the last stored timestep to repeat process
    
    # write out the current atmospheric state to a file 
    atm_snapshot = merge(external_storage, unflatten_atm(n, actively_solved; GV.num_layers))
    
    plot_atm(atm_snapshot, results_dir*sim_folder_name*"/atm_peek_$(plotnum).png", abs_tol_for_plot; t="$(round(t, digits=rounding_digits))", globvars...)
    write_atmosphere(atm_snapshot, results_dir*sim_folder_name*"/atm_dt=$(round(t, digits=rounding_digits)).h5"; GV.alt, GV.num_layers)
    global plotnum += 1
end

function chemJmat(n_active_longlived, n_active_shortlived, n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; # E FIX ATTEMPT
                  check_eigen=false, globvars...) 
    #=
    Collects coordinate tuples of (I, J, V) [row index, column index, value] for a sparse matrix
    representing the chemical jacobian of the atmospheric system. 

    Input:
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
    Output:
        sparse matrix representing the chemical jacobian 
    =#              

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :activellsp, :activeslsp, :inactivesp, :num_layers, :H2Oi, :HDOi, :upper_lower_bdy_i])

    nmat_llsp = reshape(n_active_longlived, (length(GV.activellsp), GV.num_layers))
    nmat_slsp = reshape(n_active_shortlived, (length(GV.activeslsp), GV.num_layers))
    nmat_inactive = reshape(n_inactive, (length(GV.inactivesp), GV.num_layers))
    
    # For storing the jacobian indices and values
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = ftype_chem[]

    # tc___ are the coordinate tuples containing (I, J, V) to be used to fill a sparse matrix.
    argvec = [nmat_llsp[:, 1];  # active_longlived;
              nmat_llsp[:, 2]; # active_longlived_above;
              fill(1.0, length(GV.activellsp)); # active_longlived_below;
              nmat_slsp[:, 1]; # active_shortlived;
              nmat_inactive[:,1];  # inactive_species;
              Jrates[:,1]; # Jratelist;
              GV.Tn[1]; GV.Ti[1]; GV.Te[1];  #:Tn; :Ti; :Te;
              M[1]; E[1];  # E FIX ATTEMPT
              tup[:,1]; tlower[:,1]; # local_transport_rates
              tdown[:,2]; tlower[:,2]]
    argvec = convert(Array{ftype_chem}, argvec)

    (tclocal, tcupper, tclower) = chemJmat_local(argvec...) 


    # add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])

    # and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(GV.activellsp))
    append!(chemJval, tcupper[3])

    for ialt in 2:(GV.num_layers-1)
        argvec = [nmat_llsp[:, ialt];
                  nmat_llsp[:, ialt+1];
                  nmat_llsp[:, ialt-1];
                  nmat_slsp[:, ialt];
                  nmat_inactive[:, ialt];
                  Jrates[:, ialt];
                  GV.Tn[ialt]; GV.Ti[ialt]; GV.Te[ialt];
                  M[ialt]; E[ialt];  # E FIX ATTEMPT
                  tup[:, ialt];
                  tdown[:, ialt];
                  tdown[:, ialt+1];
                  tup[:, ialt-1]]
        argvec = convert(Array{ftype_chem}, argvec)        

        (tclocal, tcupper, tclower) = chemJmat_local(argvec...)

        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(GV.activellsp))
        append!(chemJj, tclocal[2].+(ialt-1)*length(GV.activellsp))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(GV.activellsp))
        append!(chemJj, tcupper[2].+(ialt  )*length(GV.activellsp))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(GV.activellsp))
        append!(chemJj, tclower[2].+(ialt-2)*length(GV.activellsp))
        append!(chemJval, tclower[3])
    end

    argvec = [nmat_llsp[:,end];
              fill(1.0, length(GV.activellsp));
              nmat_llsp[:,end-1];
              nmat_slsp[:, end];
              nmat_inactive[:,end];
              Jrates[:,end];
              GV.Tn[end]; GV.Ti[end]; GV.Te[end];
              M[end]; E[end]; # E FIX ATTEMPT
              tupper[:,1]; tdown[:,end];
              tupper[:,2]; tup[:,end-1]]
    argvec = convert(Array{ftype_chem}, argvec)
    
    (tclocal, tcupper, tclower) = chemJmat_local(argvec...)

    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(GV.num_layers-1)*length(GV.activellsp))
    append!(chemJj, tclocal[2].+(GV.num_layers-1)*length(GV.activellsp))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(GV.num_layers-1)*length(GV.activellsp))
    append!(chemJj, tclower[2].+(GV.num_layers-2)*length(GV.activellsp))
    append!(chemJval, tclower[3])

    # fix water below whatever we set as upper/lower atmosphere boundary.
    # This only runs if water is designated as an active species; if it's in inactive_species, this won't run,
    # When it is active, this finds all the H2O and HDO indices for the lower atmosphere. 
    # It's like above where we add (ialt-1)*length(active_species), but this way it's outside the loop.
    if in(:H2O, GV.activellsp) && in(:HDO, GV.activellsp)
        H2Opositions = GV.H2Oi .+ length(GV.activellsp)*collect(0:GV.upper_lower_bdy_i-1)
        HDOpositions = GV.HDOi .+ length(GV.activellsp)*collect(0:GV.upper_lower_bdy_i-1)
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

function ratefn(n_active_longlived, n_active_shortlived, n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; 
                globvars...)
    #=
    at each altitude, get the appropriate group of concentrations, coefficients, and rates to pass to ratefn_local.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :activellsp, :activeslsp, :inactivesp, :num_layers, :H2Oi, :HDOi, :upper_lower_bdy_i])

    nmat_llsp = reshape(n_active_longlived, (length(GV.activellsp), GV.num_layers))
    nmat_slsp = reshape(n_active_shortlived, (length(GV.activeslsp), GV.num_layers))
    nmat_inactive = reshape(n_inactive, (length(GV.inactivesp), GV.num_layers))
    returnrates = zeros(size(nmat_llsp))
    
    # fill the first altitude entry with information for all species
    argvec = [nmat_llsp[:,1]; # densities for active_longlived;
              nmat_llsp[:,2]; # active_longlived_above;
              fill(1.0, length(GV.activellsp)); #active_longlived_below;
              nmat_slsp[:, 1]; # active_shortlived;
              nmat_inactive[:,1];  # inactive_species;
              Jrates[:,1];  # Jratelist;
              GV.Tn[1]; GV.Ti[1]; GV.Te[1];  # :Tn; :Ti; :Te;
              M[1]; E[1];  # E FIX ATTEMPT
              tup[:,1]; tlower[:,1]; tdown[:,2]; tlower[:,2]]
    argvec = convert(Array{ftype_chem}, argvec)
    
    returnrates[:,1] .= ratefn_local(argvec...) # local_transport_rates

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(GV.num_layers-1)
        argvec = [nmat_llsp[:, ialt]; # active_longlived;
                  nmat_llsp[:, ialt+1];
                  nmat_llsp[:, ialt-1];
                  nmat_slsp[:, ialt]; # active_shortlived;
                  nmat_inactive[:,ialt];
                  Jrates[:,ialt];
                  GV.Tn[ialt]; GV.Ti[ialt]; GV.Te[ialt];
                  M[ialt]; E[ialt];  # E FIX ATTEMPT
                  tup[:, ialt];
                  tdown[:, ialt];
                  tdown[:, ialt+1];
                  tup[:, ialt-1]]
        argvec = convert(Array{ftype_chem}, argvec)
        
        returnrates[:,ialt] .= ratefn_local(argvec...)
    end

    # fill in the last level of altitude
    argvec = [nmat_llsp[:, end];
              fill(1.0, length(GV.activellsp));
              nmat_llsp[:, end-1];
              nmat_slsp[:, end]; # active_shortlived;
              nmat_inactive[:,end];
              Jrates[:,end];
              GV.Tn[end]; GV.Ti[end]; GV.Te[end];
              M[end]; E[end]; # E FIX ATTEMPT
              tupper[:,1];
              tdown[:,end];
              tupper[:,2];
              tup[:,end-1]]
    argvec = convert(Array{ftype_chem}, argvec)
    returnrates[:,end] .= ratefn_local(argvec...)


    # NEW: Overwrite the entries for water in the lower atmosphere with 0s so that it will behave as fixed.
    # Only runs when water is in the active_species list. If neutrals are set to inactive, it will be taken care of already.
    if in(:H2O, GV.activellsp) && in(:HDO, GV.activellsp)
        returnrates[GV.H2Oi, 1:GV.upper_lower_bdy_i] .= 0
        returnrates[GV.HDOi, 1:GV.upper_lower_bdy_i] .= 0
    end

    return [returnrates...;]
end

function compile_ncur_all(n_long, n_short, n_inactive; globvars...)
    #=
    While the simulation runs, "n", the vector passed to the solver, only contains densities
    for long-lived, active species. Every time the atmospheric state changes, the transport coefficients
    and Jrates must be updated, but those all depend on the densities of ALL species. It's easiest for 
    the functions updating those things to pull from one atmospheric state dictionary, 
    so this function combines disparate density vectors back into one dictionary.

    n_long: active, long-lived species densities
    n_short: same but for short-lived species
    n_inactive: inactive species densities (truly, these never change)
    activellsp, activeslsp, inactivesp: lists of species symbols for longlived, shortlived, inactive.

    Returns: atmospheric state dictionary of species densities only (no Jrates).
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:activellsp, :activeslsp, :inactivesp, :num_layers])


    n_cur_active_long = unflatten_atm(n_long, GV.activellsp; num_layers=GV.num_layers)
    n_cur_active_short = unflatten_atm(n_short, GV.activeslsp; num_layers=GV.num_layers)
    n_cur_inactive = unflatten_atm(n_inactive, GV.inactivesp; num_layers=GV.num_layers)

    n_cur_all = Dict(vcat([k=>n_cur_active_long[k] for k in keys(n_cur_active_long)],
                          [k=>n_cur_active_short[k] for k in keys(n_cur_active_short)],
                          [k=>n_cur_inactive[k] for k in keys(n_cur_inactive)]))
    
    return n_cur_all
end

function update_Jrates!(n_cur_densities::Dict{Symbol, Array{ftype_ncur, 1}}; globvars...)
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Jratelist, :crosssection, :num_layers, :absorber, :dz])

    # Initialize an array, length=number of active layers
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, GV.num_layers)
    for i in range(1, length=GV.num_layers)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in GV.Jratelist
        species = GV.absorber[jspecies]

        jcolumn = convert(Float64, 0.)
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituent
            jcolumn += convert(Float64, n_cur_densities[species][ialt])*GV.dz

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
            BLAS.axpy!(nlambda, jcolumn, GV.crosssection[jspecies][ialt+1], 1, solarabs[ialt], 1)
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
    for j in GV.Jratelist
        n_cur_densities[j] = zeros(GV.num_layers) # TODO: Can probably change size(n_cur_densities[:CO2]) to num_layers.
        for ialt in [1:nalt;]
            n_cur_densities[j][ialt] = ftype_ncur(BLAS.dot(nlambda, solarabs[ialt], 1, GV.crosssection[j][ialt+1], 1))
        end
    end
end

# Julia DifferentialEquations.jl solver =======================================
# TODO: pub globvars in parameters.
function make_jacobian(n, p, t)
    #=
    Constructs the chemical jacobian in the normal way, including stuff to calculate parameters for chemJmat.
    =#

    # Unpack the parameters ---------------------------------------------------------------
    n_inactive, inactivesp, activesp, activellsp, activeslsp, Tn, Ti, Te, Tp, D_arr, = p 

    # get the concentrations of species assumed to be in photochemical equilibrium. 
    n_short = flatten_atm(external_storage, activeslsp; num_layers)  # retrieve the shortlived species from their storage and flatten them

    # Update Jrates
    n_cur_all = compile_ncur_all(n, n_short, n_inactive; activellsp, activeslsp, inactivesp, num_layers)

    update_Jrates!(n_cur_all; Jratelist, crosssection, num_layers, absorber, dz)
    # copy all the Jrates into an external dictionary for storage
    for jr in Jratelist                # time for this is ~0.000005 s
        global external_storage[jr] = n_cur_all[jr]
    end

    # Retrieve Jrates 
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:num_layers])

    # and update the shortlived species with the new Jrates - assuming not needed to be done in this function
    # n_short_updated = set_concentrations!(external_storage, n, n_short, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)

    tlower, tup, tdown, tupper = update_transport_coefficients(transport_species, # Species for which to update coefficients so it's not a mistake to pass it twice.
                                                               n_cur_all, activellsp, D_arr;
                                                               Tn, Tp, Hs_dict, bcdict=speciesbclist, 
                                                               all_species, neutral_species, transport_species, molmass, n_alt_index, 
                                                               polarizability, alt, num_layers, n_all_layers, dz, T_for_diff=Tprof_for_diffusion, q)
    
    return chemJmat(n, n_short, n_inactive, Jrates, tup, tdown, tlower, tupper; activellsp, activeslsp, inactivesp, Tn, Ti, Te, num_layers, H2Oi, HDOi, upper_lower_bdy_i)
end

function jacobian_wrapper(J, n, p, t)
    #=
    A wrapper since the Jacobian is sometimes generated outside of the normal simulation run
    =#

    J .= make_jacobian(n, p, t)
end

# TODO: globvars in parameters
function PnL_eqn(dndt, n, p, t)
    #=
    This is the primary function that is solved by the solver. For our purposes that means it calls ratefn().

    dndt is written as such because the output of this function is silently multiplied by the timestep within the solver's code.
    n is a vector of active species densities by altitude.
    =#

    # Unpack the parameters needed to call ratefn 
    n_inactive, inactivesp, activesp, activellsp, activeslsp, Tn, Ti, Te, Tp, D_arr = p 

    # Every time we increase dt by a factor of 10, print a progress report to the console,
    # plot the atmosphere, and save the atmosphere to an .h5 file. 
    if t / timestorage >= 10
        record_atmospheric_state(t, n, activellsp; neutral_species, ion_species, plot_grid, speciescolor, speciesstyle, zmax, alt, num_layers)
    end

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, activeslsp; num_layers)

    # Retrieve the Jrates
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:num_layers])

    # set the concentrations of species assumed to be in photochemical equilibrium. 
    n_short_updated = set_concentrations!(external_storage, n, n_short, n_inactive, Jrates; activellsp, activeslsp, inactivesp, Tn, Ti, Te, num_layers)
    # for i in 1:2
    #     n_short_updated = set_concentrations!(external_storage, n, n_short_updated, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)
    # end

    # Get the updated transport coefficients, taking into account short-lived species update
    updated_ncur_all = compile_ncur_all(n, n_short_updated, n_inactive; activellsp, activeslsp, inactivesp, num_layers)
    tlower, tup, tdown, tupper = update_transport_coefficients(transport_species, updated_ncur_all, activellsp, D_arr; 
                                                               Tn, Tp, Hs_dict, bcdict=speciesbclist,
                                                               all_species, neutral_species, transport_species, 
                                                               molmass, alt, n_alt_index, polarizability, num_layers, dz, T_for_diff=Tprof_for_diffusion, n_all_layers, q)

    dndt .= ratefn(n, n_short_updated, inactive, Jrates, tup, tdown, tlower, tupper; activellsp, activeslsp, inactivesp, Tn, Ti, Te, num_layers, H2Oi, HDOi, upper_lower_bdy_i)
end

function set_concentrations!(external_storage, n_active_long, n_active_short, n_inactive, Jrates, M, E; globvars...) # E FIX ATTEMPT
    #=
    at each altitude, sets the concentrations for short-lived species assumed to be in photochemical equilibrium
    and sends them back into the storage dictionary, external_storage

    external_storage: dictionary storing densities for short-lived and inactive species, as well as Jrates.
    n_active_short, n_active_long, n_inactive: density of short-lived, long-lived, and inactive species
    activeslsp, activellsp, inactivesp:: list of short- and long-lived species names
    Jrates: Jrates for each species, for a particular altitude
    Tn, Ti, Te: temperature arrays
    E: electrons 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:activellsp, :activeslsp, :inactivesp, :Tn, :Ti, :Te, :num_layers])
    
    # rows = species, columns = altitudes. 
    nmat_shortlived = reshape(n_active_short, (length(GV.activeslsp), GV.num_layers))
    
    # auxiliary information that is needed. 
    nmat_longlived = reshape(n_active_long, (length(GV.activellsp), GV.num_layers))
    nmat_inactive = reshape(n_inactive, (length(GV.inactivesp), GV.num_layers))
    
    # storage for the updated concentrations
    new_densities = zeros(size(nmat_shortlived))
    # dist_zero = zeros(length(nmat_shortlived)) # TODO: Figure out how to make useful

    # fill the first altitude entry with information for all species   
    argvec = [nmat_shortlived[:,1]; nmat_longlived[:, 1]; nmat_inactive[:, 1]; Jrates[:, 1]; GV.Tn[1]; GV.Ti[1]; GV.Te[1]; M[1]; E[1]] # E FIX ATTEMPT
    argvec = convert(Array{ftype_chem}, argvec)
    new_densities[:,1] .= set_concentrations_local(argvec...)
    # dist_zero[1] = check_zero_distance([nmat_shortlived[:,1]; nmat_longlived[:, 1]; nmat_inactive[:, 1]; Jrates[:, 1]; Tn[1]; Ti[1]; Te[1]]...)

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(GV.num_layers-1)
        argvec = [nmat_shortlived[:,ialt]; nmat_longlived[:, ialt]; nmat_inactive[:, ialt]; Jrates[:, ialt]; GV.Tn[ialt]; GV.Ti[ialt]; GV.Te[ialt]; M[ialt]; E[ialt]] # E FIX ATTEMPT
        argvec = convert(Array{ftype_chem}, argvec)
        new_densities[:,ialt] .= set_concentrations_local(argvec...)
        # dist_zero[ialt] = check_zero_distance([nmat_shortlived[:,ialt]; nmat_longlived[:, ialt]; nmat_inactive[:, ialt]; Jrates[:, ialt]; Tn[ialt]; Ti[ialt]; Te[ialt]]...)
    end

    # fill in the last level of altitude
    argvec = [nmat_shortlived[:, end]; nmat_longlived[:, end]; nmat_inactive[:, end]; Jrates[:, end]; GV.Tn[end]; GV.Ti[end]; GV.Te[end]; M[end]; E[end]] # E FIX ATTEMPT
    argvec = convert(Array{ftype_chem}, argvec)
    new_densities[:,end] .= set_concentrations_local(argvec...)
    # dist_zero[end] = check_zero_distance([nmat_shortlived[:, end]; nmat_longlived[:, end]; nmat_inactive[:, end]; Jrates[:, end]; Tn[end]; Ti[end]; Te[end]]...)
    
    # write out the new densities for shortlived species to the external storage
    for (s, ssp) in enumerate(GV.activeslsp)
        external_storage[ssp] .= new_densities[s, :]
    end
    
    return vec(new_densities)

    # Look at distance from zero
    # if any(x->x>100, dist_zero)
    #     println("elements >100 from zero: $(dist_zero[findall(x->x>100, dist_zero)])")
    # end
end

# Custom Gear Solver ===========================================================

function get_rates_and_jacobian(n, p, t; globvars...)
    #=
    Three major tasks:
    1. Prints progress updates to the console and saves present state to files every time
       dt has increased by a factor of 10.
    2. Updates and saves Jrates
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :neutral_species, :ion_species, 
                                   :inactivesp, :active_species, :activellsp, :activeslsp, :transport_species,
                                   :n_inactive, 
                                   :Tn, :Ti, :Te, :Tp, 
                                   
                                   :plot_grid, :speciescolor, :speciesstyle, 
                                   :zmax, :alt, :num_layers, :non_bdy_layers, :dz,
                                   :absorber, :Jratelist, :crosssection, 
                                   
                                   
                                   :Hs_dict, :bcdict, :H2Oi, :HDOi, :upper_lower_bdy_i, 
                                   :molmass, :n_alt_index, :polarizability, :n_all_layers, :Tprof_for_diffusion, :q])

    # Unpack the parameters ---------------------------------------------------------------
    #=n_inactive, inactivesp, activesp, activellsp, activeslsp, Tn, Ti, Te, Tp,=# D_arr, M, E = p # E FIX ATTEMPT

    if t / timestorage >= 10
        record_atmospheric_state(t, n, GV.activellsp; globvars...)#neutral_species, ion_species, plot_grid, speciescolor, speciesstyle, zmax, alt, num_layers)
    end

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, GV.activeslsp; GV.num_layers)

    # Update Jrates
    n_cur_all = compile_ncur_all(n, n_short, GV.n_inactive; GV.activellsp, GV.activeslsp, GV.inactivesp, GV.num_layers)

    update_Jrates!(n_cur_all; GV.Jratelist, GV.crosssection, GV.num_layers, GV.absorber, GV.dz)
    # copy all the Jrates into an external dictionary for storage
    for jr in GV.Jratelist                # time for this is ~0.000005 s
        global external_storage[jr] = n_cur_all[jr]
    end

    # Retrieve Jrates 
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in GV.Jratelist, ialt in 1:GV.num_layers])

    # set the concentrations of species assumed to be in photochemical equilibrium. 
    n_short_updated = set_concentrations!(external_storage, n, n_short, GV.n_inactive, Jrates, M, E; GV.activellsp, GV.activeslsp, GV.inactivesp, GV.Tn, GV.Ti, GV.Te, GV.num_layers)# E FIX ATTEMPT
    # for i in 1:2
    #     n_short_updated = set_concentrations!(external_storage, n, n_short_updated, n_inactive, activellsp, activeslsp, inactivesp, Jrates, Tn, Ti, Te)
    # end

    # Get the updated transport coefficients, taking into account short-lived species update
    updated_ncur_all = compile_ncur_all(n, n_short_updated, GV.n_inactive; GV.activellsp, GV.activeslsp, GV.inactivesp, GV.num_layers)
    tlower, tup, tdown, tupper = update_transport_coefficients(GV.transport_species, updated_ncur_all, GV.activellsp, D_arr; globvars...)
                                                            #    GV.Tn, GV.Tp, GV.Hs_dict, bcdict=speciesbclist,
                                                            #    all_species, neutral_species, transport_species, 
                                                            #    molmass, n_alt_index, polarizability, alt, num_layers, n_all_layers, dz, Tprof_for_diffusion)

    return (ratefn(n, n_short_updated, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; # E FIX ATTEMPT
                   globvars...),#activellsp, activeslsp, inactivesp, Tn, Ti, Te, num_layers, H2Oi, HDOi, upper_lower_bdy_i),
            chemJmat(n, n_short, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, M, E;  # E FIX ATTEMPT
                     globvars...) )#activellsp, activeslsp, inactivesp, Tn, Ti, Te, num_layers, H2Oi, HDOi, upper_lower_bdy_i) ) # dt unused in chemjmat
end

function solve_sparse(A, b)
    LU = ilu(A, τ = 0.1) # get incomplete LU preconditioner
    x = bicgstabl(A, b, 2, Pl = LU, reltol=1e-25)
    return x
end

struct TooManyIterationsException <: Exception end

function next_timestep(nstart, params, t, dt; reltol=1e-2, abstol=1e-12, globvars...)
    #=
    Mike's notes with my formatting: 
    Moves to the next timestep using Newton's method on the linearized coupled transport and chemical reaction network.

    Input:
        nstart: concentrations before the timestep
        params: parameters needed by ratefn and chemical jacobian at each timestep
        t: current time (only used to print output)
        dt: timestep to be taken
        reltol, abstol: relative and absolute tolerances on solution
    Output: 
        nthis: concentrations after the timestep
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :neutral_species, :ion_species, 
                                   :inactivesp, :active_species, :activellsp, :activeslsp, 
                                   :n_inactive, 
                                   :Tn, :Ti, :Te, :Tp, 
                                   :plot_grid, :speciescolor, :speciesstyle, 
                                   :zmax, :alt, :num_layers, :non_bdy_layers, :dz,
                                   :absorber, :Jratelist, :crosssection, 
                                   :electron_val, 
                                   :Hs_dict, :bcdict, :H2Oi, :HDOi, :upper_lower_bdy_i, 
                                   :molmass, :n_alt_index, :polarizability, :n_all_layers, :Tprof_for_diffusion, :q])
    
    # absolute and relative tolerance on rate update
    f_abstol = 1e-2 #sqrt(abstol)
    f_reltol = reltol
    
    # assume concentrations after timestep are similar to those before
    nthis = deepcopy(nstart)
    dndt = 0.0*deepcopy(nthis)
    
    converged = false
    iter = 0
    while !converged
        # println("iter = $iter")
        if iter>20
            println("Recording last known good atmospheric state")
            record_atmospheric_state(t, nthis, GV.activellsp; globvars...)
            throw(TooManyIterationsException)
        end
        
        nold = deepcopy(nthis)
        
        # we perform an update using newton's method on
        #     f = nthis - nstart - dt*dndt
        
        # to do this we need the rates, dndt, and the jacobian of fval,
        #     d(fval)/d(nthis) = I - dt*chemJ
        dndt, chemJ = get_rates_and_jacobian(nthis, params, t; globvars...)

        # we want fval = 0, corresponding to a good update. chemical_jacobian.pdf eqn 3.21 
        # if it isn't 0, then there's a big discrepancy between the present value and the past value + rate of change, which is bad
        fval = nthis - nstart - dt*dndt  # 

        # construct the update matrix according to Jacobson Equation 12.77. h=dt, β=1
        identity = sparse(I, length(nthis), length(nthis))
        updatemat = identity - dt*chemJ  

        # now we can update using newton's method, chemical_jacobian.pdf equation 3.25
        # n_i+1 = n_i - f(n_i)/f'(n_i)
        # Here, updatemat = ∂f/∂n = I - dt*J  = f'(n_i),  fval = f(n^i)
        nthis = nthis - solve_sparse(updatemat, fval)

        # restrict values to be positive definite
        nthis[nthis .< 0.] .= 0.
        
        # check for convergence
        check_n_abserr = (nthis .<= abstol) .|| (nold .<= abstol) .|| (abs.(nthis - nold) .<= abstol)
        n_relerr = abs.(nthis-nold)./nold
        println("max(n_relerr) = $(max(n_relerr[.!check_n_abserr]...))")
        check_n_relerr = n_relerr .< reltol

        check_f_abserr = (abs.(fval) .<= f_abstol)
        if length(fval[.!check_n_abserr]) > 0
            println("max(fval) = $(max((abs.(fval[.!check_n_abserr]))...))")
        else
            println("max(f_relerr) = 0.0")
        end
        # println("max(fval) = $(max(fval...))\n") # OLD 

        # The following: fval / (nstart+ dt*dndt) = nthis/(nstart+ dt*dndt) - 1, and the division term should be ~=1, so you can check whether this val < rel tol. 
        f_relerr = abs.(fval./(nthis)) # f_relerr = abs.(fval./(nstart + dt*dndt)) # OLD
        if length(f_relerr[.!check_n_abserr .&& .!check_f_abserr]) > 0
            println("max(f_relerr) = $(max(f_relerr[.!check_n_abserr .&& .!check_f_abserr]...))\n")
        else
            println("max(f_relerr) = 0.0\n")
        end
        #= OLD
        if length(f_relerr[.!check_f_abserr]) > 0
            println("max(f_relerr) = $(max(f_relerr[.!check_f_abserr]...))\n")
        else
            println("max(f_relerr) = 0.0\n")
        end 
        =#
     
        check_f_relerr = f_relerr .< f_reltol
        
        #= old
        # Suggest: Remove [.!check_n_abserr] and see what happens, put it back in if breaks.
        # Suggest: && all(check_f_relerr .|| check_f_abserr) 
        # converged = all(check_n_relerr .|| check_n_abserr) && all(check_f_abserr#=[.!check_n_abserr]=# .|| check_f_relerr) # **
        old =# 
        converged = (all(check_n_relerr .|| check_n_abserr)
                     && all(check_f_abserr[.!check_n_abserr] .|| check_f_relerr[.!check_n_abserr]))
        
        iter += 1
    end

    return nthis
end

function update(n_current::Dict{Symbol, Array{ftype_ncur, 1}}, t, dt; abstol=1e-12, reltol=1e-2, globvars...)
    # update n_current using the coupled reaction network, moving to
    # the next timestep

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :neutral_species, :ion_species, 
                                   :inactivesp, :active_species, :activellsp, :activeslsp, 
                                   :n_inactive, 
                                   :Tn, :Ti, :Te, :Tp, 
                                   :plot_grid, :speciescolor, :speciesstyle, 
                                   :zmax, :alt, :num_layers, :non_bdy_layers, :dz,
                                   :absorber, :Jratelist, :crosssection, 
                                   :electron_val, 
                                   :Hs_dict, :bcdict, :H2Oi, :HDOi, :upper_lower_bdy_i, 
                                   :molmass, :n_alt_index, :polarizability, :n_all_layers, :Tprof_for_diffusion, 
                                   :Dcoef_arr_template, :q])
        
    # println("n_current[:O] = $(n_current[:O])")

    M = sum([n_current[sp] for sp in GV.all_species])
    if electron_val=="constant" # E FIX ATTEMPT 
        E = [1e5 for i in GV.non_bdy_layers]
    elseif electron_val=="quasineutral"
        E = sum([n_current[sp] for sp in GV.ion_species])
    else
        throw("Unhandled electron profile specification: $(GV.electron_val)")
    end

    # global params for simulation
    params = [
              #GV.inactive,
              #GV.inactive_species,
              ##GV.active_species,
              #GV.active_longlived, GV.active_shortlived,
              # GV.Tn, GV.Ti, GV.Te, GV.Tp,
              GV.Dcoef_arr_template, M, E] # E FIX ATTEMPT

    # get current long-lived species concentrations
    nstart = flatten_atm(n_current, GV.activellsp; GV.num_layers)

    # update to next timestep
    nend = next_timestep(nstart, params, t, dt; abstol=abstol, reltol=reltol, globvars...)
    #    println("max(nend-nstart) = $(max((nend-nstart)...))")

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, GV.activeslsp; GV.num_layers)  

    n_current = compile_ncur_all(nend, n_short, GV.n_inactive; GV.activellsp, GV.activeslsp, GV.inactivesp, GV.num_layers)
    # ensure Jrates are included in n_current
    update_Jrates!(n_current; GV.Jratelist, GV.crosssection, GV.num_layers, GV.absorber, GV.dz)

    return n_current
end

function converge(n_current::Dict{Symbol, Array{ftype_ncur, 1}}, log_t_start, log_t_end;#, n_steps=1000;
                  abstol=1e-12, reltol=1e-2, globvars...)
    #= 
    update in logarithmiclly spaced timesteps until convergence, returning converged atmosphere 
    =#
    
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :neutral_species, :ion_species, 
                                    :inactivesp, :active_species, :activellsp, :activeslsp, 
                                    :n_inactive, 
                                    :Tn, :Ti, :Te, :Tp, 
                                    :plot_grid, :speciescolor, :speciesstyle, 
                                    :zmax, :alt, :num_layers, :non_bdy_layers, :dz,
                                    :absorber, :Jratelist, :crosssection, 
                                    :electron_val, :gearmode, :n_steps, 
                                    :Hs_dict, :bcdict, :H2Oi, :HDOi, :upper_lower_bdy_i, 
                                    :molmass, :n_alt_index, :polarizability, :n_all_layers, :Tprof_for_diffusion, 
                                    :Dcoef_arr_template, :q])
    
    if GV.gearmode=="static"
        println("Using static timesteps")
        # DUMB but perhaps fast:
        log_time_steps = 10. .^(range(log_t_start, stop=log_t_end, length=GV.n_steps))

        total_time = 0.0
        for dt in log_time_steps
            total_time += dt
            println("t  = $total_time")
            println("dt = $dt")
            n_old = deepcopy(n_current)
            n_current = update(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
            println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in all_species]...))\n\n")
        end
    elseif GV.gearmode=="dynamic"
        println("Using dynamically adjusting timesteps")
        # FANCY: Modifies timesteps when it gets stuck due to large timesteps being > rel error.
        dt = 10.0^log_t_start
        
        total_time = 0.0
        goodsteps = 0
        goodstep_limit = 1
        goodfactor = 1.1
        failfactor = 10.0

        write_to_log(logfile, ["DYNAMIC STEP PARAMETERS", "Timesteps increase on success by: $((goodfactor-1)*100)%", "On fail, timesteps divided by: $(failfactor)"])
        while dt < 10.0^log_t_end
            total_time += dt
            println("t  = $total_time")
            println("dt = $dt")
            n_old = deepcopy(n_current)
            try
                n_current = update(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
                goodsteps += 1
                if goodsteps >= goodstep_limit
                    dt *= goodfactor
                    goodsteps = 0
                end
            catch e
                if e == TooManyIterationsException
                    total_time -= dt
                    dt = dt/failfactor
                    goodsteps = 0
                else
                    rethrow(e)
                end
            end
            println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in all_species]...))\n\n")
        end
    end
    
    return n_current
end

# **************************************************************************** #
#                                                                              #
#                          INITIAL MAIN MODEL SETUP                            #
#                                                                              #
# **************************************************************************** #

#                              Files and folders                                #
#===============================================================================#

create_folder(sim_folder_name, results_dir)
const logfile = results_dir*sim_folder_name*"/simulation_params_"*sim_folder_name*".txt"
const chemrxns_logfile = results_dir*sim_folder_name*"/included_chem_rxns.txt"

#                          Load reaction network                                #
#===============================================================================#
reaction_network, n_inds, i_inds = load_reaction_network(reaction_network_spreadsheet; 
                                                        ions_on=ions_included, get_inds=true, Jratelist, absorber, photolysis_products, all_species)
reactionnetwork_source = "spreadsheet"
write_to_log(chemrxns_logfile, "SPREADSHEET INDICES OF NEUTRAL REACTIONS", mode="w")
write_to_log(chemrxns_logfile, n_inds, mode="a")
if ions_included
    write_to_log(chemrxns_logfile, "SPREADSHEET INDICES OF ION REACTIONS", mode="a")
    write_to_log(chemrxns_logfile, i_inds, mode="a")
end

#                          Change the vertical extent                           #
#===============================================================================#
if make_new_alt_grid==true
    throw("The code for extending the altitude grid needs to be redone.")
    # const alt = convert(Array, (0:2e5:200e5))
    # n_current = get_ncurrent(initial_atm_file)

    # new_zmax = parse(Int64, input("Enter the new top of the atmosphere in km: "))
    # extra_entries = Int64((new_zmax - (zmax / 1e5))/(dz/1e5))

    # # Extend the grid
    # for (k,v) in zip(keys(n_current), values(n_current))
    #    append!(v, fill(v[end], extra_entries))  # repeats the last value in the array for the upper atmo as an initial value.
    # end

    # const alt = convert(Array, (0:dz:new_zmax*1e5))

    # const max_alt = new_zmax*1e5
end

#                           Load starting atmosphere                            #
#===============================================================================#
n_current = get_ncurrent(initial_atm_file)


#                       Establish new species profiles                          #
#===============================================================================#

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

#           Define storage for species/Jrates not solved for actively           #
#===============================================================================#
# Short lived species, when defined, are solved for assuming photochemical equilibrium
# outside of the primary ODE solver. Inactive species never change during simulation.
# Jrates must be stored here because they have to be updated alongside evolution
# of the atmospheric densities--the solver doesn't handle their values currently.
const external_storage = Dict([j=>n_current[j] for j in union(short_lived_species, inactive_species, Jratelist)])
const inactive = flatten_atm(n_current, inactive_species; num_layers)


#                          Set up the water profile                             #
#===============================================================================#

println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up the water profile...")

setup_water_profile!(n_current; all_species, num_layers, non_bdy_layers, DH, alt, plot_grid, 
                                hygropause_alt, H2O_excess, HDO_excess, ealt=excess_peak_alt, H2Osat, water_mixing_ratio, n_alt_index, results_dir, sim_folder_name)

# Calculate precipitable microns, including boundary layers (assumed same as nearest bulk layer)
H2Oprum = precip_microns(:H2O, [n_current[:H2O][1]; n_current[:H2O]; n_current[:H2O][end]]; molmass)#, all_species, dz, alt, molmass)
HDOprum = precip_microns(:HDO, [n_current[:HDO][1]; n_current[:HDO]; n_current[:HDO][end]]; molmass)#, all_species, dz, alt, molmass)

# **************************************************************************** #
#                                                                              #
#                   SETUP CHEMISTRY AND TRANSPORT NETWORK                      #
#                                                                              #
# **************************************************************************** #

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

#                               Transport network                               #
#===============================================================================#

const upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
                    Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
                    for s in transport_species];

const downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
                      Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
                      for s in transport_species];

const local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in transport_species]
                                [Symbol("t"*string(s)*"_down") for s in transport_species]
                                [Symbol("t"*string(s)*"_above_down") for s in transport_species]
                                [Symbol("t"*string(s)*"_below_up") for s in transport_species]]...;];

const transportnet = [[upeqns...;]; [downeqns...;]];

# define names for all the species active in the coupled rates:
const active_longlived_above = [Symbol(string(s)*"_above") for s in active_longlived];
const active_longlived_below = [Symbol(string(s)*"_below") for s in active_longlived];

#                               Chemical jacobian                               #
#===============================================================================#
# Create symbolic expressions for the chemical jacobian at a local layer with influence from that same layer, 
# the one above, and the one below
const chemJ_local = chemical_jacobian(active_longlived, active_longlived; diff_wrt_e=ediff, diff_wrt_m=mdiff, ion_species, chem_species, transport_species, chemnet=reaction_network, transportnet);
const chemJ_above = chemical_jacobian(active_longlived, active_longlived_above; diff_wrt_e=ediff, diff_wrt_m=mdiff, ion_species, chem_species, transport_species, chemnet=reaction_network, transportnet);
const chemJ_below = chemical_jacobian(active_longlived, active_longlived_below; diff_wrt_e=ediff, diff_wrt_m=mdiff, ion_species, chem_species, transport_species, chemnet=reaction_network, transportnet);

#                     Photochemical equilibrium setup                           #
#===============================================================================#

# TODO: Make this not run if Gear is the problem type. Have to figure out which of these things we would need instead. 

function setup_photochemical_equilibrium()
    #=
    Output
        active_longlived_species_rates: 
        short_lived_density_eqn:
        shortlived_density_inputs:
        equilibrium_eqn_terms: 
    =#


    # ------------------ Long-lived species expression array ------------------------ #

    # An array to store the rate equations for active, long-lived species, which are 
    # solved for in the production and loss equation.
    # each row is for each species; each column is for chemical production, chemical loss, 
    # transport production, transport loss, in that order.
    active_longlived_species_rates = Array{Array{Expr}}(undef, length(active_longlived), 4)
    for (i, sp) in enumerate(active_longlived)
        active_longlived_species_rates[i, :] .= getrate(sp, sepvecs=true; chemnet=reaction_network, transportnet, chem_species, transport_species, )
    end

    # ------------------ Short-lived species expression array ----------------------- #

    # TODO: Try to set these matrices to constants


    # Similarly, this array stores expressions for the concentrations of
    # active, short-lived species, which are assumed to be in photochemical equilibrium. 
    # Each row corresponds to a different chemical species. Instead of rates, the expressions
    # in the columns are the solution to the equation P(n_s) - L(n_s) = 0.
    # If the loss term L is linear in n_s, the solution is n_s = P/Lcoef. (Lcoef=L with one n_s factored out).
    # If it's quadratic in n_s, the solution is a quadratic equation for which two solutions
    # are possible. 
    # The columns are as follows:
    #
    # P-L=0 is:           Column 1                Column 2 
    # Linear             :(P/Lcoef)                 :(0)
    # Quadratic   :((-b+sqrt(b^2-4ac)/2a)   :((-b-sqrt(b^2-4ac)/2a)
    # 
    short_lived_density_eqn = Array{Expr}(undef, length(short_lived_species), 2) 

    # However, sometimes both P and L will be < machine epsilon, causing calculation problems,
    # and more rarely, the determinant for the quadratic solution may be negative, 
    # giving a complex solution. The following array stores P, Lcoef, and the determinant
    # separately so we can easily do a check on them for these problems.
    # The columns are as follows:
    #
    # P-L=0 is:       Column 1            Column 2  
    # Linear            :(P)                :(L)          
    # Quadratic     :(b^2 - 4ac)            :(NaN)   
    # 
    shortlived_density_inputs = Array{Expr}(undef, length(short_lived_species), 2)


    # This array just stores the same P - L = 0, but in the format P - nLcoef = 0 for linear 
    # and an^2 + bn + c = 0 for quadratic. This allows us to check how good a job the 
    # densities solved for do in actually getting the system in equilibrium.
    # The reason it's possible for P - nLcoef != 0 or an^2 + bn + c != 0 is that 
    # we solve for the density of each species one-by-one rather than solving for all densities
    # simultaneously as if it were a vector system. This may eventually need to be changed. 
    # It may also not be possible to solve it as a vector system and satisfy the constraints... 
    # which would just be the entire original problem that we were trying to avoid (i.e. can't 
    # satisfy density equations for neutrals and ions at the same time on the same timescales) 
    # by assuming photochemical equilibrium. 
    equilibrium_eqn_terms = Array{Expr}(undef, length(short_lived_species), 1)


    for (i, sp) in enumerate(active_shortlived)
        # Removes from consideration any reaction where a species appears on both sides of the equation (as an observer)
        ret = rxns_where_species_is_observer(sp, reaction_network)
        if ret == nothing
            chemnet = reaction_network
        else
            chemnet = filter(x->!in(x, ret), reaction_network)        
        end
        
        # Get the species production rate and loss rate by chemistry. These are obtained as vectors of reaction vectors
        # in the form [[:R1, :R2], [:P1, :P2], :(rate)]
        chem_prod_rate = production_rate(sp, chemnet, return_peqn_unmapped=true)
        chem_loss_rate = loss_rate(sp, chemnet, return_leqn_unmapped=true)
        
        if linear_in_species_density(sp, chem_loss_rate)
            # factors out the n_s so we can do n_s = P/L and converts to a big expression
            Lcoef_val = make_net_change_expr(loss_coef(chem_loss_rate, sp)) 
            P_val = make_net_change_expr(chem_prod_rate) # convert production to a big expression

            # Fill in the solutions array
            short_lived_density_eqn[i, 1] = :($P_val / $Lcoef_val)
            short_lived_density_eqn[i, 2] = :(0+0) # no second solution for linear.

            # Fill in the array that stores the separate components so they can be checked
            shortlived_density_inputs[i, 1] = :($P_val)
            shortlived_density_inputs[i, 2] = :($Lcoef_val)

            # Fill in the array that stores the entire expression, to compare to zero
            equilibrium_eqn_terms[i, 1] = :($P_val - $sp*($Lcoef_val))
        else # if it's quadratic in the species in question
            # Get the quadratic coefficients A, B, C for P - L = A(n^2_s) + B(n_s) + C = 0
            println("Note: $(sp) is not linear in density")
            qc = construct_quadratic(sp, chem_prod_rate, chem_loss_rate)

            # Fill in the solutions array with the quadratic formula
            short_lived_density_eqn[i, 1] = :((-$(qc["B"]) + sqrt($(qc["B"])^2 - 4*$(qc["A"])*$(qc["C"])))/(2*$(qc["A"])) )
            short_lived_density_eqn[i, 2] = :((-$(qc["B"]) - sqrt($(qc["B"])^2 - 4*$(qc["A"])*$(qc["C"])))/(2*$(qc["A"])) )

            # Fill in the array that stores the separate components so they can be checked
            shortlived_density_inputs[i, 1] = :(sqrt($(qc["B"])^2 - 4*$(qc["A"])*$(qc["C"])))
            shortlived_density_inputs[i, 2] = :(NaN)

            # Populate the array that lets us check if the densities give us 0
            equilibrium_eqn_terms[i, 1] = :($(qc["A"])*(($sp)^2) + $(qc["B"])*($sp) + $(qc["C"]))
        end
    end

    return active_longlived_species_rates, short_lived_density_eqn, shortlived_density_inputs, equilibrium_eqn_terms
end

const active_longlived_species_rates, short_lived_density_eqn, shortlived_density_inputs, equilibrium_eqn_terms = setup_photochemical_equilibrium()

# **************************************************************************** #
#                                                                              #
#                               METAPROGRAMMING                                #
#                                                                              #
# **************************************************************************** #

#          Arguments and expressions for metaprogramming functions              #
#===============================================================================#
const ratefn_arglist = [active_longlived; active_longlived_above; active_longlived_below; active_shortlived; inactive_species; Jratelist; 
                        :Tn; :Ti; :Te; :M; :E; local_transport_rates]; # E FIX ATTEMPT
const ratefn_arglist_typed = [:($s::ftype_chem) for s in ratefn_arglist];
const set_concentration_arglist = [active_shortlived; active_longlived; inactive_species; Jratelist; :Tn; :Ti; :Te; :M; :E]; # E FIX ATTEMPT
const set_concentration_arglist_typed = [:($s::ftype_chem) for s in set_concentration_arglist];

# Expressions for more accurate values of M (total third-bodies) and E (total electrons)
# const Mexpr = Expr(:call, :+, all_species...) 
# const Eexpr = Expr(:call, :+, ion_species...) # E FIX ATTEMPT

#                           Metaprogramming functions                           #
#===============================================================================#
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
        # is used. 
        # M = $Mexpr
        # E = $Eexpr # E FIX ATTEMPT

        # stack overflow - answer (Cite)
        # create a result array for evaluating the production and loss expressions
        result = map(_ -> ftype_chem[], $active_longlived_species_rates) 

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

        rates_local = net_chem_change .+ net_trans_change
        
        return convert(Array{ftype_ncur}, rates_local)

    end
end

@eval begin
    function check_zero_distance($(set_concentration_arglist_typed...))

        # M = $Mexpr
        # E = $Eexpr
        
        # Make a result array in which to store evaluated expressions 
        result = map(_ -> 0., $equilibrium_eqn_terms)

        # then we evaluate the expressions in the unrolled loop as below 
        $( (quote
                result[$i] = $(expr)
            end 
            for (i, expr) in enumerate(equilibrium_eqn_terms)                                                                                                                 
            )...  
         )
        
        # Now return the distance from zero--this is a single value. 
        return sqrt(sum(result .^ 2))
    end
end

@eval begin
    function set_concentrations_local($(set_concentration_arglist_typed...))
        #=
        Calculates the best possible solution for the new density of each short-lived species.
        =#

        # M = $Mexpr
        # E = $Eexpr  # E FIX ATTEMPT
         
        # Make a result array in which to store evaluated expressions 
        density_result = map(_ -> 0., $short_lived_density_eqn)  # will hold evaluated expressions.
        evaluated_inputs = map(_ -> 0., $shortlived_density_inputs) 
        
        # Evaluate the inputs 
        $( (quote 
               evaluated_inputs[$i] = $(expr)
            end
            for (i, expr) in enumerate(shortlived_density_inputs)
            )...
        )
        

        # 1) Check if P and Lcoef are both < machine epsilon (eps(Float64)). If yes, then we set col 1 in short_lived_density_eqn to 0.
        # values below ε will cause math problems.
        c1 = evaluated_inputs[:, 1] .< eps()
        c2 = evaluated_inputs[:, 2] .< eps()
        density_result[c1 .& c2, 1] .= 0.

        # 2) Check that the determinant b^2-4ac of a quadratic equation is positive. If negative, then we throw an error
        #    because I don't know what to do in that situation.
        quad_rows = findall(x->isnan(x), evaluated_inputs[:, 2])
        if any(x->x<0, evaluated_inputs[quad_rows, 1])
            throw("There is a complex solution for n in some species but I have not actually handled this error in any sort of way")
        end
        
        # Evaluate the new densities n in photochemical equilibrium
        $( (quote
                density_result[$i] = $(expr)
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
        best_solutions = choose_solutions(density_result, prev_densities)

        return convert(Vector{ftype_ncur}, best_solutions)
    end
end

@eval begin
    function chemJmat_local($(ratefn_arglist_typed...))
        #=
        Generates a matrix of I, J, and V values for a sparse matrix, for the 
        local layer, for influences from the layer above, and influences 
        from the layer below. 
        =#

        # M = $Mexpr  # E FIX ATTEMPT
        # E = $Eexpr

        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = convert(Array{ftype_ncur}, $(Expr(:vcat, chemJ_local[3]...)))

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = convert(Array{ftype_ncur}, $(Expr(:vcat, chemJ_above[3]...)))

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = convert(Array{ftype_ncur}, $(Expr(:vcat, chemJ_below[3]...)))

        # return the actual values of I, J, and V (indices and numerical value):
        return ((localchemJi, localchemJj, localchemJval),
                (abovechemJi, abovechemJj, abovechemJval),
                (belowchemJi, belowchemJj, belowchemJval))
    end
end

# **************************************************************************** #
#                                                                              #
#                        PHOTOCHEMICAL CROSS SECTIONS                          #
#                                                                              #
# **************************************************************************** #
println("$(Dates.format(now(), "(HH:MM:SS)")) Populating cross section dictionary...")

const crosssection = populate_xsect_dict(photochem_data_files; ion_xsects=ions_included, Tn=Tn_arr, n_all_layers, Jratelist)

# **************************************************************************** #
#                                                                              #
#                                SOLAR INPUT                                   #
#                                                                              #
# **************************************************************************** #
solarflux = readdlm(code_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]
solarflux[:,2] = solarflux[:,2] * cosd(SZA)  # Adjust the flux according to specified SZA

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

# **************************************************************************** #
#                                                                              #
#                                 LOGGING                                      #
#                                                                              #
# **************************************************************************** #

println("$(Dates.format(now(), "(HH:MM:SS)")) Creating the simulation log file...")

# Boundary condition write out messages
bc_log_msgs = ["\nBOUNDARY CONDITIONS",
               "n: number density at surface, f: flux at top, v: velocity at top"]
for k2 in keys(speciesbclist)
    push!(bc_log_msgs, "$(string(k2)): $(join([join(speciesbclist[k2][1, :], "="), join(speciesbclist[k2][2, :], "=")], ", "))") 
end

# Crosssection log messages
xsect_log_msgs = ["\nCROSS SECTIONS"]
for k in keys(photochem_data_files)  # cross sections
    push!(xsect_log_msgs, "$(k):")
    for k2 in keys(photochem_data_files[k])
        push!(xsect_log_msgs, "    $(k2): $(photochem_data_files[k][k2])")
    end
end

write_me_out = ["Description: $(optional_logging_note)", 
                "Primary experiment type: $(simtype)", 
                "Ions included: $(ions_included)",
                "Surface temperature T(z=0 km) (T_surf)=$(controltemps[1]) K", 
                "Mesosphere temperature (T(z~=100 km)) (T_meso)=$(controltemps[2]) K", 
                "Exosphere temperature (T(z>250 km) (T_exo)=$(controltemps[3]) K", 
                "Initial water mixing ratio at surface: $(water_mixing_ratio)",
                "DH=$(DH)", 
                "SVP fixed: $(fix_SVP)", 
                "Non-zero initial profiles used: $(use_nonzero_initial_profiles)", 
                "Mean temperatures: T_surf=$(meanTs) K, T_meso=$(meanTm) K, T_exo=$(meanTe) K", 
                "Include d(dn_s/dt)/dM in chemical jacobian: $(mdiff)", 
                "Include d(dn_s/dt)/dE in chemical jacobian: $(ediff)", 
                "Electron profile setting: $(electron_val)",
                
                # 
                "\nWATER PROFILE INFORMATION",
                "Total H2O (pr μm): $(H2Oprum)",
                "Total HDO (pr μm): $(HDOprum)",
                "Total H2O+HDO: $(H2Oprum + HDOprum)",
                # 
                "\nSOURCE FILES",
                "Initial atmosphere state: $(initial_atm_file)",
                "Final atmosphere state: $(final_atm_file)",
                "Solar cycle status: $(solarfile)",
                "Reaction network source: $(reactionnetwork_source)",
                # 
                "\nSPECIES LISTS",
                "All species: $(join([string(i) for i in all_species], ", "))", 
                "No-chem species (density unaffected by chemistry): $(join([string(i) for i in no_chem_species], ", "))", 
                "Chem species: $(join([string(i) for i in chem_species], ", "))", 
                "No-transport species (density unaffected by transport): $(join([string(i) for i in no_transport_species], ", "))", 
                "Transport species: $(join([string(i) for i in transport_species], ", "))", 
                "Active long-lived species: $(join(sort([string(i) for i in active_longlived]), ", "))", 
                "Active short-lived species: $(join(sort([string(i) for i in active_shortlived]), ", "))",
                #
                bc_log_msgs..., 
                #
                xsect_log_msgs...,
                #
                ] 
                
write_to_log(logfile, write_me_out, mode="w")

# **************************************************************************** #
#                                                                              #
#                       FINAL SETUP TO CONVERGE!                               #
#                                                                              #
# **************************************************************************** #

# Uncomment this line if you'd like to add an extra parcel to some species. You must specify the species.
# n_current[:D] = map(x->1e5*exp(-((x-184)/20)^2), non_bdy_layers/1e5) + n_current[:D]

# write initial atmospheric state ==============================================
write_atmosphere(n_current, results_dir*sim_folder_name*"/initial_state.h5"; alt, num_layers)

# Plot initial temperature and water profiles ==================================
plot_temp_prof(Tn_arr; savepath=results_dir*sim_folder_name, Tprof_2=Ti_arr, Tprof_3=Te_arr, alt)

# Absolute tolerance
if problem_type == "Gear"
    const abs_tol = 1e-10#1e-12 # absolute tolerance in ppm, used by Gear solver # NOTE: I think this is actually #/cm³ not ppm, because n_i+1 - n_i is compared against it.--Eryn
    const abs_tol_for_plot = fill(abs_tol, length(n_tot(n_current; all_species)))
else
    # absolute tolerance relative to total atmosphere density, used by DifferentialEquations.jl solvers
    const abs_tol = 1e-12 .* [[n_tot(n_current, a; n_alt_index, all_species) for sp in active_longlived, a in non_bdy_layers]...] 
    const abs_tol_for_plot = 1e-12 .* n_tot(n_current; n_alt_index, all_species) # calculates 1 ppt of the total density at each altitude.
end
    
# Plot initial atmosphere condition  ===========================================
println("$(Dates.format(now(), "(HH:MM:SS)")) Plotting the initial condition")
plot_atm(n_current, results_dir*sim_folder_name*"/initial_atmosphere.png", abs_tol_for_plot, t="initial state"; 
         neutral_species, ion_species, plot_grid, speciescolor, speciesstyle, zmax) #=plot_grid, speciescolor, speciesstyle, zmax,=#

# Create a list to keep track of stiffness ratio ===============================
# const stiffness = []

# Simulation time range and when to save a snapshot 
const mindt = dt_min_and_max[converge_which][1]
const maxdt = dt_min_and_max[converge_which][2]
times_to_save = [10.0^t for t in mindt:maxdt]
timestorage = times_to_save[1] / 10  # This variable lets us periodically print the timestep being worked on for tracking purposes
plotnum = 1 # To order the plots for each timestep correctly in the folder so it's easy to page through them.

# Record setup time
t5 = time()
write_to_log(logfile, "\n$(Dates.format(now(), "(HH:MM:SS)")) Setup time $(format_sec_or_min(t5-t4))", mode="a")

# **************************************************************************** #
#                                                                              #
#             First compile and call of Jacobian when using Double64           #
#                                                                              #
# **************************************************************************** #

# This takes a long time (~30 min) when running with Double64, move the first call outside evolve_atmosphere
if ftype_ncur==Double64
    println("$(Dates.format(now(), "(HH:MM:SS)")) Compiling and calling the chemical jacobian outside evolve_atmosphere (this will take ~45 min)...")
    write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Started first chemical jacobian compile")

    # Set up the initial state and check for any problems 
    M = sum([n_current[sp] for sp in all_species])
    if electron_val=="constant" # E FIX ATTEMPT 
        E = [1e5 for i in non_bdy_layers]
    elseif electron_val=="quasineutral"
        E = sum([n_current[sp] for sp in ion_species])
    else
        throw("Unhandled electron profile specification: $(elecval)")
    end

    nstart = flatten_atm(n_current, active_longlived; num_layers)
    find_nonfinites(nstart, collec_name="nstart")

    # Set up parameters    
    Dcoef_arr_template = zeros(size(Tn_arr))  # For making diffusion coefficient calculation go faster 
    params = [#inactive, inactive_species, active_species, active_longlived, active_shortlived, Tn_arr, Ti_arr, Te_arr, Tplasma_arr, 
              Dcoef_arr_template, M, E] # E FIX ATTEMPT
    params_exjac = deepcopy(params)  # I think this is so the Dcoef doesn't get filled in with the wrong info?

    # Call the expensive rate function
    t_before_jac = time()
    dndt, example_jacobian = get_rates_and_jacobian(nstart, params_exjac, 0.0)
    t_after_jac = time()

    write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Finished first chemical jacobian compile")
    println("$(Dates.format(now(), "(HH:MM:SS)")) ...finished.\nFirst jacobian compile+call took $(format_sec_or_min(t_after_jac-t_before_jac))\n")

    find_nonfinites(dndt, collec_name="example_rates")
    find_nonfinites(example_jacobian, collec_name="example_jacobian")    
end

println("Time to beginning convergence is $(format_sec_or_min(time()-t1))\n\n")
# Usually I split everthing below this out into a separate file and
# load the above into an interactive session -Mike

# **************************************************************************** #
#                                                                              #
#                         CONVERGE THE ATMOSPHERE                              #
#                                                                              #
# **************************************************************************** #

ti = time()
println("$(Dates.format(now(), "(HH:MM:SS)")) Beginning convergence")

Dcoef_arr_template = zeros(size(Tn_arr))
atm_soln = evolve_atmosphere(n_current, mindt, maxdt; t_to_save=times_to_save, abstol=abs_tol, reltol=rel_tol, 
                             # glob vars from here.  
                             all_species, neutral_species, ion_species, transport_species, # species lists 
                             inactivesp=inactive_species, active_species, activellsp=active_longlived, activeslsp=active_shortlived,  # species lists 
                             n_inactive=inactive, # inactive densities 
                             Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, Tprof_for_diffusion, # temperatures
                             alt, zmax, num_layers, non_bdy_layers, dz, n_all_layers, plot_grid, #anything to do with the alt grid
                             speciescolor, speciesstyle, # stuff to plot nicely 
                             absorber, Jratelist, crosssection, # Photolysis stuff
                             electron_val, gearmode=gear_timestep_type, n_steps, # simulation parameters that need to get used or logged
                             bcdict=speciesbclist, H2Oi, HDOi, upper_lower_bdy_i, # Things concerning transport 
                             Hs_dict, molmass, n_alt_index, polarizability, Dcoef_arr_template, q)# Basic species characteristics


tf = time() 

write_to_log(logfile, "Simulation active convergence runtime $(format_sec_or_min(tf-ti))", mode="a")

println("$(Dates.format(now(), "(HH:MM:SS)")) Simulation active convergence runtime $((tf-ti)/60) minutes")

# **************************************************************************** #
#                                                                              #
#                      WRITE OUT THE SIMULATION RESULTS                        #
#                                                                              #
# **************************************************************************** #

if problem_type == "SS"
    # Update short-lived species one more time
    n_short = flatten_atm(external_storage, active_shortlived; num_layers)
    Jrates = deepcopy(Float64[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:num_layers])
    set_concentrations!(external_storage, atm_soln.u, n_short, inactive, 
                        active_longlived, active_shortlived, inactive_species, Jrates, Tn_arr, Ti_arr, Te_arr)
    nc_all = merge(external_storage, unflatten_atm(atm_soln.u, active_longlived; num_layers))

    # Make final atmosphere plot
    plot_atm(nc_all, [neutral_species, ion_species], results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state", plot_grid, 
                      speciescolor, speciesstyle, zmax, abs_tol_for_plot)

    # Write out final atmosphere
    write_atmosphere(nc_all, results_dir*sim_folder_name*"/"*final_atm_file; alt, num_layers)   # write out the final state to a specially named file
elseif problem_type == "ODE"

    L = length(atm_soln.u)
    i = 1
    # Write all states to individual files.
    for (timestep, atm_state) in zip(atm_soln.t, atm_soln.u)

        if i != L 
            # NOTE: This will eventually result in the final Jrates being written out at every timestep.
            # Currently there's no workaround for this and you just have to remember NOT TO TRUST Jrates
            # at any timestep except the very last. 
            # TODO: Fix this so we just don't write Jrates in these iterations...
            local nc_all = merge(external_storage, unflatten_atm(atm_state, active_longlived; num_layers))
            write_atmosphere(nc_all, results_dir*sim_folder_name*"/atm_state_t_$(timestep).h5"; alt, num_layers) 
        elseif i == L
            # Update short-lived species one more time
            local n_short = flatten_atm(external_storage, active_shortlived; num_layers)
            local Jrates = deepcopy(Float64[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:num_layers])
            set_concentrations!(external_storage, atm_state, n_short, inactive, Jrates; activellsp=active_longlived, activeslsp=active_shortlived, 
                               inactivesp=inactive_species, Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, num_layers)
            local nc_all = merge(external_storage, unflatten_atm(atm_state, active_longlived; num_layers))

            # Make final atmosphere plot
            plot_atm(nc_all, results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state", abs_tol_for_plot; neutral_species, ion_species, 
                     plot_grid, speciescolor, speciesstyle, zmax)

            # Write out final atmosphere
            write_atmosphere(nc_all, results_dir*sim_folder_name*"/"*final_atm_file; alt, num_layers)   # write out the final state to a specially named file
        end
        global i += 1 
    end
elseif problem_type == "Gear"
    nc_all = atm_soln # this solver returns a standard dictionary rather than a matrix like the other solvers

    plot_atm(nc_all, results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state", abs_tol_for_plot; 
             neutral_species, ion_species, plot_grid, speciescolor, speciesstyle, zmax)

    # Write out final atmosphere
    write_atmosphere(nc_all, results_dir*sim_folder_name*"/"*final_atm_file; alt, num_layers)   # write out the final state to a specially named file
else
    throw("Invalid problem_type")
end 

t7 = time()

println("Wrote all states to files")
write_to_log(logfile, "Simulation total runtime $(format_sec_or_min(t7-t1))", mode="a")

# write_to_log(logfile, "Stiffness ratio evolution:")
# write_to_log(logfile, "$(stiffness)")
# write_to_log(logfile, "\n\n")
# write_to_log(logfile, "Total runtime $(format_sec_or_min(t9-t1))")
