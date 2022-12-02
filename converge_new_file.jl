###############################################################################
# converge_new_file.jl
# TYPE: (1) Model files - required
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere.
# 
# Eryn Cangi
# Created 2018
# Last edited: October 2022
# Currently tested for Julia: 1.7.1
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
# using Logging: global_logger
# using TerminalLoggers: TerminalLogger
# global_logger(TerminalLogger())

using Revise
using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using DelimitedFiles
using SparseArrays
using LinearAlgebra
photochemistry_source_dir = "$(@__DIR__)/Photochemistry/src/"
push!(LOAD_PATH, photochemistry_source_dir)
using Photochemistry  # custom module
using DifferentialEquations
using IterativeSolvers
using IncompleteLU
using GeneralizedGenerated: mk_function
using Dates
using DataFrames
using XLSX

t2 = time()
println("$(Dates.format(now(), "(HH:MM:SS)")) Loaded modules in $(format_sec_or_min(t2-t1))")

# **************************************************************************** #
#                                                                              #
#                              LOAD PARAMETERS                                 #
#                                                                              #
# **************************************************************************** #

include("CONSTANTS.jl")
include("CUSTOMIZATIONS.jl")

paramfile = get_paramfile(code_dir; use_default=true) # Set use_default to false to select from available files
include(paramfile)

t3 = time()
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

function evolve_atmosphere(atm_init::Dict{Symbol, Array{ftype_ncur, 1}}, log_t_start, log_t_end; t_to_save=[], abstol=1e-12, reltol=1e-6, globvars...)
    #=
    Sets up the initial conditions for the simulation and calls the ODE solver. 

    Input:
        atm_init: dictionary of species densities by altitude for the current state.
        log_t_start: Starting time will be 10^log_t_start.
        log_t_end: Similar, but end time.
    optional input:
        t_to_save: timesteps at which to save a snapshot of the atmosphere; currently not working yet
        abstol: absolute tolerance for simulation. Should usually be passed in.
        reltol: relative tolerance
    Output:
        sol: If using DifferentialEquations.jl, this is a solver object, the solution of the atmospheric system at the end of the simulation time.
             If using the Gear solver, this is a normal atmospheric state dictionary.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, 
                                   :collision_xsect, :crosssection, :Dcoef_arr_template, :dt_decr_factor, :dt_incr_factor, :dz, 
                                   :e_profile_type, :error_checking_scheme, :timestep_type, :H2Oi, :HDOi, 
                                   :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs, 
                                   :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Hs_dict, 
                                   :inactive_species, :ion_species, :Jratelist, :molmass, :n_all_layers, :n_alt_index, :n_inactive, :n_steps, 
                                   :neutral_species, :non_bdy_layers, :num_layers, :plot_grid, :polarizability, :q, :reaction_network, :solarflux,
                                   :speciesbclist, :speciescolor, :speciesstyle, :Te, :Ti, :Tn, :Tp, :Tprof_for_diffusion, :transport_species, 
                                   :upper_lower_bdy_i, :zmax])
        

    println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up initial state")

    # Define timespan - mostly for logging but also for use in SS and ODE 
    tspan = (10.0^log_t_start, 10.0^log_t_end)
    
    # Set up the initial state and check for any problems 
    nstart = flatten_atm(atm_init, GV.active_longlived; GV.num_layers)
    find_nonfinites(nstart, collec_name="nstart")

    # Set up parameters
    M = n_tot(n_current; GV.all_species)
    E = electron_density(n_current; GV.e_profile_type, GV.non_bdy_layers, GV.ion_species)
    params_Gear = [GV.Dcoef_arr_template, M, E]
    params_J = [globvars, GV.Dcoef_arr_template, M, E] # kwargs can't be passed to the julia ODE solver functions 
    params_exjac = deepcopy(params_Gear)  # make sure not to have a pointer problem

    # Set up an example Jacobian and report any lack of sparsity
    println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up example jacobian")
    dndt, example_jacobian = get_rates_and_jacobian(nstart, params_exjac, 0.0; globvars...)
    find_nonfinites(example_jacobian, collec_name="example_jacobian")
    sparsity = round(length(example_jacobian.nzval)*100/(size(example_jacobian)[1] * size(example_jacobian)[2]), digits=2)
    if sparsity > 1
        println("Warning! Sparsity of the jacobian is rather high: $(sparsity)%")
    end

    # Now define the problem function to be solved
    if problem_type != "Gear"
        println("$(Dates.format(now(), "(HH:MM:SS)")) Defining ODE function")
        odefunc = ODEFunction(PnL_eqn, jac=jacobian_wrapper, jac_prototype=example_jacobian)
    end

    # Log solver options and run the solver
    push!(PARAMETERS_SOLVER, ("JAC_NONZERO", sparsity))
    push!(PARAMETERS_SOLVER, ("JAC_LEN", length(nstart)^2))
    push!(PARAMETERS_SOLVER, ("SYS_SIZE", length(nstart)))
    push!(PARAMETERS_GEN, ("T_I", tspan[1]))
    push!(PARAMETERS_GEN, ("T_F", tspan[2]))
    
    
    write_time_stuff = ["\nTIMING",
                        "Start time: $(Dates.format(now(), "yyyy-mm-dd at HH:MM:SS"))"]
    
    julia_DE_algorithm = "QNDF"
    if problem_type == "SS"
        push!(PARAMETERS_SOLVER, ("ALGORITHM", "SS: $(julia_DE_algorithm)"))
        push!(PARAMETERS_SOLVER, ("ABSTOL", abstol))
        push!(PARAMETERS_SOLVER, ("RELTOL", reltol))
        write_to_log(logfile, write_time_stuff, mode="a")
        probSS = SteadyStateProblem{isinplace}(odefunc, nstart, params_J)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Calling the solver")
        sol = solve(probSS, DynamicSS(QNDF(), abstol=abstol, reltol=reltol), saveat=t_to_save, progress=true, save_everystep=false, dt=10.0^log_t_start,
                    abstol=abstol, reltol=reltol, isoutofdomain=(u,p,t)->any(x->x<0,u))
    elseif problem_type == "ODE"
        push!(PARAMETERS_SOLVER, ("ALGORITHM", "ODE: $(julia_DE_algorithm)"))
        push!(PARAMETERS_SOLVER, ("ABSTOL", abstol))
        push!(PARAMETERS_SOLVER, ("RELTOL", reltol))
        write_to_log(logfile, write_time_stuff, mode="a")
        probODE = ODEProblem(odefunc, nstart, tspan, params_J)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Calling the solver")
        sol = solve(probODE, QNDF(), saveat=t_to_save, progress=true, save_everystep=false, dt=10.0^log_t_start,
                    abstol=abstol, reltol=reltol, isoutofdomain=(u,p,t)->any(x->x<0,u))
    elseif problem_type == "Gear"
        push!(PARAMETERS_SOLVER, ("ALGORITHM", problem_type))
        push!(PARAMETERS_SOLVER, ("TIMESTEP_MANAGEMENT", timestep_type))
        push!(PARAMETERS_SOLVER, ("N_STEPS", n_steps))
        push!(PARAMETERS_SOLVER, ("DT_INCR", dt_incr_factor))
        push!(PARAMETERS_SOLVER, ("DT_DECR", dt_decr_factor))
        push!(PARAMETERS_SOLVER, ("ERROR_SCHEME", error_checking_scheme))
        push!(PARAMETERS_SOLVER, ("ABSTOL", abstol))
        push!(PARAMETERS_SOLVER, ("RELTOL", reltol))
        write_to_log(logfile, write_time_stuff, mode="a")
        println("$(Dates.format(now(), "(HH:MM:SS)")) Calling the solver")
        sol, sim_time = converge(atm_init, log_t_start, log_t_end; abstol=abstol, reltol=reltol, globvars...) # this is where the action happens
    else
        throw("Invalid problem_type")
    end 

    return sol, sim_time
end

function record_atmospheric_state(t, n, actively_solved, E_prof; opt="", globvars...)
    #=
    Input:
        t: current simulation time
        n: current densities of active longlived species 
        actively_solved: species which are actively solved for. Required to compile the full atmospheric state.
        E_prof: Electron profile at the present time, needed to plot the whole atmosphere.
        Optional arguments:
            opt: Optional string to append to each filename
    Output:
        plot of the atmospheric densities at the present time
        .h5 file containing the present atmospheric state 
        Also prints a progress alert to the terminal and writes out the time to a logfile. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:hrshortcode, :neutral_species, :ion_species, :plot_grid, :rshortcode, :speciescolor, :speciesstyle, :zmax, :alt, :num_layers])

    # This is just to change how many decimal places to include depending if t >= 1.
    rounding_digits = t <= 1 ? Int64(ceil(abs(log10(t)))) : 0 

    progress_alert = "$(Dates.format(now(), "(HH:MM:SS)")) reached timestep $(t), total elapsed time=$(format_sec_or_min(time() - ti))"
    write_to_log(logfile, progress_alert)
    println(progress_alert)
    
    # write out the current atmospheric state to a file and plot it
    atm_snapshot = merge(external_storage, unflatten_atm(n, actively_solved; GV.num_layers))
    plot_atm(atm_snapshot, results_dir*sim_folder_name*"/atm_peek_$(plotnum)$(opt).png", abs_tol_for_plot, E_prof; t="$(round(t, digits=rounding_digits))", globvars...)
    write_atmosphere(atm_snapshot, results_dir*sim_folder_name*"/atm_state_$(lpad(plotnum,2,"0"))$(opt).h5"; t=round(t, digits=rounding_digits), globvars...)

    # Turn this on if you'd like to take a peek at the Jrates
    # for Jspc in values(GV.neutral_species)
    #     plot_Jrates(Jspc, atm_snapshot, results_dir*sim_folder_name*"/Jrateplots/"; filenameext="$(plotnum)", globvars...)               
    # end 

    global plotnum += 1
end

function chemJmat(n_active_longlived, n_active_shortlived, n_inactive, Jrates, tup, tdown, tlower, tupper, M, E;
                  check_eigen=false, globvars...)
    #=
    Collects coordinate tuples of (I, J, V) [row index, column index, value] for a sparse matrix
    representing the chemical jacobian of the atmospheric system. 

    Input:
        n_active_longlived: The atmospheric densities array, but flattened, in the form 
                            [n_CO(z=0), n_CO2(z=0)...n_N2Dpl(z=0), n_CO(z=2)...n_N2Dpl(z=250)], for active and chemically long-lived species.
                            I am NOT saying that CO is the first species in the order. Just describing how it goes.         
        n_active_shortlived: active shortlived species densities, necessary to to calculations for longlived species.
        n_inactive: A flattened array of the atmospheric densities of any inactive species, same format as nthis. Functionally constant.
        Jrates: Flattened array of Jrates, same format as nthis.
        tup, tdown: Transport coefficients
        tlower, tupper: Transport coefficients
        M: Total density by altitude for entire atmosphere
        E: Electron profile at the present time
    optional input:
        check_eigen: Will check the eigenvalues of the jacobian for non-real or real/positive values if true. Not currently used
    Output:
        sparse matrix representing the chemical jacobian 
    =#              

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:active_longlived, :active_shortlived, :H2Oi, :HDOi, :inactive_species, :num_layers, :Tn, :Ti, :Te, :upper_lower_bdy_i])

    nmat_llsp = reshape(n_active_longlived, (length(GV.active_longlived), GV.num_layers))
    nmat_slsp = reshape(n_active_shortlived, (length(GV.active_shortlived), GV.num_layers))
    nmat_inactive = reshape(n_inactive, (length(GV.inactive_species), GV.num_layers))
    
    # For storing the jacobian indices and values
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = ftype_chem[]

    # tc___ are the coordinate tuples containing (I, J, V) to be used to fill a sparse matrix.
    argvec = [nmat_llsp[:, 1];                        # active_longlived;
              nmat_llsp[:, 2];                        # active_longlived_above;
              fill(1.0, length(GV.active_longlived)); # active_longlived_below;
              nmat_slsp[:, 1];                        # active_shortlived;
              nmat_inactive[:,1];                     # inactive_species;
              Jrates[:,1];                            # Jratelist;
              GV.Tn[1]; GV.Ti[1]; GV.Te[1];           #:Tn; :Ti; :Te;
              M[1]; E[1];                             # total density and electrons,
              tup[:,1]; tlower[:,1];                  # local_transport_rates
              tdown[:,2]; tlower[:,2]]
    argvec = convert(Array{ftype_chem}, argvec)

    (tclocal, tcupper, tclower) = chemJmat_local(argvec...) 

    # add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])

    # and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(GV.active_longlived))
    append!(chemJval, tcupper[3])

    for ialt in 2:(GV.num_layers-1)
        argvec = [nmat_llsp[:, ialt];
                  nmat_llsp[:, ialt+1];
                  nmat_llsp[:, ialt-1];
                  nmat_slsp[:, ialt];
                  nmat_inactive[:, ialt];
                  Jrates[:, ialt];
                  GV.Tn[ialt]; GV.Ti[ialt]; GV.Te[ialt];
                  M[ialt]; E[ialt];
                  tup[:, ialt];
                  tdown[:, ialt];
                  tdown[:, ialt+1];
                  tup[:, ialt-1]]
        argvec = convert(Array{ftype_chem}, argvec)        

        (tclocal, tcupper, tclower) = chemJmat_local(argvec...)

        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(GV.active_longlived))
        append!(chemJj, tclocal[2].+(ialt-1)*length(GV.active_longlived))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(GV.active_longlived))
        append!(chemJj, tcupper[2].+(ialt  )*length(GV.active_longlived))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(GV.active_longlived))
        append!(chemJj, tclower[2].+(ialt-2)*length(GV.active_longlived))
        append!(chemJval, tclower[3])
    end

    argvec = [nmat_llsp[:,end];
              fill(1.0, length(GV.active_longlived));
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
    append!(chemJi, tclocal[1].+(GV.num_layers-1)*length(GV.active_longlived))
    append!(chemJj, tclocal[2].+(GV.num_layers-1)*length(GV.active_longlived))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(GV.num_layers-1)*length(GV.active_longlived))
    append!(chemJj, tclower[2].+(GV.num_layers-2)*length(GV.active_longlived))
    append!(chemJval, tclower[3])

    # fix water below whatever we set as upper/lower atmosphere boundary.
    # This only runs if water is designated as an active species; if it's in inactive_species, this won't run,
    # When it is active, this finds all the H2O and HDO indices for the lower atmosphere. 
    # It's like above where we add (ialt-1)*length(active_species), but this way it's outside the loop.
    if in(:H2O, GV.active_longlived) && in(:HDO, GV.active_longlived)
        H2Opositions = GV.H2Oi .+ length(GV.active_longlived)*collect(0:GV.upper_lower_bdy_i-1)
        HDOpositions = GV.HDOi .+ length(GV.active_longlived)*collect(0:GV.upper_lower_bdy_i-1)
        water_positions = sort(union(H2Opositions, HDOpositions))

        i_remove = findall(x->in(x, water_positions), chemJi)
        j_remove = findall(x->in(x, water_positions), chemJj) # these are removed because if a species is inert, a derivative with respect to it is a derivative of a constant 
        remove_these = sort(union(i_remove, j_remove)) # This makes a set, since it describes the locations where the H2O and HDO indices are.
                                                       # Kinda confusing since we're talking about indices of indices.
        chemJval[remove_these] .= 0 
    end

    # Uncomment the following to check the eigenvalues of the jacobian. Requires a global variable called stiffness.
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

    Arguments are basically the same as chemJmat.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:active_longlived, :active_shortlived, :H2Oi, :HDOi, :inactive_species, :num_layers, :Tn, :Ti, :Te, :upper_lower_bdy_i])

    nmat_llsp = reshape(n_active_longlived, (length(GV.active_longlived), GV.num_layers))
    nmat_slsp = reshape(n_active_shortlived, (length(GV.active_shortlived), GV.num_layers))
    nmat_inactive = reshape(n_inactive, (length(GV.inactive_species), GV.num_layers))
    returnrates = zeros(size(nmat_llsp))
    
    # fill the first altitude entry with information for all species
    argvec = [nmat_llsp[:,1];                         # densities for active_longlived;
              nmat_llsp[:,2];                         # active_longlived_above;
              fill(1.0, length(GV.active_longlived)); # active_longlived_below;
              nmat_slsp[:, 1];                        # active_shortlived;
              nmat_inactive[:,1];                     # inactive_species;
              Jrates[:,1];                            # Jratelist;
              GV.Tn[1]; GV.Ti[1]; GV.Te[1];           # :Tn; :Ti; :Te;
              M[1]; E[1];  
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
                  M[ialt]; E[ialt]; 
                  tup[:, ialt];
                  tdown[:, ialt];
                  tdown[:, ialt+1];
                  tup[:, ialt-1]]
        argvec = convert(Array{ftype_chem}, argvec)
        
        returnrates[:,ialt] .= ratefn_local(argvec...)
    end

    # fill in the last level of altitude
    argvec = [nmat_llsp[:, end];
              fill(1.0, length(GV.active_longlived));
              nmat_llsp[:, end-1];
              nmat_slsp[:, end]; # active_shortlived;
              nmat_inactive[:,end];
              Jrates[:,end];
              GV.Tn[end]; GV.Ti[end]; GV.Te[end];
              M[end]; E[end];
              tupper[:,1];
              tdown[:,end];
              tupper[:,2];
              tup[:,end-1]]
    argvec = convert(Array{ftype_chem}, argvec)
    returnrates[:,end] .= ratefn_local(argvec...)

    # NEW: Overwrite the entries for water in the lower atmosphere with 0s so that it will behave as fixed.
    # Only runs when water is in the active_species list. If neutrals are set to inactive, it will be taken care of already.
    if in(:H2O, GV.active_longlived) && in(:HDO, GV.active_longlived)
        returnrates[GV.H2Oi, 1:GV.upper_lower_bdy_i] .= 0
        returnrates[GV.HDOi, 1:GV.upper_lower_bdy_i] .= 0
    end

    return [returnrates...;]
end

function update_Jrates!(n_cur_densities::Dict{Symbol, Array{ftype_ncur, 1}}; globvars...)
    #=
    this function updates the photolysis rates stored in n_cur_densities to
    reflect the altitude distribution of absorbing species.

    Input:
        n_cur_densities: The present atmospheric state. This will be updated to include Jrates by this function.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:absorber, :dz, :crosssection, :Jratelist, :num_layers, :solarflux])

    # Initialize an array, length=number of active layers
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, GV.num_layers)
    for i in range(1, length=GV.num_layers)
        solarabs[i] = zeros(Float64, 2000)#2000)
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
        solarabs[ialt] = GV.solarflux[:,2] .* exp.(-solarabs[ialt])
    end

    # You can uncomment these to plot the extinction at each atmospheric level, but you have to feed it a specific Jrate
    # and also I can't promise that it's currently working.
    # plot_extinction(fill(fill(0.,size(solarflux, 1)), num_layers+2); xsect_info=[GV.crosssection[:JCO2toOplpCplpO], :JCO2toOplpCplpO], GV.zmax)
    # plot_extinction(fill(fill(0.,size(solarflux, 1)), num_layers+2); xsect_info=[GV.crosssection[:JCO2toCOplpOpl], :JCO2toCOplpOpl], GV.zmax)
    # plot_extinction(fill(fill(0.,size(solarflux, 1)), num_layers+2); xsect_info=[GV.crosssection[:JH2OtoHpOH], :JH2OtoHpOH], GV.zmax)

    # each species absorbs according to its cross section at each
    # altitude times the actinic flux.
    # BLAS.dot includes an integration (sum) across wavelengths, i.e:
    # (a·b) = aa + ab + ab + bb etc that kind of thing
    for j in GV.Jratelist
        n_cur_densities[j] = zeros(GV.num_layers)
        for ialt in [1:nalt;]
            n_cur_densities[j][ialt] = ftype_ncur(BLAS.dot(nlambda, solarabs[ialt], 1, GV.crosssection[j][ialt+1], 1))
        end
    end
end

# Julia DifferentialEquations.jl solver =======================================
function make_jacobian(n, p, t)
    #=
    Constructs the chemical jacobian in the normal way, including stuff to calculate parameters for chemJmat.
    ! MAY NOT BE UP TO DATE AND WORKING !
    =#

    # Unpack the parameters ---------------------------------------------------------------
    globvars, D_arr, M, E = p 
    
    GV = values(globvars)
    # @assert all(x->x in keys(globvars), [#:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, 
    #                                    # :bcdict, :collision_xsect, :crosssection, :dz, :H2Oi, :HDOi, :Hs_dict,  
    #                                    # :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs, 
    #                                    # :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs,
    #                                    # :inactive_species, :ion_species,  :Jratelist,
    #                                    # :molmass, :neutral_species, :non_bdy_layers, :num_layers, :n_all_layers, :n_alt_index, :n_inactive, 
    #                                    # :plot_grid, :polarizability, :q, :speciescolor, :speciesstyle, :Tn, :Ti, :Te, :Tp, :Tprof_for_diffusion, 
    #                                    # :transport_species, :upper_lower_bdy_i, :zmax
    #                                ])

    # get the concentrations of species assumed to be in photochemical equilibrium. 
    n_short = flatten_atm(external_storage, GV.active_shortlived; GV.num_layers)  # retrieve the shortlived species from their storage and flatten them

    # Update Jrates
    n_cur_all = compile_ncur_all(n, n_short, GV.n_inactive; GV.active_longlived, GV.active_shortlived, GV.inactive_species, GV.num_layers)

    update_Jrates!(n_cur_all; GV.Jratelist, GV.crosssection, GV.num_layers, GV.absorber, GV.dz, GV.solarflux)
    # copy all the Jrates into an external dictionary for storage
    for jr in GV.Jratelist                # time for this is ~0.000005 s
        global external_storage[jr] = n_cur_all[jr]
    end

    # Retrieve Jrates 
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in GV.Jratelist, ialt in 1:GV.num_layers])

    # and update the shortlived species with the new Jrates - assuming not needed to be done in this function
    # n_short_updated = set_concentrations!(external_storage, n, n_short, n_inactive, active_longlived, active_shortlived, inactive_species, Jrates, Tn, Ti, Te)

    tlower, tup, tdown, tupper = update_transport_coefficients(GV.transport_species, # Species for which to update coefficients so it's not a mistake to pass it twice.
                                                               n_cur_all, D_arr, M; calc_nonthermal=nontherm, globvars...)
                                                               # Tn, Tp, Hs_dict, bcdict=speciesbclist, 
                                                               # all_species, neutral_species, transport_species, molmass, n_alt_index, 
                                                               # polarizability, alt, num_layers, n_all_layers, dz, T_for_diff=Tprof_for_diffusion, q)
    
    return chemJmat(n, n_short, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; globvars...)
                    # active_longlived, active_shortlived, inactive_species, Tn, Ti, Te, num_layers, H2Oi, HDOi, upper_lower_bdy_i)
end

function jacobian_wrapper(J, n, p, t)
    #=
    A wrapper since the Jacobian is sometimes generated outside of the normal simulation run
    =#

    J .= make_jacobian(n, p, t)
end

function PnL_eqn(dndt, n, p, t)
    #=
    This is the primary function that is solved by the solver. For our purposes that means it calls ratefn().

    dndt is written as such because the output of this function is silently multiplied by the timestep within the solver's code.
    n is a vector of active species densities by altitude.
    ! MAY NOT BE UP TO DATE AND WORKING !
    =#

    # Unpack the parameters needed to call ratefn 
    # n_inactive, inactive_species, activesp, active_longlived, active_shortlived, Tn, Ti, Te, Tp, D_arr = p 
    globvars, D_arr, M, E = p 

    GV = values(globvars)
    # @assert all(x->x in keys(globvars), [#:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, 
    #                                # :bcdict, :collision_xsect, :crosssection, :dz, :H2Oi, :HDOi, :Hs_dict,  
    #                                # :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs, 
    #                                # :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs,
    #                                # :inactive_species, :ion_species,  :Jratelist,
    #                                # :molmass, :neutral_species, :non_bdy_layers, :num_layers, :n_all_layers, :n_alt_index, :n_inactive, 
    #                                # :plot_grid, :polarizability, :q, :speciescolor, :speciesstyle, :Tn, :Ti, :Te, :Tp, :Tprof_for_diffusion, 
    #                                # :transport_species, :upper_lower_bdy_i, :zmax
    #                                ])

    # Every time we increase dt by a factor of 10, print a progress report to the console,
    # plot the atmosphere, and save the atmosphere to an .h5 file. 
    if t > checkpoint
        record_atmospheric_state(t, n, GV.active_longlived, E; opt=GV.opt, globvars...)
        global t_i += 1
        global checkpoint = times_to_save[t_i]
    end

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, GV.active_shortlived; GV.num_layers)

    # Retrieve the Jrates
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in GV.Jratelist, ialt in 1:GV.num_layers])

    # set the concentrations of species assumed to be in photochemical equilibrium. 
    n_short_updated = set_concentrations!(external_storage, n, n_short, GV.n_inactive, Jrates, M, E; globvars...)
                                          #active_longlived, active_shortlived, inactive_species, Tn, Ti, Te, num_layers)
    # for i in 1:2
    #     n_short_updated = set_concentrations!(external_storage, n, n_short_updated, n_inactive, active_longlived, active_shortlived, inactive_species, Jrates, Tn, Ti, Te)
    # end

    # Get the updated transport coefficients, taking into account short-lived species update
    updated_ncur_all = compile_ncur_all(n, n_short_updated, GV.n_inactive; GV...)#active_longlived, active_shortlived, inactive_species, num_layers)
    tlower, tup, tdown, tupper = update_transport_coefficients(GV.transport_species, updated_ncur_all, D_arr, M; 
                                                               calc_nonthermal=nontherm, globvars...)

    dndt .= ratefn(n, n_short_updated, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; globvars...)
end

function set_concentrations!(external_storage, n_active_long, n_active_short, n_inactive, Jrates, M, E; globvars...) 
    #=
    at each altitude, sets the concentrations for short-lived species assumed to be in photochemical equilibrium
    and sends them back into the storage dictionary, external_storage

    Inputs:
        external_storage: dictionary storing densities for short-lived and inactive species, as well as Jrates.
        n_active_short, n_active_long, n_inactive: density of short-lived, long-lived, and inactive species
        active_shortlived, active_longlived, inactive_species: list of short- and long-lived species names
        Jrates: Jrates for each species, for a particular altitude
        M: Total atmospheric density
        E: electron density profile
    Outputs:
        Updates the contents of external_storage.

    Note: dist_zero array and calculations were used to determine how far the vector of density solutions are from "zero" i.e. a perfect solution.
    Those lines are probably out of date...
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:active_longlived, :active_shortlived, :inactive_species, :Tn, :Ti, :Te, :num_layers])
    
    # rows = species, columns = altitudes. 
    nmat_shortlived = reshape(n_active_short, (length(GV.active_shortlived), GV.num_layers))
    
    # auxiliary information that is needed. 
    nmat_longlived = reshape(n_active_long, (length(GV.active_longlived), GV.num_layers))
    nmat_inactive = reshape(n_inactive, (length(GV.inactive_species), GV.num_layers))
    
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
    for (s, ssp) in enumerate(GV.active_shortlived)
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
    1. Updates short lived species when DifferentialEquations.jl is used
    2. Updates and saves Jrates
    3. Updates short lived species again
    4. Calculates transport coefficients
    5. Collects rates (dn/dt) and the chemical jacobian.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, 
                                   :collision_xsect, :crosssection, :dz, :H2Oi, :HDOi, :Hs_dict,  
                                   :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs, 
                                   :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs,
                                   :inactive_species, :ion_species,  :Jratelist,
                                   :molmass, :neutral_species, :non_bdy_layers, :num_layers, :n_all_layers, :n_alt_index, :n_inactive, 
                                   :plot_grid, :polarizability, :q, :reaction_network, :solarflux, :speciesbclist, :speciescolor, :speciesstyle, :Tn, :Ti, :Te, :Tp, :Tprof_for_diffusion, 
                                   :transport_species, :upper_lower_bdy_i, :zmax])
    # Unpack the parameters ---------------------------------------------------------------
    D_arr, M, E = p 

    # Records the atmospheric state whenever the total time elapsed has progressed 10x. Does not land on perfect powers of 10.
    if t >= checkpoint  # (checkpoint - (sol_in_sec/2)) <= t <= (checkpoint + (sol_in_sec/2))  #
        record_atmospheric_state(t, n, GV.active_longlived, E; opt=GV.opt, globvars...)
        global t_i += 1

        if t_i <= length(times_to_save)
            global checkpoint = times_to_save[t_i]
        else
            global checkpoint = times_to_save[end] * 10 # This will update the checkpoint so the simulation doesn't stop at every dt in the last t. 
        end
        println("New checkpoint = $(checkpoint)")
    end

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, GV.active_shortlived; GV.num_layers)
    # Turn on following 3 lines for debugging
    # if nans_present(n_short)
    #     throw("n_short has nans")
    # end

    # Update Jrates
    n_cur_all = compile_ncur_all(n, n_short, GV.n_inactive; GV.active_longlived, GV.active_shortlived, GV.inactive_species, GV.num_layers)
    # Turn on the following 5 lines for debugging - if the active species may have nans.
    # for k in keys(n_cur_all)
    #     if nans_present(n_cur_all[k])
    #         throw("n_cur_all has nans")
    #     end
    # end

    update_Jrates!(n_cur_all; GV.Jratelist, GV.crosssection, GV.num_layers, GV.absorber, GV.dz, GV.solarflux)
    # copy all the Jrates into an external dictionary for storage
    for jr in GV.Jratelist                # time for this is ~0.000005 s
        global external_storage[jr] = n_cur_all[jr]
    end

    # Retrieve Jrates 
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in GV.Jratelist, ialt in 1:GV.num_layers])

    # set the concentrations of species assumed to be in photochemical equilibrium. 
    n_short_updated = set_concentrations!(external_storage, n, n_short, GV.n_inactive, Jrates, M, E; 
                                          GV.active_longlived, GV.active_shortlived, GV.inactive_species, GV.Tn, GV.Ti, GV.Te, GV.num_layers)
    # Debugging
    # if nans_present(n_short_updated)
    #     throw("n_short_updated has nans")
    # end

    # Reconstruct the dictionary that holds densities
    updated_ncur_all = compile_ncur_all(n, n_short_updated, GV.n_inactive; GV.active_longlived, GV.active_shortlived, GV.inactive_species, GV.num_layers)

    # Get the updated transport coefficients, taking into account short-lived species update
    tlower, tup, tdown, tupper = update_transport_coefficients(GV.transport_species, updated_ncur_all, D_arr, M; 
                                                               calc_nonthermal=nontherm, results_dir, sim_folder_name, 
                                                               Jratedict=Dict([j=>n_cur_all[j] for j in GV.Jratelist]), # Needed for nonthermal BCs
                                                               globvars...)

    # Turn on the following block for more nan debugging
    # if nans_present(tlower) || nans_present(tupper) || nans_present(tup) || nans_present(tdown)
    #     if nans_present(tlower) 
    #         println("tlower has nans: ", tlower)
    #     end
    #     if nans_present(tupper)
    #         println("tupper has nans: ", tupper)
    #     end
    #     if nans_present(tup)
    #         println("tup has nans: ", tup)
    #     end
    #     if nans_present(tdown)
    #         println("tdown has nans: ", tdown)
    #     end
    #     throw("transport coefficients have nans")
    # end

    return (ratefn(n, n_short_updated, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; globvars...),# E FIX ATTEMPT   
            chemJmat(n, n_short, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, M, E; globvars...) ) # E FIX ATTEMPT
end

function solve_sparse(A, b)
    #=
    A faster way to solve the Jacobian
    =#
    LU = ilu(A, τ = 0.1) # get incomplete LU preconditioner
    x = bicgstabl(A, b, 2, Pl = LU, reltol=1e-25)
    return x
end

struct TooManyIterationsException <: Exception end

function next_timestep(nstart, params, t, dt; reltol=1e-2, abstol=1e-12, verbose=false, globvars...)
    #=
    Moves to the next timestep using Newton's method on the linearized coupled transport and chemical reaction network.

    Input:
        nstart: concentrations before the timestep
        params: parameters needed by ratefn and chemical jacobian at each timestep. 3rd entry should be electron profile. 
        t: current time (only used to print output)
        dt: timestep to be taken
        Optional inputs:
            reltol, abstol: relative and absolute tolerances on solution
    Output: 
        nthis: concentrations after the timestep
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, 
                                   :collision_xsect, :crosssection, :dz, :e_profile_type, :error_checking_scheme, 
                                   :H2Oi, :HDOi, :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs, 
                                   :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Hs_dict, 
                                   :inactive_species, :ion_species,  :Jratelist, :molmass, 
                                   :n_all_layers, :n_alt_index, :n_inactive, :neutral_species, :non_bdy_layers, :num_layers, 
                                   :plot_grid, :polarizability, :q, :reaction_network,  :speciesbclist, :speciescolor, :speciesstyle, :Te, :Ti, :Tn, :Tp, 
                                   :Tprof_for_diffusion, :upper_lower_bdy_i, :zmax])

    # absolute and relative tolerance on rate update
    f_abstol = 1e-2
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
            record_atmospheric_state(t, nthis, GV.active_longlived, params[3]; opt=GV.opt, globvars...)
            write_to_log(logfile, ["Too many iterations exception reached at t=$(t), dt=$(dt)"])
            throw(TooManyIterationsException)
        end
        
        nold = deepcopy(nthis)
        
        # we perform an update using newton's method on
        #     f = nthis - nstart - dt*dndt
        
        # to do this we need the rates, dndt, and the jacobian of fval,
        #     d(fval)/d(nthis) = I - dt*chemJ
        dndt, chemJ = get_rates_and_jacobian(nthis, params, t; globvars...)

        # println("dndt: $(dndt)")

        if GV.error_checking_scheme == "old" # ==========================================================
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

            # density absoluite error
            check_n_abserr = (nthis .<= abstol) .|| (nold .<= abstol) .|| (abs.(nthis - nold) .<= abstol)

            # density relative error
            n_relerr = abs.(nthis-nold)./nold
            # println("max(n_relerr) = $(max(n_relerr[.!check_n_abserr]...))")
            check_n_relerr = n_relerr .< reltol

            # # fval absolute error
            # check_f_abserr = (abs.(fval) .<= f_abstol)
            # println("max(fval) = $(max(fval...))\n")
            # f_relerr = abs.(fval./(dt*dndt)) # maybe should measure relative to nstart
            # if length(f_relerr[.!check_f_abserr]) > 0
            #     println("max(f_relerr) = $(max(f_relerr[.!check_f_abserr]...))\n")
            # else
            #     println("max(f_relerr) = 0.0\n")
            # end
            # # fval relative error
            # check_f_relerr = f_relerr .< f_reltol
            
            # Suggest: && all(check_f_relerr .|| check_f_abserr). put [.!check_n_abserr] back if this breaks it
            # Final convergence check
            converged = all(check_n_relerr .|| check_n_abserr) #&& all(check_f_abserr#=[.!check_n_abserr]=# .|| check_f_relerr) # **
        elseif GV.error_checking_scheme == "new" # =======================================================
            # fval - now with dt divided out to keep things stable at long timescales
            fval = (nthis - nstart)/dt - dndt
            # println("newton update (nthis-nstart)/dt: $((nthis-nstart)/dt)")
            identity = sparse(I, length(nthis), length(nthis))
            updatemat = identity/dt - chemJ  
            nthis = nthis - updatemat \ fval # solve_sparse(updatemat, fval)
            nthis[nthis .< 0.] .= 0.

            # density absolute error
            check_n_abserr = (nthis .<= abstol) .|| (nold .<= abstol) .|| (abs.(nthis - nold) .<= abstol)

            # density relative error
            n_relerr = abs.(nthis-nold)./nold
            if verbose==true
                println("max(n_relerr) = $(max(n_relerr[.!check_n_abserr]...))")
            end

            # TROUBLESHOOTING - what's going on at the index where the maximum occurs?
            imax_nrel = argmax(n_relerr[.!check_n_abserr])
            # println("nthis[imax]: $(nthis[.!check_n_abserr][imax_nrel])")
            # println("nold[imax]: $(nold[.!check_n_abserr][imax_nrel])")
            # println("dndt[imax]: $(dndt[.!check_n_abserr][imax_nrel])")
            
            check_n_relerr = n_relerr .< reltol

            # absolute fval error 
            check_f_abserr = (abs.(fval) .<= f_abstol)
            if verbose==true
                if length(fval[.!check_n_abserr]) > 0
                    println("max(fval) = $(max((abs.(fval[.!check_n_abserr]))...))")
                else
                    println("max(fval) = 0.0")
                end
            end

            # relative fval error
            # The following: fval / (nstart+ dt*dndt) = nthis/(nstart+ dt*dndt) - 1, and the division term should be ~=1, so you can check whether this val < rel tol. 
            f_relerr = abs.(fval./(nthis)) #abs.(fval ./ dndt) # NEW when dt is divided out # 
            if verbose==true
                if length(f_relerr[.!check_n_abserr .&& .!check_f_abserr]) > 0
                    println("max(f_relerr) = $(max(f_relerr[.!check_n_abserr .&& .!check_f_abserr]...))")
                else
                    println("max(f_relerr) = 0.0")
                end     
            end
            check_f_relerr = f_relerr .< f_reltol
            
            # Debugging:
            # println("density error: $(all(check_n_relerr .|| check_n_abserr))")
            # println("fval error: $(all(check_f_abserr[.!check_n_abserr] .|| check_f_relerr[.!check_n_abserr]))")
            # println()

            # Final convergence check
            # Turning on the f error check seems to cause it to fail very very early in the run --6 April 
            converged = (all(check_n_relerr .|| check_n_abserr)) #&& all(check_f_abserr[.!check_n_abserr] .|| check_f_relerr[.!check_n_abserr]))
        end # new =====================================================================================

        iter += 1
    end

    return nthis
end

function update!(n_current::Dict{Symbol, Array{ftype_ncur, 1}}, t, dt; abstol=1e-12, reltol=1e-2, globvars...)
    #= 
    update n_current using the coupled reaction network, moving to the next timestep.
    Input:
        n_current: Present atmospheric state
        t: total time elapsed
        dt: timestep
    Output:
        n_current but updated
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, :crosssection, :Dcoef_arr_template, :dz, 
                                   :e_profile_type, :error_checking_scheme, :H2Oi, :HDOi, :Hs_dict, :inactive_species, :ion_species, :Jratelist, :molmass, 
                                   :n_all_layers, :n_alt_index, :n_inactive, :neutral_species, :non_bdy_layers, :num_layers, 
                                   :plot_grid, :polarizability, :q, :reaction_network, :solarflux, :speciesbclist, :speciescolor, :speciesstyle, :Te, :Ti, :Tn, :Tp, :Tprof_for_diffusion, 
                                   :upper_lower_bdy_i, :zmax,])
        

    M = n_tot(n_current; GV.all_species) 
    E = electron_density(n_current; GV.e_profile_type, GV.non_bdy_layers, GV.ion_species)

    # global params for simulation
    params = [GV.Dcoef_arr_template, M, E] 

    # get current long-lived species concentrations
    nstart = flatten_atm(n_current, GV.active_longlived; GV.num_layers)

    # update to next timestep
    nend = next_timestep(nstart, params, t, dt; abstol=abstol, reltol=reltol, globvars...)
    #    println("max(nend-nstart) = $(max((nend-nstart)...))")

    # retrieve the shortlived species from their storage and flatten them
    n_short = flatten_atm(external_storage, GV.active_shortlived; GV.num_layers)  

    n_current = compile_ncur_all(nend, n_short, GV.n_inactive; GV.active_longlived, GV.active_shortlived, GV.inactive_species, GV.num_layers)

    # ensure Jrates are included in n_current
    update_Jrates!(n_current; GV.Jratelist, GV.crosssection, GV.num_layers, GV.absorber, GV.dz, GV.solarflux)

    return n_current
end

function converge(n_current::Dict{Symbol, Array{ftype_ncur, 1}}, log_t_start, log_t_end; verbose=false, t_to_save=nothing,
                  abstol=1e-12, reltol=1e-2, globvars...)
    #= 
    Calls update! in logarithmiclly spaced timesteps until convergence, returning converged atmosphere 

    Input:
        n_current: Starting atmospheric state
        log_t_start: power of 10 for starting timetsep size
        log_t_end; power of 10 for largest possible timestep size 
    Output: 
        n_current, but fully converged (hopefully).
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:absorber, :active_species, :active_longlived, :active_shortlived, :all_species, :alt, :crosssection, 
                                   :Dcoef_arr_template, :dt_decr_factor, :dt_incr_factor, :dz, :e_profile_type, :error_checking_scheme, 
                                   :H2Oi, :HDOi, :Hs_dict, :inactive_species, :ion_species, :Jratelist, :molmass, 
                                   :n_all_layers, :n_alt_index, :n_inactive, :n_steps, :neutral_species, :non_bdy_layers, :num_layers, 
                                   :plot_grid, :polarizability, :q, :reaction_network, :solarflux, :speciesbclist, :speciescolor, :speciesstyle, 
                                   :Te, :Ti, :Tn, :Tp, :timestep_type, :Tprof_for_diffusion, :upper_lower_bdy_i, :zmax])
    
    # A combination of log timesteps when simulation time is low, and linear after
    if GV.timestep_type=="log-linear"
        println("Using a combo of log and linear timesteps")
        log_timesteps = 10. .^(range(log_t_start, stop=log_t_end, length=GV.n_steps));
        linear_timestep = sol_in_sec*7

        # Set up the trackers
        total_time = 0
        ts_i = 1
        dt = log_timesteps[ts_i]
        iters = 0

        # Do log steps for the first day of simulation time while the system is still quite stiff
        while (total_time+dt <= sol_in_sec)
            iters += 1 

            # println("t  = $(total_time)")
            # println("dt = $(dt)")

            # Update total time and model 
            total_time += dt # update total time
            # n_old = deepcopy(n_current)
            n_current = update!(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
            # if verbose==true
            #     println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in GV.all_species]...))\n\n")
            # end

            # Increase timestep
            ts_i += 1
            dt = log_timesteps[ts_i]

            # Check to see if timestep increase will put us over the one-day barrier
            if (total_time+dt > sol_in_sec)
                dt = sol_in_sec - total_time
                total_time += dt
                # n_old = deepcopy(n_current)
                n_current = update!(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
                # if verbose==true
                #     println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in GV.all_species]...))\n\n")
                # end
            end
        end

        write_to_log(logfile, ["Final time after log while loop: $(total_time) with $(iters) iters\n"], mode="a")

        # Reset iters
        iters=0

        # Now do the human-scale timesteps 
        while (total_time < season_length_in_sec)
            iters += 1
            println("$(total_time) + $(linear_timestep) = new total time $(total_time+linear_timestep)")
            total_time += linear_timestep
            # n_old = deepcopy(n_current)
            n_current = update!(n_current, total_time, linear_timestep; abstol=abstol, reltol=reltol, globvars...)
            # if verbose==true
            #     println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in GV.all_species]...))\n\n")
            # end
            
            # Check if the next step will put us past the end
            if (total_time+linear_timestep > season_length_in_sec)

                dt = season_length_in_sec - total_time
                println("Adjusting dt for the last timestep. dt=$(dt)")
                total_time += dt
                # n_old = deepcopy(n_current)
                n_current = update!(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
                # if verbose==true
                #     println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in GV.all_species]...))\n\n")
                # end
            end
        end

        write_to_log(logfile, ["Final time after linear timesteps: $(total_time) with $(iters) iters"], mode="a")

    # Static timesteps, logarithmically spaced
    elseif GV.timestep_type=="static-log"
        log_time_steps = 10. .^(range(log_t_start, stop=log_t_end, length=GV.n_steps))

        total_time = 0.0
        for dt in log_time_steps
            total_time += dt
            # if verbose==true
            #     println("t  = $total_time")
            #     println("dt = $dt")
            #     println()
            # end
            # n_old = deepcopy(n_current)
            n_current = update!(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
            # if verbose==true
            #     println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in GV.all_species]...))\n\n")
            # end
        end
    # Dynamically adjusting timesteps, logarithmic
    elseif GV.timestep_type=="dynamic-log"
        println("Using dynamically adjusting timesteps")
        # FANCY: Modifies timesteps when it gets stuck due to large timesteps being > rel error.
        dt = 10.0^log_t_start
        
        total_time = 0.0
        goodsteps = 0
        goodstep_limit = 1 # Need at least this many successful iterations before increasing timestep
       
        # the first clause before the & will help prevent the simulation stalling out if it has to reduce dt too much.
        while (dt < 10.0^log_t_end) | (total_time <= season_length_in_sec)
            
            # n_old = deepcopy(n_current)
            # n_temp = deepcopy(n_current)

            # SPECIAL BLOCK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # if seasonal_cycle == true
            #     println("Alert: Using the special block for seasonal cycle, forcing certain timesteps")
            #     if ((total_time + dt) > checkpoint) | (checkpoint == season_length_in_sec)
            #         println("total_time+dt = $(total_time+dt). Entering special block to try and meet $(checkpoint)")
            #         temp_dt = checkpoint - total_time
            #         n_temp = update!(n_temp, total_time+temp_dt, temp_dt; abstol=abstol, reltol=reltol, globvars...)
            #     end 

            #     # This should take care of the case where we are trying to collect the season end
            #     if (checkpoint == season_length_in_sec)
            #         temp_total_time = total_time + temp_dt 
            #         println("Entering special block to try and meet $(checkpoint)")
            #         temp_dt = checkpoint - temp_total_time
            #         n_temp = update!(n_temp, temp_total_time+temp_dt, temp_dt; abstol=abstol, reltol=reltol, globvars...)
            #     end 
            # else
            # END SPECIAL BLOCK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            total_time += dt
        
            try
                n_current = update!(n_current, total_time, dt; abstol=abstol, reltol=reltol, globvars...)
                goodsteps += 1
                if goodsteps >= goodstep_limit
                    dt *= GV.dt_incr_factor
                    goodsteps = 0
                end
            catch e
                if e == TooManyIterationsException
                    total_time -= dt
                    dt = dt/GV.dt_decr_factor
                    goodsteps = 0
                else
                    rethrow(e)
                end
            end
                # println("max(n_current-n_old) = $(max([max((abs.(n_current[sp] - n_old[sp]))...) for sp in GV.all_species]...))\n\n")
        end
    end

    if (log10(dt) / log10(total_time) > 1e-6)
        println("Quitting because the timestep is 1 ppm of the total time or smaller, so we'll never advance further")
    end
    
    return n_current, total_time
end


# **************************************************************************** #
#                                                                              #
#                          INITIAL MAIN MODEL SETUP                            #
#                                                                              #
# **************************************************************************** #

#                              Files and folders                                #
#===============================================================================#

create_folder(sim_folder_name, results_dir)
# create_folder("Jrateplots", results_dir*sim_folder_name*"/")
const logfile = results_dir*sim_folder_name*"/log_"*sim_folder_name*".txt"

#                          Load reaction network                                #
#===============================================================================#
#=
Important note: the hot networks returned for H, D, H2, HD implicitly contain
reations that only PRODUCE the given species, because they have been marked
as such in the reaction spreadsheet. 
=#
println("$(Dates.format(now(), "(HH:MM:SS)")) Loading reaction network")
reaction_network, hot_H_network, hot_D_network, 
hot_H2_network, hot_HD_network = load_reaction_network(reaction_network_spreadsheet; saveloc=results_dir*sim_folder_name*"/active_rxns.xlsx",
                                                                       write_rxns=true, get_hot_rxns=ions_included, ions_on=ions_included, all_species)

#              Create evaluatable rate coefficients by reaction                 #
#===============================================================================#
# TODO: These need to ignore the Jrates. 
const hot_H_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_H_network]);
const hot_D_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_D_network]);
const hot_H2_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_H2_network]);
const hot_HD_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hot_HD_network]);

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

#                        Load starting atmosphere                               #
#===============================================================================#
println("$(Dates.format(now(), "(HH:MM:SS)")) Loading atmosphere")
n_current = get_ncurrent(initial_atm_file)

#                       Establish new species profiles                          #
#===============================================================================#

if adding_new_species==true
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
end


#                        Initialize electron profile                            #
#===============================================================================#
E = electron_density(n_current; e_profile_type, non_bdy_layers, ion_species)

#                          Set up the water profile                             #
#===============================================================================#
# If you want to completely wipe out the water profile and install the initial one
# (i.e. when running a stand alone simulation)
if reinitialize_water_profile
    println("Initializing the water profile anew (reinitialize_water_profile=true)")
    # hygropause_alt is an optional argument. If using, must be a unit of length in cm i.e. 40e5 = 40 km.
    println("$(Dates.format(now(), "(HH:MM:SS)")) Setting up the water profile...")
    setup_water_profile!(n_current; dust_storm_on=dust_storm_on, tanh_prof=water_case, all_species, num_layers, non_bdy_layers, DH, alt, plot_grid, # hygropause_alt, 
                                    H2O_excess, HDO_excess, ealt=excess_peak_alt, H2Osat, water_mixing_ratio, n_alt_index, results_dir, 
                                    sim_folder_name, upper_lower_bdy_i)
end

# If you want to just modify the water profile, i.e. when running several simulations
# in succession to simulate seasons: 
if update_water_profile     
    fdict = Dict("high"=>50, "low"=>0.05)
    if water_case!="standard"
        new_frac_H2O = (n_current[:H2O] ./ n_tot(n_current; n_alt_index, all_species)) .* water_tanh_prof(non_bdy_layers./1e5; f=fdict[water_case])
        new_frac_HDO = (n_current[:HDO] ./ n_tot(n_current; n_alt_index, all_species))  .* water_tanh_prof(non_bdy_layers./1e5; f=fdict[water_case])
        # Update the densities, changing only below the fixed point
        n_current[:H2O][1:upper_lower_bdy_i] = new_frac_H2O[1:upper_lower_bdy_i] .* n_current[:H2O][1:upper_lower_bdy_i]
        n_current[:HDO][1:upper_lower_bdy_i] = new_frac_HDO[1:upper_lower_bdy_i] .* n_current[:HDO][1:upper_lower_bdy_i]
        plot_water_profile(new_frac_H2O, new_frac_HDO, n_current[:H2O], n_current[:HDO], results_dir*sim_folder_name; plot_grid) 
    end
end

#           Define storage for species/Jrates not solved for actively           #
#===============================================================================#
# Short lived species, when defined, are solved for assuming photochemical equilibrium
# outside of the primary ODE solver. Inactive species never change during simulation.
# Jrates must be stored here because they have to be updated alongside evolution
# of the atmospheric densities--the solver doesn't handle their values currently.
const external_storage = Dict([j=>n_current[j] for j in union(short_lived_species, inactive_species, Jratelist)])
const n_inactive = flatten_atm(n_current, inactive_species; num_layers)

# Calculate precipitable microns, including boundary layers (assumed same as nearest bulk layer)
H2Oprum = precip_microns(:H2O, [n_current[:H2O][1]; n_current[:H2O]; n_current[:H2O][end]]; molmass)
HDOprum = precip_microns(:HDO, [n_current[:HDO][1]; n_current[:HDO]; n_current[:HDO][end]]; molmass)

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
        active_longlived_species_rates: rate expressions for species which are actively solved for.
        short_lived_density_eqn: array containing expressions for density of short-lived species, i.e. n = P/L, 
                                 also includes expression when a species is quadratic (like OH + OH)
        shortlived_density_inputs: Array containing P, L, and the determinant for the short lived density
                                   expressions. Needed because if P and L are both < machine epsilon, problems occur
        equilibrium_eqn_terms: Another storage array, see below for explanation I'm tired
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

# This can be used if you want to calculate M, E at computation
# const Mexpr = Expr(:call, :+, all_species...) 
# const Eexpr = Expr(:call, :+, ion_species...)
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
        
        # This section calculates the net change but arranges entries by size, so we don't have floating point errors.
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

        # M = $Mexpr 
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

const crosssection = populate_xsect_dict(photochem_data_files; ion_xsects=ions_included, Tn=Tn_arr, n_all_layers)

# **************************************************************************** #
#                                                                              #
#                                SOLAR INPUT                                   #
#                                                                              #
# **************************************************************************** #
solarflux = readdlm(code_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:] # 2000
solarflux[:,2] = solarflux[:,2] * cosd(SZA)  # Adjust the flux according to specified SZA

lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

if !(setdiff(solarflux[:,1],lambdas)==[])
    throw("Solar flux wavelengths don't match cross section wavelengths!")
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
bc_type = Dict("n"=>"density", "f"=>"thermal flux", "v"=>"velocity", "ntf"=>"nonthermal flux")
for k in keys(speciesbclist)
    for k2 in keys(speciesbclist[k])
        push!(PARAMETERS_BCS, ("$(string(k))", "$(bc_type[string(k2)])", "$(speciesbclist[k][k2][1])", "$(speciesbclist[k][k2][2])")) 
    end
end

# Crosssection log messages
for k in keys(photochem_data_files)  # cross sections
    for k2 in keys(photochem_data_files[k])
        push!(PARAMETERS_XSECTS, ("$(k)", "$(k2)", photochem_data_files[k][k2]))
    end
end

push!(PARAMETERS_CONDITIONS, ("TOTAL_H2O", H2Oprum, "pr micrometers"))
push!(PARAMETERS_CONDITIONS, ("TOTAL_HDO", HDOprum, "pr micrometers"))
push!(PARAMETERS_CONDITIONS, ("TOTAL_WATER", H2Oprum+HDOprum, "pr micrometers"))
    
write_to_log(logfile, ["Description: $(optional_logging_note)"], mode="w")

# **************************************************************************** #
#                                                                              #
#                       FINAL SETUP TO CONVERGE!                               #
#                                                                              #
# **************************************************************************** #

# Uncomment this line if you'd like to add an extra parcel to some species. You must specify the species.
# n_current[:D] = map(x->1e5*exp(-((x-184)/20)^2), non_bdy_layers/1e5) + n_current[:D]

# write initial atmospheric state ==============================================
write_atmosphere(n_current, results_dir*sim_folder_name*"/initial_atmosphere.h5"; alt, num_layers, hrshortcode, rshortcode)

# Plot initial temperature and water profiles ==================================
plot_temp_prof(Tn_arr; savepath=results_dir*sim_folder_name, Tprof_2=Ti_arr, Tprof_3=Te_arr, alt)

# Absolute tolerance
if problem_type == "Gear"
    const atol = 1e-12 # absolute tolerance in ppm, used by Gear solver # NOTE: I think this is actually #/cm³ not ppm, because n_i+1 - n_i is compared against it.--Eryn
    const abs_tol_for_plot = fill(atol, length(n_tot(n_current; all_species)))
else
    # absolute tolerance relative to total atmosphere density, used by DifferentialEquations.jl solvers
    const atol = 1e-12 .* [[n_tot(n_current, a; n_alt_index, all_species) for sp in active_longlived, a in non_bdy_layers]...] 
    const abs_tol_for_plot = 1e-12 .* n_tot(n_current; n_alt_index, all_species) # calculates 1 ppt of the total density at each altitude.
end
    
# Plot initial atmosphere condition  ===========================================
println("$(Dates.format(now(), "(HH:MM:SS)")) Plotting the initial condition")
plot_atm(n_current, results_dir*sim_folder_name*"/initial_atmosphere.png", abs_tol_for_plot, E; # initial electron profile
         t="initial state", neutral_species, ion_species, plot_grid, speciescolor, speciesstyle, zmax, hrshortcode, rshortcode) 

# Create a list to keep track of stiffness ratio ===============================
# const stiffness = [] # Tirn this on if you are trying to check Jacobian eigenvalues.

# Simulation time range and when to save a snapshot 
const mindt = dt_min_and_max[converge_which][1]
const maxdt = dt_min_and_max[converge_which][2]
# times to save for cycling:
if seasonal_cycle==true
    const times_to_save = [1., 60., 3600., sol_in_sec/2] # save at 1 sec, 1 min, 1 hour, 1/2 day
    append!(times_to_save, collect(sol_in_sec:sol_in_sec*7:season_length_in_sec)) # save every 7 mars days
    push!(times_to_save, times_to_save[end]+sol_in_sec*( (season_length_in_sec - times_to_save[end])/sol_in_sec) ) # The last save step is less than a week forward so we don't go over the timeframe.
    println("I will try to save the following times: $(times_to_save)")
else
    const times_to_save = [10.0^t for t in mindt:maxdt]
end
t_i = 1  # index for iterating through times_to_save. 
checkpoint = times_to_save[t_i]  # This variable lets us save the atmospheric state at roughly every t = power of 10 with the Gear method, 
                                 # or periodically save the state using the Julia solvers in case they don't reach the end.
plotnum = 1 # To order the plots for each timestep correctly in the folder so it's easy to page through them.

# Record setup time
t5 = time()
write_to_log(logfile, ["\nRequesting timesteps $(join(times_to_save, ", "))", "\n$(Dates.format(now(), "(HH:MM:SS)")) Setup time $(format_sec_or_min(t5-t4))"], mode="a")

# **************************************************************************** #
#                                                                              #
#             First compile and call of Jacobian when using Double64           #
#                                                                              #
# **************************************************************************** #

# This code needs to run once outside evolve_atmosphere because it takes around 30-40 
# minutes when running with Double64
if ftype_ncur==Double64
    println("$(Dates.format(now(), "(HH:MM:SS)")) Compiling and calling the chemical jacobian outside evolve_atmosphere (this will take ~45 min)...")
    write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Started first chemical jacobian compile")

    # Set up the initial state and check for any problems 
    M = n_tot(n_current; all_species)
    E = electron_density(n_current; e_profile_type, non_bdy_layers, ion_species)

    nstart = flatten_atm(n_current, active_longlived; num_layers)
    find_nonfinites(nstart, collec_name="nstart")

    # Set up parameters    
    Dcoef_arr_template = zeros(size(Tn_arr))  # For making diffusion coefficient calculation go faster 
    params = [#inactive, inactive_species, active_species, active_longlived, active_shortlived, Tn_arr, Ti_arr, Te_arr, Tplasma_arr, 
              Dcoef_arr_template, M, E] # E FIX ATTEMPT
    params_exjac = deepcopy(params)  # I think this is so the Dcoef doesn't get filled in with the wrong info?

    # Call the expensive rate function
    t_before_jac = time()
    dndt, example_jacobian = get_rates_and_jacobian(nstart, params_exjac, 0.0) # TODO: Fill in global variables 
    t_after_jac = time()

    write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Finished first chemical jacobian compile")
    println("$(Dates.format(now(), "(HH:MM:SS)")) ...finished.\nFirst jacobian compile+call took $(format_sec_or_min(t_after_jac-t_before_jac))\n")

    find_nonfinites(dndt, collec_name="example_rates")
    find_nonfinites(example_jacobian, collec_name="example_jacobian")    
end

println("Time to beginning convergence is $(format_sec_or_min(time()-t1))\n\n")

# **************************************************************************** #
#                                                                              #
#                         CONVERGE THE ATMOSPHERE                              #
#                                                                              #
# **************************************************************************** #

ti = time()
println("$(Dates.format(now(), "(HH:MM:SS)")) Beginning convergence")

Dcoef_arr_template = zeros(size(Tn_arr)) # initialize diffusion coefficient array

atm_soln = Dict()

try
    global atm_soln, sim_time = evolve_atmosphere(n_current, mindt, maxdt; t_to_save=times_to_save, abstol=atol, reltol=rel_tol, 
                                 # glob vars from here.  
                                 absorber, active_species, active_longlived, active_shortlived, all_species, alt, 
                                 collision_xsect, crosssection, Dcoef_arr_template, dt_incr_factor, dt_decr_factor, dz, 
                                 e_profile_type, error_checking_scheme, timestep_type, H2Oi, HDOi, 
                                 hot_H_network, hot_H_rc_funcs, hot_D_network, hot_D_rc_funcs, hot_H2_network, hot_H2_rc_funcs, hot_HD_network, hot_HD_rc_funcs,
                                 hrshortcode, Hs_dict,
                                 ion_species, inactive_species, Jratelist, molmass, 
                                 neutral_species, non_bdy_layers, num_layers, n_all_layers, n_alt_index, n_inactive, n_steps, 
                                 polarizability, plot_grid, q, reaction_network, rshortcode, solarflux, speciesbclist, speciescolor, speciesstyle, 
                                 Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, Tprof_for_diffusion, transport_species, opt="",
                                 upper_lower_bdy_i, zmax)
catch y
    XLSX.writetable("$(results_dir)$(sim_folder_name)/PARAMETERS.xlsx", "General"=>PARAMETERS_GEN, "AtmosphericConditions"=>PARAMETERS_CONDITIONS, "SpeciesLists"=>PARAMETERS_SPLISTS,
                    "Solver"=>PARAMETERS_SOLVER, "Crosssections"=>PARAMETERS_XSECTS, "BoundaryConditions"=>PARAMETERS_BCS)
    write_to_log(logfile, "Terminated before completion at $(format_sec_or_min(time()-ti))", mode="a")
    throw("ERROR: Simulation terminated before completion with exception:")
end

tf = time() 

write_to_log(logfile, "Finished!\nSimulation active convergence runtime $(format_sec_or_min(tf-ti))", mode="a")

println("$(Dates.format(now(), "(HH:MM:SS)")) Simulation active convergence runtime $((tf-ti)/60) minutes")

# **************************************************************************** #
#                                                                              #
#                      WRITE OUT THE SIMULATION RESULTS                        #
#                                                                              #
# **************************************************************************** #

# ! NO GUARANTEE THAT SS AND ODE SOLVERS WORK RIGHT NOW !
if problem_type == "SS"
    # Update short-lived species one more time
    println("One last update of short-lived species")
    n_short = flatten_atm(external_storage, active_shortlived; num_layers)
    Jrates = deepcopy(Float64[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:num_layers])
    set_concentrations!(external_storage, atm_soln.u, n_short, inactive, 
                        active_longlived, active_shortlived, inactive_species, Jrates, Tn_arr, Ti_arr, Te_arr)
    nc_all = merge(external_storage, unflatten_atm(atm_soln.u, active_longlived; num_layers))

    println("Plotting final atmosphere, writing out state")
    # Make final atmosphere plot
    plot_atm(nc_all, [neutral_species, ion_species], results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state", plot_grid, 
                      speciescolor, speciesstyle, zmax, abs_tol_for_plot)


    write_final_state(nc_all, results_dir, sim_folder_name, final_atm_file; alt, num_layers, hrshortcode, Jratedict=Jrates, rshortcode, external_storage)
    write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Making production/loss plots", mode="a")
    println("Making production/loss plots (this tends to take several minutes)")
    plot_production_and_loss(nc_all, results_dir, sim_folder_name; nonthermal=nontherm, all_species, alt, chem_species, collision_xsect, 
                              dz, hot_D_rc_funcs, hot_H_rc_funcs, hot_H2_rc_funcs, hot_HD_rc_funcs, Hs_dict, 
                              hot_H_network, hot_D_network, hot_H2_network, hot_HD_network, hrshortcode, ion_species, Jratedict,
                              molmass, neutral_species, non_bdy_layers, num_layers, n_all_layers, n_alt_index, polarizability, 
                              plot_grid, q, rshortcode, reaction_network, speciesbclist, Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, 
                              Tprof_for_Hs, Tprof_for_diffusion, transport_species, upper_lower_bdy_i, upper_lower_bdy, zmax)
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
            write_atmosphere(nc_all, results_dir*sim_folder_name*"/atm_state_t_$(timestep).h5"; alt, num_layers, hrshortcode, rshortcode) 
        elseif i == L
            # Update short-lived species one more time
            println("One last update of short-lived species")
            local n_short = flatten_atm(external_storage, active_shortlived; num_layers)
            local Jrates = deepcopy(Float64[external_storage[jr][ialt] for jr in Jratelist, ialt in 1:num_layers])
            set_concentrations!(external_storage, atm_state, n_short, inactive, Jrates; active_longlived, active_shortlived, 
                               inactive_species, Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, num_layers)
            local nc_all = merge(external_storage, unflatten_atm(atm_state, active_longlived; num_layers))

            # Make final atmosphere plot
            println("Plotting final atmosphere, writing out state")
            plot_atm(nc_all, results_dir*sim_folder_name*"/final_atmosphere.png", t="final converged state", abs_tol_for_plot; neutral_species, ion_species, 
                     plot_grid, speciescolor, speciesstyle, zmax)

            write_final_state(nc_all, results_dir, sim_folder_name, final_atm_file; alt, num_layers, hrshortcode, Jratedict=Jrates, rshortcode, external_storage)
            write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Making production/loss plots", mode="a")
            println("Making production/loss plots (this tends to take several minutes)")
            plot_production_and_loss(nc_all, results_dir, sim_folder_name; nonthermal=nontherm, all_species, alt, chem_species, collision_xsect, 
                                      dz, hot_D_rc_funcs, hot_H_rc_funcs, hot_H2_rc_funcs, hot_HD_rc_funcs, Hs_dict, 
                                      hot_H_network, hot_D_network, hot_H2_network, hot_HD_network, hrshortcode, ion_species, Jratedict,
                                      molmass, neutral_species, non_bdy_layers, num_layers, n_all_layers, n_alt_index, polarizability, 
                                      plot_grid, q, rshortcode, reaction_network, speciesbclist, Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, 
                                      Tprof_for_Hs, Tprof_for_diffusion, transport_species, upper_lower_bdy_i, upper_lower_bdy, zmax)

        end
        global i += 1 
    end
elseif problem_type == "Gear"
    # Plot the final atmospheric state
    println("Plotting final atmosphere, writing out state")
    plot_atm(atm_soln, results_dir*sim_folder_name*"/final_atmosphere.png", abs_tol_for_plot, sum([atm_soln[i] for i in ion_species]); 
             t="final converged state, total time = $(sim_time)", neutral_species, ion_species, plot_grid, speciescolor, speciesstyle, zmax, hrshortcode, rshortcode)

    # Collect the J rates
    Jratedict = Dict([j=>external_storage[j] for j in keys(external_storage) if occursin("J", string(j))])
    # Write out the final state to a unique file for easy finding
    write_final_state(atm_soln, results_dir, sim_folder_name, final_atm_file; alt, num_layers, hrshortcode, Jratedict, rshortcode, external_storage)
    
    write_to_log(logfile, "$(Dates.format(now(), "(HH:MM:SS)")) Making production/loss plots", mode="a")
    println("$(Dates.format(now(), "(HH:MM:SS)")) Making production/loss plots (this tends to take several minutes)")
    # make production and loss plots
    plot_production_and_loss(atm_soln, results_dir, sim_folder_name; nonthermal=nontherm, all_species, alt, chem_species, collision_xsect, 
                              dz, hot_D_rc_funcs, hot_H_rc_funcs, hot_H2_rc_funcs, hot_HD_rc_funcs, Hs_dict, 
                              hot_H_network, hot_D_network, hot_H2_network, hot_HD_network, hrshortcode, ion_species, Jratedict,
                              molmass, neutral_species, non_bdy_layers, num_layers, n_all_layers, n_alt_index, polarizability, 
                              plot_grid, q, rshortcode, reaction_network, speciesbclist, Tn=Tn_arr, Ti=Ti_arr, Te=Te_arr, Tp=Tplasma_arr, 
                              Tprof_for_Hs, Tprof_for_diffusion, transport_species, upper_lower_bdy_i, upper_lower_bdy, zmax)

else
    throw("Invalid problem_type")
end 

t7 = time()

# Write out the parameters as the final step 

XLSX.writetable("$(results_dir)$(sim_folder_name)/PARAMETERS.xlsx", "General"=>PARAMETERS_GEN, "AtmosphericConditions"=>PARAMETERS_CONDITIONS, "SpeciesLists"=>PARAMETERS_SPLISTS,
                "Solver"=>PARAMETERS_SOLVER, "Crosssections"=>PARAMETERS_XSECTS, "BoundaryConditions"=>PARAMETERS_BCS)
println("Saved parameter spreadsheet")
write_to_log(logfile, "Simulation total runtime $(format_sec_or_min(t7-t1))", mode="a")
println("Simulation finished!")

