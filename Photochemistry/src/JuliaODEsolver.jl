# **************************************************************************** #
#                                                                              #
#                   Function which utilize the Julia ODE solvers               #
#                                                                              #
# **************************************************************************** #

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
    n_short = flatten_atm(external_storage, GV.active_shortlived, n_horiz; GV.num_layers)  # retrieve the shortlived species from their storage and flatten them

    # Update Jrates
    n_cur_all = compile_ncur_all(n, n_short, GV.n_inactive; GV.active_longlived, GV.active_shortlived, GV.inactive_species, GV.num_layers)

    update_Jrates!(n_cur_all, n_horiz; GV.Jratelist, GV.crosssection, GV.num_layers, GV.absorber, GV.dz, GV.solarflux, enable_horiz_transport=GV.enable_horiz_transport)
    # copy all the Jrates into an external dictionary for storage
    for jr in GV.Jratelist                # time for this is ~0.000005 s
        global external_storage[jr] = n_cur_all[jr]
    end

    # Retrieve Jrates 
    Jrates = deepcopy(ftype_ncur[external_storage[jr][ialt] for jr in GV.Jratelist, ialt in 1:GV.num_layers])

    # and update the shortlived species with the new Jrates - assuming not needed to be done in this function
    # n_short_updated = set_concentrations!(external_storage, n, n_short, n_inactive, active_longlived, active_shortlived, inactive_species, Jrates, Tn, Ti, Te)

    tlower, tup, tdown, tupper = update_transport_coefficients(GV.transport_species, # Species for which to update coefficients so it's not a mistake to pass it twice.
                                                               n_cur_all, D_arr, M, n_horiz; calc_nonthermal=nontherm, globvars...)
                                                               # Tn, Tp, Hs_dict, bcdict=speciesbclist, 
                                                               # all_species, neutral_species, transport_species, molmass, n_alt_index, 
                                                               # polarizability, alt, num_layers, n_all_layers, dz, T_for_diff=Tprof_for_diffusion, q)

    tbackedge, tforwards, tbackwards, tfrontedge =
        update_horiz_transport_coefficients(
            GV.transport_species, n_cur_all, D_arr, M, n_horiz;
            calc_nonthermal=nontherm, globvars...
        )
    
    return chemJmat(n, n_short, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, tforwards, tbackwards, tfrontedge, tbackedge, M, E; globvars...)
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
    n_short = flatten_atm(external_storage, GV.active_shortlived, n_horiz; GV.num_layers)

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
    tlower, tup, tdown, tupper = update_transport_coefficients(GV.transport_species, updated_ncur_all, D_arr, M, n_horiz;
                                                               calc_nonthermal=nontherm, globvars...)

    tbackedge, tforwards, tbackwards, tfrontedge =
        update_horiz_transport_coefficients(
            GV.transport_species, updated_ncur_all, D_arr, M, n_horiz;
            calc_nonthermal=nontherm, globvars...
        )

    dndt .= ratefn(n, n_short_updated, GV.n_inactive, Jrates, tup, tdown, tlower, tupper, tforwards, tbackwards, tfrontedge, tbackedge, M, E; globvars...)
end
