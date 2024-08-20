# **************************************************************************** #
#                                                                              #
#          Non-core functions for analyzing the chemistry and transport        #
#                                                                              #
# These functions are generally used after the model has run and you want to   #
# Do some analysis on the output. I don't believe any of them are called       #
# inside of the core of the model (Core.jl). If so, they should be in there.   #
# If such functions are found, please open an issue on Github.                 #
# **************************************************************************** #

#===============================================================================#
#                           Chemistry functions                                 #
#===============================================================================#

function chemical_lifetime(s::Symbol, atmdict; globvars...)
    #=
    Calculates chemical lifetime of a molecule s in the atmosphere atmdict. 
    
    Good for comparing with the results of diffusion_timescale.
    =#
    GV = values(globvars)
    required = [:all_species, :Jratelist, :n_alt_index, :ion_species, :num_layers, :reaction_network, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    loss_all_rxns, ratecoefs = get_volume_rates(s, atmdict; species_role="reactant", which="all", remove_sp_density=true, 
                                               GV.all_species, GV.ion_species, GV.num_layers, GV.reaction_network, 
                                               Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1])

    total_loss_by_alt = zeros(size(GV.Tn[2:end-1]))

    for k in keys(loss_all_rxns)
        total_loss_by_alt += loss_all_rxns[k]
    end

    return chem_lt = 1 ./ total_loss_by_alt
end

function get_column_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; which="all", sp2=nothing, role="product", startalt_i=1, returntype="df", globvars...)
    #=
    Input:
        sp: species for which to search for reactions
        atmdict: the present atmospheric state to calculate on
        Tn, Ti, Te: Arrays of the temperature profiles including boundary layers
        bcdict: Boundary conditions dictionary specified in parameters file
        which: whether to do photochemistry, or just bimolecular reactions. "all", "Jrates" or "krates"
        sp2: optional second species to include, i.e. usually sp's ion.
        role: "product" or "reactant" only.
        startalt_i: Index of the altitude at which to start. This lets you only calculate column rate down to a certain altitude, which is
                    useful, for example, for water, which is fixed at 80 km so we don't care what produces/consumes it,
                    or if you want to know the column rate only above a certain altitude.
        returntype: whether to return a dataframe ("df") or a dictionary "dict"
    Output:
        sorted: Total column rates for all reactions of species sp. Sorted, in order of largest rate to smallest. NOT a dictionary.
                sorted[1] is the top production mechanism, e.g.
    =#
    GV = values(globvars)
    required = [:Tn, :Ti, :Te, :all_species, :ion_species, :reaction_network, :num_layers, :dz]
    check_requirements(keys(GV), required)
    
    rxd, coefs = get_volume_rates(sp, atmdict; species_role=role, which=which, globvars...)
                                   
    # Make the column rates dictionary for production
    columnrate = Dict()
    for k in keys(rxd)
        columnrate[k] = sum(rxd[k][startalt_i:end] .* GV.dz)
    end
    
    # Optionally one can specify a second species to include in the sorted result, i.e. a species' ion.
    if sp2 != nothing
        rxd2, coefs2 = get_volume_rates(sp2, atmdict; species_role=role, which=which, globvars...)

        columnrate2 = Dict()

        for k in keys(rxd2)
            columnrate2[k] = sum(rxd2[k][startalt_i:end] .* GV.dz)
        end

        colrate_dict = merge(columnrate, columnrate2)
    else
        colrate_dict = columnrate
    end
    
    sorted = sort(collect(colrate_dict), by=x->x[2], rev=true)

    if returntype=="df"
        return DataFrame([[names(DataFrame(sorted))]; collect.(eachrow(DataFrame(sorted)))], [:Reaction; :ColumnRate])
    else 
        return sorted 
    end
end

function get_volume_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; species_role="both", which="all", remove_sp_density=false, globvars...)
    #=
    Input:
        sp: Species name
        atmdict: Present atmospheric state dictionary
        Tn, Ti, Te: temperature arrays
        species_role: whether to look for the species as a "reactant", "product", or "both".  If it has a value, so must species.
        which: "all", "Jrates", "krates". Whether to fill the dictionary with all reactions, only photochemistry/photoionization 
               (Jrates) or only chemistry (krates).
       remove_sp_density: if set to true, the density of sp will be removed from the calculation k[sp][B][C].... Useful for chemical lifetimes.
    Output: 
        rxn_dat: Evaluated rates, i.e. k[A][B], units #/cm^3/s for bimolecular rxns
        rate_coefs: Evaluated rate coefficients for each reaction 
    =#

    GV = values(globvars)
    required = [:all_species, :ion_species, :num_layers, :reaction_network, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    # Make sure temperatures are correct format
    @assert length(GV.Tn)==GV.num_layers 
    @assert length(GV.Ti)==GV.num_layers
    @assert length(GV.Te)==GV.num_layers

    # Fill in the rate x density dictionary ------------------------------------------------------------------------------
    rxn_dat =  Dict{String, Array{ftype_ncur, 1}}()
    rate_coefs = Dict{String, Array{ftype_ncur, 1}}()

    filtered_rxn_list = filter_network(sp, which, species_role; GV.reaction_network)

    for rxn in filtered_rxn_list
        # get the reactants and products in string form for use in plot labels
        rxn_str = format_chemistry_string(rxn[1], rxn[2])

        # Fill in rate coefficient * species density for all reactions
        if typeof(rxn[3]) == Symbol # for photodissociation
            if remove_sp_density==false
                rxn_dat[rxn_str] = atmdict[rxn[1][1]] .* atmdict[rxn[3]]
            else 
                rxn_dat[rxn_str] = 1 .* atmdict[rxn[3]]  # this will functionally be the same as the rate coefficient for photodissociation.
            end
            rate_coefs[rxn_str] = atmdict[rxn[3]]
        else                        # bi- and ter-molecular chemistry
            remove_me = remove_sp_density==true ? sp : nothing
            density_prod = reactant_density_product(atmdict, rxn[1]; removed_sp=remove_me, globvars...)
            thisrate = typeof(rxn[3]) != Expr ? :($rxn[3] + 0) : rxn[3]
            rate_coef = eval_rate_coef(atmdict, thisrate; globvars...)

            rxn_dat[rxn_str] = density_prod .* rate_coef # This is k * [R1] * [R2] where [] is density of a reactant. 
            if typeof(rate_coef) == Float64
                rate_coef = rate_coef * ones(GV.num_layers)
            end
            rate_coefs[rxn_str] = rate_coef
        end
    end

    return rxn_dat, rate_coefs
end

function get_volume_rates(sp::Symbol, source_rxn::Vector{Any}, source_rxn_rc_func, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Mtot; globvars...)
    #=
    Override to call for a single reaction. Useful for doing non-thermal flux boundary conditions.
    Input:
        sp: Species name
        source_rxn: chemical reaction for which to get the volume rate 
        atmdict: Present atmospheric state dictionary
        Mtot: total atmospheric density
      Output: 
        vol_rates: Evaluated rates, e.g. k[A][B] [#/cm^3/s] for bimolecular rxns, for the whole atmosphere.
    =#

    GV = values(globvars)
    required = [:all_species, :ion_species, :Jratedict, :num_layers, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    # Make sure temperatures are correct format
    @assert length(GV.Tn)==GV.num_layers 
    @assert length(GV.Ti)==GV.num_layers
    @assert length(GV.Te)==GV.num_layers

    # Fill in the rate x density dictionary ------------------------------------------------------------------------------
    if typeof(source_rxn[3]) == Symbol # for photodissociation
        # Look for density of dissociating molecule in atmdict, but rate in Jratedict. This is like this because
        # of the wonky way the code is written to allow for photochemical equilibrium as an option, which 
        # requires the Jrates be stored in an external dictionary because they can't go through the Julia solvers. 
        # Honestly I could probably rewrite everything so that photochem eq is possible with the Gear solver and ditch
        # the Julia solvers entirely but I like having the option and would rather get my PhD and get a pay raise
        # println(keys(GV.Jratedict))
        vol_rates = atmdict[source_rxn[1][1]] .* GV.Jratedict[source_rxn[3]] 
    else                        # bi- and ter-molecular chemistry
        rate_coef = source_rxn_rc_func(GV.Tn, GV.Ti, GV.Te, Mtot)
        vol_rates = reactant_density_product(atmdict, source_rxn[1]; globvars...) .* rate_coef # This is k * [R1] * [R2] where [] is density of a reactant. 
    end
    return vol_rates
end

function make_chemjac_key(fn, fpath, list1, list2) 
    #=
    This somewhat superfluous function makes a key to the chemical jacobian,
    telling which index corresponds to which species. But really it just gives the 
    indices of the entries in all_species, because that's how the jacobian is ordered,
    but this function is written agnostically so that could technically change and this
    function would still work.

    fn: filename to save the key to
    fpath: where to save fn
    list1: jacobian row indices
    list2: jacobian col indices
    =#
    dircontents = readdir(fpath)
    if !(fn in dircontents)
        println("Creating the chemical jacobian row/column key")
        f = open(fpath*"/"*fn, "w")
        write(f, "Chemical jacobian rows (i):\n")
        write(f, "$([i for i in 1:length(list1)])\n")
        write(f, "$(list1)\n\n")
        write(f, "Chemical jacobian cols (j):\n")
        write(f, "$([j for j in 1:length(list2)])\n")
        write(f, "$(list2)\n\n")
        close(f)
    end
end

function reactant_density_product(atmdict::Dict{Symbol, Vector{ftype_ncur}}, reactants; removed_sp=nothing, globvars...)
    #=
    Calculates the product of all reactant densities for a chemical reaction for the whole atmosphere, 
    i.e. for A + B --> C + D, return n_A * n_B.

    Input:
        atmdict: the atmospheric state dictionary
        reactants: a list of reactant symbols.
    Output: 
        density_product: returns n_A * n_B for all altitudes for the reaction A + B --> ...
    =#

    GV = values(globvars)
    required =  [:all_species, :ion_species, :num_layers]
    check_requirements(keys(GV), required)

    if removed_sp != nothing # remove reactant if requested - useful for calculating chemical lifetimes
        deleteat!(reactants, findfirst(x->x==removed_sp, reactants))
    end

    density_product = ones(GV.num_layers)
    for r in reactants
        if r != :M && r != :E
            # species densities by altitude
            density_product .*= atmdict[r]  # multiply by each reactant density
        elseif r == :M
            density_product .*= sum([atmdict[sp] for sp in GV.all_species]) 
        elseif r == :E
            density_product .*= sum([atmdict[sp] for sp in GV.ion_species])
        else
            throw("Got an unknown symbol in a reaction rate: $(r)")
        end
    end

    return density_product 
end 

function volume_rate_wrapper(sp, source_rxns, source_rxn_rc_funcs, atmdict, Mtot; returntype="array", globvars...)
    #=
    Gets altitude-dependent volume production or loss of species sp due to reactions in source_rxns.
    Does NOT care if it is production or loss. This is mainly just a convenient wrapper to get_volume_production
    to return the information in a variety of different useful formats.
    
    Input
        sp: species
        source_rxns: reaction network--should ALREADY be filtered to be either production or loss reactions for sp.
        source_rxn_rc_funcs: Evalutable functions for each reaction. 
        atmdict: present atmospheric state dictionary
        Mtot: Atmospheric density at all altitudes
    Output: 
        array of production or loss by altitude (rows) and reaction  (columns)
    =#
    
    GV = values(globvars)
    required = [:all_species, :alt, :collision_xsect, :ion_species, :Jratedict, :molmass, :non_bdy_layers, :num_layers,  
                                   :n_alt_index, :Tn, :Ti, :Te, :dz, :zmax]
    check_requirements(keys(GV), required)

    rates = Array{ftype_ncur}(undef, GV.num_layers, length(source_rxns))
    
    i=1
    for source_rxn in source_rxns
        rates[:, i] = get_volume_rates(sp, source_rxn, source_rxn_rc_funcs[source_rxn], atmdict, Mtot; globvars..., 
                                                  Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1])
        i += 1
    end

    # Returns an array where rows represent altitudes and columns are reactions.
    if returntype=="by rxn"
        return sum(rates, dims=1)
    elseif returntype=="by alt"
        return sum(rates, dims=2)
    elseif returntype=="array"
        return rates
    elseif returntype=="df" # Useful if you want to look at the arrays yourself.
        ratesdf = DataFrame(rates, vec([format_chemistry_string(r[1], r[2]) for r in source_rxns]))
        return ratesdf
    end
end

#===============================================================================#
#                   Transport and escape functions                              #
#===============================================================================#

function diffusion_timescale(s::Symbol, T_arr::Array, atmdict; globvars...)
    #=
    Inputs:
        s: species symbol
        T_arr: temperature array for the given species 
        atmdict: atmospheric state dict
        Dcoef_template: array of 0s to use in Dcoef!
    Output: Molecular and eddy diffusion timescale (s) by alt
    =#
    
    GV = values(globvars)
    required = [:all_species, :alt, :molmass, :n_alt_index, :neutral_species, :polarizability, :planet, :q, :speciesbclist, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)

    # Get diffusion coefficient array template
    Dcoef_template = zeros(size(T_arr)) 

    # Other stuff
    Hs = scaleH(GV.alt, s, T_arr; GV.molmass)
    ncur_with_bdys =  ncur_with_boundary_layers(atmdict; GV.all_species, GV.n_alt_index)
    
    # Molecular diffusion timescale: H_s^2 / D, scale height over diffusion constant
    D = Dcoef!(Dcoef_template, T_arr, s, ncur_with_bdys; globvars...)
    molec_or_ambi_timescale = (Hs .^ 2) ./ D
   
    # Eddy timescale... this was in here only as scale H... 
    K = Keddy(alt, n_tot(ncur_with_bdys; GV.all_species, GV.molmass); GV.planet) 
    eddy_timescale = (Hs .^ 2) ./ K

    # Combined timescale?!??
    combined_timescale = (Hs .^ 2) ./ (K .+ D)

    return molec_or_ambi_timescale, eddy_timescale, combined_timescale
end

function final_escape(thefolder, thefile; globvars...)
    #=
    thefolder: Folder in which an atmosphere file lives
    thefile: the file containing an atmosphere for which you'd like to calculate the final escape fluxes of H and D.
    =#
    
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
                :q, :molmass, :polarizability, :collision_xsect,
                # From CUSTOMIZATIONS.jl
                :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, 
                # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)
    
    # First load the atmosphere and associated variables.
    atmdict = get_ncurrent(thefolder*thefile);

    vardict = load_from_paramlog(thefolder; globvars...);
    
    # Get Jrate list 
    Jratelist = format_Jrates(thefolder*"active_rxns.xlsx", GV.all_species, "Jratelist"; hot_atoms=true, ions_on=true)[1];
    Jratedict = Dict([j=>atmdict[j] for j in Jratelist])
    
    # Make a dataframe to store things
    escdf = DataFrame("EscapeType"=>["Thermal", "Nonthermal", "Total"], 
                      "H"=>[0, 0, 0], "D"=>[0, 0, 0], "H2"=>[0, 0, 0], "HD"=>[0, 0, 0])

    # Now collect non-thermal and thermal fluxes for each species. 
    for s in ["H", "D", "H2", "HD"]
        nonthermal_esc, thermal_esc = get_transport_PandL_rate(Symbol(s), atmdict; returnfluxes=true, Jratedict, zmax=GV.alt[end],
                                                               hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, hot_HD_network=GV.hHDnet,
                                                               hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, hot_HD_rc_funcs=GV.hHDrc, 
                                                               Hs_dict=vardict["Hs_dict"], ion_species=vardict["ion_species"], neutral_species=vardict["neutral_species"],
                                                               speciesbclist=vardict["speciesbclist"],
                                                               Tprof_for_Hs=vardict["Tprof_for_Hs"], Tprof_for_diffusion=vardict["Tprof_for_diffusion"], 
                                                               transport_species=vardict["transport_species"], 
                                                               Tn=vardict["Tn_arr"], Ti=vardict["Ti_arr"], Te=vardict["Te_arr"], Tp=vardict["Tplasma_arr"],
                                                               globvars...)

        escdf.:($s) = [thermal_esc, nonthermal_esc, thermal_esc+nonthermal_esc]
    end
    
    # Calculate total H atoms lost
    escdf."TotalHAtomsLost" = sum(eachcol(escdf[!, Not([:EscapeType, :H2, :D])])) .+ 2 .* escdf[:, :H2]
    # Calculate total H atoms lost
    escdf."TotalDAtomsLost" = sum(eachcol(escdf[!, Not([:EscapeType, :H2, :H, :TotalHAtomsLost])])) # adds up 1 * D and 1 * HD
    # Calculate total atoms lost
    escdf."TotalHnDAtomsLost" = sum(eachcol(escdf[!, Not([:EscapeType, :D, :H, :H2, :HD])]))
    
    return escdf
end

function fractionation_factor(esc_df, h2o_0, hdo_0; ftype="total")
    #=
    Calculates fractionation factor if given the H2O and HDO at the bottom of the atmosphere, as well as a dataframe
    of escape rates as generated by final_escape.
    =#
    flux_t_D = df_lookup(esc_df, "EscapeType", "Thermal", "D")[1] + df_lookup(esc_df, "EscapeType", "Thermal", "HD")[1]
    flux_t_H = df_lookup(esc_df, "EscapeType", "Thermal", "H")[1] + 2*df_lookup(esc_df, "EscapeType", "Thermal", "H2")[1] + df_lookup(esc_df, "EscapeType", "Thermal", "HD")[1]

    flux_nt_D = df_lookup(esc_df, "EscapeType", "Nonthermal", "D")[1] + df_lookup(esc_df, "EscapeType", "Nonthermal", "HD")[1]
    flux_nt_H = df_lookup(esc_df, "EscapeType", "Nonthermal", "H")[1] + 2*df_lookup(esc_df, "EscapeType", "Nonthermal", "H2")[1] + df_lookup(esc_df, "EscapeType", "Nonthermal", "HD")[1]

    if ftype=="thermal"
        flux_nt_D = 0
        flux_nt_H = 0
    elseif ftype=="nonthermal"
        flux_t_D = 0
        flux_t_H = 0
    end
    
    return f = ((flux_t_D + flux_nt_D) / (flux_t_H + flux_nt_H)) / (hdo_0 / (2 * h2o_0))
end

function get_transport_PandL_rate(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; returnfluxes=false, nonthermal=true, globvars...)
    #=
    Input:
        sp: species for which to return the transport production and loss
        atmdict: species number density by altitude
        returnfluxes: whether to return fluxes (thermal and nonthermal) instead of production/loss
        nonthermal: whether to consider nonthermal escape 
    Output
        Array of production and loss (#/cm³/s) at each atmospheric layer boundary.
        i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
        and the end of the array is the boundary at 249 km.
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :dz, :Hs_dict, :molmass,  :n_alt_index,
                :neutral_species, :num_layers, :polarizability, :q, :speciesbclist, :Te, :Ti, :Tn, :Tp, 
                :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)

    if nonthermal
        required = [:hot_H_network, :hot_D_network, :hot_H_rc_funcs, :hot_D_rc_funcs, 
                    :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Jratedict]
        check_requirements(keys(GV), required)
    end

    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(GV.Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, D_arr; globvars...) 
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict; globvars...)

    # For the bulk layers only to make the loops below more comprehendable: 
    fluxcoefs_bulk_layers = Dict([s=>fluxcoefs_all[s][2:end-1, :] for s in keys(fluxcoefs_all)])

    bc_dict = boundaryconditions(fluxcoefs_all, atmdict, sum([atmdict[sp] for sp in GV.all_species]); nonthermal=nonthermal, globvars...)

    # each element in thesebcs has the format [downward, upward]
    thesebcs = bc_dict[sp]

    # Fill array 
    transport_PL = fill(convert(ftype_ncur, NaN), GV.num_layers)

    # These are the derivatives, which should be what we want (check math)
    transport_PL[1] = ((atmdict[sp][2]*fluxcoefs_bulk_layers[sp][2, 1]  # in from layer above
                        -atmdict[sp][1]*fluxcoefs_bulk_layers[sp][1, 2]) # out to layer above
                    +(-atmdict[sp][1]*thesebcs[1, 1] # out to boundary layer
                      +thesebcs[1, 2])) # in from the boundary layer
    for ialt in 2:length(transport_PL) - 1
        transport_PL[ialt] = ((atmdict[sp][ialt+1]*fluxcoefs_bulk_layers[sp][ialt+1, 1]  # coming in from above
                               -atmdict[sp][ialt]*fluxcoefs_bulk_layers[sp][ialt, 2])    # leaving out to above layer
                             +(-atmdict[sp][ialt]*fluxcoefs_bulk_layers[sp][ialt, 1]     # leaving to the layer below
                               +atmdict[sp][ialt-1]*fluxcoefs_bulk_layers[sp][ialt-1, 2]))  # coming in from below
    end
    transport_PL[end] = ((thesebcs[2, 2] # in from upper boundary layer - (non-thermal loss from flux bc)
                          - atmdict[sp][end]*thesebcs[2, 1]) # (#/cm³) * (#/s) out to space from upper bdy (thermal loss from velocity bc)
                        + (-atmdict[sp][end]*fluxcoefs_bulk_layers[sp][end, 1] # leaving out to layer below
                           +atmdict[sp][end-1]*fluxcoefs_bulk_layers[sp][end-1, 2])) # coming in to top layer from layer below

    # Use these for a sanity check if you like. 
    # println("Activity in the top layer for sp $(sp) AS FLUX:")
    # println("Flux calculated from flux bc. for H and D, this should be the nonthermal flux: $(thesebcs[2, 2]*GV.dz)")
    # println("Calculated flux from velocity bc. For H and D this should be thermal escape: $(atmdict[sp][end]*thesebcs[2, 1]*GV.dz)")
    # println("Down to layer below: $(-atmdict[sp][end]*fluxcoefs_all[sp][end, 1]*GV.dz)")
    # println("In from layer below: $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2]*GV.dz)")

    if returnfluxes
        tflux = atmdict[sp][end]*thesebcs[2, 1]*GV.dz
        if nonthermal
            ntflux = thesebcs[2, 2]*GV.dz
            if sp in [:H, :D, :H2, :HD]
                ntflux = ntflux < 0 ? abs(ntflux) : throw("I somehow got a positive nonthermal flux, meaning it's going INTO the atmosphere? for $(sp)")
            else 
                ntflux = 0 
            end
            return ntflux, tflux
        else 
            return tflux 
        end
    else 
        return transport_PL
    end
end

function get_directional_fluxes(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; nonthermal=true, return_up_n_down=false, globvars...)
    #=
    Returns the flux up and down from each atmospheric cell. 

    Input:
        sp: species for which to return the transport production and loss
        atmdict: species number density by altitude
        returnfluxes: whether to return fluxes (thermal and nonthermal) instead of production/loss
        nonthermal: whether to consider nonthermal escape 
    Output
        Array of production and loss (#/cm³/s) at each atmospheric layer boundary.
        i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
        and the end of the array is the boundary at 249 km.
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :dz, :Hs_dict, :molmass, :n_alt_index,
               :neutral_species, :num_layers, :polarizability, :q, :speciesbclist, :Te, :Ti, :Tn, :Tp, 
               :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species, :use_ambipolar, :use_molec_diff]
    check_requirements(keys(GV), required)

    if nonthermal
        required = [:hot_H_network, :hot_D_network, :hot_H_rc_funcs, :hot_D_rc_funcs, 
                    :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Jratedict]
        check_requirements(keys(GV), required)
    end

    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(GV.Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, D_arr; globvars...) 
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict; globvars...)

    # For the bulk layers only to make the loops below more comprehendable: 
    fluxcoefs_bulk_layers = Dict([s=>fluxcoefs_all[s][2:end-1, :] for s in keys(fluxcoefs_all)])

    bc_dict = boundaryconditions(fluxcoefs_all, atmdict, sum([atmdict[sp] for sp in GV.all_species]); nonthermal=nonthermal, globvars...)

    # each element in thesebcs has the format [downward, upward]
    thesebcs = bc_dict[sp]

    # Fill array 
    flux = fill(convert(ftype_ncur, NaN), (length(GV.alt), 1)) # will store positive and negative values at each alt, with + meaning up, - meaning down.
    up = fill(convert(ftype_ncur, NaN), (length(GV.alt), 1)) # Bethan asked for just up fluxes on time so now I also track them separately. 
    down = fill(convert(ftype_ncur, NaN), (length(GV.alt), 1))
    flux[1] = NaN # 0 alt
    up[1] = NaN
    down[1] = NaN

    # Lower boundary 
    up2 = thesebcs[1, 2] # in from the boundary layer (upwards)
    down2 = atmdict[sp][1]*thesebcs[1, 1] # out to boundary layer (downwards)
    flux[2] = up2 - down2
    up[2] = up2
    down[2] = down2

    # println("Alt $(alt[2]): up $(up2*dz), down $(down2*dz), net $(up2*dz-down2*dz)")

    for ialt in 3:length(flux) - 2
        up_i = atmdict[sp][ialt-1]*fluxcoefs_bulk_layers[sp][ialt-1, 2] # in from layer below (upwards)
        down_i = atmdict[sp][ialt]*fluxcoefs_bulk_layers[sp][ialt, 1]     # out to the layer below [downwards]
        # println("Alt $(GV.alt[ialt]): up $(up*dz), down $(down*dz), net $(up*dz-down*dz)")
        flux[ialt] = up_i - down_i
        up[ialt] = up_i
        down[ialt] = down_i
    end

    up_penult = atmdict[sp][end]*thesebcs[2, 1] # (#/cm³) * (#/s) out to space from upper bdy (thermal loss from velocity bc) (upwards)
    down_penult = thesebcs[2, 2] # in from upper boundary layer (non-thermal loss from flux bc) (downwards) (has to be negative to come out positive)
    flux[end-1] = up_penult - down_penult
    up[end-1] = up_penult
    down[end-1] = down_penult
    # println("Alt $(alt[end-1]): up $(up_penult*dz), down $(down_penult*dz), net $(up_penult*dz-down_penult*dz)")

    flux[end] = NaN
    up[end] = NaN 
    down[end] = NaN

    # mult by dz so the array is truly a flux
    flux = flux .* GV.dz
    up = up .* GV.dz 
    down = down .* GV.dz

    # Use these for a sanity check if you like. 
    # println("Activity in the top layer for sp $(sp) AS FLUX:")
    # println("Flux calculated from flux bc. for H and D, this should be the nonthermal flux: $(thesebcs[2, 2]*GV.dz)")
    # println("Calculated flux from velocity bc. For H and D this should be thermal escape: $(atmdict[sp][end]*thesebcs[2, 1]*GV.dz)")
    # println("Down to layer below: $(-atmdict[sp][end]*fluxcoefs_all[sp][end, 1]*GV.dz)")
    # println("In from layer below: $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2]*GV.dz)")
    # println("net in the second to last layer $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2]*GV.dz -atmdict[sp][end]*fluxcoefs_all[sp][end, 1]*GV.dz)")

    if !return_up_n_down
        return flux
    else 
        return flux, up, down
    end
end

function flux_pos_and_neg(fluxarr) 
    #=
    Input:
        fluxarr: the output of function get_flux. 
    Outputs: 
        This generates two arrays, one with the positive flux
        and one with the negative flux, but all values are positive. This is just so 
        you can easily plot flux on a log axis with different markers for positive and negative.
    =#
    pos = []
    abs_val_neg = []

    for f in fluxarr
        if f > 0
            append!(pos, f)
            append!(abs_val_neg, NaN)
        else
            append!(abs_val_neg, abs(f))
            append!(pos, NaN)
        end
    end
    return pos, abs_val_neg
end

# Note: These functions are probably misleading. They were used to create plots that never made it 
# to publication.
function limiting_flux(sp, atmdict, T_arr; treat_H_as_rare=false, full_equation=true, globvars...)
    #=
    Calculate the limiting upward flux (Hunten, 1973; Zahnle, 2008). 
    Inputs:
        sp: A species that is traveling upwards
        atmdict: present atmospheric state
        T_arr: Array of neutral temperatures
    Output:
        Φ, limiting flux for a hydrostatic atmosphere
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :dz, :non_bdy_layers, :molmass, :n_alt_index, :n_all_layers]
    check_requirements(keys(GV), required)
    
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s
    if treat_H_as_rare==true
        if sp==:H
            thedensity = atmdict[:D]
        elseif sp==:H2 
            thedensity = atmdict[:HD]
        end
    elseif treat_H_as_rare==false
        thedensity = atmdict[sp]
    end

    if length(T_arr)==length(GV.alt)
        T_arr = T_arr[2:end-1]
    end

    Ha = scaleH(atmdict, T_arr; ignore=[sp], globvars..., alt=GV.non_bdy_layers)
    bi = binary_dcoeff_inCO2(sp, T_arr) # AT^s

    if full_equation
        dTdz = zeros(GV.n_all_layers)
        dTdz = dTdz[2:end] = @. (T_arr[2:end] - T_arr[1:end-1]) / GV.dz # make the temp gradient
        print(dTdz)
        fi = thedensity ./ n_tot(atmdict; ignore=[sp], globvars...)
        ma = meanmass(atmdict; ignore=[sp], globvars...) 

        return @. ((bi*fi)/(1+fi)) * ( mH*(ma - GV.molmass[sp]) * (g/(kB*T_arr)) - (thermaldiff(sp)/T_arr) * dTdz[1:end-1])
    else
        D = Dcoef_neutrals(non_bdy_layers, sp, bi, atmdict; globvars...)    
        return (D .* atmdict[sp] ./ Ha) .* (1 .- GV.molmass[sp] ./ meanmass(atmdict; ignore=[sp], globvars...))
    end
end

function limiting_flow_velocity(sp, atmdict, T_arr; globvars...)
    #=
    Calculate the limiting upward flux (Hunten, 1973; Zahnle, 2008). 
    Inputs:
        sp: A species that is traveling upwards
        atmdict: present atmospheric state
        T_arr: Array of neutral temperatures
    Output:
        Φ, limiting flux for a hydrostatic atmosphere
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :n_alt_index, :non_bdy_layers, :molmass]
    check_requirements(keys(GV), required)
    
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s
    Ha = scaleH(atmdict, T_arr[2:end-1]; ignore=[sp], globvars..., alt=GV.non_bdy_layers)
    Hi = scaleH(GV.non_bdy_layers, sp, T_arr[2:end-1]; GV.molmass)
    bi = binary_dcoeff_inCO2(sp, T_arr[2:end-1]) # AT^s
    na = n_tot(atmdict; ignore=[sp], GV.all_species)

    return @. (bi / na) * (1/Ha - 1/Hi)
end

function limiting_flux_molef(sp, atmdict, T_arr; globvars...)
    #=
    Roger requested the limiting flux in in mole fraction. This is actually the same result as above. But this way we're sure
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :non_bdy_layers, :molmass, :n_alt_index]
    check_requirements(keys(GV), required)

    avogadro = 6.022e23

    X = (atmdict[sp] ./ avogadro) ./ (n_tot(atmdict; globvars...) ./ avogadro)
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s

    Ha = scaleH(atmdict, T_arr; globvars..., alt=GV.non_bdy_layers)
    bi = binary_dcoeff_inCO2(sp, T_arr)
    Hi = scaleH(non_bdy_layers, sp, T_arr; globvars...)

    return bi .* X .* (1 ./ Ha - 1 ./ Hi), X
end
