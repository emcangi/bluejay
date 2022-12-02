# **************************************************************************** #
#                                                                              #
#                     Boundary conditions and flux functions                   #
#                                                                              #
# **************************************************************************** #

function boundaryconditions(fluxcoef_dict, atmdict, M; nonthermal=true, globvars...)
    #= 
    Inputs:
        fluxcoef_dict: a dictionary containing the K and D flux coefficients for every species throughout
                       the atmosphere. Format species=>Array(length(all_species), length(alt)).
                       Because it has length alt, it ALREADY INCLUDES boundary layer info in the 
                       1st and last elements. 2nd and penultimate elements are for the edge bulk layers.
        atmdict: Atmospheric state dictionary, required for the nonthermal escape boundary condition.
        M: total atmospheric density, required for the nonthermal escape boundary condition.
    Outputs:
        boundary conditions for species in a 2 x 2 matrix, format:
        [n_1 -> n_0, n_0 -> n_1;      
         n_(nl) -> n_(nl+1), n_(nl+1) -> n_(n_l)]

        where n_0 is the boundary layer from [-1 km, 1 km], n_1 is the first bulk layer from [1 km, 3 km],
        n_(nl) is the topmost bulk layer, and n_(nl+1) is the top boundary layer.

        Form of the output is:

         Surface [↓, ↑;     [density-dependent, density-independent;    [#/s, #/cm³/s.;
         Top      ↑, ↓]      density-dependent, density-independent]     #/s, #/cm³/s]

        Each row has two elements:
            1st element: n_bulk  -> NULL (depends on species concentration in bulk layer)
            2nd element: NULL -> n_bulk (independent of species concentration in bulk layer)

            note, technically, these are chemical equations.

        More specifically, when the return value of this function is used in other functions, the first
        element in each row will eventually be multiplied by a density taken from the atmospheric 
        state dictionary, and the second element will be used as-is. That way, eventually the total
        change recorded in other functions is always #/cm³/s. 
    =#
    
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :speciesbclist, :dz])
    
    bc_dict = Dict{Symbol, Array{ftype_ncur}}([s=>[0 0; 0 0] for s in GV.all_species])

    for sp in keys(GV.speciesbclist)
        try 
            global these_bcs = GV.speciesbclist[sp]
        catch KeyError
            println("No entry $(sp) in bcdict")
            continue
        end
 
        # DENSITY
        try 
            n_lower = [fluxcoef_dict[sp][2, :][1], fluxcoef_dict[sp][1, :][2]*these_bcs["n"][1]]
            try
                @assert all(x->!isnan(x), n_lower)
                bc_dict[sp][1, :] .+= n_lower
            catch y
                if !isa(y, AssertionError)
                    throw("Unhandled exception in lower density bc: $(y)")
                end
            end
            try 
                n_upper = [fluxcoef_dict[sp][end-1, :][2], fluxcoef_dict[sp][end, :][1]*these_bcs["n"][2]]
                @assert all(x->!isnan(x), n_upper)
                bc_dict[sp][2, :] .+= n_upper
            catch y
                if !isa(y, AssertionError)
                    throw("Unhandled exception in upper density bc: $(y)")
                end
            end
        catch y
            if !isa(y, KeyError)
                throw("Unhandled exception in density bcs for $(sp): $(y)")
            end
        end
    
        # FLUX 
        try 
            f_lower = [0, -these_bcs["f"][1]/GV.dz]
            try        
                @assert all(x->!isnan(x), f_lower)
                bc_dict[sp][1, :] .+= f_lower
            catch y
                if !isa(y, AssertionError)
                    throw("Unhandled exception in lower flux bc: $(y)")
                end
            end
            try 
                f_upper = [0, -these_bcs["f"][2]/GV.dz]
                @assert all(x->!isnan(x), f_upper)
                bc_dict[sp][2, :] .+= f_upper
            catch y
                if !isa(y, AssertionError)
                    throw("Unhandled exception in upper flux bc: $(y)")
                end
            end
        catch y
            if !isa(y, KeyError)
                throw("Unhandled exception in flux bcs for $(sp)")
            end
        end
    
        # VELOCITY
        try 
            v_lower = [these_bcs["v"][1]/GV.dz, 0]

            try
                @assert all(x->!isnan(x), v_lower)
                bc_dict[sp][1, :] .+= v_lower
            catch y
                if !isa(y, AssertionError)
                    throw("Unhandled exception in lower velocity bc: $(y)")
                end
            end

            try 
                v_upper = [these_bcs["v"][2]/GV.dz, 0]
                @assert all(x->!isnan(x), v_upper)
                bc_dict[sp][2, :] .+= v_upper
            catch y
                if !isa(y, AssertionError)
                    throw("Unhandled exception in lower velocity bc: $(y)")
                end
            end
        catch y
            if !isa(y, KeyError)
                throw("Unhandled exception in velocity bcs for $(sp)")
            end
        end 
    end
    
    # SPECIAL CASE: add on the non-thermal escape for H and D. 
    if nonthermal
        @assert all(x->x in keys(GV), [:hot_H_network, :hot_D_network, :hot_H_rc_funcs, :hot_D_rc_funcs, 
                                       :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Jratedict])
        prod_hotH = escaping_hot_atom_production(:H, GV.hot_H_network, GV.hot_H_rc_funcs, atmdict, M; globvars...)
        prod_hotD = escaping_hot_atom_production(:D, GV.hot_D_network, GV.hot_D_rc_funcs, atmdict, M; globvars...)
        prod_hotH2 = escaping_hot_atom_production(:H2, GV.hot_H2_network, GV.hot_H2_rc_funcs, atmdict, M; globvars...)
        prod_hotHD = escaping_hot_atom_production(:HD, GV.hot_HD_network, GV.hot_HD_rc_funcs, atmdict, M; globvars...)

        # DIAGNOSTIC: produced hot H
        # if :results_dir in keys(GV)
        #     fig, ax = subplots()
        #     plot(prod_hotH, GV.plot_grid)
        #     xlabel("production rate")
        #     ylabel("altitude")
        #     xscale("log")
        #     xlim(left=1e-10)
        #     savefig(GV.results_dir*GV.sim_folder_name*"/prod_hotH.png")
        #     close(fig)
        # end

        bc_dict[:H][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_H_network, prod_hotH; returntype="number", globvars...)]
        bc_dict[:D][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_D_network, prod_hotD; returntype="number", globvars...)]
        bc_dict[:H2][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_H2_network, prod_hotH2; returntype="number", globvars...)]
        bc_dict[:HD][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_HD_network, prod_hotHD; returntype="number", globvars...)]
    end 
    return bc_dict
end

# Thermal escape: 
function effusion_velocity(Texo, m; globvars...)
    #=
    Returns effusion velocity for a species in cm/s

    Inputs:
        Texo: temperature of the exobase (upper boundary) in K
        m: mass of one molecule of species in amu
        zmax: max altitude in cm
    Outputs:
        v: effusion velocity for species of mass m 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:zmax])
    
    # lambda is the Jeans parameter (Gronoff 2020), basically the ratio of the 
    # escape velocity GmM/z to the thermal energy, kT.
    lambda = (m*mH*bigG*marsM)/(kB*Texo*(radiusM+GV.zmax))
    vth = sqrt(2*kB*Texo/(m*mH))
    v = exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)

    return v
end

# Nonthermal escape functions: 
function escape_probability(sp, atmdict; globvars...)::Array
    #=
    Returns an exponential profile of escape probability by altitude that accounts for collisions with the background 
    atmosphere. from Bethan Gregory, A and a for H. Could be redone for D, possibly.
    Input
        sp: species escaping (H or D, generally)
        atmdict: Atmospheric state dictionary
    Output
        Array by altitude of escape probabilities for hot atoms. 0-1.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :collision_xsect])
    
    A = 0.916 # escape probability at altitude where above column = 0, for high energy particles. upper limit
    a = 0.039 # how "transparent" the atmosphere is to an escaping atom. smaller for higher energy so this is for an upper limit.
    
    return A .* exp.(-a .* GV.collision_xsect[sp] .* column_density_above(n_tot(atmdict; GV.all_species))) 
end

function escaping_hot_atom_production(sp, source_rxns, source_rxn_rc_funcs, atmdict, Mtot; returntype="array", globvars...)
    #=
    Solves the equation k[R1][R2] * P to get the total volume escape of hot atoms of species sp
    from the exobase region where P is the escape probability.
    
    Input
        sp: species
        source_rxns: reaction network that will cause hot atoms to be produced
        atmdict: present atmospheric state dictionary
        Mtot: total atmospheric density array
    Output: 
        array of production by altitude (rows) and reaction  (columns)
    =#
    
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :collision_xsect, :ion_species, :Jratedict, :molmass, :non_bdy_layers, :num_layers,  
                                   :n_alt_index, :Tn, :Ti, :Te, :dz, :zmax])

    produced_hot = volume_rate_wrapper(sp, source_rxns, source_rxn_rc_funcs, atmdict, Mtot; returntype="array", globvars...) 

    # Returns an array where rows represent altitudes and columns are reactions. Multiplies each vertical profile (each column) by escape_probability. 
    if returntype=="array" # Used within the code to easily calculate the total flux later on. 
        return produced_hot .* escape_probability(sp, atmdict; globvars...)
    elseif returntype=="df" # Useful if you want to look at the arrays yourself.
        return DataFrame(produced_hot .* escape_probability(sp, atmdict; globvars...), vec([format_chemistry_string(r[1], r[2]) for r in source_rxns]))
    end
end

function nonthermal_escape_flux(source_rxn_network, prod_rates_by_alt; verbose=false, returntype="dataframe", globvars...) 
    #=
    Given a matrix where each column is a vertical profile of the production rates (#/cm³/s) of escaping hot atoms
    from a single reaction, this calculates the total effective flux by doing a simple sum * dz. 
    Then it puts the flux into a dataframe so it is easier to sort, and it returns either the sorted dataframe or 
    the collapsed sum (i.e. the total flux across all the given reactions). Species is not passed because this depends
    on prod_rates_by_alt which does have species implicit in it.
    
    Input:
        source_rxn_network: List of chemical reactions producing either H or D that are hot and liable to escape
        prod_rates_by_alt: Matrix where each column is a vertical profile of the production rates (#/cm³/s) of escaping hot atoms
                            from a single reaction
    Output:
        Dataframe of total fluxes by reaction, sorted by dominance, or a simple number which is the total flux
        for whatever species is calculated for.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:dz])

    # Get a vector of strings that contain the chemical reactions
    rxn_strings = vec([format_chemistry_string(r[1], r[2]) for r in source_rxn_network])
    
    # Calculate the column rate for each reaction 
    sum_all_alts = sum(prod_rates_by_alt .* GV.dz, dims=1)
    
    # Convert to a dataframe because it's convenient to sort
    df = DataFrame("Rxn"=>rxn_strings, "Value"=>vec(sum_all_alts))
    sorted_total_esc_by_rxn = sort(df, [:Value], rev=true)

    if verbose
        println("Total hot atoms from all reactions: $(sum(sorted_total_esc_by_rxn.Value))" )
    end
    
    if returntype=="dataframe"
        return sorted_total_esc_by_rxn
    elseif returntype=="number"
        return sum(sorted_total_esc_by_rxn.Value)
    end
end

# Other fluxes and such:
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

function get_flux(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; nonthermal=true, globvars...)
    #=
    NEW VERSION : THIS IS THE BETTER VERSION NOW! But only for fluxes.
    
    Input:
        atmdict: Array; species number density by altitude
        sp: Symbol
        Tn, Ti, Te, Tp: Temperature arrays (neutral, ion, electron, plasma)
        bcdict: the boundary condition dictionary.

    Output: 
        Array of flux values (#/cm²/s) at each atmospheric layer boundary.
        i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
        and the end of the array is the boundary at 249 km.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :speciesbclist, :dz, :Hs_dict, :molmass, :neutral_species, :num_layers, :n_all_layers, :n_alt_index, 
                                    :polarizability, :q, :Tn, :Ti, :Te, :Tp, :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species])
    
    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(GV.Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, D_arr; globvars...)
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict; globvars...)
    bc_dict = boundaryconditions(fluxcoefs_all, atmdict, sum([atmdict[sp] for sp in GV.all_species]); nonthermal=nonthermal, globvars...)

    # each element in bulk_layer_coefs has the format [downward flow (i to i-1), upward flow (i to i+1)].  units 1/s
    bulk_layer_coefs = fluxcoefs_all[sp][2:end-1, :]

    bcs = bc_dict[sp]
    
    net_bulk_flow = fill(convert(ftype_ncur, NaN), GV.n_all_layers-1)  # units #/cm^3/s; tracks the cell boundaries, of which there are length(alt)-1

    # We will calculate the net flux across each boundary, with sign indicating direction of travel.
    # Units for net bulk flow are always: #/cm³/s. 
    # NOTE: This might not actually represent the flow correctly, because I was assuming 
    # that the 1st bc was into the layer, and the 2nd was out, but it's actually just about in/dependence on density.
    net_bulk_flow[1] = (bcs[1, 2]                  # increase of the lowest atmospheric layer's density. 0 unless the species has a density or flux condition
                       - atmdict[sp][1]*bcs[1, 1]) # lowest atmospheric layer --> surface ("depositional" term). UNITS: #/cm³/s. 
                        
    for ialt in 2:GV.num_layers  # now iterate through every cell boundary within the atmosphere. boundaries at 3 km, 5...247. 123 elements.
        # UNITS for both of these terms:  #/cm³/s. 
        net_bulk_flow[ialt] = (atmdict[sp][ialt-1]*bulk_layer_coefs[ialt-1, 2]   # coming up from below: cell i-1 to cell i. Should be positive * positive
                              - atmdict[sp][ialt]*bulk_layer_coefs[ialt, 1])     # leaving to the layer below: downwards: cell i to cell i-1
    end

    # now the top boundary - between 124th atmospheric cell (alt = 249 km)
    net_bulk_flow[end] = (atmdict[sp][end]*bcs[2, 1] # into exosphere from the cell. UNITS: #/cm³/s. 
                         - bcs[2, 2]) # into top layer from exosphere. negative because the value in bcs is negative. do not question this. UNITS: #/cm³/s. 
                
    return net_bulk_flow .* GV.dz # now it is a flux. hurrah.
end

function limiting_flux(sp, atmdict, T_arr; globvars...)
    #=
    Calculate the limiting upward flux (Hunten, 1974; Zahnle, 2008). 
    Inputs:
        sp: A species that is traveling upwards
        atmdict: present atmospheric state
        T_arr: Array of neutral temperatures
    Output:
        Φ, limiting flux for a hydrostatic atmosphere
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :non_bdy_layers, :molmass, :n_alt_index])
    
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s
    fi = atmdict[sp] ./ n_tot(atmdict; globvars...)
    Ha = scaleH(atmdict, T_arr; globvars..., alt=GV.non_bdy_layers)
    bi = binary_dcoeff_inCO2(sp, T_arr)

    mass_ratio = GV.molmass[sp] / meanmass(atmdict; globvars...) 

    if (all(m->m<0.1, mass_ratio)) & (all(f->f<0.0001, fi)) # Light minor species approximation
        return bi .* fi ./ Ha
    else # Any species
        D = Dcoef_neutrals(non_bdy_layers, sp, bi, atmdict; globvars...)    
        return (D .* atmdict[sp] ./ Ha) .* (1 .- GV.molmass[sp] ./ meanmass(atmdict; globvars...))
    end
end

function limiting_flux_molef(sp, atmdict, T_arr; globvars...)
    #=
    Roger requested the limiting flux in in mole fraction. This is actually the same result as above. But this way we're sure
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :non_bdy_layers, :molmass, :n_alt_index])

    avogadro = 6.022e23

    X = (atmdict[sp] ./ avogadro) ./ (n_tot(atmdict; globvars...) ./ avogadro)
    # Calculate some common things: mixing ratio, scale height, binary diffusion coefficient AT^s

    Ha = scaleH(atmdict, T_arr; globvars..., alt=GV.non_bdy_layers)
    bi = binary_dcoeff_inCO2(sp, T_arr)
    Hi = scaleH(non_bdy_layers, sp, T_arr; globvars...)

    return bi .* X .* (1 ./ Ha - 1 ./ Hi), X
end

# **************************************************************************** #
#                                                                              #
#                        Vertical Transport Functions                          #
#                                                                              #
# **************************************************************************** #
#=

TRANSPORT DESCRIPTION 

at each level of the atmosphere, density can be transferred up or
down with a specified rate coefficient.

                         | n_i+1
      ^  tspecies_i_up   v tspecies_i+1_down
  n_i |
      v tspecies_i_down  ^ tspecies_i-1_up
                         | n_i-1

the flux at each cell boundary is the sum of the upward flux from
the cell below and the downward flux of the cell above. These fluxes
are determined using flux coefficients that come from the diffusion
equation. Care must be taken at the upper and lower boundary so that
tspecies_top_up and tspecies_bottom_down properly reflect the
boundary conditions of the atmosphere.

This is handled in the code with the population of the appropriate
reactions, with variable rate coeffecients that are populated
between each timestep (similar to the way photolysis rates are
included). We need to make reactions at each interior altitude
level:
         n_i -> n_i+1  tspecies_i_up
         n_i -> n_i-1  tspecies_i_down

At the upper and lower boundary we omit the species on the RHS, so
that these reactions are potentially non-conservative:

        n_top    -> NULL  tspecies_top_up
        n_bottom -> NULL  tspecies_bottom_down

These coefficients then describe the diffusion velocity at the top
and bottom of the atmosphere.
=#

function binary_dcoeff_inCO2(sp, T)
    #=
    Calculate the bindary diffusion coefficient for species sp.

    Currently, this is set up to only work for diffusion through CO2 since that's the Mars atm.
    Could be extended to be for any gas, but that will require some work.
    =#
    return diffparams(sp)[1] .* 1e17 .* T .^ (diffparams(sp)[2])
end

function Dcoef_neutrals(z, sp::Symbol, b, atmdict::Dict{Symbol, Vector{ftype_ncur}}; globvars...)
    #=
    Calculate the basic diffusion coefficient, AT^s/n.
    Inputs:
        z: An altitude or array of altitudes in cm.
        sp: Species for which to calculate.
        T: Temperature in K or array of temperatures.
        atmdict: Present atmospheric state.
    Outputs: 
        D: Diffusion coefficient AT^s/n

    Usable at either a specific altitude or all altitudes (array format). z and T must be same type.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :n_alt_index])

    if (typeof(z)==Float64) & (typeof(b)==Float64)
        return b ./ n_tot(atmdict, z; GV.all_species, GV.n_alt_index)
    else 
        return b ./ n_tot(atmdict; GV.all_species, GV.n_alt_index)
    end
end

function Dcoef!(D_arr, T_arr, sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; globvars...) 
    #=
    Calculates the molecular diffusion coefficient for an atmospheric layer.
    For neutrals, returns D = AT^s/n, from Banks and Kockarts Aeronomy, part B, pg 41, eqn 
    15.30 and table 15.2 footnote.
    For ions, it returns ambipolar diffusion coefficients according to Krasnopolsky 2002 and 
    Schunk & Nagy equation 5.55. Yes, the equation is done in one line and it's ugly, but it works.
    Units: cm/s

    Inputs:
        D_arr: the container for the diffusion coefficients for ONE species.
        T_arr: temperature (K). Neutral temp for neutrals, plasma temp for ions.
        sp: whichever species we are calculating for
        atmdict: state of the atmosphere; should include boundary layers, i.e. be the result of calling atmdict_with_boundary_layers.
        bcdict: Boundary conditions dictionary specified in parameters file
    Outputs:
        D_arr: An array of the diffusion coefficients by altitude for species
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :molmass, :neutral_species, :n_alt_index, :polarizability, :q, :speciesbclist])
   
    # Calculate as if it was a neutral - not using function above because this is faster than going into 
    # the function and using an if/else block since we know we'll always have vectors in this case.
    D_arr[:] .= (binary_dcoeff_inCO2(sp, T_arr)) ./ n_tot(atmdict; GV.all_species, GV.n_alt_index)

    # If an ion, overwrite with the ambipolar diffusion
    if charge_type(sp) == "ion"
        # a place to store the density array and nu_in
        species_density = zeros(size(T_arr))
        sum_nu_in = zeros(size(T_arr))

        # mi = GV.molmass[sp] .* mH
        # create the sum of nu_in. Note that this depends on density, but we only have density for the real layers,
        # so we have to assume the density at the boundary layers is the same as at the real layers.
        for n in GV.neutral_species
            species_density = atmdict[n]

            # This sets the species density to a boundary condition if it exists. 
            if haskey(GV.speciesbclist, n)
                if haskey(GV.speciesbclist[n], "n") 
                    if !isnan(GV.speciesbclist[n]["n"][1])
                        species_density[1] = GV.speciesbclist[n]["n"][1]
                    end
                    if !isnan(GV.speciesbclist[n]["n"][2]) # currently this should never apply.
                        species_density[end] = GV.speciesbclist[n]["n"][2]
                    end
                end
            end
            
            sum_nu_in .+= 2 .* pi .* (((GV.polarizability[n] .* GV.q .^ 2) ./ reduced_mass(GV.molmass[sp], GV.molmass[n])) .^ 0.5) .* species_density

        end
        
        D_arr .= (kB .* T_arr) ./ (GV.molmass[sp] .* mH .* sum_nu_in)

    end
    return D_arr
end

#=
molecular diffusion parameters. value[1] = A, value[2] = s in the equation
D = AT^s / n given by Banks & Kockarts Aeronomy part B eqn. 15.30 and Hunten
1973, Table 1.

molecular diffusion is different only for small molecules and atoms
(H, D, HD, and H2), otherwise all species share the same values (Krasnopolsky
1993 <- Hunten 1973; Kras cites Banks & Kockarts, but this is actually an
incorrect citation.)

D and HD params are estimated by using Banks & Kockarts eqn 15.29 (the
coefficient on T^0.5/n) to calculate A for H and H2 and D and HD, then using
(A_D/A_H)_{b+k} = (A_D/A_H)_{hunten} to determine (A_D)_{hunten} since Hunten
1973 values are HALF what the calculation provides.

s, the power of T, is not calculated because there is no instruction on how to
do so and is assumed the same as for the hydrogenated species.
=#  
diffparams(s) = get(Dict(:H=>[8.4, 0.597], :H2=>[2.23, 0.75],
                             :D=>[5.98, 0.597], :HD=>[1.84, 0.75],
                             :Hpl=>[8.4, 0.597], :H2pl=>[2.23, 0.75],
                             :Dpl=>[5.98, 0.597], :HDpl=>[1.84, 0.75]),
                        s,[1.0, 0.75])

function fluxcoefs(sp::Symbol, Kv, Dv, H0v; globvars...)
    #= 
    base function to generate flux coefficients of the transport network. 
    
    For all the arrays, length = num_layers 

    Inputs:
        sp: species symbol 
        z: altitude array in cm.
        dz: altitude layer thickness ("resolution")
        for all the following, length = num_layers 
        Kv: eddy diffusion coefficient
        Dv: molecular diffusion coefficient
        Tv_n: neutral temperature
        Tv_p: plasma temperature
        Hsv: scale height by species
        H0v: mean atmospheric scale height
    Outputs:
        Arrays of coefficients (units 1/s) at each atmospheric layer for downward and upward flux.
        Note that even though it's defined as being between a layer and the one above or below, the value is 
        evaluated at the center of the layer 

    v just means "vector"
    u refers to "upper" (a given layer coupled to the one above)
    l refers to "lower" (a given layer coupled to the one below)
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Tp, :Hs_dict, :n_all_layers, :dz])

    # Initialize arrays for downward (i to i-1) and upward (i to i+1) coefficients
    Dl = zeros(GV.n_all_layers)
    Kl = zeros(GV.n_all_layers)
    Tl_n = zeros(GV.n_all_layers)
    Tl_p = zeros(GV.n_all_layers)
    dTdzl_n = zeros(GV.n_all_layers)
    dTdzl_p = zeros(GV.n_all_layers)
    Hsl = zeros(GV.n_all_layers)
    H0l = zeros(GV.n_all_layers)

    Du = zeros(GV.n_all_layers)
    Ku = zeros(GV.n_all_layers)
    Tu_n = zeros(GV.n_all_layers)
    Tu_p = zeros(GV.n_all_layers)
    dTdzu_n = zeros(GV.n_all_layers)
    dTdzu_p = zeros(GV.n_all_layers)
    Hsu = zeros(GV.n_all_layers)
    H0u = zeros(GV.n_all_layers)

    # Calculate the coefficients between this layer and the lower layer. 
    Dl[2:end] = @. (Dv[sp][1:end-1] + Dv[sp][2:end]) /  2.0
    Kl[2:end] = @. (Kv[1:end-1] + Kv[2:end]) / 2.0
    Tl_n[2:end] = @. (GV.Tn[1:end-1] + GV.Tn[2:end]) / 2.0
    Tl_p[2:end] = @. (GV.Tp[1:end-1] + GV.Tp[2:end]) / 2.0
    dTdzl_n[2:end] = @. (GV.Tn[2:end] - GV.Tn[1:end-1]) / GV.dz
    dTdzl_p[2:end] = @. (GV.Tp[2:end] - GV.Tp[1:end-1]) / GV.dz
    Hsl[2:end] = @. (GV.Hs_dict[sp][1:end-1] + GV.Hs_dict[sp][2:end]) / 2.0
    H0l[2:end] = @. (H0v[charge_type(sp)][1:end-1] + H0v[charge_type(sp)][2:end]) / 2.0

    # Handle the lower boundary layer:
    Dl[1] = @. (1 + Dv[sp][1]) /  2.0
    Kl[1] = @. (1 + Kv[1]) / 2.0
    Tl_n[1] = @. (1 + GV.Tn[1]) / 2.0
    Tl_p[1] = @. (1 + GV.Tp[1]) / 2.0
    dTdzl_n[1] = @. (GV.Tn[1] - 1) / GV.dz
    dTdzl_p[1] = @. (GV.Tp[1] - 1) / GV.dz
    Hsl[1] = @. (1 + GV.Hs_dict[sp][1]) / 2.0
    H0l[1] = @. (1 + H0v[charge_type(sp)][1]) / 2.0

    # Now the coefficients between this layer and upper layer
    Du[1:end-1] = @. (Dv[sp][1:end-1] + Dv[sp][2:end]) /  2.0
    Ku[1:end-1] = @. (Kv[1:end-1] + Kv[2:end]) / 2.0
    Tu_n[1:end-1] = @. (GV.Tn[1:end-1] + GV.Tn[2:end]) / 2.0
    Tu_p[1:end-1] = @. (GV.Tp[1:end-1] + GV.Tp[2:end]) / 2.0
    dTdzu_n[1:end-1] = @. (GV.Tn[2:end] - GV.Tn[1:end-1]) / GV.dz
    dTdzu_p[1:end-1] = @. (GV.Tp[2:end] - GV.Tp[1:end-1]) / GV.dz
    Hsu[1:end-1] = @. (GV.Hs_dict[sp][1:end-1] + GV.Hs_dict[sp][2:end]) / 2.0
    H0u[1:end-1] = @. (H0v[charge_type(sp)][1:end-1] + H0v[charge_type(sp)][2:end]) / 2.0

    # Handle upper boundary layer:
    Du[end] = @. (Dv[sp][end] + 1) /  2.0
    Ku[end] = @. (Kv[end] + 1) / 2.0
    Tu_n[end] = @. (GV.Tn[end] + 1) / 2.0
    Tu_p[end] = @. (GV.Tp[end] + 1) / 2.0
    dTdzu_n[end] = @. (1 - GV.Tn[end]) / GV.dz
    dTdzu_p[end] = @. (1 - GV.Tp[end]) / GV.dz
    Hsu[end] = @. (GV.Hs_dict[sp][end] + 1) / 2.0
    H0u[end] = @. (H0v[charge_type(sp)][end] + 1) / 2.0


    # two flux terms: eddy diffusion and gravity/thermal diffusion.
    # these are found in line 5 of Mike's transport_as_chemistry.pdf:
    # sumeddy = (D+K)/(Δz²), gravthermal = ☐/(2Δz), where ☐ = {D(1/H + 1+(α/T)(dT/dz)) + K(1/H_H + (1/T)(dT/dz))}
    sumeddyl = @. (Dl+Kl)/GV.dz/GV.dz
    if charge_type(sp) == "neutral"
        gravthermall = @. (Dl*((1/Hsl) + ((1+thermaldiff(sp))/Tl_n)*dTdzl_n) +
                        Kl*((1/H0l) + (1/Tl_n)*dTdzl_n))/(2*GV.dz)
    elseif charge_type(sp) == "ion"
        gravthermall = @. (Dl*((1/Hsl) + ((1+thermaldiff(sp))/Tl_p)*dTdzl_p) +
                        Kl*((1/H0l) + (1/Tl_n)*dTdzl_n))/(2*GV.dz)
    elseif charge_type(sp) == "electron"
        throw("Electrons not handled as individual species")
    end

    sumeddyu = @. (Du+Ku)/GV.dz/GV.dz  # this is the line where we divide by cm^2
    if charge_type(sp) == "neutral"
        gravthermalu = @. (Du*((1/Hsu) + ((1 + thermaldiff(sp))/Tu_n)*dTdzu_n) +
                        Ku*((1/H0u) + (1/Tu_n)*dTdzu_n))/(2*GV.dz)
    elseif charge_type(sp) == "ion"
        gravthermalu = @. (Du*((1/Hsu) + ((1 + thermaldiff(sp))/Tu_p)*dTdzu_p) +
                        Ku*((1/H0u) + (1/Tu_n)*dTdzu_n))/(2*GV.dz)
    elseif charge_type(sp) == "electron"
        throw("Electrons not handled as individual species")
    end
    
    # this results in the following coupling coefficients; sumeddy + gravthermal = (D+K)/(Δz²) + ☐/(2Δz), units 1/s <-----_!!!!! important
    # first return is this term between layer i and i-1 for whole atmosphere.
    # second return is between layer i and i+1
    return sumeddyl .+ gravthermall,  # down
            sumeddyu .- gravthermalu # up; negative because gravity points down. I think that's why.
end

function fluxcoefs(species_list::Vector, K, D, H0; globvars...) 
    #=
    New optimized version of fluxcoefs that calls the lower level version of fluxcoefs,
    producing a dictionary that contains both up and down flux coefficients for each layer of
    the atmosphere including boundary layers. Created to optimize calls to this function
    during the solution of the production and loss equation.

    Here, D and Hs depend on the current atmospheric densities, and need to be pre-calculated
    within the upper level function which calls this one.
    The parameters below which vary by species are dictionaries, and those that are arrays
    don't depend on the species. All profiles are by altitude. All lengths are the same 
    as for the alt variable (full altitude grid including boundary layers).
    
    Inputs:
        species_list: Species for which to generate transport coefficients. This allows the code to only do it for
                transport species during the main simulation run, and for all species when trying to plot 
                rate balances after the run.
        T_neutral: 1D neutral temperature profile
        T_plasma: the same, but for the plasma temperature
        K: Array; 1D eddy diffusion profile by altitude for current atmospheric state
        D: Dictionary (key=species); 1D molecular diffusion profiles for current atmospheric state
        H0: Dictionary (key="neutral" or "ion"); 1D mean atmospheric scale height profiles for each type
        Hs: Dictionary (key=species); 1D species scale height profiles
    Outputs:
        fluxcoef_dict: dictionary of flux coefficients of the form [flux down, flux up] by altitude 

    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Tp, :Hs_dict, :n_all_layers, :dz])
    
    # the return dictionary: Each species has 2 entries for every layer of the atmosphere.
    fluxcoef_dict = Dict{Symbol, Array{ftype_ncur}}([s=>fill(0., GV.n_all_layers, 2) for s in species_list])

    for s in species_list
        layer_below_coefs, layer_above_coefs = fluxcoefs(s, K, D, H0; globvars...)
        fluxcoef_dict[s][:, 1] .= layer_below_coefs
        fluxcoef_dict[s][:, 2] .= layer_above_coefs
    end

    return fluxcoef_dict
end

function get_transport_PandL_rate(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; returnfluxes=false, nonthermal=true, globvars...)
    #=
    Input:
        sp: species for which to return the transport production and loss
        atmdict: species number density by altitude
        Tn, Ti, Te, Tp: Temperature arrays
        bcdict: Boundary conditions dictionary specified in parameters file
    Output
        Array of production and loss (#/cm³/s) at each atmospheric layer boundary.
        i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
        and the end of the array is the boundary at 249 km.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :dz, :Hs_dict, :molmass, :n_all_layers, :n_alt_index, 
                                   :neutral_species, :num_layers, :polarizability, :q, :speciesbclist, :Te, :Ti, :Tn, :Tp, 
                                   :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species])

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

function Keddy(z::Vector, nt::Vector)
    #=
    Input:
        z: Altitude in cm
        nt: Total atmospheric density
    Ouptut:
        k: eddy diffusion coefficients at all altitudes.
    =#

    k = zeros(size(z))
    upperatm = findall(i->i .> 60e5, z)
    k[findall(i->i .<= 60e5, z)] .= 10. ^ 6
    k[upperatm] .= 2e13 ./ sqrt.(nt[upperatm])

    return k
end

# thermal diffusion factors
thermaldiff(sp) = get(Dict(:H=>-0.25, :H2=>-0.25, :D=>-0.25, :HD=>-0.25,
                                :He=>-0.25, 
                                :Hpl=>-0.25, :H2pl=>-0.25, :Dpl=>-0.25, :HDpl=>-0.25,
                                :Hepl=>-0.25), sp, 0)

function update_diffusion_and_scaleH(species_list, atmdict::Dict{Symbol, Vector{ftype_ncur}}, D_coefs; globvars...) 
    #=
    Input:
        atmdict: Atmospheric state dictionary without boundary layers
        Tn: Neutral temperature profile
        Tp: Plasma temperature profile
        D_coefs: placeholder array for diffusion coefficients to speed up performance.
        bcdict: boundary condition dictionary
        species_list: Species for which to generate molecular diffusion coefficients. This allows the code to only do it for
                 transport species during the main simulation run, and for all species when trying to plot 
                 rate balances after the run.
    Output:
        K: Vector of eddy diffusion coefficient by altitude. Independent of species.
        Dcoefs: Dictionary of molecular diffusion by altitude. Keys are species: species=>[D by altitude] 
        H0: Dictionary of mean atmospheric scale height by altitude. Keys are "neutral" and "ion". 
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :speciesbclist, :molmass, :neutral_species, :n_alt_index, :polarizability, :q,
                                   :Tn, :Tp, :Tprof_for_diffusion])

    ncur_with_bdys = ncur_with_boundary_layers(atmdict; GV.n_alt_index, GV.all_species)
    
    K = Keddy(GV.alt, n_tot(ncur_with_bdys; GV.all_species, GV.n_alt_index))
    H0_dict = Dict{String, Vector{ftype_ncur}}("neutral"=>scaleH(ncur_with_bdys, GV.Tn; globvars...),
                                               "ion"=>scaleH(ncur_with_bdys, GV.Tp; globvars...))
    
    # Molecular diffusion is only needed for transport species, though.  
    Dcoef_dict = Dict{Symbol, Vector{ftype_ncur}}([s=>deepcopy(Dcoef!(D_coefs, GV.Tprof_for_diffusion[charge_type(s)], s, ncur_with_bdys; globvars...)) for s in species_list])

    return K, H0_dict, Dcoef_dict
end

function update_transport_coefficients(species_list, atmdict::Dict{Symbol, Vector{ftype_ncur}}, D_coefs, M; 
                                       calc_nonthermal=true, globvars...) 
    #=
    Input:
        species_list: Species which will have transport coefficients updated
        atmdict: Atmospheric state dictionary for bulk layers
        tspecies: species which need transport coefficients calculated. May vary depending on sim parameters.
        Tn, Tp: Neutrals and plasma temperatures
        D_coefs: placeholder for molecular diffusion coefficients
        bcdict: Dictionary of boundary conditions, needed for pretty much everything.
        species_list: Species for which to generate molecular diffusion coefficients. This allows the code to only do it for
                 transport species during the main simulation run, and for all species when trying to plot 
                 rate balances after the run.

    Return: 
        Transport coefficients for all atmospheric layers, units 1/s

        tlower: transport coefficients at the lower boundary layer. shape: # TODO
        tup: upward coefficients at inner bulk layers. shape: #TODO 
        tdown: downward coefficients at inner bulk layers. shape: #TODO 
        tupper: at the upper boundary layer. shape: #TODO 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :speciesbclist, :dz, :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs, 
                                   :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs,  
                                   :Hs_dict, :ion_species, :molmass, :neutral_species, :non_bdy_layers, :num_layers, :n_all_layers, :n_alt_index, 
                                   :polarizability, :q, :Tn, :Ti, :Te, :Tp, :Tprof_for_diffusion, :transport_species, :zmax])
    
    # Update the diffusion coefficients and scale heights
    K_eddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(species_list, atmdict, D_coefs; globvars...)

    # Get flux coefficients
    fluxcoefs_all = fluxcoefs(species_list, K_eddy_arr, Dcoef_dict, H0_dict; globvars...)
    
    # Transport coefficients, non-boundary layers
    tup = fill(-999., length(GV.transport_species), GV.num_layers)
    tdown = fill(-999., length(GV.transport_species), GV.num_layers)
    for (i, s) in enumerate(GV.transport_species)
        tup[i, :] .= fluxcoefs_all[s][2:end-1, 2]
        tdown[i, :] .= fluxcoefs_all[s][2:end-1, 1]
    end

    bc_dict = boundaryconditions(fluxcoefs_all, atmdict, M; nonthermal=calc_nonthermal, globvars...)

    # transport coefficients for boundary layers
    tlower = permutedims(reduce(hcat, [bc_dict[sp][1,:] for sp in GV.transport_species]))
    tupper = permutedims(reduce(hcat, [bc_dict[sp][2,:] for sp in GV.transport_species]))

    return tlower, tup, tdown, tupper
end
