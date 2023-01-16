# **************************************************************************** #
#                                                                              #
#         Functions for atmospheric dictionary object or its constituents      #
#                                                                              #
# **************************************************************************** #

function atm_dict_to_matrix(atmdict::Dict{Symbol, Vector{ftype_ncur}}, species_list)
    #=
    Converts atmospheric state dictionary atmdict to a matrix,
    such that rows correspond to species in the order listed in species_list
    and columns correspond to altitudes in the order lowest-->highest.
    =#
    
    num_alts = length(atmdict[collect(keys(atmdict))[1]])
    n_mat = zeros(length(species_list), num_alts)
    
    for i in 1:length(species_list)
        n_mat[i, :] = atmdict[species_list[i]]
    end
    
    return n_mat
end

function atm_matrix_to_dict(n_matrix, species_list)
    #=
    Input:
        n_matrix: matrix of the atmospheric state
    Output:
        dictionary for only the species in species_list
    =#
    atmdict = Dict{Symbol, Vector{ftype_ncur}}([species_list[k]=>n_matrix[k, :] for k in 1:length(species_list)])
    
    return atmdict
end

function electron_density(atmdict; globvars...)
    #=
    Calculate the electron profile for the current atmospheric state. Usually used in 
    set up.
    Inputs:
        atmdict: Current atmospheric state
    Outputs:
        electron density array by altitude
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:e_profile_type, :ion_species, :non_bdy_layers])

    if GV.e_profile_type=="constant"
        E = [1e5 for i in GV.non_bdy_layers]
    elseif GV.e_profile_type=="quasineutral"
        E = sum([atmdict[sp] for sp in GV.ion_species])
    elseif GV.e_profile_type=="none"  # For neutrals-only simulation but without changing how E is passed to other functions. 
        E = [0. for i in GV.non_bdy_layers]
    else
        throw("Unhandled electron profile specification: $(e_profile_type)")
    end
    return E
end

function column_density(n::Vector; start_alt=1)
    #=
    Returns column density above a given atmospheric layer. 

    Input
        n: species number density (#/cm³) by altitude
    Output
        Column density (#/cm²)
    =#
    return sum(n[start_alt:end] .* dz)
end

function column_density_above(n_tot_by_alt::Vector)
    #=
    Returns an array where entries are the total integrated column density above
    that level of the atmosphere. e.g. the value at the topmost altitude is 
    called 0 since we assume anything beyond that level can escape. 
    
    n_tot_by_alt: Total atmospheric density at each altitude layer.
    =#
    col_above = zeros(size(n_tot_by_alt))

    for i in 1:num_layers
        col_above[i] = column_density(n_tot_by_alt, start_alt=i+1)
    end

    return col_above
end

function compile_ncur_all(n_long, n_short, n_inactive; globvars...)
    #=
    While the simulation runs, "n", the vector passed to the solver, only contains densities
    for long-lived, active species. Every time the atmospheric state changes, the transport coefficients
    and Jrates must be updated, but those all depend on the densities of ALL species. It's easiest for 
    the functions updating those things to pull from one atmospheric state dictionary, 
    so this function combines disparate density vectors back into one dictionary.

    Input:
        n_long: active, long-lived species densities
        n_short: same but for short-lived species
        n_inactive: inactive species densities (truly, these never change)
    Output:
        atmospheric state dictionary of species densities only (no Jrates).
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:active_longlived, :active_shortlived, :inactive_species, :num_layers])

    n_cur_active_long = unflatten_atm(n_long, GV.active_longlived; num_layers=GV.num_layers)
    n_cur_active_short = unflatten_atm(n_short, GV.active_shortlived; num_layers=GV.num_layers)
    n_cur_inactive = unflatten_atm(n_inactive, GV.inactive_species; num_layers=GV.num_layers)

    n_cur_all = Dict(vcat([k=>n_cur_active_long[k] for k in keys(n_cur_active_long)],
                          [k=>n_cur_active_short[k] for k in keys(n_cur_active_short)],
                          [k=>n_cur_inactive[k] for k in keys(n_cur_inactive)]))
    
    return n_cur_all
end

function find_exobase(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; returntype="index", verbose=false, globvars...)
    #=
    Finds the exobase altitude, where mean free path is equal to a scale height.

    Inputs:
        s: species 
        atmdict: Atmospheric state dictionary
        returntype: whether to return the "altitude" in km or the "index" in the n_alt_index dictionary. 
    Output:
        Altitude of exobase in cm
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:non_bdy_layers, :all_species, :Tn, :molmass, :alt, :collision_xsect, :n_alt_index, :zmax])

    H_s = scaleH(GV.non_bdy_layers, sp, GV.Tn[2:end-1]; globvars...)
    mfp_sp = 1 ./ (GV.collision_xsect[sp] .* n_tot(atmdict; GV.all_species, GV.n_alt_index))
    exobase_alt = findfirst(mfp_sp .>= H_s)

    if typeof(exobase_alt)==Nothing # If no exobase is found, use the top of the atmosphere.
        if verbose
            println("Warning: No exobase found for species $(sp); assuming top of atmosphere, but this is not guaranteed to be true.")
        end
        returnme = Dict("altitude"=>GV.zmax, "index"=>GV.n_alt_index[GV.zmax])
    else
        returnme = Dict("altitude"=>GV.alt[exobase_alt], "index"=>exobase_alt)
    end
    return returnme[returntype]
end

function flatten_atm(atmdict::Dict{Symbol, Vector{ftype_ncur}}, species_list; globvars...) 
    #=
    Input:
        atmdict: atmospheric densities by altitude
        species_list: Included species which will have profiles flattened
    Output:
        Vector of form [n_sp1(z=0), n_sp2(z=0)...n_sp1(z=zmax)...n_spN(z=zmax)]
    
    This function is the reverse of unflatten_atm. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:num_layers])

    return deepcopy(ftype_ncur[[atmdict[sp][ialt] for sp in species_list, ialt in 1:GV.num_layers]...])
end

function get_deuterated(sp_list; exclude=[:O1D, :Nup2D])
    #=
    Just returns the same list but containing only deuterated species.
    Inputs:
        sp_list: List to search for the D-bearing species.
        Optional:
            exclude: a list of species names that may be identified as D-bearing but really are not 
                     (usually because their name involves a D that represents an excited state).
    Output:
        The same list, but with only the D-bearing species.

    =#
    return [s for s in setdiff(sp_list, exclude) if occursin('D', string(s))];
end

function meanmass(atmdict::Dict{Symbol, Vector{ftype_ncur}}, z; globvars...)
    #= 
    find the mean molecular mass at a given altitude z

    atmdict: species number density by altitude
    z: float; altitude in atmosphere in cm

    return: mean molecular mass in amu
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:all_species, :molmass, :n_alt_index])

    thisaltindex = GV.n_alt_index[z]
    c = [atmdict[sp][thisaltindex] for sp in GV.all_species]
    m = [GV.molmass[sp] for sp in GV.all_species]
    return sum(c.*m)/sum(c)
end

function meanmass(atmdict::Dict{Symbol, Vector{ftype_ncur}}; globvars...)
    #= 
    Override for vector form. Calculates mean molecular mass at all atmospheric layers.

    atmdict: Array; species number density by altitude
    returns: mean molecular mass in amu for all atmospheric layers.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:all_species, :molmass, :n_alt_index])

    # Gets the atmosphere as a matrix with rows = altitudes and cols = species
    # so we can do matrix multiplication.
    n_mat = transpose(atm_dict_to_matrix(atmdict, GV.all_species))

    m = [GV.molmass[sp] for sp in GV.all_species] # this will always be 1D

    weighted_mm = zeros(size(n_mat)[1]) # This will store the result

    # Multiply densities of each species by appropriate molecular mass 
    mul!(weighted_mm, n_mat, m)

    return weighted_mm ./ n_tot(atmdict; GV.all_species, GV.n_alt_index)
end

function ncur_with_boundary_layers(atmdict_no_bdys::Dict{Symbol, Vector{ftype_ncur}}; globvars...)
    #=
    Here's a weird one. The atmospheric density matrix stores values for each
    species (column of the matrix) at each altitude (row of the matrix) of the atmosphere. 
    This means there are num_layers (124 for a 0-250 km grid) rows. But, a lot of our
    functions like Keddy, Dcoef, and so on need to be evaluated at alts 0 through 250 km, or at
    num_layers+2 (126 for 0-250 km grid) altitudes, AND they depend on the atmospheric density dictionary,
    which is only defined for num_layers altitudes. That means that to calculate Keddy etc. as
    vectors, we need to "duplicate" the first and last entries in the atmospheric density matrix,
    i.e. we put a "fake" density in at alt = 0 and alt = 250 km which is the same as the 
    densities at alt = 2 km and alt = 248 km, respectively. 
    
    In the non-vectorized version, Keddy etc. get the atmospheric density at the boundary layers
    by using n_alt_index, which uses clamp, and maps alt=0 => i=1, alt=2 => i=1, alt=4 => i=2...
    alt=248 => i=124, alt=250 => i=124. 
    
    So this does the same thing but as an array.
    
    atmdict_no_bdys: This is the atmospheric state DICTIONARY without the boundary layers,
                       i.e. only num_layers rows.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:n_alt_index, :all_species])

    # This gets a sorted list of the clamped indices, so it's [1, 1, 2, 3...end-1, end, end].
    clamped_n_alt_index = sort(collect(values(GV.n_alt_index)))
    
    atmdict_with_bdy_layers = Dict{Symbol, Vector{ftype_ncur}}()
    
    # Fill the dictionary with the profile. This duplicates the lowest and highest altitude values.
    for i in 1:length(GV.all_species)
        atmdict_with_bdy_layers[GV.all_species[i]] = atmdict_no_bdys[GV.all_species[i]][clamped_n_alt_index]
    end
    return atmdict_with_bdy_layers
end

function n_tot(atmdict::Dict{Symbol, Vector{ftype_ncur}}, z; globvars...)
    #= 
    Calculates total atmospheric density at altitude z.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        z: altitude, in cm
    Output: 
        Density of the atmosphere at altitude z
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:n_alt_index, :all_species])

    thisaltindex = GV.n_alt_index[z]
    return sum( [atmdict[s][thisaltindex] for s in GV.all_species] )
end

function n_tot(atmdict::Dict{Symbol, Vector{ftype_ncur}}; globvars...)
    #= 
    Override to calculate total atmospheric density at all altitudes.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
    Output: 
        Density of the atmosphere at all non-boundary layer altitudes.

    This function is agnostic as to the number of atmospheric layers. it collects it directly from atmdict.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:all_species])

    ndensities = zeros(length(GV.all_species), length(atmdict[collect(keys(atmdict))[1]]))
    for i in 1:length(GV.all_species)
        ndensities[i, :] = atmdict[GV.all_species[i]]
    end

    # returns the sum over all species at each altitude as a vector.
    return vec(sum(ndensities, dims=1)) 
end

function reduced_mass(mA, mB)
    #=
    Returns reduced mass.
    Input:
        mA, mB: species masses in AMU
        Uses global variable mH which is mass of hydrogen in GRAMS.
    Output:
        reduced mass in grams.
    =#
    try
        @assert floor(log10(mH)) == -24
    catch AssertionError
        throw("mH is somehow set to the wrong units")
    end

    return ((1/(mA*mH)) + (1/(mB*mH)))^(-1)
end

function scaleH(z, sp::Symbol, T; globvars...)
    #=
    Input:
        z: Altitudes in cm
        sp: Speciecs to calculate for
        T: temperature array for this species
    Output: 
        species-specific scale height at all altitudes (in cm)
    =#  
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:molmass])

    return @. kB*T/(GV.molmass[sp]*mH*marsM*bigG)*(((z+radiusM))^2)
end

function scaleH(atmdict::Dict{Symbol, Vector{ftype_ncur}}, T::Vector; globvars...)
    #= 
    Input:
        atmdict: Present atmospheric state dictionary
        T: temperature array for the neutral atmosphere
    Output:
        Mean atmospheric scale height at all altitudes (in cm)
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :molmass, :n_alt_index])

    mm_vec = meanmass(atmdict; globvars...)#, all_species, mmass) # vector version.
    return @. kB*T/(mm_vec*mH*marsM*bigG)*(((GV.alt+radiusM))^2)
end

function unflatten_atm(n_vec, species_list; globvars...)
    #=
    Input:
        n_vec: flattened density vector for the species in species_list: [n_sp1(z=0), n_sp2(z=0)...n_sp1(z=250)...n_spN(z=250)] 
    Output:
        dictionary of atmospheric densities by altitude with species as keys 

    This function is the reverse of flatten_atm.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:num_layers])

    n_matrix = reshape(n_vec, (length(species_list), GV.num_layers))

    return atm_matrix_to_dict(n_matrix, species_list)
end