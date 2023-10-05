# **************************************************************************** #
#                                                                              #
#         Functions for atmospheric dictionary object or its constituents      #
#                                                                              #
# **************************************************************************** #

function column_density(n::Vector; start_alt=1, end_alt=9999)
    #=
    Returns column density of n above ONE atmospheric layer defined by start_alt.

    Input
        n: species number density (#/cm³) by altitude
    Optional input:
        start_alt: index of starting altitude
        end_alt: index of end altitude. If 9999, end_alt will be set to the length of n.
    Output
        Column density (#/cm²)
    =#

    end_alt = end_alt==9999 ? length(n) : end_alt 
    return sum(n[start_alt:end_alt] .* dz)
end

function column_density_above(n_tot_by_alt::Vector)
    #=
    Returns an array where entries are the total integrated column density above
    that level of the atmosphere. e.g. the value at the topmost altitude is 
    called 0 since we assume anything beyond that level can escape. 

    This is NOT redundant with column_density.
    
    n_tot_by_alt: Total atmospheric density at each altitude layer.
    =#
    col_above = zeros(size(n_tot_by_alt))

    for i in 1:num_layers
        col_above[i] = column_density(n_tot_by_alt; start_alt=i+1)
    end

    return col_above
end

function column_density_species(atmdict, sp; start_alt=0., end_alt=250e5, globvars...)
    #=
    Returns the column density of species sp in atmosphere atmdict between the two altitudes (inclusive).
    =#
    GV = values(globvars)
    required = [:n_alt_index, :dz]
    check_requirements(keys(GV), required)

    return column_density(atmdict[sp]; start_alt=n_alt_index[start_alt], end_alt=n_alt_index[end_alt])
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
    required = [:e_profile_type, :ion_species, :non_bdy_layers]
    check_requirements(keys(GV), required)

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
    required =  [:non_bdy_layers, :all_species, :Tn, :molmass, :alt, :collision_xsect, :n_alt_index, :zmax]
    check_requirements(keys(GV), required)

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

function meanmass(atmdict::Dict{Symbol, Vector{ftype_ncur}}, z; globvars...)
    #= 
    find the mean molecular mass at a given altitude z

    atmdict: species number density by altitude
    z: float; altitude in atmosphere in cm

    return: mean molecular mass in amu
    =#
    GV = values(globvars)
    required =  [:all_species, :molmass, :n_alt_index]
    check_requirements(keys(GV), required)

    counted_species = setdiff(GV.all_species, ignore)

    thisaltindex = GV.n_alt_index[z]
    c = [atmdict[sp][thisaltindex] for sp in counted_species]
    m = [GV.molmass[sp] for sp in counted_species]
    return sum(c.*m)/sum(c)
end

function meanmass(atmdict::Dict{Symbol, Vector{ftype_ncur}}; ignore=[], globvars...)
    #= 
    Override for vector form. Calculates mean molecular mass at all atmospheric layers.

    Inputs:
        atmdict: Array; species number density by altitude
        ignore: Set; contains symbols representing species to ignore in the calculation

    Outputs:
        returns: mean molecular mass in amu for all atmospheric layers.
    =#

    GV = values(globvars)
    required = [:all_species, :molmass, :n_alt_index]
    check_requirements(keys(GV), required)

    counted_species = setdiff(GV.all_species, ignore)

    # Delete ignored species from the dictionary since we have to transform it
    trimmed_atmdict = deepcopy(atmdict)
    for isp in ignore
        delete!(trimmed_atmdict, isp)
    end

    # Gets the atmosphere as a matrix with rows = altitudes and cols = species
    # so we can do matrix multiplication.
    n_mat = transpose(atm_dict_to_matrix(trimmed_atmdict, counted_species))

    m = [GV.molmass[sp] for sp in counted_species] # this will always be 1D

    weighted_mm = zeros(size(n_mat)[1]) # This will store the result

    # Multiply densities of each species by appropriate molecular mass 
    mul!(weighted_mm, n_mat, m)

    return weighted_mm ./ n_tot(trimmed_atmdict; all_species=counted_species, GV.n_alt_index)
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
    required =  [:n_alt_index, :all_species]
    check_requirements(keys(GV), required)

    # This gets a sorted list of the clamped indices, so it's [1, 1, 2, 3...end-1, end, end].
    clamped_n_alt_index = sort(collect(values(GV.n_alt_index)))
    
    atmdict_with_bdy_layers = Dict{Symbol, Vector{ftype_ncur}}()
    
    # Fill the dictionary with the profile. This duplicates the lowest and highest altitude values.
    for i in 1:length(GV.all_species)
        atmdict_with_bdy_layers[GV.all_species[i]] = atmdict_no_bdys[GV.all_species[i]][clamped_n_alt_index]
    end
    return atmdict_with_bdy_layers
end

function n_tot(atmdict::Dict{Symbol, Vector{ftype_ncur}}, z; ignore=[], globvars...)
    #= 
    Calculates total atmospheric density at altitude z.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        z: altitude, in cm
        ignore: Set; contains symbols representing species to ignore in the calculation
    Output: 
        Density of the atmosphere at altitude z
    =#
    GV = values(globvars)
    required = [:n_alt_index, :all_species]
    check_requirements(keys(GV), required)

    counted_species = setdiff(GV.all_species, ignore)

    thisaltindex = GV.n_alt_index[z]
    return sum( [atmdict[s][thisaltindex] for s in counted_species] )
end

function n_tot(atmdict::Dict{Symbol, Vector{ftype_ncur}}; ignore=[], globvars...)
    #= 
    Override to calculate total atmospheric density at all altitudes.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        ignore: Set; contains symbols representing species to ignore in the calculation
    Output: 
        Density of the atmosphere at all non-boundary layer altitudes.

    This function is agnostic as to the number of atmospheric layers. it collects it directly from atmdict.
    =#
    GV = values(globvars)
    required =  [:all_species]
    check_requirements(keys(GV), required)

    counted_species = setdiff(GV.all_species, ignore)
    ndensities = zeros(length(counted_species), length(atmdict[collect(keys(atmdict))[1]]))

    for i in 1:length(counted_species)
        ndensities[i, :] = atmdict[counted_species[i]]
    end

    # returns the sum over all species at each altitude as a vector.
    return vec(sum(ndensities, dims=1)) 
end

function optical_depth(n_cur_densities; globvars...)
    #=
    Given the current state (atmdict), this populates solarabs, a 1D array of 1D arrays 
    (which is annoying, but required for using BLAS.axpy! for some inscrutable reason) 
    with the optical depth of the atmosphere. The shape of solar abs is 124 elements, each 
    its own array of 2000 elements. 
    =#
    
    GV = values(globvars)
    required = [:num_layers, :Jratelist, :absorber, :crosssection, :dz]
    check_requirements(keys(GV), required)
    
    nlambda = 2000
    
    # Initialize the solar absorption array with 0s for all wavelengths.
    solarabs = Array{Array{Float64}}(undef, GV.num_layers)
    for i in range(1, length=GV.num_layers)
        solarabs[i] = zeros(Float64, nlambda)
    end
    
    for jspecies in GV.Jratelist
        species = GV.absorber[jspecies]

        jcolumn = convert(Float64, 0.)

        for ialt in [GV.num_layers:-1:1;]
            #get the (overhead) vertical column of the absorbing constituent
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
    return solarabs
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
        ignore: Set; contains symbols representing species to ignore in the calculation
    Output: 
        species-specific scale height at all altitudes (in cm)
    =#  
    GV = values(globvars)
    required = [:molmass]
    check_requirements(keys(GV), required)

    return @. kB*T/(GV.molmass[sp]*mH*marsM*bigG)*(((z+radiusM))^2)
end

function scaleH(atmdict::Dict{Symbol, Vector{ftype_ncur}}, T::Vector; ignore=[], globvars...)
    #= 
    Input:
        atmdict: Present atmospheric state dictionary
        T: temperature array for the neutral atmosphere
        ignore: Set; contains symbols representing species to ignore in the calculation
    Output:
        Mean atmospheric scale height at all altitudes (in cm)
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :molmass, :n_alt_index]
    check_requirements(keys(GV), required)

    counted_species = setdiff(GV.all_species, ignore)

    mm_vec = meanmass(atmdict; ignore=ignore, globvars...)#, all_species, mmass) # vector version.
    return @. kB*T/(mm_vec*mH*marsM*bigG)*(((GV.alt+radiusM))^2)
end

