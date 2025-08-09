# **************************************************************************** #
#                                                                              #
#                             PHOTOCHEMICAL CORE                               #
#                                                                              #
# This file forms the logical core of the photochemical model. Some things     #
# not included in this file are:                                               #
# 1. Anything related to file input/output/writing is in FileIO.jl.            #
# 2. Everything related to photochemical cross sections is in Crosssections.jl.#
# 3. Everything related to building the chemical network object from the       #
#    provided spreadsheet is in ReactionNetwork.jl.                            #
# 4. Really small functions necessary to the model but that are not as         #
#    "complicated" as these core functions are in BasicUtilities.jl.           #
#                                                                              #
# **************************************************************************** #

# Special Exception for this model.
struct TooManyIterationsException <: Exception end

#===============================================================================#
#                      Atmospheric attribute calculations                       #
#===============================================================================#


function column_density(n::Vector; start_alt=1, end_alt=9999, globvars...)
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
    GV = values(globvars)
    required =  [:dz]
    check_requirements(keys(GV), required)

    end_alt = end_alt==9999 ? length(n) : end_alt 
    return sum(n[start_alt:end_alt] .* GV.dz)
end

function column_density_above(n_tot_by_alt::Vector; globvars...)
    #=
    Returns an array where entries are the total integrated column density above
    that level of the atmosphere. e.g. the value at the topmost altitude is 
    called 0 since we assume anything beyond that level can escape. 

    This is NOT redundant with column_density.
    
    n_tot_by_alt: Total atmospheric density at each altitude layer.
    =#
    GV = values(globvars)
    required =  [:dz, :num_layers]
    check_requirements(keys(GV), required)

    col_above = zeros(size(n_tot_by_alt))

    for i in 1:GV.num_layers
        col_above[i] = column_density(n_tot_by_alt; start_alt=i+1, globvars...)
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

    return column_density(atmdict[sp]; start_alt=GV.n_alt_index[start_alt], end_alt=GV.n_alt_index[end_alt], globvars...)
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

    # Determine how many horizontal columns are present.
    n_horiz = get(GV, :n_horiz, 1)

    if GV.e_profile_type == "constant"
        # Fill with a constant electron density profile for every column.
        template = fill(1e5, length(GV.non_bdy_layers))
        return [copy(template) for _ in 1:n_horiz]

    elseif GV.e_profile_type == "quasineutral"
        # Multi column atmospheres store species as vectors of altitude arrays.
        # Handle both the legacy single column format and the new multicolumn
        # format.
        if atmdict[GV.ion_species[1]] isa Vector
            return [sum(atmdict[sp][ihoriz] for sp in GV.ion_species)
                    for ihoriz in 1:n_horiz]
        else
            return [sum(atmdict[sp] for sp in GV.ion_species)]
        end

    elseif GV.e_profile_type == "none"
        template = zeros(length(GV.non_bdy_layers))
        return [copy(template) for _ in 1:n_horiz]

    else
        throw("Unhandled electron profile specification: $(GV.e_profile_type)")
    end
end

function find_exobase(sp::Symbol, atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}};
                      returntype="index", verbose=false, globvars...)
    #=
    Finds the exobase altitude, where mean free path is equal to a scale height.
    ONLY VALID FOR H AND D (because collision_xsect only includes xsects for those, and they're H and D on O)

    Inputs:
        sp: species 
        atmdict: Atmospheric state dictionary
        returntype: whether to return the "altitude" in km or the "index" in the n_alt_index dictionary. 
    Output:
        Array of exobase altitudes or indices for each horizontal column.
    =#

    if !(sp in [:H, :D])
        throw("find_exobase is only defined for species H and D at this time.")
    end

    GV = values(globvars)
    required = [:all_species, :alt, :collision_xsect, :M_P, :molmass, :non_bdy_layers, 
                :n_alt_index, :R_P, :Tn, :zmax, :n_horiz]
    check_requirements(keys(GV), required)

    exobase_altitudes = fill(NaN, GV.n_horiz)
    exobase_indices = fill(0, GV.n_horiz)

    # Loop over each horizontal column to find individual exobase altitudes
    for ihoriz in 1:GV.n_horiz
        Tn_col = GV.Tn[ihoriz, :]

        # Calculate scale height using temperature for this column
        # H_s = scaleH(GV.non_bdy_layers, sp, Tn_col[ihoriz, 2:end-1]; globvars...)
        H_s = scaleH(GV.non_bdy_layers, sp, Tn_col[2:end-1]; globvars...)

        # Mean free path calculation for each vertical column separately
        mfp_sp = 1 ./ (GV.collision_xsect[sp] .* n_tot(atmdict, ihoriz; GV.all_species, GV.n_alt_index))

        exobase_alt = findfirst(mfp_sp .>= H_s)

        if exobase_alt === nothing # If no exobase is found, use the top of the atmosphere.
            if verbose
                println("Warning: No exobase found for species $(sp) at column $(ihoriz); assuming top of atmosphere, but this is not guaranteed to be true.")
            end
            exobase_altitudes[ihoriz] = GV.zmax
            exobase_indices[ihoriz] = GV.n_alt_index[GV.zmax]
        else
            exobase_altitudes[ihoriz] = GV.alt[exobase_alt]
            exobase_indices[ihoriz] = exobase_alt
        end
    end

    # Choose the return type: altitude or index
    if returntype=="altitude"
        return exobase_altitudes
    elseif returntype=="index"
        return exobase_indices
    else
        throw("Invalid returntype specified. Choose 'altitude' or 'index'")
    end
end

function meanmass(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, n_horiz::Int64, ihoriz::Int64; ignore=[], globvars...)
    #= 
    Override for vector form. Calculates mean molecular mass at all atmospheric layers.

    Inputs:
        atmdict: Array; species number density by altitude
        n_horiz: Integer; number of vertical columns in the model
        ihoriz: Integer; vertical column index
        ignore: Set; contains symbols representing species to ignore in the calculation

    Outputs:
        returns: mean molecular mass in amu for all atmospheric layers.
    =#

    GV = values(globvars)
    required = [:all_species, :molmass, :n_alt_index]
    check_requirements(keys(GV), required)

    # counted_species = setdiff(GV.all_species, ignore)
    counted_species = [s for s in GV.all_species if haskey(atmdict, s) && !(s in ignore)]

    # Delete ignored species from the dictionary since we have to transform it
    trimmed_atmdict = deepcopy(atmdict)
    for isp in ignore
        delete!(trimmed_atmdict, isp)
    end

    # Gets the atmosphere as a matrix with rows = altitudes, cols = species, and the third dimension as vertical columns
    # so we can do matrix multiplication. Only counted species are included so
    # that ignored species do not contribute to the mean.
    n_mat = permutedims(atm_dict_to_matrix(trimmed_atmdict, counted_species, n_horiz), (2, 1, 3))

    m = [GV.molmass[sp] for sp in counted_species] # this will always be 1D

    weighted_mm = zeros(size(n_mat)[1]) # This will store the result

    # Multiply the density of each species in the requested column by its
    # molecular mass.  The third dimension of `n_mat` corresponds to the
    # horizontal column number.
    mul!(weighted_mm, n_mat[:, :, ihoriz], m)

    return weighted_mm ./ n_tot(trimmed_atmdict, ihoriz; all_species=counted_species, GV.n_alt_index)
end

function n_tot(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, z, ihoriz::Int64; ignore=[], globvars...)
    #= 
    Calculates total atmospheric density at altitude z.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        z: altitude, in cm
        ignore: Set; contains symbols representing species to ignore in the calculation
        ihoriz: vertical column index
    Output: 
        Density of the atmosphere at altitude z
    =#
    GV = values(globvars)
    required = [:n_alt_index, :all_species]
    check_requirements(keys(GV), required)

    # counted_species = setdiff(GV.all_species, ignore)
    counted_species = [s for s in GV.all_species if haskey(atmdict, s) && !(s in ignore)]

    thisaltindex = GV.n_alt_index[z]
    # Sum the densities of all counted species at the specified altitude and
    # horizontal column.
    return sum(atmdict[s][ihoriz][thisaltindex] for s in counted_species)
end

function n_tot(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, ihoriz::Int64; ignore=[], globvars...)
    #= 
    Override to calculate total atmospheric density at all altitudes.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        ignore: Set; contains symbols representing species to ignore in the calculation
        ihoriz: vertical column index
    Output: 
        Density of the atmosphere at all non-boundary layer altitudes.

    This function is agnostic as to the number of atmospheric layers. it collects it directly from atmdict.
    =#
    GV = values(globvars)
    required =  [:all_species]
    check_requirements(keys(GV), required)
    
    # counted_species = setdiff(GV.all_species, ignore)
    counted_species = [s for s in GV.all_species if haskey(atmdict, s) && !(s in ignore)]
    # allocate an array to gather density profiles for this vertical column
    if isempty(counted_species)
        return zeros(length(atmdict[collect(keys(atmdict))[1]][ihoriz]))
    end
    ndensities = zeros(length(counted_species), length(atmdict[collect(keys(atmdict))[1]][ihoriz]))

    for i in 1:length(counted_species)
        # copy the density profile for each species at this column
        ndensities[i, :] = atmdict[counted_species[i]][ihoriz]
    end

    # returns the sum over all species at each altitude as a vector.
    return vec(sum(ndensities, dims=1))
end

function optical_depth(n_cur_densities; n_horiz::Int64, globvars...)
    #=
    Given the current state (atmdict), this populates solarabs, a 2D array structured as [n_horiz, num_layers],
    with each element containing an array of 2000 elements (one per wavelength).

    The optical depth is calculated independently for each horizontal column.
    =#

    GV = values(globvars)
    required = [:num_layers, :Jratelist, :absorber, :crosssection, :dz]
    check_requirements(keys(GV), required)

    nlambda = 2000

    # Initialize the solar absorption array for all wavelengths and horizontal columns.
    solarabs = [[zeros(Float64, nlambda) for ialt in 1:GV.num_layers] for ihoriz in 1:n_horiz]

    for jspecies in GV.Jratelist
        species = GV.absorber[jspecies]

        # Skip species absent from either the current atmosphere or the cross-section dictionary
        if !(haskey(n_cur_densities, species) && haskey(GV.crosssection, jspecies))
            continue
        end

        for ihoriz in 1:n_horiz
            jcolumn = convert(Float64, 0.)

            for ialt in GV.num_layers:-1:1
                # Vertical column of the absorbing constituent for this horizontal column
                jcolumn += convert(Float64, n_cur_densities[species][ihoriz][ialt]) * GV.dz

                # Add total extinction to solarabs for this horizontal column and altitude
                # multiplies air column density (N, #/cm^2) at all wavelengths by crosssection (σ)
                # to get optical depth (τ). This is an override of axpy! to use the
                # full arguments. For the equation Y' = alpha*X + Y:
                # ARG 1: n (length of arrays in ARGS 3, 5)
                # ARG 2: alpha, a scalar.
                # ARG 3: X, an array of length n.
                # ARG 4: the increment of the index values of X, maybe?
                # ARG 5: Y, an array of length n
                # ARG 6: increment of index values of Y, maybe?

                # updated optical_depth to read column-specific cross sections when computing extinction
                BLAS.axpy!(nlambda, jcolumn, GV.crosssection[jspecies][ihoriz][ialt+1], 1, solarabs[ihoriz][ialt], 1)
            end
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
    required = [:molmass, :M_P, :R_P]
    check_requirements(keys(GV), required)

    return @. kB*T/(GV.molmass[sp]*mH*GV.M_P*bigG)*(((z+GV.R_P))^2)
end

function scaleH(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, T::Vector, n_horiz::Int64, ihoriz::Int64; ignore=[], globvars...)
    #= 
    Input:
        atmdict: Present atmospheric state dictionary
        T: temperature array for the neutral atmosphere
        n_horiz: number of vertical columns in the model
        ihoriz: vertical column index
        ignore: Set; contains symbols representing species to ignore in the calculation
    Output:
        Mean atmospheric scale height at all altitudes (in cm)
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :M_P, :molmass, :n_alt_index, :R_P]
    check_requirements(keys(GV), required)

    counted_species = setdiff(GV.all_species, ignore)

    mm_vec = meanmass(atmdict, n_horiz, ihoriz; ignore=ignore, globvars...)
    return @. kB*T/(mm_vec*mH*GV.M_P*bigG)*(((GV.alt+GV.R_P))^2)
end

function scaleH_lowerboundary(z::Float64, T::Float64; globvars...)
    #= 
    Input:
        z: Some altitude (in cm)
        atmdict: atmospheric state dictionary, needed to calculate over other species
        T: Temperature at altitude z
    Output:
        Special case of mean atmospheric scale height at one altitude in cm, relevant for
        lower boundary condition. Here we will assume the mean mass at the lower boundary is constant,
        since we are assuming a density boundary condition for CO2, this will be a pretty ok assumption.
        By using this function, we don't have to make the lower bc for velocities a function of 
        the present atmosphere like we do for non-thermal escape (which would require a big rewrite of 
        the boundaryconditions() velocity section). 
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:M_P, :molmass, :R_P, :zmin])
    @assert z==GV.zmin

    mm = GV.molmass[:CO2]
    return kB*T/(mm*mH*GV.M_P*bigG)*(((z+GV.R_P))^2)
end 

# Subsection - functions that manipulate the atmospheric dictionary/matrix object.
#---------------------------------------------------------------------------------#

function atm_dict_to_matrix(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, species_list, n_horiz::Int64)
    #=
    Converts atmospheric state dictionary atmdict to a matrix,
    such that rows correspond to species in the order listed in species_list
    and columns correspond to altitudes in the order lowest-->highest
    and the third dimension corresponds to the vertical column in the order lowest-->highest.
    =#
    
    num_alts = length(atmdict[collect(keys(atmdict))[1]][1])
    n_mat = zeros(length(species_list), num_alts, n_horiz)
    
    for i in 1:length(species_list)
        for ihoriz in 1:n_horiz
            n_mat[i, :, ihoriz] = atmdict[species_list[i]][ihoriz]
	end
    end
    
    return n_mat
end

function atm_matrix_to_dict(n_matrix, species_list, n_horiz::Int64)
    #=
    Input:
        n_matrix: matrix of the atmospheric state
    Output:
        dictionary for only the species in species_list
    =#
    atmdict = Dict{Symbol, Vector{Array{ftype_ncur}}}([species_list[k]=>[n_matrix[k, :, ihoriz] for ihoriz in 1:n_horiz] for k in 1:length(species_list)])
    
    return atmdict
end

function compile_ncur_all(n_long, n_horiz::Int64, n_short, n_inactive; globvars...)
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
    required = [:active_longlived, :active_shortlived, :inactive_species, :num_layers]
    check_requirements(keys(GV), required)

    n_cur_active_long = unflatten_atm(n_long, GV.active_longlived, n_horiz; num_layers=GV.num_layers)
    n_cur_active_short = unflatten_atm(n_short, GV.active_shortlived, n_horiz; num_layers=GV.num_layers)
    n_cur_inactive = unflatten_atm(n_inactive, GV.inactive_species, n_horiz; num_layers=GV.num_layers)

    n_cur_all = Dict(vcat([k=>n_cur_active_long[k] for k in keys(n_cur_active_long)],
                          [k=>n_cur_active_short[k] for k in keys(n_cur_active_short)],
                          [k=>n_cur_inactive[k] for k in keys(n_cur_inactive)]))
    
    return n_cur_all
end

function flatten_atm(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, species_list, n_horiz::Int64; globvars...) 
    #=
    Input:
        atmdict: atmospheric densities by altitude
        species_list: Included species which will have profiles flattened
    Output:
        Vector of form [n_sp1(h=1,z=0), n_sp2(h=1,z=0)...n_sp1(h=1,z=zmax)...n_spN(h=1,z=zmax),n_sp1(h=2,z=0), n_sp2(h=2,z=0)...n_sp1(h=2,z=zmax)...n_spN(h=2,z=zmax)]
    
    This function is the reverse of unflatten_atm. 
    =#

    GV = values(globvars)
    required =  [:num_layers]
    check_requirements(keys(GV), required)

    # Construct a matrix in the (species, altitude, column) layout and then
    # flatten it with Julia's column-major ordering so that unflatten_atm can
    # correctly reshape it.
    return vec(atm_dict_to_matrix(atmdict, species_list, n_horiz))
end

function ncur_with_boundary_layers(atmdict_no_bdys::Dict{Symbol, Vector{Array{ftype_ncur}}}, n_horiz::Int64; globvars...)
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
    
    atmdict_with_bdy_layers = Dict{Symbol, Vector{Array{ftype_ncur}}}()
    
    # Fill the dictionary with the profile. This duplicates the lowest and highest altitude values.
    for i in 1:length(GV.all_species)
        atmdict_with_bdy_layers[GV.all_species[i]] = ([atmdict_no_bdys[GV.all_species[i]][ihoriz][clamped_n_alt_index] for ihoriz in 1:n_horiz])
    end
    return atmdict_with_bdy_layers
end

function unflatten_atm(n_vec, species_list, n_horiz::Int64; globvars...)
    #=
    Input:
        n_vec: flattened density vector for the species in species_list: [n_sp1(h=1,z=0), n_sp2(h=1,z=0)...n_sp1(h=1,z=zmax)...n_spN(h=1,z=zmax),n_sp1(h=2,z=0), n_sp2(h=2,z=0)...n_sp1(h=2,z=zmax)...n_spN(h=2,z=zmax)]

    Output:
        dictionary of atmospheric densities by altitude with species as keys 

    This function is the reverse of flatten_atm.
    =#
    GV = values(globvars)
    required =  [:num_layers]
    check_requirements(keys(GV), required)

    n_matrix = reshape(n_vec, (length(species_list), GV.num_layers, n_horiz))

    return atm_matrix_to_dict(n_matrix, species_list, n_horiz)
end

#===============================================================================#
#                             Chemistry functions                               #
#===============================================================================#

function calculate_stiffness(J)
    #=
    Input:
        J: a jacobian matrix
    Output:
        r: stiffness r = max(|Re(λ)|) / min(|Re(λ)), where λ is the matrix of eigenvalues of J.
    =#

    if typeof(J) == SparseMatrixCSC{ftype_chem, Int64}
        J = Array(J)
    end

    r = maximum(abs.(real(eigvals(J)))) / minimum(abs.(real(eigvals(J))))
    # if r == Inf
    #     println("Eigenvalues:")
    #     println(eigvals(J))
    #     throw("Error: Infinite stiffness")
    # end
    return r
end   

function charge_type(sp::Symbol)
    #=
    Returns a string representing the type of sp, i.e. ion, neutral, or electron
    =#
    if occursin("pl", String(sp))
        return "ion"
    elseif sp==:E
        return "electron"
    else
        return "neutral"
    end
end

function check_jacobian_eigenvalues(J, path)
    #=
    Check a jacobian matrix to see if it has complex eigenvalues.
    Per Jacob 2003, Models of Atmospheric Transport and Chemistry,
    all stable atmospheric chemistry models should result in real and 
    negative eigenvalues of jacobians. 

    Input:
        J: a Jacobian matrix, sparse or normal.
    Output:
        Print statements about eigenvalue types.
    =#

    # Warning: This tends to take a long time.
    if typeof(J) == SparseMatrixCSC{ftype_chem, Int64}
        J_nonsparse = Array(J)
    end

    if any(i->isnan(i), eigvals(J_nonsparse)) || any(i->isinf(i), eigvals(J_nonsparse))
        throw("ValueError: Jacobian eigenvalues have inf values: $(any(i->isinf(i), eigvals(J_nonsparse))); NaN values: $(any(i->isnan(i), eigvals(J_nonsparse)))")
    end

    if all(i->typeof(i) != ComplexF64, eigvals(J_nonsparse)) # all eigenvalues are real
        if all(i->i<0, eigvals(J_nonsparse)) # all eigenvalues are real and negative
            return 0
        else
            println("Warning: Some Jacobian eigenvalues are real and positive. Solution will grow without bound")
        end

    elseif any(i->typeof(i) == ComplexF64, eigvals(J_nonsparse))  # complex eigenvalues are present
        # f = open(path*"/jacobian_eigenvalues.txt", "w")
        println("Warning: Some Jacobian eigenvalues are complex. Solution behavior is unpredictable and may oscillate (possibly forever?).")

        # # find the indices of the eigenvalues that are complex
        # i = findall(i->imag(i)!=0, eigvals(J))

        # # get the associated eigenvectors - they are in the ith column of the eigenvector matrix.
        # problem_eigenvecs = eigvecs(J)[:, i]


        # if any(i->real(i)>0, eigvals(J)[imag(eigvals(J)) .== 0])
        #     println("Warning: Some Jacobian eigenvalues are real and positive. Solution will grow without bound")
        #     # find the indices of the eigenvalues that are real and positive
        #     i = findall(i->i>0, eigvals(J))

        #     # get the associated eigenvectors - they are in the ith column of the eigenvector matrix.
        #     problem_eigenvecs = eigvecs(J)[:, i]
            
        # end
        # close(f)
    end
end

function chemical_jacobian(specieslist, dspecieslist; diff_wrt_e=true, diff_wrt_m=true, transportnet_horiz, globvars...)
    #= 
    Compute the symbolic chemical jacobian of a supplied chemnet and transportnet
    for the specified specieslist. 

    Input:
        chemnet: chemical reactions
        transportnet: transport equations
        specieslist: list of species to calculate for; i.e. the species whose equation is differentiated
        dspecieslist: species with respect to which the equations are differentiated
        chem_species: list of chemistry species
        transportsp: list of transport species
        ionsp: the list of ions, needed to add in the terms arising from dependence on E, 
               since density of electrons is just the sum of all the ion densities.
    Output:
        three arrays suitable for constructing a sparse matrix: 
            i: row indices
            j: column indices
            v: values to place at (i, j)
    =#

    GV = values(globvars)
    required = [:chem_species, :transport_species, :chemnet, :transportnet]
    check_requirements(keys(GV), required)

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)  # this is the active species. 
    ndspecies = length(dspecieslist)  # this is the species with respect to which we differentiate

    if diff_wrt_e==true
        required = [:ion_species]
        check_requirements(keys(GV), required)
        ion_cols = indexin(GV.ion_species, specieslist) # for filling in the derivatives wrt electrons 
        # println("The ion columns are $(ion_cols), expect to see new terms in those columns") # For testing
    end

    for i in 1:nspecies # for each species
        ispecies = specieslist[i]

        # get the production and loss equations
        peqn = []
        leqn = []
        if issubset([ispecies], GV.chem_species)
            peqn = [peqn; production_equations(ispecies, GV.chemnet)] 
            leqn = [leqn; loss_equations(ispecies, GV.chemnet)]
        end
        if issubset([ispecies],GV.transport_species)
            peqn = [peqn; production_equations(ispecies, GV.transportnet)]
            leqn = [leqn; loss_equations(ispecies, GV.transportnet)]
	    peqn = [peqn; production_equations(ispecies, transportnet_horiz)]
	    leqn = [leqn; loss_equations(ispecies, transportnet_horiz)]
	end

        # Account for e's
        if diff_wrt_e
            ppos_electrons = map(x->deletefirst(peqn[x[1]], :E), getpos(peqn, :E))
            lpos_electrons = map(x->deletefirst(leqn[x[1]], :E), getpos(leqn, :E)) # differentiate wrt e term using deletefirst
            # Following line will keep track of other reactants in the case you have molecules and an electron reacting. 
            # Presently this is not used as there are no reactions like that. 
            # other_reactants = [i[1:end-1] for i in lpos_electrons] 
        end

        # Account for dependence on M
        if diff_wrt_m
            ppos_M = map(x->deletefirst(peqn[x[1]], :M), getpos(peqn, :M))
            lpos_M = map(x->deletefirst(leqn[x[1]], :M), getpos(leqn, :M))
        end

        for j in 1:ndspecies # now take the derivative with respect to the other species
            jspecies = dspecieslist[j]
            #= find the places where the production rates depend on
            jspecies, and return the list rates with the first
            occurrance of jspecies deleted. (Note: this seamlessly
            deals with multiple copies of a species on either side of
            an equation, because it is found twice wherever it lives) =#

            ppos = map(x->deletefirst(peqn[x[1]],jspecies), getpos(peqn, jspecies))
            lpos = map(x->deletefirst(leqn[x[1]],jspecies), getpos(leqn, jspecies))

            # Add in dependence on M. Since M is a sum of all species, this has to be added to every column. 
            if diff_wrt_m
                append!(ppos, ppos_M)
                append!(lpos, lpos_M)
            end

            if diff_wrt_e
                if j in ion_cols
                    append!(ppos, ppos_electrons)
                    append!(lpos, lpos_electrons)
                end
            end

            if length(ppos)+length(lpos)>0 #if there is a dependence
                #make note of where this dependency exists 
                append!(ivec,[i])
                append!(jvec,[j])

                #= smash the production and loss rates together,
                multiplying for each distinct equation, adding
                together the production and loss seperately, and
                subtracting loss from production. =#
                if length(ppos)==0
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($lval))
                elseif length(lpos)==0
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    tval = :(+($pval))
                else
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($pval,$lval))
                end
                # attach the symbolic expression to the return values
                append!(tvec,[tval])

            end
        end
    end

    # println(ivec, jvec, tvec) # TEST

    return (ivec, jvec, tvec)
end

# function eval_rate_coef(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, krate::Expr, ihoriz::Int64; globvars...)
function eval_rate_coef(
    atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}},
    krate::Expr,
    ihoriz::Int64;
    globvars...
)
    #=
    Evaluates a chemical reaction rate coefficient, krate, for all levels of the atmosphere. 

    Input:
        atmdict: the atmospheric state dictionary
        krate: rate coefficient for a single reaction
	ihoriz: Vertical column index
        tn, _i, _e: temperature profiles for neutrals, ions, electrons
    Output:
        rate_coefficient: evaluated rate coefficient at all atmospheric layers
    =#
    GV = values(globvars)
    required = [:Tn, :Ti, :Te, :all_species]
    check_requirements(keys(GV), required)

    # Set stuff up
    eval_k = mk_function(:((Tn, Ti, Te, M) -> $krate))

    # Grab the temperature profiles for this column and drop the boundary layers
    Tn_data = GV.Tn
    Ti_data = GV.Ti
    Te_data = GV.Te

    # When get_volume_rates is called for a single column, Tn, Ti and Te may
    # already be 1-D vectors. Support both cases gracefully.
    if ndims(Tn_data) == 2
        Tn_col = Tn_data[ihoriz, 2:end-1]
        Ti_col = Ti_data[ihoriz, 2:end-1]
        Te_col = Te_data[ihoriz, 2:end-1]
    else
        # The arrays already exclude boundary layers
        Tn_col = Tn_data
        Ti_col = Ti_data
        Te_col = Te_data
    end

    return eval_k(
        Tn_col,
        Ti_col,
        Te_col,
        sum([atmdict[sp][ihoriz] for sp in GV.all_species]),
    )
end 

# function getrate(sp::Symbol; chemistry_on=true, transport_on=true, sepvecs=false, globvars...)
function getrate(sp::Symbol; chemistry_on=true, transport_on=true, sepvecs=false,
                  transportnet_horiz=[], globvars...)
    #=
    Creates a symbolic expression for the rate at which a given species is
    either produced or lost due to chemical reactions or transport.

    Input:
        sp: species for which to get the rate 
        chemnet: chemistry reaction array
        transportnet: vertical transport network array
        transportnet_horiz: horizontal transport network array
        chem_species: species with active chemistry
        transport_species: species which transport
        chemistry_on: set to false to disallow chemical changes to species
        transport_on: set to false to disallow transport of a species
        sepvecs: Allows this function to return a vector of expressions for chemical production and loss 
                 and transport production and loss
    Output: either
        rate: the final value of dn/dt for sp from all processes, or
        chemprod_rate, chemloss_rate, transprod_rate, transloss_rate: dn/dt for sp due to these processes, calculated separately.
    =#

    GV = values(globvars)
    required =  [:chemnet, :transportnet, :chem_species, :transport_species]
    check_requirements(keys(GV), required)


    # This block will return the total net change for the species, P - L.
    if sepvecs == false
        rate = :(0.0)
        if issubset([sp], GV.chem_species) && chemistry_on
            # rate = :($rate 
            #          + $(production_rate(sp, GV.chemnet, sepvecs=sepvecs)) 
            #          - $(loss_rate(sp, GV.chemnet, sepvecs=sepvecs)) 
            rate = :($rate
                     + $(production_rate(sp, GV.chemnet, sepvecs=sepvecs))
                     - $(loss_rate(sp, GV.chemnet, sepvecs=sepvecs))
                    )
        end
        if issubset([sp], GV.transport_species) && transport_on
            # rate = :($rate 
            #          + $(production_rate(sp, GV.transportnet, sepvecs=sepvecs)) 
            #          - $(loss_rate(sp, GV.transportnet, sepvecs=sepvecs))
            rate = :($rate
                     + $(production_rate(sp, GV.transportnet, sepvecs=sepvecs))
                     - $(loss_rate(sp, GV.transportnet, sepvecs=sepvecs))
                     + $(production_rate(sp, transportnet_horiz, sepvecs=sepvecs))
                     - $(loss_rate(sp, transportnet_horiz, sepvecs=sepvecs))
                    )
        end
        return rate
    else  # if we want a vector of expressions for each production and loss (4 terms, 2 each for chemistry and transport)
        if issubset([sp], GV.chem_species) && chemistry_on
            chemprod_rate = production_rate(sp, GV.chemnet, sepvecs=sepvecs)
            chemloss_rate = loss_rate(sp, GV.chemnet, sepvecs=sepvecs)
        else
            chemprod_rate = [:(0.0 + 0.0)]  # Doing it this way because it's the easiest way to make a vector of one expression that's just 0
            chemloss_rate = [:(0.0 + 0.0)]
        end
        
        if issubset([sp], GV.transport_species) && transport_on
            # transprod_rate = production_rate(sp, GV.transportnet, sepvecs=sepvecs)
            # transloss_rate = loss_rate(sp, GV.transportnet, sepvecs=sepvecs)
            transprod_rate = vcat(production_rate(sp, GV.transportnet, sepvecs=sepvecs),
                                  production_rate(sp, transportnet_horiz, sepvecs=sepvecs))
            transloss_rate = vcat(loss_rate(sp, GV.transportnet, sepvecs=sepvecs),
                                  loss_rate(sp, transportnet_horiz, sepvecs=sepvecs))
        else
            transprod_rate = [:(0.0 + 0.0)]
            transloss_rate = [:(0.0 + 0.0)]

        end

        return chemprod_rate, chemloss_rate, transprod_rate, transloss_rate
    end
end

function loss_equations(sp::Symbol, network)
    #=  
    Input:
        sp: Species for which to construct a loss equation
        network: The type of loss process to consider, i.e. either a chemical reaction network or a transport network
    Output: 
        losseqns: loss equations and relevant rate coefficients for species sp.
                  the form is an array where each entry is of the form [reactants..., rate].
                  For example, [[:O2, :JO2toOpO], [:O1D, :O2, :k]] are two possible entries 
                  for loss of O2 (k will be some more complicated expression, but not important for this example).

    Automatically accounts for cases where a species occurs twice on the LHS by
    reporting those reactions twice.
    =#

    # Identify all positions of the species within the network.
    # The format is [R, side, el] where R = index of reaction vector in network,
    # side = 1 (reactant) or 2 (product) and finally
    # el = index of species position in either reactant or product vector.
    speciespos = getpos(network, sp)

    # This collects only R (see above comment) for any reaction here species is on LHS:
    lhspos = map(x->x[1], map(x->speciespos[x], findall(x->x[2]==1, speciespos)))
    # and RHS:
    rhspos = map(x->x[1], map(x->speciespos[x], findall(x->x[2]==2, speciespos)))

    # Ignore reactions where species occurs on both sides of the equation as an observer.
    # Since we are counting loss equations here, only need to remove it from the LHS list.
    for i in intersect(lhspos, rhspos)
        lhspos = deletefirst(lhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    # format is a vector of vectors of symbols.
    losseqns = map(x->vcat(Any[network[x][1]...,network[x][3]]), lhspos)
end

function loss_rate(sp::Symbol, network; return_leqn_unmapped=false, sepvecs=false) 
    #=  
    Input:
        sp: Species for which to calculate the loss rate
        network: either chemical reaction network or transport network
        return_leqn_unmapped: If true, will return leqn
        sepvecs: If true, a vector of expressions will be returned. If false,
                 lval will be returned.
    Output: either:
        leqn: Vector of vectors that contain reactants and a rate coefficient
        lval: symbolic expression for the summed loss rate of species in network, 
              i.e. :(H .* O .* eval(k for O + H))
    =#
    leqn = loss_equations(sp, network) # select relevant reactions from network
    if isempty(leqn)
        throw("NetworkError: $(sp) is missing a loss pathway in the chemical reaction network.")
    end

    if return_leqn_unmapped # gives vector of vectors like [:CO2, :O, :(k)]
        return leqn 
    end

    if sepvecs # returns a vector of expressions
        return map(x->:(*($(x...))), leqn)
    else  # returns one massive expression
        lval = make_net_change_expr(leqn)
        return lval
    end
end

function make_net_change_expr(network_vectors) 
    #=
    Input:
        network_vectors: a subset of the chemical or transport network, pre-selected for a 
                         specific species and either production or loss. format is a vector 
                         of vectors of symbols; each vector of symbols is a specific 
                         reaction rate or transport terms.

    Output: A single, massive expression that gives the net rate of change due to the 
            equations in network_vectors by multiplying the densities of each
            reactant with the reaction rate coefficient (chemistry) or transport
            coefficient (transport) and summing the products.

    I.e. to create a term for the rate of change of the CO2 density due to the
    reactions [[[:CO, :O], [:CO2], :(k1)], [[:CO2], [:CO, :O], :(J1)]], this function 
    would generate the expression :(CO * O * k1 + CO2 * J1). 
    =#
    net_change = :(+($(map(x->:(*($(x...))), network_vectors)...)))
    return net_change
end

function production_equations(sp::Symbol, network) 
    #=  
    Input:
        sp: Species for which to construct a production equation
        network: The type of production process to consider, i.e. either a chemical reaction network or a transport network
    Output:
        prodeqns: production equations and relevant rate coefficients for species sp.
                  the form is an array where each entry is of the form [reactants..., rate].
                  For example, [[:O2, :JO2toOpO], [:O1D, :O2, :k]] are two possible entries 
                  for production of O2 (k will be some more complicated expression, but not important for this example).

    Automatically accounts for cases where a species occurs twice on the RHS by
    reporting those reactions twice.
    =#

    speciespos = getpos(network, sp) 
    lhspos = map(x->x[1], map(x->speciespos[x], findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1], map(x->speciespos[x], findall(x->x[2]==2, speciespos)))

    for i in intersect(rhspos, lhspos)
        rhspos = deletefirst(rhspos, i)
    end

    prodeqns = map(x->vcat(Any[network[x][1]...,network[x][3]]), rhspos)

    return prodeqns
end

function production_rate(sp::Symbol, network; return_peqn_unmapped=false, sepvecs=false) 
    #= 
    Same as loss_rate but for production.
    =#

    peqn = production_equations(sp, network)
    if isempty(peqn)
        throw("NetworkError: $(sp) is missing a production pathway in the chemical reaction network.")
    end


    if return_peqn_unmapped
        return peqn 
    end

    if sepvecs
        return map(x->:(*($(x...))), peqn)  
    else
        pval = make_net_change_expr(peqn)
        return pval
    end
end

function solve_sparse(A, b)
    #=
    A faster way to solve the Jacobian
    =#
    LU = ilu(A, τ = 0.1) # get incomplete LU preconditioner
    x = bicgstabl(A, b, 2, Pl = LU, reltol=1e-25)
    return x
end

function subtract_difflength(a::Array, b::Array)
    #=
    A very specialized function that accepts two vectors, a and b, sorted
    by value (largest first), of differing lengths. It will subtract b from a
    elementwise up to the last index where they are equal, and then add any 
    extra values in a, and subtract any extra values in b.

    Used exclusively in ratefn_local. 

    a: production rates 
    b: loss rates

    both a and b should be positive for the signs to work!
    =#

    shared_size = min(size(a), size(b))[1]

    extra_a = 0
    extra_b = 0
    if shared_size < length(a)
        extra_a += sum(a[shared_size+1:end])
    end

    if shared_size < length(b)
        extra_b += sum(b[shared_size+1:end])
    end

    return sum(a[1:shared_size] .- b[1:shared_size]) + extra_a - extra_b
end

function update_Jrates!(n_cur_densities::Dict{Symbol, Vector{Array{ftype_ncur}}},
                        n_horiz::Int64; nlambda=2000, globvars...)
    #=
    Updates photolysis rates (Jrates) in n_cur_densities for each altitude and horizontal column
    considering altitude distribution of absorbing species.

    Input:
        n_cur_densities: The present atmospheric state. This will be updated to include Jrates by this function.
    n_cur_densities structure: 
        Dict{Symbol, Vector{Array{ftype_ncur}}} (species -> [horizontal columns][altitudes])
    =#

    GV = values(globvars)
    required = [:absorber, :dz, :crosssection, :Jratelist, :num_layers, :solarflux, :enable_horiz_transport]
    check_requirements(keys(GV), required)

    # Calculate optical depth (now 2-D)
    solarabs = optical_depth(n_cur_densities; n_horiz=n_horiz, globvars...)
    # solarabs now records the total optical depth of the atmosphere at each wavelength and altitude

     # Determine whether solar flux is provided per column.  If a single array is
    # given, replicate it for all columns to maintain backwards compatibility.
    flux_per_column = if GV.solarflux isa Vector
        GV.solarflux
    else
        [GV.solarflux for _ in 1:n_horiz]
    end

    # Actinic flux at each wavelength and each horizontal column is the solar
    # flux for that column diminished by the integrated optical depth above that
    # layer.
    for ihoriz in 1:n_horiz
        col_flux = flux_per_column[ihoriz][:, 2]
        for ialt in 1:GV.num_layers
            solarabs[ihoriz][ialt] .= col_flux .* exp.(-solarabs[ihoriz][ialt])
        end
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

    # Initialize and calculate Jrates independently for each horizontal column
    for j in GV.Jratelist
        # Skip J-rates lacking cross sections or whose absorber species is absent
        if !(haskey(GV.crosssection, j) && haskey(n_cur_densities, GV.absorber[j]))
            continue
        end

        n_cur_densities[j] = [zeros(Float64, GV.num_layers) for ihoriz in 1:n_horiz]

        for ihoriz in 1:n_horiz
            for ialt in 1:GV.num_layers
                n_cur_densities[j][ihoriz][ialt] = ftype_ncur(
                    BLAS.dot(nlambda, solarabs[ihoriz][ialt], 1, GV.crosssection[j][ihoriz][ialt+1], 1) # updated update_Jrates! so Jrates integrate column-specific cross sections
                )
            end
        end
    end

    # # Verify identical Jrates when horizontal transport is disabled
    # if n_horiz > 1 && !GV.enable_horiz_transport
    #     for j in GV.Jratelist
    #         if haskey(n_cur_densities, j)
    #             first = n_cur_densities[j][1]
    #             for col in n_cur_densities[j][2:end]
    #                 @assert isapprox(col, first; rtol=1e-8, atol=0.0) "Jrate mismatch across columns for $(j)"
    #             end
    #         end
    #     end
    # end
end

#===============================================================================#
#                             Escape functions                                  #
#===============================================================================#

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
    required =  [:zmax, :M_P, :R_P]
    check_requirements(keys(GV), required)
    
    # lambda is the Jeans parameter (Gronoff 2020), basically the ratio of the 
    # escape velocity GmM/z to the thermal energy, kT.
    lambda = (m*mH*bigG*GV.M_P)/(kB*Texo*(GV.R_P+GV.zmax))
    vth = sqrt(2*kB*Texo/(m*mH))
    v = exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)

    return v
end

# Nonthermal escape functions: 
function escape_probability(sp, atmdict, ihoriz; globvars...)::Array
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
    required = [:all_species, :collision_xsect, :dz, :planet]
    check_requirements(keys(GV), required)
    
    # Parameters determined through Bethan's Monte Carlo model. 
    # [1] = A = escape probability at altitude where above column = 0, for high energy particles. upper limit
    # [2] = a = how "transparent" the atmosphere is to an escaping atom. smaller for higher energy so this is for an upper limit.
    params = Dict("Mars"=>[0.916, 0.039],
                  "Venus"=>[0.868, 0.058]
                 )[GV.planet]

    totdens = n_tot(atmdict, ihoriz; GV.all_species, GV.dz)
    return params[1] .* exp.(-params[2] .* GV.collision_xsect[sp] .* column_density_above(totdens; globvars...))
end

function escaping_hot_atom_production(sp, source_rxns, source_rxn_rc_funcs, atmdict, Mtot, ihoriz; returntype="array", globvars...)
    #=
    Solves the equation k[R1][R2] * P to get the total volume escape of hot atoms of species sp
    from the exobase region where P is the escape probability.
    
    Input
        sp: species
        source_rxns: reaction network that will cause hot atoms to be produced
        atmdict: present atmospheric state dictionary
        Mtot: total atmospheric density array
	ihoriz: Vertical column index
    Output: 
        array of production by altitude (rows) and reaction  (columns)
    =#
    
    GV = values(globvars)
    required = [:all_species, :alt, :collision_xsect, :ion_species, :Jratedict, :molmass, :non_bdy_layers, :num_layers,  
                :n_alt_index, :Tn, :Ti, :Te, :dz]
    check_requirements(keys(GV), required)

    produced_hot = volume_rate_wrapper(sp, source_rxns, source_rxn_rc_funcs, atmdict, Mtot, ihoriz; returntype="array", zmax=GV.alt[end], globvars...) 

    # Returns an array where rows represent altitudes and columns are reactions. Multiplies each vertical profile (each column) by escape_probability. 
    if returntype=="array" # Used within the code to easily calculate the total flux later on.
        return produced_hot .* escape_probability(sp, atmdict, ihoriz; globvars...)
    elseif returntype=="df" # Useful if you want to look at the arrays yourself.
        return DataFrame(produced_hot .* escape_probability(sp, atmdict, ihoriz; globvars...),
                         vec([format_chemistry_string(r[1], r[2]) for r in source_rxns]))
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
    required = [:dz]
    check_requirements(keys(GV), required)

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

#===============================================================================#
#                        Temperature profile functions                          #
#===============================================================================#

# TO DO: Ideally the temperature profile would be generic, but currently that's hard to do because
# we are loading stuff from a file for Venus. 
function T_Mars(Tsurf, Tmeso, Texo; lapserate=-1.4e-5, z_meso_top=108e5, weird_Tn_param=8, globvars...)
    #= 
    Input:
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
    Opt inputs:
        lapserate: adiabatic lapse rate for the lower atmosphere. 1.4e-5 from Zahnle+2008, accounts for dusty atmo.
        z_meso_top: height of the top of the mesosphere (sometimes "tropopause")
        weird_Tn_param: part of function defining neutral temp in upper atmo. See Krasnopolsky 2002, 2006, 2010.
    Output: 
        Arrays of temperatures in K for neutrals, ions, electrons
    =#
    GV = values(globvars)
    required = [:alt]
    check_requirements(keys(GV), required)

    # Altitudes at which various transitions occur -------------------------------
    z_meso_bottom = GV.alt[searchsortednearest(GV.alt, (Tmeso-Tsurf)/(lapserate))]
    
    # These are the altitudes at which we "stitch" together the profiles 
    # from fitting the tanh profile in Ergun+2015,2021 to Hanley+2021 DD8
    # ion and electron profiles, and the somewhat arbitary profiles defined for
    # the region roughly between z_meso_top and the bottom of the fitted profiles.
    z_stitch_electrons = 142e5
    z_stitch_ions = 164e5
    
    # Various indices to define lower atmo, mesosphere, and atmo -----------------
    i_lower = findall(z->z < z_meso_bottom, GV.alt)
    i_meso = findall(z->z_meso_bottom <= z <= z_meso_top, GV.alt)
    i_upper = findall(z->z > z_meso_top, GV.alt)
    i_meso_top = findfirst(z->z==z_meso_top, GV.alt)
    i_stitch_elec = findfirst(z->z==z_stitch_electrons, GV.alt)
    if i_stitch_elec === nothing
        i_stitch_elec = length(GV.alt)
    end
    i_stitch_ions = findfirst(z->z==z_stitch_ions, GV.alt)
    if i_stitch_ions === nothing
        i_stitch_ions = length(GV.alt)
    end

    function NEUTRALS()
        function upper_atmo_neutrals(z_arr)
            @. return Texo - (Texo - Tmeso)*exp(-((z_arr - z_meso_top)^2)/(weird_Tn_param*1e10*Texo))
        end
        Tn = zeros(size(GV.alt))

        if Tsurf==Tmeso==Texo # isothermal case
            Tn .= Tsurf
        else
            Tn[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
            Tn[i_meso] .= Tmeso
            Tn[i_upper] .= upper_atmo_neutrals(GV.alt[i_upper])
        end

        return Tn 
    end 

    function ELECTRONS() 
        function upper_atmo_electrons(z_arr, TH, TL, z0, H0)
            #=
            Functional form from Ergun+ 2015 and 2021, but fit to data for electrons and ions
            in Hanley+ 2021, DD8 data.
            =#

            @. return ((TH + TL) / 2) + ((TH - TL) / 2) * tanh(((z_arr / 1e5) - z0) / H0)
        end
        Te = zeros(size(GV.alt))

        if Tsurf==Tmeso==Texo # isothermal case
            Te .= Tsurf
        else
            Te[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
            Te[i_meso] .= Tmeso

            # This region connects the upper atmosphere with the isothermal mesosphere
            Te[i_meso_top+1:i_stitch_elec] .= upper_atmo_electrons(GV.alt[i_meso_top+1:i_stitch_elec], -1289.05806755, 469.31681082, 72.24740123, -50.84113252)

            # Don't let profile get lower than specified meso temperature
            Te[findall(t->t < Tmeso, Te)] .= Tmeso

            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            Te[i_stitch_elec:end] .= upper_atmo_electrons(GV.alt[i_stitch_elec:end], 1409.23363494, 292.20319103, 191.39012079, 36.64138724)
        end 

        return Te
    end

    function IONS()
        function upper_atmo_ions(z_arr)
            #=
            fit to Gwen's DD8 profile, SZA 40, Hanley+2021, with values M = 47/13, B = -3480/13
            =#
            # New values for fit to SZA 60,Hanley+2022:
            M = 3.40157034 
            B = -286.48716122
            @. return M*(z_arr/1e5) + B
        end
        
        function meso_ions_byeye(z_arr)
            # This is completely made up! Not fit to any data!
            # It is only designed to make a smooth curve between the upper atmospheric temperatures,
            # which WERE fit to data, and the mesosphere, where we demand the ions thermalize.

            # Values used for Gwen's DD8 profile, SZA 40:  170/49 * (z/1e5) -11990/49
            @. return 136/49 * (z_arr/1e5) -9000/49 # These values are for the new fit to SZA 60, Hanley+2022. (11/2/22)
        end

        Ti = zeros(size(GV.alt))

        if Tsurf==Tmeso==Texo # isothermal case
            Ti .= Tsurf
        else
            Ti[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
            Ti[i_meso] .= Tmeso

            # try as an average of neutrals and electrons. There is no real physical reason for this.
            Ti[i_meso_top:i_stitch_ions] = meso_ions_byeye(GV.alt[i_meso_top:i_stitch_ions])

            # Don't let profile get lower than specified meso temperature
            Ti[findall(t->t < Tmeso, Ti)] .= Tmeso

            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            Ti[i_stitch_ions:end] .= upper_atmo_ions(GV.alt[i_stitch_ions:end])
        end 
        
        return Ti
    end 

    return Dict("neutrals"=>NEUTRALS(), "ions"=>IONS(), "electrons"=>ELECTRONS())
end 

function T_Venus(Tsurf::Float64, Tmeso::Float64, Texo::Float64, file_for_interp; z_meso_top=80e5, lapserate=-8e-5, weird_Tn_param=8, globvars...)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    
    Uses the typical temperature structure for neutrals, but interpolates temperatures
    for ions and electrons above 108 km according to the profiles in Fox & Sung 2001.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:alt])
    
    # Subroutines -------------------------------------------------------------------------------------

    function upperatmo_i_or_e(z, particle_type; Ti_array=Ti_interped, Te_array=Te_interped, select_alts=new_a)
        #=
        Finds the index for altitude z within new_a
        =#
        i = something(findfirst(isequal(z), select_alts), 0)
        returnme = Dict("electron"=>Te_array[i], "ion"=>Ti_array[i])
        return returnme[particle_type]
    end

    function NEUTRALS(new_a)
        Tn = zeros(size(GV.alt))

        Tn[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Tn[i_meso] .= Tmeso

        # upper atmo
        Tfile = readdlm(file_for_interp, ',')
        T_n = Tfile[2,:]
        alt_i = Tfile[1,:] .* 1e5
        interp_ion = LinearInterpolation(alt_i, T_n)
        Tn_interped = [interp_ion(a) for a in new_a];
        Tn[i_upper] .= Tn_interped # upper_atmo_neutrals(GV.alt[i_upper])

        return Tn
    end 

    function ELECTRONS(new_a; spc="electron") 
        Te = zeros(size(GV.alt))

        Te[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Te[i_meso] .= Tmeso

        # Upper atmo
        Tfile = readdlm(file_for_interp, ',')
        T_e = Tfile[4,:]
        alt_e = Tfile[1,:] .* 1e5
        interp_elec = LinearInterpolation(alt_e, T_e)
        Te_interped = [interp_elec(a) for a in new_a];
        Te[i_upper] .= Te_interped

        return Te
    end

    function IONS(new_a; spc="ion") 
        Ti = zeros(size(GV.alt))

        Ti[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Ti[i_meso] .= Tmeso

        # Upper atmo
        Tfile = readdlm(file_for_interp, ',')
        T_i = Tfile[3,:]
        alt_i = Tfile[1,:] .* 1e5
        interp_ion = LinearInterpolation(alt_i, T_i)
        Ti_interped = [interp_ion(a) for a in new_a];
        Ti[i_upper] .= Ti_interped

        return Ti
    end

    # Define mesosphere scope.
    z_meso_bottom = GV.alt[searchsortednearest(GV.alt, (Tmeso-Tsurf)/(lapserate))]

    # Various indices to define lower atmo, mesosphere, and atmo -----------------
    i_lower = findall(z->z < z_meso_bottom, GV.alt)
    i_meso = findall(z->z_meso_bottom <= z <= z_meso_top, GV.alt)
    i_upper = findall(z->z > z_meso_top, GV.alt)
    # i_meso_top = findfirst(z->z==z_meso_top, GV.alt)

    # For interpolating upper atmo temperatures from Fox & Sung 2001 - only from 90 km up
    interp_alts = collect(90e5:GV.alt[2]-GV.alt[1]:GV.alt[end])

    return Dict("neutrals"=>NEUTRALS(interp_alts), "ions"=>IONS(interp_alts), "electrons"=>ELECTRONS(interp_alts))
end


#===============================================================================#
#                   Transport and boundary condition functions                  #
#===============================================================================#

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
    Calculate the bindary diffusion coefficient for species sp, b = AT^s.

    Currently, this is set up to only work for diffusion through CO2 since that's the Mars atm.
    Could be extended to be for any gas, but that will require some work.
    =#
    
    A = diffparams(sp)[1] .* 1e17 # Empirical parameter (determined it by experiment)
    s = (diffparams(sp)[2]) # Empirical parameter (det. by exp.)
    b = A .* T .^ s # T = temperature
    return b
end

function boundaryconditions(fluxcoef_dict, atmdict, M, n_horiz::Int64; nonthermal=true, globvars...)
    #=
    Inputs:
        fluxcoef_dict: a dictionary containing the K and D flux coefficients for every species throughout
                       the atmosphere. Format species=>Array(length(all_species), length(alt)).
                       Because it has length alt, it ALREADY INCLUDES boundary layer info in the 
                       1st and last elements. 2nd and penultimate elements are for the edge bulk layers.
        atmdict: Atmospheric state dictionary, required for the nonthermal escape boundary condition.
        M: total atmospheric density, required for the nonthermal escape boundary condition.
	n_horiz: Number of vertical columns in simulation
    Outputs:
        boundary conditions for species in a 2 x 2 matrix, format:
        [n_1 -> n_0, n_0 -> n_1;      
         n_(nl) -> n_(nl+1), n_(nl+1) -> n_(n_l)] for each vertical column

        where n_0 is the boundary layer from [-1 km, 1 km], n_1 is the first bulk layer from [1 km, 3 km],
        n_(nl) is the topmost bulk layer, and n_(nl+1) is the top boundary layer.

        Form of the output, for each species for each vertical column is:

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

        FROM VENUS VERSION, MIKE:
        Sign convention: The density-dependent terms (bc_dict[sp][ihoriz][:, 1]) are multiplied by -1 when the  
                         transport rates are computed in get_transport_PandL_rate. Density independent 
                         terms (bc_dict[sp][ihoriz][:, 2]) are not.
    =#
    
    GV = values(globvars)
    required = [:all_species, :speciesbclist, :dz, :planet]
    check_requirements(keys(GV), required)
    
    bc_dict = Dict{Symbol, Vector{Array{ftype_ncur}}}([s=>[[0 0; 0 0] for ihoriz in 1:n_horiz] for s in GV.all_species])

    for sp in keys(GV.speciesbclist)
        try 
            global these_bcs = GV.speciesbclist[sp]
        catch KeyError
            println("No entry $(sp) in bcdict")
            continue
        end

        for ihoriz in 1:n_horiz
 
            # DENSITY
            try 
                # lower boundary...
                if GV.planet=="Mars"
                    n_lower = [fluxcoef_dict[sp][ihoriz][2, :][1], fluxcoef_dict[sp][ihoriz][1, :][2]*these_bcs["n"][ihoriz][1]]
                elseif GV.planet=="Venus"
                    # get the eddy+molecular mixing velocities at the lower boundary of the atmosphere
                    v_lower_boundary_up = fluxcoef_dict[sp][ihoriz][1, # lower boundary cell, outside atmosphere
                                                                    2] # upward mixing velocity
                    v_lower_boundary_dn = fluxcoef_dict[sp][ihoriz][2, # bottom cell of atmosphere
                                                                    1] # downward mixing velocity

                    n_lower = [v_lower_boundary_dn, v_lower_boundary_up*these_bcs["n"][ihoriz][1]]
                end
                try
                    # TODO for Venus only, this note from Mike originally: 
                    # throw an error if density boundary condition
                    # is specified simultaneous with any flux or velocity condition
                    @assert all(x->!isnan(x), n_lower)
                    bc_dict[sp][ihoriz][1, :] .+= n_lower
                catch y
                    if !isa(y, AssertionError)
                        throw("Unhandled exception in lower density bc: $(y)")
                    end
                end

                # upper boundary...
                try 
                    if GV.planet=="Mars"
                        n_upper = [fluxcoef_dict[sp][ihoriz][end-1, :][2], fluxcoef_dict[sp][ihoriz][end, :][1]*these_bcs["n"][ihoriz][2]]
                    elseif GV.planet=="Venus"
                        # get the eddy+molecular mixing velocities at the upper boundary of the atmosphere
                        v_upper_boundary_up = fluxcoef_dict[sp][ihoriz][end-1, # top cell of atmosphere
                                                                        2]     # upward mixing velocity
                        v_upper_boundary_dn = fluxcoef_dict[sp][ihoriz][end, # upper boundary cell, outside atmosphere
                                                                        1]   # downward mixing velocity

                        n_upper = [v_upper_boundary_up, v_upper_boundary_dn*these_bcs["n"][ihoriz][2]]
                    
                    end

                    # TODO for Venus only, this note from Mike originally: 
                    # throw an error if density boundary condition
                    # is specified simultaneous with any flux or velocity condition
                    @assert all(x->!isnan(x), n_upper)
                    bc_dict[sp][ihoriz][2, :] .+= n_upper
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
                # lower boundary...
                if GV.planet=="Mars"
                    f_lower = [0, -these_bcs["f"][ihoriz][1]/GV.dz]
                elseif GV.planet=="Venus"
                    f_lower = [0, these_bcs["f"][ihoriz][1]/GV.dz]
                    #             ^ no (-) sign, negative flux at lower boundary represents loss to surface
                end
                try        
                    @assert all(x->!isnan(x), f_lower)
                    bc_dict[sp][ihoriz][1, :] .+= f_lower
                catch y
                    if !isa(y, AssertionError)
                        throw("Unhandled exception in lower flux bc: $(y)")
                    end
                end
                try 
                    # upper boundary...
                    f_upper = [0, -these_bcs["f"][ihoriz][2]/GV.dz]
                    #             ^ (-) sign needed so that positive flux at upper boundary represents loss to space
                    #             (see "Sign convention" note above)
                    @assert all(x->!isnan(x), f_upper)
                    bc_dict[sp][ihoriz][2, :] .+= f_upper
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
                # lower boundary...
                if GV.planet=="Mars" 
                    v_lower = [these_bcs["v"][ihoriz][1]/GV.dz, 0]
                elseif GV.planet=="Venus"
                    v_lower = [-these_bcs["v"][ihoriz][1]/GV.dz, 0]
                #          ^ (-) sign needed so that negative velocity at lower boundary represents loss to surface
                #          (see "Sign convention" note above)
                end

                try
                    @assert all(x->!isnan(x), v_lower)
                    bc_dict[sp][ihoriz][1, :] .+= v_lower
                catch y
                    if !isa(y, AssertionError)
                        throw("Unhandled exception in lower velocity bc: $(y)")
                    end
                end

                try 
                    # upper boundary...
                    v_upper = [these_bcs["v"][ihoriz][2]/GV.dz, 0]
                    #          ^ no (-) sign needed,  positive velocity at upper boundary represents loss to space
                    #          (see "Sign convention" note above)
                    @assert all(x->!isnan(x), v_upper)
                    bc_dict[sp][ihoriz][2, :] .+= v_upper
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
    end
    
    # SPECIAL CASE: add on the non-thermal escape for H and D. 
    if nonthermal
        required = [:hot_H_network, :hot_D_network, :hot_H_rc_funcs, :hot_D_rc_funcs, 
                    :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Jratedict]
        check_requirements(keys(GV), required)

        for ihoriz in 1:n_horiz
	    if :H in keys(atmdict) && :D in keys(atmdict) && :H2 in keys(atmdict) && :HD in keys(atmdict) 
               prod_hotH = escaping_hot_atom_production(:H, GV.hot_H_network, GV.hot_H_rc_funcs, atmdict, M, ihoriz; globvars...)
               prod_hotD = escaping_hot_atom_production(:D, GV.hot_D_network, GV.hot_D_rc_funcs, atmdict, M, ihoriz; globvars...)
               prod_hotH2 = escaping_hot_atom_production(:H2, GV.hot_H2_network, GV.hot_H2_rc_funcs, atmdict, M, ihoriz; globvars...)
               prod_hotHD = escaping_hot_atom_production(:HD, GV.hot_HD_network, GV.hot_HD_rc_funcs, atmdict, M, ihoriz; globvars...)

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

                bc_dict[:H][ihoriz][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_H_network, prod_hotH; returntype="number", globvars...)]
            	bc_dict[:D][ihoriz][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_D_network, prod_hotD; returntype="number", globvars...)]
            	bc_dict[:H2][ihoriz][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_H2_network, prod_hotH2; returntype="number", globvars...)]
            	bc_dict[:HD][ihoriz][2, :] .+= [0, -(1/GV.dz)*nonthermal_escape_flux(GV.hot_HD_network, prod_hotHD; returntype="number", globvars...)]
	 #   else
	 #       println("At least one of H, D, H2 and HD is not in the model.",'\n')
            end
	end
    end
    return bc_dict
end

# function boundaryconditions_horiz(globvars)
function boundaryconditions_horiz(
    atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}},
    horiz_wind_v::Vector{Vector{Float64}},
    n_horiz::Int64;
    cyclic::Bool=true,
    globvars...
)
    #= 
    Inputs:
        M: total atmospheric density, required for the nonthermal escape boundary condition.
    Outputs:
        boundary conditions for species in a 2 x 2 matrix, format:  
        [n_1 -> n_0, n_0 -> n_1;      
         n_(nhoriz) -> n_(nhoriz+1), n_(nhoriz+1) -> n_(nhoriz)] for each bulk layer altitude bin

        where n_0 is outside the model behind the back edge, n_1 is the first vertical column (ihoriz=1),
        n_(nhoriz) is the front-most verical column (ihoriz=n_horiz), and n_(nhoriz+1) is outside the model in front of the front edge.

        Form of the output, for each species for each vertical column is:

         Back edge   [←, →;     [density-dependent, density-independent;    [#/s, #/cm³/s.;
         Front edge   →, ←]      density-dependent, density-independent]     #/s, #/cm³/s]
         where ← is backwards from higher to lower values of ihoriz; → is forwards from lower to higher

        Each row has two elements:
            1st element: inside model  -> outside model (depends on species concentration in cell)
            2nd element: outside model -> inside moel (independent of species concentration in cell)

            note, technically, these are chemical equations.

        More specifically, when the return value of this function is used in other functions, the first
        element in each row will eventually be multiplied by a density taken from the atmospheric 
        state dictionary, and the second element will be used as-is. That way, eventually the total
        change recorded in other functions is always #/cm³/s. 

        FROM VENUS VERSION, MIKE:
        Sign convention: The density-dependent terms (bc_dict[sp][ialt][:, 1]) are multiplied by -1 when the 
                         transport rates are computed in get_transport_PandL_rate. Density independent
                         terms (bc_dict[sp][ialt][:, 2]) are not.

        However, please note that the model is currently set up to use zero flux edge boundary conditions only. The above comments have been left for future development and flexibility.
        By default the dictionary ``speciesbclist_horiz`` in ``MODEL_SETUP.jl`` supplies
        zero-flux edge conditions for all species.  You may override these values
        to impose custom influxes or outfluxes at either edge.  The sign
        convention is the same as for vertical boundary conditions: positive
        numbers inject material and negative numbers remove it.
    =#
    
    GV = values(globvars)
    required = [:all_species, :speciesbclist_horiz, :dx, :planet]
    check_requirements(keys(GV), required)
    
    # bc_dict_horiz = Dict{Symbol, Vector{Array{ftype_ncur}}}([s=>[[0.0 0.0; 0.0 0.0] for ialt in 1:GV.num_layers] for s in GV.all_species])
    bc_dict_horiz = Dict{Symbol, Vector{Array{ftype_ncur}}}(
        [s => [fill(0.0, 2, 2) for ialt in 1:GV.num_layers] for s in GV.all_species]
    )

    for sp in keys(GV.speciesbclist_horiz)
        these_bcs_horiz = GV.speciesbclist_horiz[sp]

        for ialt in 1:GV.num_layers
            if cyclic
                # Periodic domain: no exchange with the environment
                bc_dict_horiz[sp][ialt] .= 0.0
            else
                back_flux  = these_bcs_horiz["f"][1][ialt]
                front_flux = these_bcs_horiz["f"][2][ialt]
                # if GV.planet == "Mars"
                #     f_backedge = [0, -back_flux / GV.dx]
                # elseif GV.planet == "Venus"
                #     f_backedge = [0, back_flux / GV.dx]
                # end
                f_backedge  = [0, -back_flux / GV.dx]
                f_frontedge = [0,  front_flux / GV.dx]
                bc_dict_horiz[sp][ialt][1, :] .+= f_backedge
                bc_dict_horiz[sp][ialt][2, :] .+= f_frontedge
            end
        end
    end
    

    return bc_dict_horiz
end

function Dcoef_neutrals(z, sp::Symbol, b, atmdict::Dict{Symbol, Vector{ftype_ncur}}; ihoriz=1, globvars...)
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
    required = [:all_species, :n_alt_index]
    check_requirements(keys(GV), required)

    if (typeof(z)==Float64) & (typeof(b)==Float64)
        return b ./ n_tot(atmdict, z, ihoriz; GV.all_species, GV.n_alt_index)
    else
        return b ./ n_tot(atmdict, ihoriz; GV.all_species, GV.n_alt_index)
    end
end

function Dcoef!(D_arr, T_arr_2D, sp::Symbol, atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, n_horiz::Int64; globvars...)
    #=
    Calculates the molecular diffusion coefficient for an atmospheric layer.
    For neutrals, returns D = AT^s/n, from Banks and Kockarts Aeronomy, part B, pg 41, eqn 
    15.30 and table 15.2 footnote.
    For ions, it returns ambipolar diffusion coefficients according to Krasnopolsky 2002 and 
    Schunk & Nagy equation 5.55. Yes, the equation is done in one line and it's ugly, but it works.
    Units: cm^2/s

    Inputs:
        D_arr: the container for the diffusion coefficients for ONE species.
               D_arr[ihoriz] is a 1D array with altitude for column ihoriz.
        T_arr_2D: 2D temperature (K), shape (n_horiz, num_layers+2).
                  Neutral temp if sp is neutral; plasma temp if sp is ion.
        sp: whichever species we are calculating for
        atmdict: state of the atmosphere; should include boundary layers, i.e.
                 be the result of calling atmdict_with_boundary_layers.
        n_horiz: number of vertical columns
    Outputs:
        D_arr: The same container, now filled with the diffusion coefficients by altitude
               for this species, for each column.

    =#

    GV = values(globvars)
    required = [
        :all_species, :molmass, :neutral_species, :n_alt_index, :polarizability,
        :q, :speciesbclist, :use_ambipolar, :use_molec_diff
    ]
    check_requirements(keys(GV), required)

    # Loop over each vertical column
    for ihoriz in 1:n_horiz

        # Extract this column's temperature profile
        local T_col = T_arr_2D[ihoriz, :]

        # If using molecular diffusion, compute D
        if GV.use_molec_diff == true
            # Calculate as if it was a neutral - not using function above because this is faster than going into 
            # the function and using an if/else block since we know we'll always have vectors in this case.
            # This equation is: D = b/n
            D_arr[ihoriz] .= (binary_dcoeff_inCO2(sp, T_col)) ./ n_tot(atmdict, ihoriz; GV.all_species, GV.n_alt_index)
        else
            D_arr[ihoriz] .= 0 
        end

        # If ambipolar diffusion is turned on, and this is an ion, overwrite with ambipolar D
        if GV.use_ambipolar==true
            if charge_type(sp) == "ion"
                # Prepare arrays to store the sum over collisions
                # (We build them for each column to match T_col's shape.)
                species_density = zeros(size(T_col))
                sum_nu_in      = zeros(size(T_col))

                # mi = GV.molmass[sp] .* mH
                # create the sum of nu_in. Note that this depends on density, but we only have density for the real layers,
                # so we have to assume the density at the boundary layers is the same as at the real layers.
                for n in GV.neutral_species
                    # Add the neutral's density profile for this column
                    species_density .= atmdict[n][ihoriz]

                    # Optionally override with boundary condition if specified
                    if haskey(GV.speciesbclist, n)
                        if haskey(GV.speciesbclist[n], "n")
                            # lower boundary
                            if !isnan(GV.speciesbclist[n]["n"][ihoriz][1])
                                species_density[1] = GV.speciesbclist[n]["n"][ihoriz][1]
                            end
                            # upper boundary 
                            if !isnan(GV.speciesbclist[n]["n"][ihoriz][2])
                                species_density[end] = GV.speciesbclist[n]["n"][ihoriz][2]
                            end
                        end
                    end

                    sum_nu_in .+= 2π .* sqrt.(
                        (GV.polarizability[n] .* (GV.q^2)) ./
                        reduced_mass(GV.molmass[sp], GV.molmass[n])
                    ) .* species_density
                end

                # Finally set D_ambipolar = kB * T / (m_sp * sum_nu_in)
                D_arr[ihoriz] .= (kB .* T_col) ./ (GV.molmass[sp] .* mH .* sum_nu_in)
            end
        end

    end # for ihoriz
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

function fluxcoefs(sp::Symbol, Kv, Dv, H0v, ihoriz::Int64; globvars...)
    #= 
    base function to generate flux coefficients of the transport network. 
    
    For all the arrays, length = num_layers. "n" typically refers to "neutrals", and "p" to plasma.

    Inputs:
        sp: species symbol 
        z: altitude array in cm.
        dz: altitude layer thickness ("resolution")
        for all the following, length = num_layers 
        Kv: eddy diffusion coefficient
        Dv: molecular diffusion coefficient
        Tv_n: neutral temperature (2D, size num_layers × n_horiz)
        Tv_p: plasma temperature (2D, size num_layers × n_horiz)
        Hsv: scale height by species
        H0v: mean atmospheric scale height
        ihoriz: vertical column index
    Outputs:
        Arrays of coefficients (units 1/s) at each atmospheric layer for downward and upward flux.
        Note that even though it's defined as being between a layer and the one above or below, the value is 
        evaluated at the center of the layer 

    v just means "vector"
    u refers to "upper" (a given layer coupled to the one above)
    l refers to "lower" (a given layer coupled to the one below)
    =#

    GV = values(globvars)
    required = [:Tn, :Tp, :Hs_dict, :n_all_layers, :dz, :planet]
    check_requirements(keys(GV), required)

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

    # Local column references
    Tn_col = GV.Tn[ihoriz, :]   # 1D array for neutrals
    Tp_col = GV.Tp[ihoriz, :]   # 1D array for plasma
    Hs_col = GV.Hs_dict[sp][ihoriz]

    # Retrieve arrays specific to this vertical column
    thisD  = Dv[sp][ihoriz]
    thisK  = Kv[ihoriz]
    thisH0 = H0v[charge_type(sp)][ihoriz]

    # Calculate the coefficients between this layer and the lower layer.
    Dl[2:end]     = @. (thisD[1:end-1] + thisD[2:end]) /  2.0
    Kl[2:end]     = @. (thisK[1:end-1] + thisK[2:end]) / 2.0
    Tl_n[2:end]   = @. (Tn_col[1:end-1] + Tn_col[2:end]) / 2.0
    Tl_p[2:end]   = @. (Tp_col[1:end-1] + Tp_col[2:end]) / 2.0
    dTdzl_n[2:end]= @. (Tn_col[2:end] - Tn_col[1:end-1]) / GV.dz
    dTdzl_p[2:end]= @. (Tp_col[2:end] - Tp_col[1:end-1]) / GV.dz
    Hsl[2:end]    = @. (Hs_col[1:end-1] + Hs_col[2:end]) / 2.0
    H0l[2:end]    = @. (thisH0[1:end-1] + thisH0[2:end]) / 2.0

    if GV.planet=="Mars"
        # Handle the lower boundary layer:
        Dl[1]      = @. (1 + thisD[1]) /  2.0
        Kl[1]      = @. (1 + thisK[1]) / 2.0
        Tl_n[1]    = @. (1 + Tn_col[1]) / 2.0
        Tl_p[1]    = @. (1 + Tp_col[1]) / 2.0
        dTdzl_n[1] = @. (Tn_col[1] - 1) / GV.dz
        dTdzl_p[1] = @. (Tp_col[1] - 1) / GV.dz
        Hsl[1]     = @. (1 + Hs_col[1]) / 2.0
        H0l[1]     = @. (1 + thisH0[1]) / 2.0
    elseif GV.planet=="Venus"
        # Downward transport away from the lower boundary layer, which is outside the model
        # These should never be used but we need to fill the array
        Dl[1] = Float64(NaN)
        Kl[1] = Float64(NaN)
        Tl_n[1] = Float64(NaN)
        Tl_p[1] = Float64(NaN)
        dTdzl_n[1] = Float64(NaN)
        dTdzl_p[1] = Float64(NaN)
        Hsl[1] = Float64(NaN)
        H0l[1] = Float64(NaN)
    end

    # Upward transport from each altitude to the cell above
    Du[1:end-1]     = @. (thisD[1:end-1] + thisD[2:end]) /  2.0
    Ku[1:end-1]     = @. (thisK[1:end-1] + thisK[2:end]) / 2.0
    Tu_n[1:end-1]   = @. (Tn_col[1:end-1] + Tn_col[2:end]) / 2.0
    Tu_p[1:end-1]   = @. (Tp_col[1:end-1] + Tp_col[2:end]) / 2.0
    dTdzu_n[1:end-1]= @. (Tn_col[2:end] - Tn_col[1:end-1]) / GV.dz
    dTdzu_p[1:end-1]= @. (Tp_col[2:end] - Tp_col[1:end-1]) / GV.dz
    Hsu[1:end-1]    = @. (Hs_col[1:end-1] + Hs_col[2:end]) / 2.0
    H0u[1:end-1]    = @. (thisH0[1:end-1] + thisH0[2:end]) / 2.0

    if GV.planet=="Mars"
        # Handle upper boundary layer:
        Du[end]      = @. (thisD[end] + 1) /  2.0
        Ku[end]      = @. (thisK[end] + 1) / 2.0
        Tu_n[end]    = @. (Tn_col[end] + 1) / 2.0
        Tu_p[end]    = @. (Tp_col[end] + 1) / 2.0
        dTdzu_n[end] = @. (1 - Tn_col[end]) / GV.dz
        dTdzu_p[end] = @. (1 - Tp_col[end]) / GV.dz
        Hsu[end]     = @. (Hs_col[end] + 1) / 2.0
        H0u[end]     = @. (thisH0[end] + 1) / 2.0
    elseif GV.planet=="Venus"
        # Upwards flux from the upper boundary layer, which is outside the model
        # These should never be used but we need to fill the array
        Du[end] = Float64(NaN)
        Ku[end] = Float64(NaN)
        Tu_n[end] = Float64(NaN)
        Tu_p[end] = Float64(NaN)
        dTdzu_n[end] = Float64(NaN)
        dTdzu_p[end] = Float64(NaN)
        Hsu[end] = Float64(NaN)
        H0u[end] = Float64(NaN)
    end


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

function fluxcoefs(species_list::Vector, K, D, H0, n_horiz::Int64; globvars...) 
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
    required = [:Tn, :Tp, :Hs_dict, :n_all_layers, :dz]
    check_requirements(keys(GV), required)
    
    # the return dictionary: Each species has 2 entries for every layer of the atmosphere.
    fluxcoef_dict = Dict{Symbol, Vector{Array{ftype_ncur}}}([s=>[fill(0., GV.n_all_layers, 2) for ihoriz in 1:n_horiz] for s in species_list])

    for s in species_list
        for ihoriz in 1:n_horiz
            layer_below_coefs, layer_above_coefs = fluxcoefs(s, K, D, H0, ihoriz; globvars...)
            fluxcoef_dict[s][ihoriz][:, 1] .= layer_below_coefs
            fluxcoef_dict[s][ihoriz][:, 2] .= layer_above_coefs
        end
    end

    return fluxcoef_dict
end

# function fluxcoefs_horiz(species_list::Vector, horiz_wind_v::Vector{Vector{Float64}}, n_horiz::Int64; globvars...) 
#     New optimized version of fluxcoefs that calls the lower level version of fluxcoefs,
#     producing a dictionary that contains both up and down flux coefficients for each layer of
#     the atmosphere including boundary layers. Created to optimize calls to this function
#     during the solution of the production and loss equation.

#     Here, D and Hs depend on the current atmospheric densities, and need to be pre-calculated
#     within the upper level function which calls this one.
#     The parameters below which vary by species are dictionaries, and those that are arrays
#     don't depend on the species. All profiles are by altitude. All lengths are the same 
#     as for the alt variable (full altitude grid including boundary layers).
    
#     Inputs:
#         species_list: Species for which to generate transport coefficients. This allows the code to only do it for
#                 transport species during the main simulation run, and for all species when trying to plot 
#                 rate balances after the run.
#         T_neutral: 1D neutral temperature profile
#         T_plasma: the same, but for the plasma temperature
#         K: Array; 1D eddy diffusion profile by altitude for current atmospheric state
#         D: Dictionary (key=species); 1D molecular diffusion profiles for current atmospheric state
#         H0: Dictionary (key="neutral" or "ion"); 1D mean atmospheric scale height profiles for each type
#         Hs: Dictionary (key=species); 1D species scale height profiles
#     Outputs:
#         fluxcoef_dict: dictionary of flux coefficients of the form [flux down, flux up] by altitude 

#     =#
"""
Compute horizontal transport coefficients for each species and column.

Each entry ``fluxcoef_dict[s][i]`` contains a ``n_all_layers x 2`` matrix of
coefficients linking column ``i`` to the column behind (column 1) and the column
in front (column 2). Diffusion coefficients are averaged between neighbouring
columns and scaled by ``dx²`` while advection uses an upwind scheme based on the
provided horizontal wind profile.  Passing `cyclic=true` wraps the indices so
that column 1 connects to column ``n_horiz`` and vice versa.
"""
function fluxcoefs_horiz(
    species_list::Vector,
    K::Vector{Vector{ftype_ncur}},
    D::Dict{Symbol, Vector{Vector{ftype_ncur}}},
    horiz_wind_v::Vector{Vector{Float64}},
    n_horiz::Int64;
    cyclic::Bool = true,
    globvars...
)
    GV = values(globvars)
    required = [:dx, :n_all_layers, :enable_horiz_transport]
    check_requirements(keys(GV), required)
    
    # the return dictionary: Each species has 2 entries for every layer of the atmosphere.
    # fluxcoef_horiz_dict = Dict{Symbol, Vector{Array{ftype_ncur}}}([s=>[fill(0., GV.n_all_layers, 2) for ihoriz in 1:n_horiz] for s in species_list])
    fluxcoef_dict = Dict{Symbol, Vector{Array{ftype_ncur}}}(
        [s => [zeros(ftype_ncur, GV.n_all_layers, 2) for _ in 1:n_horiz]
         for s in species_list],
    )

    for s in species_list
        for ihoriz in 1:n_horiz
            if cyclic
                behind_idx  = ihoriz == 1        ? n_horiz : ihoriz - 1
                infront_idx = ihoriz == n_horiz ? 1       : ihoriz + 1
            else
                behind_idx  = ihoriz - 1
                infront_idx = ihoriz + 1
            end
            for ialt in 1:GV.n_all_layers
                diff_back = 0.0
                diff_front = 0.0

                if GV.enable_horiz_transport
                    if behind_idx >= 1
                        K_back = (K[ihoriz][ialt] + K[behind_idx][ialt]) / 2
                        D_back = (D[s][ihoriz][ialt] + D[s][behind_idx][ialt]) / 2
                        diff_back = (K_back + D_back) / GV.dx^2
                    end

                    if infront_idx <= n_horiz
                        K_front = (K[ihoriz][ialt] + K[infront_idx][ialt]) / 2
                        D_front = (D[s][ihoriz][ialt] + D[s][infront_idx][ialt]) / 2
                        diff_front = (K_front + D_front) / GV.dx^2
                    end
                end

                # Calculate velocities at the interfaces for proper upwind scheme
                v_local = horiz_wind_v[ihoriz][ialt]
                v_front = infront_idx <= n_horiz ? horiz_wind_v[infront_idx][ialt] : 0.0
                v_back  = behind_idx >= 1     ? horiz_wind_v[behind_idx][ialt]  : 0.0

                # Interface velocities (average of adjacent cell velocities)
                v_interface_front = (v_local + v_front) / 2  # Velocity at interface i+½
                v_interface_back  = (v_back + v_local) / 2   # Velocity at interface i-½

                adv_front = 0.0
                adv_back  = 0.0
                if GV.enable_horiz_transport
                    # Upwind scheme: flux depends on velocity direction at the interface
                    # adv_front: flux from current column to next column (interface i+½)
                    # adv_back: flux from previous column to current column (interface i-½)
                    adv_front = (v_interface_front > 0 ? v_interface_front : 0.0) / GV.dx
                    adv_back  = (v_interface_back < 0 ? -v_interface_back : 0.0) / GV.dx
                end

                fluxcoef_dict[s][ihoriz][ialt, 1] = diff_back + adv_back
                fluxcoef_dict[s][ihoriz][ialt, 2] = diff_front + adv_front
            end
        end
    end

    # if !GV.enable_horiz_transport && all(all(v .== 0.0) for v in horiz_wind_v)
    #     maxval = maximum(abs, vcat([vec(mat) for col in values(fluxcoef_dict) for mat in col]...))
    #     println("DEBUG: horiz transport disabled, max coefficient = ", maxval)
    # end

    # return fluxcoef_dict_horiz
    return fluxcoef_dict
end

function Keddy(z::Vector, nt; ihoriz=nothing, globvars...)
    #=
    Input:
        z: Altitude in cm
        nt: Total atmospheric density.
            - For single-column case: Vector (#/cm³)
            - For multi-column, this should be a Matrix (n_horiz × #altitudes) with horizontal index first, and `ihoriz` must be specified.
    Output:
        k: eddy diffusion coefficients at all altitudes.
           Units: cm^2/s

    If running multi-column, make sure to specify the column index (ihoriz).
    This version accommodates both single-column and multi-column simulations.
    =#

    GV = values(globvars)
    required = [:planet]
    check_requirements(keys(GV), required)

    k = zeros(size(z))

    if GV.planet == "Mars"
        upperatm = findall(z .> 60e5)

        # Assign constant Keddy below 60 km altitude (lower atmosphere)
        k[z .<= 60e5] .= 10. ^ 6

        # Upper atmosphere (>60 km), altitude-dependent Keddy
        if ndims(nt) == 1  # Single-column case
            k[upperatm] .= 2e13 ./ sqrt.(nt[upperatm])
        elseif ndims(nt) == 2  # Multi-column case
            @assert !isnothing(ihoriz) "ihoriz must be specified for multi-column simulations."
            k[upperatm] .= 2e13 ./ sqrt.(nt[ihoriz, upperatm])
        else
            throw("Invalid dimensions for nt. Expecting Vector or Matrix.")
        end

    elseif GV.planet == "Venus"
        # Venus: just an altitude-dependent eddy diffusion
        if ndims(nt) == 1  # Single-column case
            k .= 8e12 .* nt .^ (-0.5)
        else  # Multi-column case
            k .= 8e12 .* (nt[ihoriz, :] .^ (-0.5))
        end
    end

    return k
end

# thermal diffusion factors
thermaldiff(sp) = get(Dict(:H=>-0.25, :H2=>-0.25, :D=>-0.25, :HD=>-0.25,
                                :He=>-0.25, 
                                :Hpl=>-0.25, :H2pl=>-0.25, :Dpl=>-0.25, :HDpl=>-0.25,
                                :Hepl=>-0.25), sp, 0)

function update_diffusion_and_scaleH(
    species_list,
    atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, 
    n_horiz::Int64;
    globvars...
) 
    #=
    Calculates:
        1) The eddy diffusion coefficient K for each column;
        2) The mean scale height H₀ for each column (both neutral and ion);
        3) The molecular (or ambipolar) diffusion coefficients, stored in Dcoef_dict.

    Inputs:
        - atmdict: Atmospheric state dictionary with boundary layers already included
        - n_horiz: Number of vertical columns
        - species_list: The species for which to generate molecular or ambipolar D.

    Outputs:
        - K: Vector of length n_horiz, each entry is a 1D array of K vs altitude
        - H0_dict: Dict("neutral"=>[H0_neutral_col1, H0_neutral_col2, ...],
                        "ion"=>[H0_ion_col1, H0_ion_col2, ...])
        - Dcoef_dict: Dict(s1 => [[D alt array col1], [D alt array col2], ...],
                            s2 => [...], ...)
    =#
    GV = values(globvars)
    required = [
        :all_species, :alt, :speciesbclist, :M_P, :molmass, :neutral_species,
        :n_alt_index, :polarizability, :planet, :R_P, :q, :Tn, :Tp,
        :Tprof_for_diffusion, :use_ambipolar, :use_molec_diff
    ]
    check_requirements(keys(GV), required)

    # Assumes atmdict already includes boundary layers
    ncur_with_bdys = ncur_with_boundary_layers(
        atmdict, n_horiz; GV.n_alt_index, GV.all_species
    )
    num_alts_with_bdy = length(ncur_with_bdys[GV.all_species[1]][1])

    # 1) Eddy diffusion, explicitly calculated for each column
    K = [
        Keddy(
            GV.alt,
            n_tot(ncur_with_bdys, ihoriz; GV.all_species, GV.n_alt_index);
            ihoriz = ihoriz,
            globvars...
        )
        for ihoriz in 1:n_horiz
    ]

    # 2) Mean atmospheric scale heights (neutral and ion), column-wise explicitly
    H0_dict = Dict{String, Vector{Vector{ftype_ncur}}}()
    H0_dict["neutral"] = [
        scaleH(ncur_with_bdys, GV.Tn[ihoriz, :], n_horiz, ihoriz; globvars...)
        for ihoriz in 1:n_horiz
    ]
    H0_dict["ion"] = [
        scaleH(ncur_with_bdys, GV.Tp[ihoriz, :], n_horiz, ihoriz; globvars...)
        for ihoriz in 1:n_horiz
    ]

    # 3) Molecular/Ambipolar diffusion, explicitly calculated per species per column
    Dcoef_dict = Dict{Symbol, Vector{Vector{ftype_ncur}}}()

    for s in species_list
        # Freshly initialized array per species and column (no reuse)
        D_coefs_for_s = [
            zeros(ftype_ncur, num_alts_with_bdy) for _ in 1:n_horiz
        ]

        # Calculate diffusion coefficient for this species
        Dcoef!(D_coefs_for_s,
               GV.Tprof_for_diffusion[charge_type(s)],  # (num_alts_with_bdy, n_horiz)
               s,
               ncur_with_bdys,
               n_horiz;
               globvars...
        )

        # Store the calculated coefficients
        Dcoef_dict[s] = deepcopy(D_coefs_for_s)
    end

    # ------------------------------------------------------------------
    # Validate shapes
    # ------------------------------------------------------------------
    expected_alt_len = GV.n_all_layers
    @assert length(K) == n_horiz
    @assert all(length(k) == expected_alt_len for k in K)
    @assert all(length(h) == expected_alt_len for h in H0_dict["neutral"])
    @assert all(length(h) == expected_alt_len for h in H0_dict["ion"])
    for s in species_list
        @assert length(Dcoef_dict[s]) == n_horiz
        @assert all(length(vec) == expected_alt_len for vec in Dcoef_dict[s])
    end

    return K, H0_dict, Dcoef_dict
end

function update_transport_coefficients(
    species_list,
    atmdict::Dict{Symbol, Vector{Array{Float64}}},
    D_coefs::Vector{Matrix{Float64}},  # Unused
    M::Matrix{Float64},  # now explicitly 2D: num_layers × n_horiz
    n_horiz::Int;
    calc_nonthermal=true,
    debug=false,
    globvars...
)
    # Call the real function that doesn't expect D_coefs
    return update_transport_coefficients(
        species_list,
        atmdict,
        M,
        n_horiz;
        calc_nonthermal=calc_nonthermal,
        debug=debug,
        globvars...
    )
end

function update_transport_coefficients(
    species_list,
    atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}},
    M::Matrix{ftype_ncur},  # now explicitly 2D: num_layers × n_horiz
    n_horiz::Int64;
    calc_nonthermal=true,
    debug=false,
    globvars...
)
    #=
    Input:
      - species_list: Species that need transport coefficients updated
      - atmdict: Current atmospheric state dictionary (WITH boundary layers) 
      - M: total atmospheric density by altitude and horizontal column (num_layers × n_horiz)
      - n_horiz: Number of columns
      - calc_nonthermal: whether to compute nonthermal BCs (like H, D escape)
    Output:
      - tlower, tup, tdown, tupper: boundary & in‐domain transport coefficients
    =#
    GV = values(globvars)
    required = [
        :all_species, :alt, :speciesbclist, :dz, :hot_H_network, :hot_H_rc_funcs, 
        :hot_D_network, :hot_D_rc_funcs, :hot_H2_network, :hot_H2_rc_funcs,
        :hot_HD_network, :hot_HD_rc_funcs, :Hs_dict, :ion_species, :M_P, :molmass,
        :neutral_species, :non_bdy_layers, :num_layers, :n_all_layers, :n_alt_index,
        :polarizability, :q, :R_P, :Tn, :Ti, :Te, :Tp, :Tprof_for_diffusion,
        :transport_species, :use_ambipolar, :use_molec_diff, :zmax
    ]
    check_requirements(keys(GV), required)

    # 1) Build the K, H0, and Dcoef_dict
    K_eddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(
        species_list, atmdict, n_horiz; globvars...
    )
    # K_eddy_arr is a Vector of length n_horiz; each is 1D array of K vs altitude
    # H0_dict = Dict("neutral" => [...], "ion"=>[...]) each an array of length n_horiz
    # Dcoef_dict[s] = [[D vs alt col1], [D vs alt col2], ...]

    # 2) Build the flux coefficients dictionary
    fluxcoefs_all = fluxcoefs(
        species_list,
        K_eddy_arr,
        Dcoef_dict,
        H0_dict,
        n_horiz; 
        globvars...
    )
    # fluxcoefs_all[s][ihoriz] is a 2D array (num_alt x 2) => e.g. [down, up] at each altitude

    # Prepare arrays to store upward/downward transport in the bulk layers
    # shape = (n_horiz, num_layers, number_of_transport_species)
    tup   = fill(-999., n_horiz, GV.num_layers, length(GV.transport_species))
    tdown = fill(-999., n_horiz, GV.num_layers, length(GV.transport_species))

    # Fill those from fluxcoefs
    for (isp, s) in enumerate(GV.transport_species)
        for ihoriz in 1:n_horiz
            # The inner bulk layers are alt indices 2:end-1 in fluxcoefs
            # with columns: col=1 => downward, col=2 => upward
            # We skip boundary layers => fluxcoefs_all[s][ihoriz][2:end-1, 1 or 2]
            tup[ihoriz, :, isp]   .= fluxcoefs_all[s][ihoriz][2:end-1, 2]
            tdown[ihoriz, :, isp] .= fluxcoefs_all[s][ihoriz][2:end-1, 1]
        end
    end

    # 3) Boundary conditions
    bc_dict = boundaryconditions(
        fluxcoefs_all,
        atmdict,
        M,
        n_horiz;
        nonthermal=calc_nonthermal,
        globvars...
    )
    # bc_dict[s][ihoriz][1, :] => lower boundary [density-dependent, density-indep]
    # bc_dict[s][ihoriz][2, :] => upper boundary

    # Convert bc_dict to tlower and tupper shape
    tlower = Vector{Array{Float64}}(undef, n_horiz)
    tupper = Vector{Array{Float64}}(undef, n_horiz)
    for ihoriz in 1:n_horiz
        # Each is shape (2, number_of_transport_species)
        # so we gather bc_dict[sp][ihoriz][1,:] for each sp in transport_species
        tlower[ihoriz] = permutedims(reduce(hcat, [
            bc_dict[sp][ihoriz][1, :] for sp in GV.transport_species
        ]))
        tupper[ihoriz] = permutedims(reduce(hcat, [
            bc_dict[sp][ihoriz][2, :] for sp in GV.transport_species
        ]))
    end

    if debug
        println("[update_transport_coefficients] vertical transport coefficients:")
        nprint = min(3, GV.num_layers)
        for ihoriz in 1:n_horiz
            println("  column ", ihoriz, " tup[1:" , nprint, "] = ",
                    tup[ihoriz, 1:nprint, :])
            println("  column ", ihoriz, " tdown[1:" , nprint, "] = ",
                    tdown[ihoriz, 1:nprint, :])
        end
        if !GV.enable_horiz_transport && n_horiz > 1
            identical = all(
                arrays_equal_with_nan(tup[ihoriz, :, :], tup[1, :, :]) &&
                arrays_equal_with_nan(tdown[ihoriz, :, :], tdown[1, :, :])
                for ihoriz in 2:n_horiz)
            println("  identical across columns? ", identical)
        end
    end

    return tlower, tup, tdown, tupper
end

function update_horiz_transport_coefficients(species_list, atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, D_coefs, M, n_horiz::Int64;
                                       calc_nonthermal=true, cyclic=true, debug=false, globvars...)
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
        cyclic: wrap the horizontal domain so that material leaving one edge
                enters from the opposite side

    Return: 
        Horizontal transport coefficients for all atmospheric layers, units 1/s

        tbackedge: horizontal transport coefficients at the back edge. shape: a vector of GV.num_layers arrays of size (2 x length(GV.transport_species))
        tforwards: forward-moving (lower to higher vertical column number) coefficients for all columns. shape: length(GV.transport_species) * GV.num_layers
        tbackwards: backward-moving (higher to lower vertical column number) coefficients for all columns. shape: length(GV.transport_species) * GV.num_layers
        tfrontedge: horizontal transport coefficients at the front edge. shape: a vector of GV.num_layers arrays of size (2 x GV.transport_species)
    =#

    GV = values(globvars)
    required = [:all_species, :alt, :speciesbclist, :dx, :hot_H_network, :hot_H_rc_funcs, :hot_D_network, :hot_D_rc_funcs,
               :hot_H2_network, :hot_H2_rc_funcs, :hot_HD_network, :hot_HD_rc_funcs, :Hs_dict,
               :ion_species, :M_P, :molmass, :neutral_species, :non_bdy_layers, :num_layers, :n_all_layers, :n_alt_index,
               :polarizability, :q, :R_P, :Tn, :Ti, :Te, :Tp, :Tprof_for_diffusion, :transport_species,
               :use_ambipolar, :use_molec_diff, :zmax, :horiz_wind_v, :enable_horiz_transport]
    check_requirements(keys(GV), required)

    # Get flux coefficients
    # Calculate diffusion coefficients and scale heights
    K_eddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(
        species_list,
        atmdict,
        n_horiz;
        globvars...
    )

    # Get horizontal flux coefficients
    fluxcoefs_horiz_all = fluxcoefs_horiz(
        species_list,
        K_eddy_arr,
        Dcoef_dict,
        GV.horiz_wind_v,
        n_horiz;
        cyclic=cyclic,
        globvars...
    )

    tforwards  = fill(0.0, n_horiz, GV.num_layers, length(species_list))
    tbackwards = fill(0.0, n_horiz, GV.num_layers, length(species_list))

    for (isp, sp) in enumerate(species_list)
        for ihoriz in 1:n_horiz
            tbackwards[ihoriz, :, isp] .= fluxcoefs_horiz_all[sp][ihoriz][2:end-1, 1]
            tforwards[ihoriz, :, isp]  .= fluxcoefs_horiz_all[sp][ihoriz][2:end-1, 2]
        end
    end

    # Validate horizontal coefficient shapes
    expected_flux_shape = (GV.n_all_layers, 2)
    @assert all(all(size(mat) == expected_flux_shape for mat in mats)
                for mats in values(fluxcoefs_horiz_all))
    expected_tb_shape = (n_horiz, GV.num_layers, length(species_list))
    @assert size(tforwards) == expected_tb_shape
    @assert size(tbackwards) == expected_tb_shape

    bc_dict_horiz = boundaryconditions_horiz(atmdict, GV.horiz_wind_v, n_horiz; cyclic=cyclic, globvars...)

    # transport coefficients for boundaries
    tbackedge = Vector{Array{Float64}}(undef, GV.num_layers)
    tfrontedge = Vector{Array{Float64}}(undef, GV.num_layers)
    for ialt in 1:GV.num_layers
        back_cols  = [bc_dict_horiz[sp][ialt][1, :] for sp in GV.transport_species]
        front_cols = [bc_dict_horiz[sp][ialt][2, :] for sp in GV.transport_species]
        tbackedge[ialt]  = permutedims(reduce(hcat, back_cols))
        tfrontedge[ialt] = permutedims(reduce(hcat, front_cols))
    end

    expected_edge_shape = (length(GV.transport_species), 2)
    @assert all(size(mat) == expected_edge_shape for mat in tbackedge) "horizontal back edge shape mismatch"
    @assert all(size(mat) == expected_edge_shape for mat in tfrontedge) "horizontal front edge shape mismatch"

    # if !GV.enable_horiz_transport && all(all(v .== 0.0) for v in GV.horiz_wind_v)
    #     max_for  = maximum(abs, tforwards)
    #     max_back = maximum(abs, tbackwards)
    #     max_backedge = maximum(abs, vcat([vec(m) for m in tbackedge]...))
    #     max_frontedge = maximum(abs, vcat([vec(m) for m in tfrontedge]...))
    #     println("DEBUG: horiz transport disabled -- max forward=", max_for,
    #             " max backward=", max_back,
    #             " back edge=", max_backedge,
    #             " front edge=", max_frontedge)
    # end    

    if debug
        println("[update_horiz_transport_coefficients] horizontal transport coefficients:")
        nprint = min(3, GV.num_layers)
        for ihoriz in 1:n_horiz
            println("  column ", ihoriz, " forwards[1:" , nprint, "] = ",
                    tforwards[ihoriz, 1:nprint, :])
            println("  column ", ihoriz, " backwards[1:" , nprint, "] = ",
                    tbackwards[ihoriz, 1:nprint, :])
        end
        if !GV.enable_horiz_transport
            max_for  = maximum(abs, tforwards)
            max_back = maximum(abs, tbackwards)
            println("  horizontal disabled -> max forward=", max_for,
                    " max backward=", max_back)
        end
    end

    return tbackedge, tforwards, tbackwards, tfrontedge
end

#===============================================================================#
#                    Water profile setup and SVP functions                      #   
#===============================================================================#

function precip_microns(sp, sp_profile; globvars...)
    #=
    Calculates precipitable microns of a species in the atmosphere.
    I guess you could use this for anything but it's only correct for H2O and HDO.

    Inputs:
        sp: Species name
        sp_mixing_ratio: mixing ratio profile of species sp
        atmdict: Present atmospheric state dictionary
    Outputs:
        Total precipitable micrometers of species sp
    =#
    GV = values(globvars)
    required =  [:molmass, :dz]
    check_requirements(keys(GV), required)

    col_abundance = column_density(sp_profile; GV.dz)
    cc_per_g = GV.molmass[sp] / GV.molmass[:H2O] # Water is 1 g/cm^3. Scale appropriately.

    #pr μm = (#/cm²) * (1 mol/molecules) * (g/1 mol) * (1 cm^3/g) * (10^4 μm/cm)
    pr_microns = col_abundance * (1/6.02e23) * (GV.molmass[sp]/1) * (cc_per_g / 1) * (1e4/1)
    return pr_microns
end

function colabund_from_prum(sp, prum; globvars...)
    #=
    Calculates precipitable microns of a species in the atmosphere.
    I guess you could use this for anything but it's only correct for H2O and HDO.

    Inputs:
        sp: Species name
        sp_mixing_ratio: mixing ratio profile of species sp
        atmdict: Present atmospheric state dictionary
    Outputs:
        Total precipitable micrometers of species sp
    =#
    GV = values(globvars)
    required =  [:molmass]
    check_requirements(keys(GV), required)

    cc_per_g = GV.molmass[sp] / GV.molmass[:H2O] # Water is 1 g/cm^3. Scale appropriately.

    #based on pr μm = (#/cm²) * (1 mol/molecules) * (g/1 mol) * (1 cm^3/g) * (10^4 μm/cm)
    col_abundance = prum * (6.02e23) * (1/GV.molmass[sp]) * (1/cc_per_g) * (1/1e4) 
    # NOTE: colabundance includes an implicit dz because when you calculate column abundance from prum you multiply by dz.
    return col_abundance
end

# 1st term is a conversion factor to convert to (#/cm^3) from Pa. Source: Marti & Mauersberger 1993
Psat(T) = (1e-6 ./ (kB_MKS .* T)) .* (10 .^ (-2663.5 ./ T .+ 12.537))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# However, this function is defined on the offchance someone studies HDO.
Psat_HDO(T) = (1e-6/(kB_MKS * T))*(10^(-2663.5/T + 12.537))

# function set_h2oinitfrac_bySVP(atmdict, h_alt; globvars...)
function set_h2oinitfrac_bySVP(atmdict, h_alt; ihoriz=1, globvars...)
    #=
    Calculates the initial fraction of H2O in the atmosphere, based on the supplied mixing ratio (global variable)
    and the saturation vapor pressure of water.
    Inputs:
        atmdict: atmospheric state
        h_alt: hygropause alt
    Output: 
        H2Oinitfrac: a mixing ratio by altitude vector.
    =#

    GV = values(globvars)
    required = [:all_species, :alt,  :num_layers, :n_alt_index,
               :H2Osat, :water_mixing_ratio]
    check_requirements(keys(GV), required)

    H2Osatfrac = GV.H2Osat ./ map(z->n_tot(atmdict, z, ihoriz; GV.all_species, GV.n_alt_index), GV.alt)  # get SVP as fraction of total atmo for this column
    # set H2O SVP fraction to minimum for all alts above first time min is reached
    H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
    H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
                   fill(minimum(H2Osatfrac), GV.num_layers-length(H2Oinitfrac))]

    # Set lower atmospheric water to be well-mixed (constant with altitude) below the hygropause
    H2Oinitfrac[findall(x->x<h_alt, GV.alt)] .= GV.water_mixing_ratio

    for i in [1:length(H2Oinitfrac);]
        H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
    end
    return H2Oinitfrac, H2Osatfrac
end

function setup_water_profile!(atmdict; constfrac=1, dust_storm_on=false, make_sat_curve=false, water_amt="standard", excess_water_in="mesosphere", 
                                       venus_special_water=false, venus_special_h2o_bot=nothing, venus_special_hdo_bot=nothing,
                                       venus_special_h2o_top=nothing, venus_special_hdo_top=nothing,
                                       showonly=false, hygropause_alt=40e5, globvars...)
    #=
    Sets up the water profile as a fraction of the initial atmosphere. 
    Input:
        atmdict: dictionary of atmospheric density profiles by altitude
        Optional:
            constfrac: Fraction of total density to use for the water profile. Currently only applied to Venus runs.
            dust_storm_on: whether to add an extra parcel of water at a certain altitude.
            tanh_prof: "low", "standard", or "high" to choose 1/10, mean, or 10x as much water in the atmosphere.
            hygropause_alt: altitude at which the water will switch from well-mixed to following the saturation vapor pressure curve.
    Output: 
        atmdict: Modified in place with the new water profile. 
    =#

    GV = values(globvars)
    check_requirements(keys(GV), [:all_species, :DH, :n_alt_index, :n_horiz, :planet, :plot_grid, :results_dir, :sim_folder_name, :enable_horiz_transport])

    # Set the initial fraction of the atmosphere for water to take up, plus the saturation fraction
    # ================================================================================================================
    # Currently this doesn't change behavior based on planet. 5/15/24
    # H2Oinitfrac, H2Osatfrac = set_h2oinitfrac_bySVP(atmdict, hygropause_alt; globvars...)

    H2Oinitfrac_all = Vector{Vector{ftype_ncur}}(undef, GV.n_horiz)
    H2Osatfrac_all = Vector{Vector{ftype_ncur}}(undef, GV.n_horiz)

    for ihoriz in 1:GV.n_horiz
        H2Oinitfrac_all[ihoriz], H2Osatfrac_all[ihoriz] =
            set_h2oinitfrac_bySVP(atmdict, hygropause_alt; ihoriz=ihoriz, globvars...)
    end

    if GV.planet=="Mars"
        required = [:alt, :H2Osat, :n_alt_index, :non_bdy_layers, :num_layers, :speciescolor, :speciesstyle, :upper_lower_bdy_i, :water_mixing_ratio,]
        check_requirements(keys(GV), required)

        # For doing high and low water cases 
        # ================================================================================================
        if (water_amt=="standard") | (excess_water_in=="loweratmo")
            println("Standard profile: water case = $(water_amt), loc = $(excess_water_in), MR = $(GV.water_mixing_ratio)")
        else # low or high in mesosphere and above - special code for paper 3
            println("$(water_amt) in $(excess_water_in)")

            toplim_dict = Dict("mesosphere"=>GV.upper_lower_bdy_i, "everywhere"=>GV.n_alt_index[GV.alt[end]])
            a = 1
            b = toplim_dict[excess_water_in]
            # H2Oinitfrac[a:b] = H2Oinitfrac[a:b] .* water_tanh_prof(GV.non_bdy_layers./1e5; z0=GV.ealt, f=GV.ffac)[a:b]
            for ihoriz in 1:GV.n_horiz
                prof = H2Oinitfrac_all[ihoriz]
                prof[a:b] .= prof[a:b] .* water_tanh_prof(GV.non_bdy_layers./1e5; z0=GV.ealt, f=GV.ffac)[a:b]
                H2Oinitfrac_all[ihoriz] = prof
            end

            # Set the upper atmo to be a constant mixing ratio, wherever the disturbance ends
            if excess_water_in=="everywhere"
                # H2Oinitfrac[GV.upper_lower_bdy_i:end] .= H2Oinitfrac[GV.upper_lower_bdy_i]
                for ihoriz in 1:GV.n_horiz
                    prof = H2Oinitfrac_all[ihoriz]
                    prof[GV.upper_lower_bdy_i:end] .= prof[GV.upper_lower_bdy_i]
                    H2Oinitfrac_all[ihoriz] = prof
                end
            end
        end

        # set the water profiles 
        # ===========================================================================================================
        for ihoriz in 1:GV.n_horiz
            atmdict[:H2O][ihoriz] = H2Oinitfrac_all[ihoriz] .* n_tot(atmdict, ihoriz; GV.n_alt_index, GV.all_species)
            atmdict[:HDO][ihoriz] = 2 * GV.DH * atmdict[:H2O][ihoriz]
        end
        HDOinitfrac_all = [atmdict[:HDO][ihoriz] ./ n_tot(atmdict, ihoriz; GV.n_alt_index, GV.all_species)
                            for ihoriz in 1:GV.n_horiz]

        # Add a gaussian parcel of water, to simulate the effect of a dust storm
        # ===========================================================================================================
        if dust_storm_on
            sigma = 12.5
            for ihoriz in 1:GV.n_horiz
                H2Oppm = 1e-6*map(z->GV.H2O_excess .* exp(-((z-GV.ealt)/sigma)^2), GV.non_bdy_layers/1e5) + H2Oinitfrac_all[ihoriz]
                HDOppm = 1e-6*map(z->GV.HDO_excess .* exp(-((z-GV.ealt)/sigma)^2), GV.non_bdy_layers/1e5) + HDOinitfrac_all[ihoriz]
                atmdict[:H2O][ihoriz][1:GV.upper_lower_bdy_i] = (H2Oppm .* n_tot(atmdict, ihoriz; GV.n_alt_index, GV.all_species))[1:GV.upper_lower_bdy_i]
                atmdict[:HDO][ihoriz][1:GV.upper_lower_bdy_i] = (HDOppm .* n_tot(atmdict, ihoriz; GV.all_species))[1:GV.upper_lower_bdy_i]
            end
        end
    elseif GV.planet=="Venus"
        ntot = n_tot(atmdict; GV.n_alt_index, GV.all_species)

        ntot_all = [n_tot(atmdict, ihoriz; GV.n_alt_index, GV.all_species) for ihoriz in 1:GV.n_horiz]

        for ihoriz in 1:GV.n_horiz
            ntot = ntot_all[ihoriz]
            atmdict[:H2O][ihoriz] = constfrac .* ntot
            atmdict[:HDO][ihoriz] = 2 * GV.DH * atmdict[:H2O][ihoriz]
        end

        # SPECIAL: Try a prescribed high water abundance in the mesosphere.
        if venus_special_water == true
            n = 11
            vmr_h2o = logrange(venus_special_h2o_bot, venus_special_h2o_top, n)
            vmr_hdo = logrange(venus_special_hdo_bot, venus_special_hdo_top, n)

            for ihoriz in 1:GV.n_horiz
                ntot = ntot_all[ihoriz]
                atmdict[:H2O][ihoriz][1:n] = vmr_h2o .* ntot[1:n]
                atmdict[:HDO][ihoriz][1:n] = vmr_hdo .* ntot[1:n]
            end
        end
    end

    # Plot the water profile 
    # ===========================================================================================================
    if make_sat_curve
        satarray = H2Osatfrac_all[1]
    else
        satarray = nothing 
    end

    plot_water_profile(atmdict, GV.results_dir*GV.sim_folder_name; watersat=satarray, H2Oinitf=H2Oinitfrac_all[1], plot_grid=GV.plot_grid, showonly=showonly, globvars...)

    # # ----------------------------------------------------------------------
    # # Consistency checks when horizontal transport is disabled
    # # ----------------------------------------------------------------------
    # if !GV.enable_horiz_transport && GV.n_horiz > 1
    #     base_H2O = atmdict[:H2O][1]
    #     base_HDO = atmdict[:HDO][1]

    #     for ihoriz in 2:GV.n_horiz
    #         @assert size(atmdict[:H2O][ihoriz]) == size(base_H2O) "H2O profile shape mismatch across columns"
    #         @assert size(atmdict[:HDO][ihoriz]) == size(base_HDO) "HDO profile shape mismatch across columns"
    #         @assert all(atmdict[:H2O][ihoriz] .== base_H2O) "H2O profiles differ across columns with horizontal transport disabled"
    #         @assert all(atmdict[:HDO][ihoriz] .== base_HDO) "HDO profiles differ across columns with horizontal transport disabled"
    #     end
    # end
end

function water_tanh_prof(z; f=10, z0=62, dz=11)
    #=
    Apply a tanh_prof to the water init fraction to add or subtract water from the atmosphere.
    =#
    return ((f .- 1)/2) * (tanh.((z .- z0) ./ dz) .+ 1) .+ 1
end

#===============================================================================#
#                       Photochemical equilibrium functions                     #   
#===============================================================================#

function choose_solutions(possible_solns, prev_densities) 
    #=
    Input:
        possible_solns: An array of possible densities for each species being solved for. 
                        If P-L=0 is linear in the species density, the solution is in 
                        column 1 of possible_solns. If quadratic, there are two solutions
                        in columns 1 and 2 and 3 possible cases:
                            - Both solutions are positive, and we use the one that minimizes the 
                              change in density from the previous state
                            - Both solutions are negative, in which case we set the density to 0
                            - There is a positive and a negative solution, in which case we use
                              the solution that minimizes the change in density from the previous
                              state and also sets the negative value to 0 if that's the chosen one.
        prev_densities: densities from the last timestep, for checking how much of a change 
                        has occurred.
    Output:
        accepted_solns: a single vector of the best possible choice of solution for each 
                        species based on these rules.
    =#

    # make an array copy to return 
    accepted_solns = deepcopy(possible_solns)

    # need to get the size of previous densities to be the same as what we will eventually subtract it from
    prev_densities_2d = hcat(prev_densities, prev_densities)

    # all indices in possible_solns that store quadratic solutions
    has_quadsoln = findall(possible_solns[:, 2] .!= 0)

    # find any quadratic equation where at least one solution is positive
    anypos = [r for r in has_quadsoln if any(possible_solns[r, :] .> 0)]

    # If there are any positive solutions, use the solution that minimizes the difference between the previous density and itself
    # and set the accepted solution to 0 if it's negative
    change = abs.(accepted_solns[anypos, :] .- prev_densities_2d[anypos, :])
    inds_of_min_change = reshape([x[2] for x in argmin(change, dims=2)], size(change)[1])
    accepted_solns[anypos, 1] .= [accepted_solns[i, j] for (i,j) in zip(anypos, inds_of_min_change)]
    accepted_solns[anypos, 2] .= NaN

    # zero out any negative solutions -- this includes quadratic solutions where both entries are negative.
    @views neg_solns = findall(accepted_solns[:, 1] .< 0)
    accepted_solns[neg_solns, 1] .= 0
    
    # The second column is now useless -- includes quadratic solutions where this second solution is negative.
    accepted_solns[:, 2] .= NaN 

    return accepted_solns[:, 1]
end

function construct_quadratic(sp::Symbol, prod_rxn_list, loss_rxn_list)
    #=
    Input:
        sp: Species with a quadratic reaction (where it reacts with itself)
        prod_rxn_list: List of reactions that produce species sp
        loss_rxn_list: List of reactions that consume species sp
    Output:
        quad_coef_dict: Expression like :((k1 + k2) * (n_sp^2) + (k3) * (n_sp))
    
    for this function, let ns represent the density term for sp, the species of interest.
    =#
    
    # Get rid of duplicates (normally present because a reaction with two instances of a species on the LHS
    # would be counted twice to account for the two species density terms). Only for loss reactions.
    # Continue letting production reactions appear twice since for those, sp is on the RHS. 
    # (easier to add a reaction in twice when figure out it's a duplicate and append "2 *")
    # Also because of the way filter! and insert! work, reaching all the way out to the global scope (??), 
    # deepcopy() has to be here to avoid adding another sp^0 term to whatever was passed in as prod_rxn_list and loss_rxn_list 
    # every time this function runs.
    loss_rxn_list = deepcopy(collect(Set(loss_rxn_list)))
    prod_rxn_list = deepcopy(prod_rxn_list) 
 
    # Reconstruct loss_rxn_list such that the first term is n_s to some power
    for r in 1:length(loss_rxn_list)
        ns_power = count(x->x==sp, loss_rxn_list[r]) # find the power of the density term

        # Replace all instances of n_s with a single term raised to ns_power
        filter!(x->x!=sp, loss_rxn_list[r]) 
        insert!(loss_rxn_list[r], 1, :($sp ^ $ns_power))
    end
 
    # Do the same for prod_rxn_list, but the power should always be 0
    for r in 1:length(prod_rxn_list)
        ns_power = count(x->x==sp, prod_rxn_list[r])

        # Replace all instances of n_s with a single term raised to ns_power
        filter!(x->x!=sp, prod_rxn_list[r]) 
        insert!(prod_rxn_list[r], 1, :($sp ^ $ns_power))
    end

    # make the list of vectors into a single expression
    quad_coef_dict = group_terms(prod_rxn_list, loss_rxn_list)
end

function group_terms(prod_rxn_arr, loss_rxn_arr)
    #=
    Input:
        prod_rxn_arr: Production reactions 
        loss_rxn_arr: Loss reactions
        Both in the form: Vector{Any}[[:(sp ^ 2), :k1], [:sp, :X, :k2]...], [[:Y, :k3]...]...]
    Output:
        terms_exprs: a dictionary of of the form Dict("A"=>:(r1*k1 + r2*k2...), "B"=>:(r3*k3 + r4*k4...))
                     for the terms A, B, C in the quadratic equation that defines a species density.
    =#
    terms_vecs = Dict("A"=>[], "B"=>[], "C"=>[])
    terms_exprs = Dict("A"=>:(), "B"=>:(), "C"=>:())
    qcoef_dict = Dict(2=>"A", 1=>"B", 0=>"C")
    
    for r in 1:length(loss_rxn_arr)

        # assign the quadratic coefficient.
        pow = (loss_rxn_arr[r][1]).args[3]
        quadratic_coef = qcoef_dict[pow]
        
        if quadratic_coef == "C"
            throw("Loss eqn has n_s^0")
        end
    
        # For each power term (i.e. (n_s)^2 or (n_s)^1), store all the coefficients
        # as a product of reactant and associated rate, i.e. :(r1*k1), :(r2*k2)
        # Note this looks a little ugly if end=2, but it works fine.
        push!(terms_vecs[quadratic_coef], :(*($(loss_rxn_arr[r][2:end]...))))
    end
    
    # Same loop, but over the production reactions. All of these should have quadratic_coef = "C"
    for r in 1:length(prod_rxn_arr)
        # assign the quadratic coefficient.
        pow = (prod_rxn_arr[r][1]).args[3]
        quadratic_coef = qcoef_dict[pow]
        
        if quadratic_coef != "C"
            throw("Production eqn has n_s power > 0")
        end
    
        push!(terms_vecs[quadratic_coef], :(*($(prod_rxn_arr[r][2:end]...))))
    end

    for k in keys(terms_vecs)
        # Collect the value vectors into sum expressions
        # again, this is ugly if length(terms[k]) = 1, but it works and keeps code simple.
        terms_exprs[k] = :(+($(terms_vecs[k]...)))
    end
    
    return terms_exprs
end

function loss_coef(leqn_vec, sp::Symbol; calc_tau_chem=false)
    #=
    Input:
        leqn: output of loss_equations (a vector of vectors of symbols)
        sp: Symbol; species for which to calculate the loss coefficient for use in 
            calculating photochemical equilibrium, n = P/L
        calc_tau_chem: if set to true, the first species density (n_s) term will be deleted
                       from all reactions, including photolysis reactions.
    
    Output:
        A vector of vectors of symbols like leqn, but with sp symbol removed from each
        sub-vector. This is basically taking (n_s*n_A*k1 + n_s*n_B*k2...) = n_s(n_A*k1 + n_B*k2...)

    FUNCTION ASSUMES THAT (P_s - L_s) = 0 FOR A GIVEN SPECIES IS LINEAR IN THE SPECIES DENSITY n_s. 
    NOT to be used if (P_s - L_s) = 0 is quadratic in n_s!!
    =#

    # Sets the number of terms a reaction must have before entries of the species density n_s
    # start getting removed. If 1, then all reactions will have n_s removed, including 
    # photolysis reactions--this is useful to calculate the chemical lifetime.
    # If 2, photolysis reactions will keep their n_s, since it's the only reactant.
    if calc_tau_chem
        min_terms = 1
    else
        min_terms = 2
    end

    for L in 1:length(leqn_vec)

        # conditions for readability
        modifiable_rxn = length(leqn_vec[L]) > min_terms ? true : false
        num_density_terms_present = count(t->t==sp, leqn_vec[L])
        
        if modifiable_rxn && num_density_terms_present == 1 
            leqn_vec[L] = deletefirst(leqn_vec[L], sp)
        else
            if num_density_terms_present > 1
                throw("Error: $(num_density_terms_present) density terms found in the reaction for $(sp). Not linear!") 
            elseif num_density_terms_present == 0
                continue
            end
        end
    end

    return leqn_vec
end

function linear_in_species_density(sp, lossnet)
    #=
    Input:
        sp: species name
        lossnet: Collection of reactions that consume sp
    Output: 
        true: all equations in lossnet are linear in sp density (only one 
              instance of sp density on the LHS)
        false: at least one equation has >=2 instances of sp density on the LHS
    =#
    for r in 1:length(lossnet)
        d = count(x->x==sp, lossnet[r])
        if d > 1
            return false
        end
    end
    return true 
end

# function setup_photochemical_equilibrium(; globvars...)
function setup_photochemical_equilibrium(; transportnet_horiz=[], globvars...)

    #=
    Output
        active_longlived_species_rates: rate expressions for species which are actively solved for.
        short_lived_density_eqn: array containing expressions for density of short-lived species, i.e. n = P/L, 
                                 also includes expression when a species is quadratic (like OH + OH)
        shortlived_density_inputs: Array containing P, L, and the determinant for the short lived density
                                   expressions. Needed because if P and L are both < machine epsilon, problems occur
        equilibrium_eqn_terms: Another storage array, see below for explanation I'm tired
    =#

    GV = values(globvars)
    required = [:active_longlived, :active_shortlived, :short_lived_species, :reaction_network, :transportnet, :chem_species, :transport_species]
    check_requirements(keys(GV), required)


    # ------------------ Long-lived species expression array ------------------------ #

    # An array to store the rate equations for active, long-lived species, which are 
    # solved for in the production and loss equation.
    # each row is for each species; each column is for chemical production, chemical loss, 
    # transport production, transport loss, in that order.
    active_longlived_species_rates = Array{Array{Expr}}(undef, length(GV.active_longlived), 4)
    for (i, sp) in enumerate(GV.active_longlived)
        # active_longlived_species_rates[i, :] .= getrate(sp, sepvecs=true; chemnet=GV.reaction_network, GV.transportnet, GV.chem_species, GV.transport_species, )
        active_longlived_species_rates[i, :] .= getrate(sp, sepvecs=true,
                                                       transportnet_horiz=transportnet_horiz;
                                                       chemnet=GV.reaction_network,
                                                       GV.transportnet,
                                                       GV.chem_species,
                                                       GV.transport_species, )

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
    short_lived_density_eqn = Array{Expr}(undef, length(GV.short_lived_species), 2) 

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
    shortlived_density_inputs = Array{Expr}(undef, length(GV.short_lived_species), 2)


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
    equilibrium_eqn_terms = Array{Expr}(undef, length(GV.short_lived_species), 1)


    for (i, sp) in enumerate(GV.active_shortlived)
        # Removes from consideration any reaction where a species appears on both sides of the equation (as an observer)
        ret = rxns_where_species_is_observer(sp, GV.reaction_network)
        if ret == nothing
            chemnet = GV.reaction_network
        else
            chemnet = filter(x->!in(x, ret), GV.reaction_network)        
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
