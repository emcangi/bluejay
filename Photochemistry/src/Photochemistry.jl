__precompile__()

module Photochemistry

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils
using GeneralizedGenerated
using DataFrames
using XLSX


export # Basic utility functions
       charge_type, create_folder, deletefirst, find_nonfinites, fluxsymbol, format_chemistry_string, format_sec_or_min, getpos, input, get_paramfile, 
       logrange, nans_present, next_in_loop, searchdir, searchsortednearest, search_subfolders, subtract_difflength, write_to_log,
       # Atmospheric basic functions
       atm_dict_to_matrix, atm_matrix_to_dict, column_density, column_density_above, find_exobase, flatten_atm, get_ncurrent, ncur_with_boundary_layers, n_tot, 
       precip_microns, setup_water_profile!, unflatten_atm, write_atmosphere,
       # Plotting functions
       get_colors, get_grad_colors, plot_atm, plot_bg, plot_extinction, plot_Jrates, plot_rxns, plot_temp_prof, plot_water_profile,                 
       # Reaction rate functions
       eval_rate_coef, get_column_rates, get_volume_rates, reactant_density_product,                                                    
       # Boundary condition functions                                                   
       boundaryconditions, effusion_velocity,
       # transport functions                                                                           
       Dcoef!, fluxcoefs, flux_param_arrays, flux_pos_and_neg, get_flux, Keddy, scaleH, update_diffusion_and_scaleH, update_transport_coefficients,                      
       # Chemistry functions
       load_reaction_network, calculate_stiffness, check_jacobian_eigenvalues, chemical_jacobian, getrate, loss_equations, 
       loss_rate, make_chemjac_key, make_net_change_expr, meanmass, production_equations, production_rate, rxns_where_species_is_observer, 
       # Photochemical equilibrium functions
       choose_solutions, construct_quadratic, group_terms, loss_coef, linear_in_species_density,
       # Photochemistry functions
       binupO2, co2xsect, h2o2xsect_l, h2o2xsect, hdo2xsect, ho2xsect_l, o2xsect, O3O1Dquantumyield, padtosolar, populate_xsect_dict, quantumyield, 
       # Temperature functions
       T_updated, T_all, Tpiecewise,                                                                                                      
       # Water profile functions  
       Psat, Psat_HDO         

code_loc = "$(@__DIR__)/../../"
println("Photochemistry.jl code_loc = $code_loc")
include(code_loc*"CUSTOMIZATIONS.jl") 
include(code_loc*"CONSTANTS.jl")
                                                                                

# **************************************************************************** #
#                                                                              #
#                          Basic utility functions                             #
#                                                                              #
# **************************************************************************** #

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

function create_folder(foldername::String, parentdir::String)
    #=
    Checks to see if foldername exists within parentdir. If not, creates it.
    =#
    println("Checking for existence of $(foldername) folder in $(parentdir)")
    dircontents = readdir(parentdir)
    if foldername in dircontents
        println("Folder $(foldername) exists")
    else
        mkdir(parentdir*foldername)
        println("Created folder ", foldername)
    end
end

function delete_old_h5file(filename::String)
    #= 
    the HDF5 module doesn't allow implicit deletion of existing names in the function
    h5write, so this function will delete previously xisting files with the name 'filename',
    i.e. in cases where a run did not complete successfully and left behind garbage. 
    =#
    current_folder = pwd()
    dircontents = readdir(current_folder)
    if filename in dircontents
        println("File $(filename) exists, deleting")
        rm(filename)
    end
end

function deletefirst(A, v)
    #=
    returns: list A with its first element equal to v removed. Used to make the derivative for 
    the chemical jacboian 
    =#
    index = something(findfirst(isequal(v), A), 0)  # this horrible syntax introduced by Julia devs
    keep = setdiff([1:length(A);],index)
    return A[keep]
end

function find_nonfinites(collection; collec_name="collection")
    #=
    Returns indices of any nonfinite values (inf or nan) in collection.
    =#
    nonfinites = findall(x->x==0, map(el->isfinite(el), collection))
    if length(nonfinites) != 0
        throw("ALERT: Found nonfinite values in $(collec_name) at indices $(nonfinites)")
    end
end

function fluxsymbol(x)
    #= 
    Converts string x to a symbol. f in the string means flux. 
    =#
    return Symbol(string("f",string(x)))
end

function format_chemistry_string(reactants::Array{Symbol}, products::Array{Symbol})
    #=
    Given arrays of the reactants and products in a reaction,
    this small function constructs a chemistry reaction string.
    =#
    return string(join(reactants, " + ")) * " --> " * string(join(products, " + "))
end

function format_sec_or_min(t) 
    #=
    Takes a time value t in seconds and formats it to read as H,M,S.
    =#

    thehour = div(t, 3600.)
    themin = div(t%3600., 60.)
    thesec = round(t%3600. - (themin*60.), digits=1)

    return "$(thehour) hours, $(themin) minutes, $(thesec) seconds"
end

function get_Jrate_symb(molecs::Array, Jrates::Array)
    #=
    Input:
        molecs: Array of reactants and products as strings, order doesn't matter, in some photodissociation/
                photoionization reaction. 
        Jrate_symbs: Array of Jrate symbols used in the model. 
    Output
        J: the Jrate associated with the list of reactants and products. 
    =#
    for J in Jrates
        if sort([m.match for m in collect(eachmatch(r"[A-I0-9K-Z]+(pl)*", string(J)))]) == sort(molecs)
            return J
            break
        end      
    end
end

function get_paramfile(working_dir)
    #=
    Helper function to let the user select a parameter file form working_dir. 
    =#
    available_paramfiles = filter(x->(occursin(".jl", x) & occursin("PARAMETERS", x)), readdir(working_dir))

    println("Available parameter files: ")
    [println("[$(i)] $(available_paramfiles[i])") for i in 1:length(available_paramfiles)]

    if length(available_paramfiles) == 1
        paramfile = available_paramfiles[1]
        println("Using the only file available, $(paramfile)")
    else 
        paramfile_selection = input("Select a parameter file or enter filename with extension: ")
        paramfile = available_paramfiles[parse(Int64, paramfile_selection)]
        println("Using parameter file $(paramfile).")
    end

    return paramfile
end

function getpos(array, test::Function, n=Any[])
    #= 
    Get the position of keys that match "test" in an array

    array: Array; any size and dimension
    test: used to test the elements in the array
    n: Array; for storing the indices (I think)

    returns: 1D array of elements in array that match test function supplied
    =#
    if !isa(array, Array)
        test(array) ? Any[n] : Any[]
    else
        vcat([ getpos(array[i], test, Any[n...,i]) for i=1:size(array)[1] ]...)
    end
end

function getpos(array, value)
    #= 
    overloading getpos for the most common use case, finding indicies
    corresponding to elements in array that match value. 
    =#
    getpos(array, x->x==value)
end

function input(prompt::String="")::String
    #=
    Prints a prompt to the terminal and accepts user input, returning it as a string.
    =#
    print(prompt)
    return chomp(readline())
end

function logrange(x1, x2, n::Int64)
    #=
    Equivalent to python's logspace, returns a generator with entries that
    are powers of 10.
    x1: start, some power of 10
    x2: end, some power of 10
    n: number of entries, spaced between the powers of x1 and x2

    e.g.: collect(logrange(1e-4, 1e2, 4)) returns 0.0001, 0.01, 1.0, 100.0

    =#
    return (10.0^y for y in range(log10(x1), log10(x2), length=n))
end 

nans_present(a) = any(x->isnan(x), a)

function next_in_loop(i::Int64, n::Int64)
    #=
    returns i+1, restricted to values within [1, n]. 
    next_in_loop(n, n) returns 1.
    =#
    return i % n + 1
end

# searches path for key
searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function searchsortednearest(a,x)
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
        return idx
    else
        return idx-1
    end
end

function search_subfolders(path::String, key; type="folders")
    #=
    path: a folder containing subfolders and files.
    key: the text pattern in folder names that you wish to find.
    type: whether to look for "folders" or "files". 

    Searches the top level subfolders within path for all folders matching a 
    certain regex given by key. Does not search files or sub-subfolders.
    =#
    folderlist = []
    filelist = []
    if type=="folders"
        for (root, dirs, files) in walkdir(path)
            if root==path
                for dir in dirs
                    push!(folderlist, joinpath(root, dir)) # path to directories
                end
            end
        end
    elseif type=="files"
        for (root, dirs, files) in walkdir(path)
            if root==path
                for fil in files
                    push!(filelist, fil) # files
                end
            end
        end
    end

    if type=="folders"
        folderlist = filter(x->occursin(key, x), folderlist)
        return folderlist
    elseif type=="files"
        filelist = filter(x->occursin(key, x), filelist)
        return filelist
    end
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

function write_to_log(full_path::String, entries; mode="a")
    #=
    full_path: full path to log file (folder structure and filename, with extension)
    entries: List of strings to write to the log file.
    mode: w or a for write or append.
    =#

    f = open(full_path, mode) 
    if isa(entries, Array)
        for e in entries 
            write(f, string(e)*"\n")
        end
    elseif isa(entries, String)
        write(f, entries*"\n")
    else
        throw("Wrong format for logging: $(typeof(entries))")
    end

    close(f)
end


# **************************************************************************** #
#                                                                              #
#                        Atmosphere - basic functions                          #
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

function find_exobase(s::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; returntype="index", globvars...)
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
    @assert all(x->x in keys(GV),  [:nbl, :all_species, :Tn, :molmass, :alt, :xsects])

    H_s = scaleH(GV.nbl, s, GV.Tn; globvars...)# GV.molmass)
    try
        mfp_s = 1 ./ (GV.xsects[s] .* n_tot(atmdict; GV.all_species, GV.n_alt_index))
    catch KeyError
        throw("You need to implement handling for collisional cross section of species heavier than D")
    end

    for i in 50:length(GV.nbl) # extremely unlikely the exobase will ever be at 100 km (i=50) in alt lol
        if mfp_s[i] < H_s[i]
            continue
        else
            # println("Exobase found for species $(s) at altitude $(alts[i]/1e5) km")
            returnme = Dict("altitude"=>GV.alt[i], "index"=>n_alt_index[GV.alt[i]])
            return returnme[returntype]
            break
        end

    end
    throw("Didn't find an exobase so the whole atmosphere is collisional I guess, this should probably be addressed")
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

function get_ncurrent(readfile::String)
    #=
    Input:
        readfile: HDF5 file containing a matrix with atmospheric density profiles
    Output:
        n_current: dictionary of atmospheric density profiles by species 
    =#

    # This accounts for the old format of some files. 
    try
        global n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
        global n_current_mat = h5read(readfile,"n_current/n_current_mat");
    catch
        global n_current_tag_list = map(Symbol, h5read(readfile,"ncur/species"))
        global n_current_mat = h5read(readfile,"ncur/ncur_mat");
    end
    
    n_current = Dict{Symbol, Array{ftype_ncur, 1}}()

    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies], length(n_current_mat[:, ispecies]))
    end
    return n_current
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

function precip_microns(sp, sp_profile; globvars...)
    #=
    Calculates precipitable microns of a species in the atmosphere.
    I guess you could use this for anything but it's probably only correct for H2O and HDO.

    Inputs:
        sp: Species name
        sp_mixing_ratio: mixing ratio profile of species sp
        atmdict: Present atmospheric state dictionary
    Outputs:
        Total precipitable micrometers of species sp
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:molmass])

    col_abundance = column_density(sp_profile)
    cc_per_g = GV.molmass[sp] / GV.molmass[:H2O] # Water is 1 g/cm^3. Scale appropriately.

    #pr μm = (#/cm²) * (1 mol/molecules) * (g/1 mol) * (1 cm^3/g) * (10^4 μm/cm)
    pr_microns = col_abundance * (1/6.02e23) * (GV.molmass[sp]/1) * (cc_per_g / 1) * (1e4/1)
    return pr_microns
end

function setup_water_profile!(atmdict; dust_storm_on=false, globvars...)
    #=
    Sets up the water profile as a fraction of the initial atmosphere. 
    Input:
        atmdict: dictionary of atmospheric density profiles by altitude
    Output: 
        atmdict: Modified in place with the new water profile. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :num_layers, :DH, :alt, :plot_grid, :n_alt_index,
                                   :hygropause_alt, :H2O_excess, :HDO_excess, :ealt, :non_bdy_layers, :H2Osat, :water_mixing_ratio,
                                   :results_dir, :sim_folder_name])

    # H2O Water Profile ================================================================================================================
    H2Osatfrac = GV.H2Osat ./ map(z->n_tot(atmdict, z; GV.all_species, GV.n_alt_index), GV.alt)  # get SVP as fraction of total atmo
    # set H2O SVP fraction to minimum for all alts above first time min is reached
    H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
    H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
                   fill(minimum(H2Osatfrac), GV.num_layers-length(H2Oinitfrac))]

    # Set lower atmospheric water to be well-mixed (constant with altitude) below the hygropause
    H2Oinitfrac[findall(x->x<GV.hygropause_alt, GV.alt)] .= GV.water_mixing_ratio

    for i in [1:length(H2Oinitfrac);]
        H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
    end

    # set the water profiles ===========================================================================================================
    atmdict[:H2O] = H2Oinitfrac.*n_tot(atmdict; GV.n_alt_index, GV.all_species)
    atmdict[:HDO] = 2 * GV.DH * atmdict[:H2O] # This should be the correct way to set the HDO profile.
    HDOinitfrac = atmdict[:HDO] ./ n_tot(atmdict; GV.n_alt_index, GV.all_species)  # Needed to make water plots.

    # Plot the water profile ===========================================================================================================
    plot_water_profile(H2Oinitfrac, HDOinitfrac, atmdict[:H2O], atmdict[:HDO], GV.results_dir*GV.sim_folder_name, watersat=H2Osatfrac, plot_grid=GV.plot_grid)

    # ADD EXCESS WATER AS FOR DUST STORMS.
    if dust_storm_on
        H2Oppm = 1e-6*map(x->GV.H2O_excess .* exp(-((x-GV.ealt)/12.5)^2), GV.non_bdy_layers/1e5) + H2Oinitfrac
        HDOppm = 1e-6*map(x->GV.HDO_excess .* exp(-((x-GV.ealt)/12.5)^2), GV.non_bdy_layers/1e5) + HDOinitfrac  # 350 ppb at 38 km (peak)
        atmdict[:H2O] = H2Oppm .* n_tot(atmdict; GV.n_alt_index, GV.all_species)
        atmdict[:HDO] = HDOppm .* n_tot(atmdict; GV.all_species)
    end
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

function write_atmosphere(atmdict::Dict{Symbol, Vector{ftype_ncur}}, filename::String; globvars...) 
    #=
    Writes out the current atmospheric state to an .h5 file

    Input: 
        atmdict: atmospheric density profile dictionary
        filename: filename to write to
    Output: none
    =# 
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:alt, :num_layers])

    sorted_keys = sort(collect(keys(atmdict)))
    atm_mat = Array{Float64}(undef, GV.num_layers, length(sorted_keys));

    for ispecies in [1:length(sorted_keys);]
        for ialt in [1:GV.num_layers;]
            atm_mat[ialt, ispecies] = convert(Float64, atmdict[sorted_keys[ispecies]][ialt])
        end
    end
    delete_old_h5file(filename)
    h5open(filename, "w") do f # this syntax is ok because we never write multiple times to a file.
        write(f, "n_current/n_current_mat", atm_mat)
        write(f, "n_current/alt", GV.alt)
        write(f, "n_current/species", map(string, sorted_keys))
    end
end

# **************************************************************************** #
#                                                                              #
#                             Plotting Functions                               #
#                                                                              #
# **************************************************************************** #

function get_colors(L::Int64, cmap; strt=0, stp=1)
    #=
    Generates some colors based on a non-gradient color map for use in plotting a 
    bunch of lines all at once.
    Input:
        L: number of colors to generate.
        cmap: color map name
    Output:
        An array of RGB values: [[R1, G1, B1], [R2, G2, B2]...]
    =#

    Cmap = get_cmap(Symbol(cmap))
    colors = [Cmap(x) for x in range(strt, stop=stp, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors))
        c[i, 1] = colors[i][1]
        c[i, 2] = colors[i][2]
        c[i, 3] = colors[i][3]
    end
    return c
end

function get_grad_colors(L::Int64, cmap; strt=0, stp=1)
    #=
    Generates some colors based on a GRADIENT color map for use in plotting a 
    bunch of lines all at once.

    Input:
        L: number of colors to generate.
        cmap: color map name
        strt and stp: By setting these to other values between 0 and 1 you can restrict 
                      the range of colors drawn from.
    Output:
        An array of RGB values: [[R1, G1, B1], [R2, G2, B2]...]

    AVAILABLE MAPS: blues, viridis, pu_or, magma, plasma, inferno
    =#

    colors = [cgrad(Symbol(cmap))[x] for x in range(strt, stop=stp, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors))
        c[i, 1] = red(colors[i])
        c[i, 2] = green(colors[i])
        c[i, 3] = blue(colors[i])
    end
    return c
end

function plot_atm(atmdict::Dict{Symbol, Vector{ftype_ncur}}, savepath::String, atol; t="", showonly=false, xlab="", xlim_1=(1e-12, 1e18), xlim_2=(1e-5, 1e5), globvars...)
    #=
    Makes a "spaghetti plot" of the species concentrations by altitude in the
    atmosphere. 

    atmdict: dictionary of vertical profiles of a plottable quantity, keys are species names.
    species_lists: a list of species lists. The first is usually neutrals and goes in column 1; the second, ions, column 2.
    savepath: path and name for saving resulting .png file
    t: title for the overall plot
    showonly: whether to just show() the plot instead of saving. If setting to true, send in junk string for savepath.
    xlab: override for x axis label
    xlim_1: override for x axis limits in column 1 or entire plot of length(species_lists)==1
    xlim_2: override for x axis limits in column 2
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:neutral_species, :plot_grid, :speciescolor, :speciesstyle, :zmax])

    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    
    # Establish logical groups for axes =======================================================

    # This dictionary was created by hand. Species may be missing. If so, an error will automatically be thrown by the code
    species_groups = Dict( # PRIMARY NEUTRALS + IONS
                            1=>[:CO2,:CO2pl,
                                :CO,:COpl,
                                :N2,:Nup2D,:N2pl,
                                :Ar,:Arpl,
                                :ArHpl,:ArDpl,
                                :HCO,:HCOpl,:DCO,:DCOpl,
                                :HOCpl,:DOCpl,
                                :HOCO,:DOCO, :HCO2pl,:DCO2pl,
                                :O1D,
                                :O,:Opl,
                                :O2,:O2pl,
                                :O3],
                            # ODD HYDROGEN + IONS
                            2=>[:H,:Hpl,:D,:Dpl,
                                :H2,:H2pl, :HD, :HDpl,
                                :H3pl,:H2Dpl, :HD2pl, 
                                :H2O,:H2Opl,:HDO,:HDOpl,
                                :H3Opl,:H2DOpl,
                                :H2O2,:HDO2,
                                :HO2,:HO2pl,:DO2,
                                :OH,:OHpl,:OD,:ODpl],
                            # NITROGEN NEUTRALS + IONS and RARE C MOLECULES
                            3=>[:C,:Cpl,:CH,:CHpl,
                                :CN,:CNpl,:HCN,:HCNpl,:HCNHpl,
                                :HNO,:HNOpl,:HN2Opl,
                                :N,:Npl,:N2O,:N2Opl,
                                :NH,:NHpl,:NH2,:NH2pl,:NH3pl,
                                :NO,:NOpl,:NO2,:NO2pl,
                                :N2Hpl,:N2Dpl]
                        );

    axes_by_sp = Dict()

    for k in keys(species_groups)
        for sp in species_groups[k]
            axes_by_sp[sp] = k
        end
    end

    # Plot neutrals and ions together =========================================================
    if haskey(GV, :ion_species)#length(species_lists)==2  # neutrals and ions 
        
        # set up the overall plot -------------------------------------------------------------
        atm_fig, atm_ax = subplots(3, 2, sharex=false, sharey=true, figsize=(14, 16))
        subplots_adjust(wspace=0, hspace=0)
        tight_layout()
                
        # only the neutral-col axes
        atm_ax[1, 1].set_title("Neutrals")
        for i in 1:3
            plot_bg(atm_ax[i, 1])
            atm_ax[i, 1].set_xlim(xlim_1[1], xlim_1[2])
            atm_ax[i, 1].fill_betweenx(GV.plot_grid,#=plotalts,=# xlim_1[1] .* ones(size(GV.plot_grid,)), x2=atol, alpha=0.1, color=medgray, zorder=10)
            atm_ax[i, 1].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end

        x_axis_label = xlab == "" ? L"Species concentration (cm$^{-3}$)" : xlab
        atm_ax[3, 1].set_xlabel(x_axis_label)
        
        # only the ion-col axes
        atm_ax[1, 2].set_title("Ions")
        for i in 1:3
            plot_bg(atm_ax[i, 2])
            atm_ax[i, 2].set_xlim(xlim_2[1], xlim_2[2])
            atm_ax[i, 2].fill_betweenx(GV.plot_grid, xlim_2[1] .* ones(size(GV.plot_grid#=plotalts=#)), x2=atol, alpha=0.1, color=medgray, zorder=10)
            atm_ax[i, 2].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end
        atm_ax[3, 2].set_xlabel(x_axis_label)
        
        # y axes labels
        for j in 1:3
            atm_ax[j, 1].set_ylabel("Altitude [km]")
        end
        
        # plot the neutrals according to logical groups -------------------------------------------------------
        for sp in GV.neutral_species #species_lists[1]
            atm_ax[axes_by_sp[sp], 1].plot(convert(Array{Float64}, atmdict[sp]), GV.plot_grid#=plotalts=#, color=get(GV.speciescolor#=scolor=#, sp, "black"),
                                           linewidth=2, label=sp, linestyle=get(GV.speciesstyle#=sstyle=#, sp, "-"), zorder=10)
        end
        
        # plot the ions according to logical groups ------------------------------------------------------------
        for sp in GV.ion_species #species_lists[2]
            atm_ax[axes_by_sp[sp], 2].plot(convert(Array{Float64}, atmdict[sp]), GV.plot_grid#=plotalts=#, color=get(GV.speciescolor#=scolor=#, sp, "black"),
                                           linewidth=2, label=sp, linestyle=get(GV.speciesstyle#=sstyle=#, sp, "-"), zorder=10)
        end

        # stuff that applies to all axes
        for a in atm_ax
            a.set_ylim(0, GV.zmax/1e5)
            a.set_xscale("log")
            handles, labels = a.get_legend_handles_labels()
            if isempty(handles) == false
                a.legend(handles, labels, fontsize=8)#bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
            end
        end

    # Plot only neutrals - to support the fractionation factor project ==========================================
    else # ion species is not defined 
        atm_fig, atm_ax = subplots(figsize=(16,6))
        tight_layout()
        for sp in GV.neutral_species #species_lists[1]
            atm_ax.plot(convert(Array{Float64}, atmdict[sp]), GV.plot_grid#=plotalts=#, color=get(GV.speciescolor, sp, "black"),
                        linewidth=2, label=sp, linestyle=get(GV.speciesstyle, sp, "-"), zorder=1)
            atm_ax.set_xlim(xlim_1[1], xlim_1[2])
            atm_ax.set_ylabel("Altitude [km]")
            # atm_ax.set_title("Neutrals")
        end
        atm_ax.tick_params(which="both", labeltop=true, top=true)
        plot_bg(atm_ax)
        atm_ax.set_ylim(0, GV.zmax/1e5)
        atm_ax.set_xscale("log")
        atm_ax.set_xlabel(x_axis_label)
        atm_ax.legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0, fontsize=16)
    end

    suptitle(t, y=1.05)

    if showonly==false  
        atm_fig.savefig(savepath, bbox_inches="tight", dpi=300)
        # show()
    else
        show()
    end
end

function plot_bg(axob; bg="#ededed")
    #=
    Make plots not look ugly. 

    intended to immitate seaborn
    =#
    axob.set_facecolor(bg)
    axob.grid(zorder=-5, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function plot_extinction(solabs, path::String, dt, iter::Int64; tauonly=false, xsect_info=nothing, solflux=nothing, globvars...)
    #=
    Used to plot the extinction in the atmosphere by wavelength and altitude

    solabs: a ROW VECTOR (sadly) of solar absorption in the atmosphere. rows = altitudes, "columns" = wavelengths,
            but being a row vector it supposedly has only one column, but each element is a row.
    path: a path at which to save the plot
    dt: timestep, used for making a unique image file name
    iter: iteration number, same purpose as dt
    tauonly: plot just tau rather than e^-tau
    xsect_info: a list containing two things:
                xsect: list of crosssections for some Jrate, shape (124, 2000)
                xsect_sp: string representation of a chemical species which is absorbing using xsect
    solflux: solar flux values of shape (1, 2000), optional
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:zmax])

    
    fig, ax = subplots()

    num_x = length(solabs[1])  # extinction is in an awkward format
    num_y = size(solabs)[1]
    X = 1:num_x  # start with 0 because we are sending it into a python plotting library, not julia 
    Y = 0:2:(GV.zmax/1e5 - 2)

    # Don't mess with this line. It is witchcraft that translates the row vector that is extinction
    # into an actual 2D matrix so that we can do normal math on it. 
    solabs = transpose(reshape(collect(Iterators.flatten(solabs[:,:])), (num_x,num_y)))

    # now convert it to actual extinction:
    if tauonly
        solabs[solabs .< 0] .= 0
        # z = log10.(solabs)
        z = solabs
        titlestr = L"\tau"
        savestr = "tau_$(dt)_$(iter).png"
        z_min = 0
        z_max = 5
        heatmap = ax.pcolor(X, Y, z, cmap="bone", vmin=z_min, vmax=z_max)
        cbarlbl = L"\tau"
    else  # this is if we want to plot the extinction.
        if xsect_info==nothing
            throw("Error! Do you want to plot tau or the extinction? Please pass either tauonly=true or xsect_info")
        end
        xsect, xsect_sp = xsect_info
        z = exp.(-solabs)
        # println("printing maximums: exp, xsect, solflux")
        # println(maximum(z))
        
        titlestr = L"e^{-\tau}\sigma_{\lambda}"
        if xsect != nothing  # multiply by the crosssection if we want to plot the Jrates
            xsect = transpose(reshape(collect(Iterators.flatten(xsect[:,:])), (num_x,num_y)))
            z = z .* xsect
            titlestr = L"e^{-\tau}\sigma_{\lambda}" * " for $(xsect_sp)"
            savestr = "extinction_$(dt)_$(iter)_$(xsect_sp).png"
        end

        if solflux != nothing
            solflux = reshape(solflux, (1, length(solflux)))
            z = solflux .* z # solflux is an an array (actually a vector, tbh) of shape (1, 2000);
                             # exp.(-solabs) is an array of shape (124, 2000). This multiplies the solar flux
                             # values by the extinction across the 2000 wavelengths. Unbelievably, this 
                             # operation works in Julia in either direction, whether you do solflux .* z or reverse.
            titlestr = L"J_{\lambda} e^{-\tau}\sigma_{\lambda}" * " for $(xsect_sp)"
        end
        heatmap = ax.pcolor(X, Y, z, cmap="bone")
        cbarlbl = "photons/sec"
    end

    cbar = fig.colorbar(heatmap, ax=ax)
    cbar.set_label(cbarlbl, rotation=270, fontsize=14, labelpad=15)
    xlabel("Wavelength (nm)")
    ylabel("Altitude (km)")
    yticks(0:50:(GV.zmax/1e5))

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(25))

    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(25))

    title(titlestr)
    ax.set_xlim(0, 300)  # turn off to see full range to 2000 nm.
    savefig(path*savestr, bbox_inches="tight")
    close(fig)
end

function plot_Jrates(sp, atmdict::Dict{Symbol, Vector{ftype_ncur}}, results_dir::String; 
                     filenameext="", source2=nothing, globvars...)                
    #=
    Plots the Jrates for each photodissociation or photoionizaiton reaction. Override for small groups of species.

    species: species for which to plot the Jrates

    TODO: Needs to be worked on especially get_volume_rates works as expected
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:rxnnet, :plot_grid, :Tn, :Ti, :Te, :all_species, :ion_species, :bcdict, :num_layers, :plot_grid ])

    # Plot setup
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

    # --------------------------------------------------------------------------------
    # make plot
    
    rxd_prod, prod_rc = get_volume_rates(sp, atmdict; which="Jrates", globvars...) # Tn, Ti, Te, GV.all_species, GV.ion_species, GV.rxnnet, numlyrs, which="Jrates")
    rxd_loss, loss_rc = get_volume_rates(sp, atmdict; which="Jrates", globvars...) # Tn, Ti, Te, GV.all_species, GV.ion_species, GV.rxnnet, numlyrs, which="Jrates")

    if source2 != nothing
        rxd_prod2, prod_rc2 = get_volume_rates(sp, atmdict; which="Jrates", globvars...) # Tn, Ti, Te, GV.all_species, GV.ion_species, source2, numlyrs, which="Jrates")
        rxd_loss2, loss_rc2 = get_volume_rates(sp, atmdict; which="Jrates", globvars...) # Tn, Ti, Te, GV.all_species, GV.ion_species, source2, numlyrs, which="Jrates")
    end 

    fig, ax = subplots(figsize=(8,6))
    plot_bg(ax)

    minx = 1e10
    maxx = 0
    
    # Collect chem production equations and total 
    if isempty(keys(rxd_prod))
        println("$(sp) has no photodissociation or photoionization reactions, no plot generated")
    else 
        for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
            lbl = "$(kv[1])"

            if source2 != nothing 
                if !all(x->x<=1e-10, abs.(kv[2] - rxd_prod2[kv[1]]))
                    ax.semilogx(kv[2] - rxd_prod2[kv[1]], GV.plot_grid, linestyle="-", linewidth=1, label=lbl)
                else
                    text(0.5, 0.5, "Everything is basically 0, nothing to plot", transform=ax.transAxes)
                end
            else
                ax.semilogx(kv[2], GV.plot_grid, linestyle="-", linewidth=1, label=lbl)
            end
            
            # set the xlimits
            if minimum(kv[2]) <= minx
                minx = minimum(kv[2])
            end
            if maximum(kv[2]) >= maxx
                maxx = maximum(kv[2])
            end
        end

        if source2 != nothing
            suptitle("J rate difference between spreadsheet and jl file for $(sp), $(filenameext)", fontsize=20)
        else
            suptitle("J rates for $(sp), $(filenameext)", fontsize=20)
        end
        ax.set_ylabel("Altitude (km)")
        ax.legend(bbox_to_anchor=(1.01, 1))
        ax.set_xlabel(L"Reaction rate (cm$^{-3}$s$^{-1}$)")
        # labels and such 
        if minx < 1e-14
            minx = 1e-14
        end
        maxx = 10^(ceil(log10(maxx)))
        ax.set_xlim([minx, maxx])
        savefig(results_dir*"J_rates_$(sp)_$(filenameext).png", bbox_inches="tight", dpi=300)
    end
end

function plot_rxns(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, results_dir::String, E_prof; 
                   shown_rxns=nothing, subfolder="", plotsfolder="", dt=nothing, num="", extra_title="", 
                   plot_timescales=false, plot_total_rate_coefs=false, showonly=false, globvars...)
    #=
    Plots the production and loss rates from chemistry, transport, and both together by altitude for a given species, sp, 
    at a given snapshot in time of the atmosphere atmdict. 

    Inputs:
        sp: species to make plot for
        atmdict: the atmospheric state dictionary
        Tn, Ti, Te: atmospheric temperature arrays
        bcdict: Boundary conditions dictionary specified in parameters file
        rxnnet: Chemical reaction network
        plot_indiv_rxns: whether to plot lines for the individual chemical reactions. if false, only the total will be plotted.
        thresh: a threshhold of reaction rate. will only plot the reactions that have values above this threshhold
        subfolder: an optional subfolder within results_dir
        plotsfolder: specific folder within which the plots will go.
        dt: timestep, optional, to be included in the plot title
        num: optional string to add to filename
        extra_title: extra text for the title, e.g. "converged"
        plot_timescales: plots the timescale over which a species population changes. 
        plot_total_rate_coefs: plots chemical produciton and loss rates only (not multiplied by density)
        showonly: if true, will show the plots instead of saving them
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:Tn, :Ti, :Te, :Tp, :bcdict, :rxnnet, 
                                    :all_species, :neutral_species, :ion_species, :transport_species, :chem_species, 
                                    :molmass, :polarizability, :Tprof_for_Hs, :Tprof_for_diffusion, :Hs_dict,
                                    :alt, :n_alt_index, :dz, :num_layers, :n_all_layers, :plot_grid, :upper_lower_bdy_i, :upper_lower_bdy, :q])

    # ================================================================================
    # Plot setup stuff
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

    if plot_timescales==true
        total_ax = 4
    else
        total_ax = 3
    end

    fig, ax = subplots(1, total_ax, sharey=true, figsize=(20,6))
    # subplot order will be chemistry, transport, total.
    subplots_adjust(wspace=0.1, bottom=0.25)
    for a in ax
        plot_bg(a)
    end

    # this will be used to determine the x limits. the three entries are the chemistry panel, transport panel, and the total panel
    # these are just starting values. 
    minx = [1e-8, 1e-8, 1e-8]  
    maxx = [1e2, 1e2, 1e2]

    # ================================================================================
    # Plot the chemical production and loss rates: Panel #1

    if shown_rxns != nothing
        shown_rxns = [format_chemistry_string(r[1], r[2]) for r in shown_rxns]
    end

    # Arrays to store the total reactions per second for this species of interest
    total_prod_rate = zeros(GV.num_layers)
    total_loss_rate = zeros(GV.num_layers)

    # Arrays to hold the total chemical production and loss 
    total_chem_prod = zeros(GV.num_layers)
    total_chem_loss = zeros(GV.num_layers)
    total_chem_prod_ratecoef = zeros(GV.num_layers)
    total_chem_loss_ratecoef = zeros(GV.num_layers)

    if sp in GV.chem_species
        # --------------------------------------------------------------------------------
        # calculate reaction rates x density for all reactions of the species at each level of the atmosphere.
        # Temperature arrays include boundary layers so that diffusion can be calculated later, so we have to 
        # index these in this way to make the evaluation of chemistry reaction rate coefficients work. 
        # Entering them separately from globvars allows us to keep passing globvars as a "packed" variable but 
        # use the most recent Tn, Ti, Te according to rightmost taking precedence.
        rxd_prod, rate_coefs_prod = get_volume_rates(sp, atmdict, E_prof; species_role="product", globvars..., Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1])
        rxd_loss, rate_coefs_loss = get_volume_rates(sp, atmdict, E_prof; species_role="reactant", globvars..., Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1])

        # Water is turned off in the lower atmosphere, so we should represent that.
        if sp in [:H2O, :HDO] 
            for (prod_k, loss_k) in zip(keys(rxd_prod), keys(rxd_loss))
                rxd_prod[prod_k][1:GV.upper_lower_bdy_i] .= NaN
                rxd_loss[loss_k][1:GV.upper_lower_bdy_i] .= NaN
            end
        end

        # set up a selection of colors and linestyles so that individual reaction rate lines are distinguishable
        cols = ["#cc156d", "#7e5fe7", "#b26145", "#257950", "#265582", "#fc5931"]  # colorgorical
        ls = ["-", "--", (0, (3, 1, 1, 1)), ":"]

        # ---------------------------------------------------------------------------------
        col_i = 1
        ls_i = 1

        # Pattern to help identify whether there is more than one product of species sp in the reaction
        pat = r"((?<=--> ).+)"

        # Chemical production - add up total, and plot individual reactions if needed
        for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
            if shown_rxns != nothing
                if kv[1] in shown_rxns 
                    ax[1].semilogx(kv[2], GV.plot_grid, linestyle=ls[ls_i], marker=9, markevery=20, color=cols[col_i], linewidth=1, label=kv[1])
                    col_i = next_in_loop(col_i, length(cols))
                    ls_i = next_in_loop(ls_i, length(ls))
                end
            end
            # Add up the total chemical production. Accounts for cases where species is produced more than once. 
            prods = split(match(pat, kv[1])[1], " + ")
            num_created = count(i->(i==string(sp)), prods)
            total_chem_prod += num_created .* kv[2]
        end
        for kv in rate_coefs_prod
            total_chem_prod_ratecoef += kv[2]
        end

        # Chemical loss 
        for kv in rxd_loss
            if shown_rxns != nothing
                if kv[1] in shown_rxns
                    ax[1].semilogx(kv[2], GV.plot_grid, linestyle=ls[ls_i], marker=8, markevery=20, color=cols[col_i], linewidth=1, label=kv[1])
                    col_i = next_in_loop(col_i, length(cols))
                    ls_i = next_in_loop(ls_i, length(ls))
                end
            end
            total_chem_loss += kv[2]
        end
        for kv in rate_coefs_loss
            total_chem_loss_ratecoef += kv[2]
        end

        # Plot the totals 
        ax[1].semilogx(total_chem_prod, GV.plot_grid, color="xkcd:forest green", linestyle=(0, (4,2)), marker=9, markevery=20, linewidth=2, label="Total chemical production", zorder=5)
        ax[1].semilogx(total_chem_loss, GV.plot_grid, color="xkcd:shamrock", linestyle=(0, (4,2)), marker=8, markevery=20, linewidth=2, label="Total chemical loss", zorder=5)

        # set the x lims for chem axis
        prod_without_nans = filter(x->!isnan(x), total_chem_prod)
        loss_without_nans = filter(x->!isnan(x), total_chem_loss)
        minx[1] = minimum( [minimum(prod_without_nans), minimum(loss_without_nans)] )
        maxx[1] = maximum( [maximum(prod_without_nans), maximum(loss_without_nans)] )
        minx[1] = 10^(floor(log10(minx[1])))
        maxx[1] = 10^(ceil(log10(maxx[1])))

        ax[1].legend(fontsize=10)

        if plot_total_rate_coefs
            ax_1_2 = ax[1].twiny()
            ax_1_2.tick_params(axis="x", labelcolor="xkcd:royal purple")
            ax_1_2.set_ylabel("Rate coefficient (cm^3/s)", color="xkcd:royal purple")
            ax_1_2.semilogx(total_chem_prod_ratecoef, GV.plot_grid, color="xkcd:royal purple", linestyle=":", linewidth=3, label="Total chemical production rate coef", zorder=6)
            ax_1_2.semilogx(total_chem_loss_ratecoef, GV.plot_grid, color="xkcd:lavender", linestyle=":", linewidth=3, label="Total chemical loss rate coef", zorder=6)
            ax_1_2.legend()
        end

    else
        println()
        minx[1] = 0
        maxx[1] = 1
        ax[1].text(0.1, 225, "Chemical production/loss is off for $(sp).")
    end

    # ================================================================================
    # Plot the transport production and loss rates: panel #2

    # Calculate the transport fluxes for the species
    plottitle_ext = "" # no extra info in the plot title if flux==false 
    if sp in GV.transport_species
        transportPL = get_transport_PandL_rate(sp, atmdict; globvars...)
        # now separate into two different arrays for ease of addition.
        production_i = transportPL .>= 0  # boolean array for where transport entries > 0 (production),
        loss_i = transportPL .< 0 # and for where transport entries < 0 (loss).
        total_transport_prod = production_i .* transportPL
        total_transport_loss = loss_i .* abs.(transportPL)

        if sp in [:H2O, :HDO] # Water is turned off in the lower atmosphere, so we should represent that.
            total_transport_prod[1:GV.upper_lower_bdy_i] .= NaN
            total_transport_loss[1:GV.upper_lower_bdy_i] .= NaN
        end

        # set the x lims for transport axis. Special because total_transport_prod, and etc are incomplete arrays.
        minx[2] = 10.0^(floor(log10(minimum(abs.(transportPL)))))
        maxx[2] = 10.0^(ceil(log10(maximum(abs.(transportPL)))))

        # Plot the transport production and loss without the boundary layers
        ax[2].scatter(total_transport_prod, GV.plot_grid, color="red", marker=9, label="Total gain this layer", zorder=4)
        ax[2].scatter(total_transport_loss, GV.plot_grid, color="blue", marker=8, label="Total loss this layer", zorder=4)
        ax[2].set_xscale("log")

        ax[2].legend(fontsize=12)

        plottitle_ext = " by chemistry & transport"
    else
        total_transport_prod = zeros(GV.num_layers)
        total_transport_loss = zeros(GV.num_layers)
        minx[2] = 0
        maxx[2] = 1
        ax[2].text(0.1, 225, "Vertical transport is off for $(sp).")
    end

    # =================================================================================
    # Total production and loss from chemistry and transport: panel #3

    total_prod_rate = total_transport_prod .+ total_chem_prod
    total_loss_rate = total_transport_loss .+ total_chem_loss

    prod_without_nans = filter(x->!isnan(x), total_prod_rate)
    loss_without_nans = filter(x->!isnan(x), total_loss_rate)

    ax[3].semilogx(total_prod_rate, GV.plot_grid, color="black", marker=9, markevery=15, linewidth=2, label="Total production", zorder=3) #linestyle=(2, (1,2)), 
    ax[3].semilogx(total_loss_rate, GV.plot_grid, color="gray", marker=8, markevery=15, linewidth=2, label="Total loss", zorder=3) #linestyle=(2, (1,2)),

    minx[3] = minimum([minimum(prod_without_nans), minimum(loss_without_nans)])
    maxx[3] = maximum([maximum(prod_without_nans), maximum(loss_without_nans)])
    minx[3] = 10^(floor(log10(minx[3])))
    maxx[3] = 10^(ceil(log10(maxx[3])))

    ax[3].legend(fontsize=12)

    # ---------------------------------------------------------------------------------
    # Plot the timescale of dn/dt if requested
    if plot_timescales==true
        timescale_of_change = atmdict[sp] ./ abs.(total_prod_rate .- total_loss_rate)
        c = "teal"
        ax[total_ax].set_xscale("log")
        ax[total_ax].plot(timescale_of_change, GV.plot_grid, color=c, zorder=5)
    end

    # Final plotting tasks ============================================================
    # check for and correct any ridiculously low limits
    for i in 1:length(minx)
        if minx[i] < 1e-12
            minx[i] = 1e-12
        end
    end

    titles = ["Chemistry", "Transport", "Total chem+trans"]
    notes = ["Chemical production/loss is off for $(sp) below $(Int64(GV.upper_lower_bdy /1e5)) km.", 
             "Vertical transport is off for $(sp) below $(Int64(GV.upper_lower_bdy /1e5)) km.", 
             "Abundance of $(sp) below $(Int64(GV.upper_lower_bdy /1e5)) km is fixed."] 
    if plot_timescales==true
        titles = ["Chemistry", "Transport", "Total chem+trans", "Timescale of dn/dt (s)"]
    end
    dtstr = dt == nothing ? "" : ", $(dt)"
    titlestr = extra_title == "" ? "" : ", "*extra_title
    for i in 1:length(ax)
        ax[i].set_title(titles[i])
        ax[i].set_xlim(minx[i], maxx[i])
        ax[i].set_ylim(0, zmax/1e5)
        if sp in [:H2O, :HDO]
            ax[i].text(2*minx[i], 5, notes[i])
        end
        
    end
    suptitle("Production & loss" * plottitle_ext * ", $(string(sp))" * dtstr * titlestr, fontsize=20)
    ax[1].set_ylabel("Altitude (km)")
    ax[1].set_xlabel("Rate ("*L"cm^{-3}s^{-1})")

    if num==""
        num=extra_title
    end

    path_folders = [results_dir[1:end-1], subfolder, plotsfolder, "chem_rates_$(sp)_$(num).png"]
    filter!(e->e≠"", path_folders)  # gets rid of empty names, in case subfolder or plotsfolder hasn't been passed in
    savepathname = join(path_folders, "/")

    if showonly
        show()
    else
        savefig(savepathname, bbox_inches="tight", dpi=300)
        close(fig)
    end
end

function plot_temp_prof(Tprof_1; lbls=["Neutrals", "Ions", "Electrons"], Tprof_2=nothing, Tprof_3=nothing, savepath=nothing, showonly=false, globvars...)
    #=
    Creates a .png image of the tepmeratures plotted by altitude in the atmosphere

    Inputs:
        Tprof_1: an array of neutral temperature by altitude
        Tprof_2: same, but for the ions
        Tprof_3: same but for the electrons
        savepath: where to save the resulting .png image
        showonly: whether to just show() the figure. If false, must be accompanied by a value for savepath.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:alt])

    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    fig, ax = subplots(figsize=(4,6))
    plot_bg(ax)

    plot(Tprof_1, GV.alt./1e5, label="Neutrals", color=medgray)

    if Tprof_2 != nothing
        ax.plot(Tprof_2, GV.alt./1e5, label="Ions", color="xkcd:bright orange")
        ax.legend(fontsize=16)
        ax.set_xscale("log")
    end
    if Tprof_3 != nothing
        ax.plot(Tprof_3, GV.alt./1e5, label="Electrons", color="cornflowerblue")
        ax.legend(fontsize=16, loc=(1.05,0.8))#"center right")
        ax.set_xscale("log")
    end

    # plot the control temps

    # surface
    ax.scatter(Tprof_1[1], 0, marker="o", color=medgray, zorder=10)
    ax.text(Tprof_1[1], 0, L"\mathrm{T}_{\mathrm{surface}}"*" = $(Int64(round(Tprof_1[1], digits=0))) K ")

    # mesosphere
    meso_ind = findfirst(x->x==minimum(Tprof_1), Tprof_1)# Int64(length(Tprof_1)/2)
    # println(meso_ind)
    ax.scatter(Tprof_1[meso_ind], GV.alt[meso_ind+5]/1e5, marker="o", color=medgray, zorder=10)
    ax.text(Tprof_1[meso_ind]+5, GV.alt[meso_ind+5]/1e5, L"\mathrm{T}_{\mathrm{meso}}"*" = $(Int64(round(Tprof_1[meso_ind], digits=0))) K ")

    # exosphere
    ax.scatter(Tprof_1[end], 250, marker="o", color=medgray, zorder=10)
    ax.text(Tprof_1[end]*1.05, 240, L"\mathrm{T}_{\mathrm{exo}}"*" = $(Int64(round(Tprof_1[end], digits=0))) K ")
    
    # final labels
    ax.set_ylabel("Altitude [km]")
    ax.set_yticks(collect(0:50:GV.alt[end]/1e5))
    ax.set_yticklabels(collect(0:50:GV.alt[end]/1e5))
    ax.set_xlabel("Temperature [K]")
    ax.set_xlim(95, 2e3)
    ax.tick_params(which="both", axis="x", top=true, labeltop=true)

    if showonly==true
        show()
    else
        try
            savefig(savepath*"/temp_profile.png", bbox_inches="tight", dpi=300) 
        catch
            println("Error: You asked to save the figure but didn't provide a savepath")
        end
    end
end

function plot_water_profile(H2Oinitf, HDOinitf, nH2O, nHDO, savepath::String; showonly=false, watersat=nothing, globvars...) 
    #=
    Plots the water profile in mixing ratio and number densities, in two panels.

    Inputs:
        H2Oinitf: initial fraction of H2O in the atmosphere
        HDOinitf: same for HDO
        nH2O: number density of H2O
        nHDO: same for HDO
        savepath: where to place the image
        showonly: show rather than save image
        watersat: optional. must be a list of the saturation fractions with HDO second, 
                  i.e. [H2Osatfrac, HDOsatfrac]
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:plot_grid])

    fig = figure(figsize=(4,6))
    ax = gca()
    plot_bg(ax)
    ax.tick_params(axis="x", which="minor", bottom=true, top=true)

    
    ax1col = "#88527F"
    
    # mixing ratio in PPM axis
    ax.semilogx(convert(Array{Float64}, H2Oinitf)/1e-6, GV.plot_grid, color=ax1col, linewidth=2)
    ax.semilogx(convert(Array{Float64}, HDOinitf)/1e-6, GV.plot_grid, color=ax1col, linestyle="--", linewidth=2)
    ax.set_xlabel("Volume Mixing Ratio [ppm]", color=ax1col)
    ax.set_ylabel("Altitude [km]")
    ax.tick_params(axis="x", labelcolor=ax1col)
    ax.set_xticks(collect(logrange(1e-4, 1e2, 4)))

    # LEgend
    L2D = PyPlot.matplotlib.lines.Line2D
    lines = [L2D([0], [0], color="black"),
             L2D([0], [0], color="black", linestyle="--")]
    ax.legend(lines, [L"H_2O", "HDO"], fontsize=12)#, loc="center left")

    ax2 = ax.twiny()
    ax2col = "#429EA6" 
    ax2.tick_params(axis="x", labelcolor=ax2col)
    ax2.semilogx(convert(Array{Float64}, nH2O), GV.plot_grid, color=ax2col, linewidth=2, label=L"H$_2$O")
    ax2.semilogx(convert(Array{Float64}, nHDO), GV.plot_grid, color=ax2col, linestyle="--", linewidth=2, label="HDO")
    ax2.set_xlabel(L"Number density [cm$^{-3}$]", color=ax2col, y=1.07)
    ax2.set_xticks(collect(logrange(1e-4, 1e16, 6)))
    for side in ["top", "bottom", "left", "right"]
        ax2.spines[side].set_visible(false)
    end

    # suptitle(L"H$_2$O and HDO vertical profiles", y=1.05)
    # save it
    if showonly==true
        show()
    else
        savefig(savepath*"/water_profiles.png", dpi=300, bbox_inches="tight")
        close(fig)
    end

    if watersat != nothing
        fig, ax = subplots(figsize=(6,9))
        plot_bg(ax)
        semilogx(convert(Array{Float64}, H2Oinitf), GV.plot_grid, color=ax1col, linewidth=3, label=L"H$_2$O initial fraction")
        semilogx(convert(Array{Float64}, watersat[2:end-1]), GV.plot_grid, color="black", alpha=0.5, linewidth=3, label=L"H$_2$O saturation")
        xlabel("Mixing ratio", fontsize=18)
        ylabel("Altitude [km]", fontsize=18)
        title(L"H$_2$O saturation fraction", fontsize=20)
        ax.tick_params("both",labelsize=16)
        legend()
        savefig(savepath*"/water_MR_and_saturation.png", dpi=300, bbox_inches="tight")
    end
end

# **************************************************************************** #
#                                                                              #
#                        Boundary condition functions                          #
#                                                                              #
# **************************************************************************** #

function boundaryconditions(fluxcoef_dict::Dict; globvars...)
    #= 
    Inputs:
        fluxcoef_dict: a dictionary containing the K and D flux coefficients for every species throughout
                       the atmosphere. Format species=>Array(length(all_species), length(alt)).
                       Because it has length alt, it ALREADY INCLUDES boundary layer info in the 
                       1st and last elements. 2nd and penultimate elements are for the edge bulk layers.
        bcdict: Boundary conditions dictionary specified in parameters file, format:
                [lower bc type, lower bc value; upper bc type, upper bc value]
    Outputs:
        boundary conditions for species in a 2 x 2 matrix, format:
        [n_1 -> n_0, n_0 -> n_1;      
         n_(nl) -> n_(nl+1), n_(nl+1) -> n_(n_l)]

         Or:

         Surface [↓, ↑;     [out, in;    [#/s, #/cm³/s.;
         Top      ↑, ↓]      out, in]     #/s, #/cm³/s]

        where n_0 is the boundary layer from [-1 km, 1 km], n_1 is the first bulk layer from [1 km, 3 km],
        n_(nl) is the topmost bulk layer, and n_(nl+1) is the top boundary layer.

        Each row has two elements:
            1st element: n_bulk  -> NULL (depends on species concentration)
            2nd element: NULL -> n_bulk (independent of species concentration)

            note, technically, these are chemical equations not indicating atmospheric layers!

        More specifically, when the return value of this function is used in other functions, the first
        element in each row will eventually be multiplied by a density taken from the atmospheric 
        state dictionary, and the second element will be used as-is. That way, eventually the total
        change recorded in other functions is always #/cm³/s. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:all_species, :transport_species, :dz, :bcdict])
    
    # This is where we will store all the boundary condition numbers that will be returned 
    bc_dict = Dict{Symbol, Array{ftype_ncur}}([s=>[0 0; 0 0] for s in GV.all_species])

    for s in GV.all_species
        # retrieves user-supplied boundary conditions. When unsupplied, falls back to flux = 0 at both boundaries.
        bcs = get(GV.bcdict, s, ["f" 0.; "f" 0.]) 
        
        if issubset([s], setdiff(GV.all_species, GV.transport_species)) # non transporting species.
            bcs = ["f" 0.; "f" 0.] # Flux = 0 at both boundaries for no_transport_species.
        end

        bcvec = ftype_ncur[0 0; 0 0]
        #= the signs on these are the opposite of what you expect because it is necessary to 
          encode loss in production equations and vice versa.
          in first element, negative means you're adding something, and vice versa.
          in second element, positive means you're adding something, and vice versa. =#

        # LOWER
        if bcs[1, 1] == "n"
            # TODO: The second index can probably replace the :. 
            bcvec[1,:]=[fluxcoef_dict[s][2, :][1],
                        fluxcoef_dict[s][1, :][2]*bcs[1,2]]
        elseif bcs[1, 1] == "f"
            bcvec[1,:] = [0.0,            # no depositional flux. Could be changed if desired. 
                          -bcs[1, 2]/GV.dz]  # boundary layer to bulk layer flux. I think the negative sign actually represents depositional flux also...
        elseif bcs[1, 1] == "v"
            bcvec[1,:] = [bcs[1, 2]/GV.dz,   # depositional velocity.
                          0.0] 
            throw("Improper lower boundary condition!")
        end

        # UPPER
        if bcs[2, 1] == "n"
            bcvec[2,:] = [fluxcoef_dict[s][end-1, :][2],         # bulk layer to boundary layer
                          fluxcoef_dict[s][end, :][1]*bcs[2,2]]  # boundary layer to bulk layer
        elseif bcs[2, 1] == "f"  
            # n_(top+1) -> n_top flux is constrained to this value (negative because of the way transport is encoded as chemistry.)
            # allows terms not proportional to density - requires another reaction.
            bcvec[2,:] = [0.0,            # bulk layer to boundary layer
                          -bcs[2, 2]/GV.dz]  # boundary layer to bulk layer. A negative sign here represents removal from the atmosphere. If it was positive it would be incoming. 
        elseif bcs[2, 1] == "v"
            bcvec[2,:] = [bcs[2, 2]/GV.dz,   # bulk layer to boundary layer: velocity of escaping particles
                          0.0]            # boundary layer to bulk layer: who cares what the velocity of incoming particles from space is?
        else
            throw("Improper upper boundary condition!")
        end
        
        bc_dict[s] .= bcvec
    end

    return bc_dict
end

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
    @assert all(x->x in keys(GV), [:all_species, :bcdict, :molmass, :polarizability, :neutral_species, :q])
   
    # Calculate as if it was a neutral
    D_arr[:] .= (diffparams(sp)[1] .* 1e17 .* T_arr .^ (diffparams(sp)[2])) ./ n_tot(atmdict; GV.all_species, GV.n_alt_index)

    # a place to store the density array
    species_density = zeros(size(T_arr))
    
    # If an ion, overwrite with the ambipolar diffusion
    if charge_type(sp) == "ion"
        sum_nu_in = zeros(size(T_arr))

        # mi = GV.molmass[sp] .* mH
        # create the sum of nu_in. Note that this depends on density, but we only have density for the real layers,
        # so we have to assume the density at the boundary layers is the same as at the real layers.
        for n in GV.neutral_species
            species_density = atmdict[n]

            # This sets the species density to a boundary condition if it exists. 
            bccheck = get(GV.bcdict, n, ["f" 0.; "f" 0.])
            if bccheck[1,1] == "n"
                species_density[1] = bccheck[1,2]
            end
            if bccheck[2,1] == "n"  # currently this should never apply.
                species_density[end] = bccheck[2,2]
            end

            # mu_in = (1 ./ mi .+ 1 ./ (GV.molmass[n] .* mH)) .^ (-1) # reduced mass in g
            sum_nu_in .+= 2 .* pi .* (((GV.polarizability[n] .* GV.q .^ 2) ./ reduced_mass(GV.molmass[sp], GV.molmass[n])) .^ 0.5) .* species_density

        end
        
        D_arr .= (kB .* T_arr) ./ (GV.molmass[sp] .* sum_nu_in)
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

function fluxcoefs(species_list::Vector{Any}, K, D, H0; globvars...) # T_neutral, T_plasma, Hs, alts, dz) 
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
    fluxcoef_dict = Dict{Symbol, Array{ftype_ncur}}([s=>fill(0., length(GV.alt), 2) for s in species_list])

    for s in species_list
        layer_below_coefs, layer_above_coefs = fluxcoefs(s, K, D, H0; globvars...) #alts, dz, T_neutral, T_plasma, Hs) 
        fluxcoef_dict[s][:, 1] .= layer_below_coefs
        fluxcoef_dict[s][:, 2] .= layer_above_coefs
    end

    return fluxcoef_dict
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

function get_flux(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, globvars...)
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
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :Tp, :bcdict, :all_species, :neutral_species, :transport_species, :molmass, :alt, :n_alt_index, :polarizability, 
                                   :num_layers, :dz, :Tprof_for_Hs, :Tprof_for_diffusion, :n_all_layers, :q])
    
    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(GV.Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, globvars...)
    Hs_dict = Dict{Symbol, Vector{ftype_ncur}}([sp=>scaleH(GV.alt, sp, GV.Tprof_for_Hs[charge_type(sp)]; globvars...) for sp in GV.all_species])
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict; globvars...)
    bc_dict = boundaryconditions(fluxcoefs_all; globvars...)

    # each element in bulk_layer_coefs has the format [downward flow (i to i-1), upward flow (i to i+1)].  units 1/s
    bulk_layer_coefs = fluxcoefs_all[sp][2:end-1]

    bcs = bc_dict[sp]
    
    net_bulk_flow = fill(convert(ftype_ncur, NaN), GV.n_all_layers-1)  # units #/cm^3/s; tracks the cell boundaries, of which there are length(alt)-1

    # We will calculate the net flux across each boundary, with sign indicating direction of travel.
    # Units for net bulk flow are always: #/cm³/s. 
    net_bulk_flow[1] = (bcs[1, 2]                  # increase of the lowest atmospheric layer's density. 0 unless the species has a density or flux condition
                       - atmdict[sp][1]*bcs[1, 1]) # lowest atmospheric layer --> surface ("depositional" term). UNITS: #/cm³/s. 
                        
    for ialt in 2:GV.num_layers  # now iterate through every cell boundary within the atmosphere. boundaries at 3 km, 5...247. 123 elements.
        # UNITS for both of these terms:  #/cm³/s. 
        net_bulk_flow[ialt] = (atmdict[sp][ialt-1]*bulk_layer_coefs[ialt-1][2]   # coming up from below: cell i-1 to cell i. Should be positive * positive
                              - atmdict[sp][ialt]*bulk_layer_coefs[ialt][1])     # leaving to the layer below: downwards: cell i to cell i-1
    end

    # now the top boundary - between 124th atmospheric cell (alt = 249 km)
    net_bulk_flow[end] = (atmdict[sp][end]*bcs[2, 1] # into exosphere from the cell. UNITS: #/cm³/s. 
                         - bcs[2, 2]) # into top layer from exosphere. negative because the value in bcs is negative. do not question this. UNITS: #/cm³/s. 
                
    return net_bulk_flow .* GV.dz # now it is a flux. hurrah.
end

function get_transport_PandL_rate(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; globvars...)
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
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :Tp, :bcdict, :all_species, :neutral_species, :transport_species, :molmass, :alt, :n_alt_index, :polarizability, 
                                   :num_layers, :dz, :Hs_dict, :Tprof_for_Hs, :Tprof_for_diffusion, :n_all_layers, :q])

    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(GV.Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(GV.all_species, atmdict, D_arr; globvars...) 
    # Hs_dict = Dict{Symbol, Vector{ftype_ncur}}([sp=>scaleH(GV.alt, sp, GV.Tprof_for_Hs[charge_type(sp)]; globvars...) for sp in all_species]) 
    fluxcoefs_all = fluxcoefs(GV.all_species, Keddy_arr, Dcoef_dict, H0_dict; globvars...)

    # For the bulk layers only to make the loops below more comprehendable: 
    fluxcoefs_bulk_layers = Dict([s=>fluxcoefs_all[s][2:end-1, :] for s in keys(fluxcoefs_all)])

    bc_dict = boundaryconditions(fluxcoefs_all; globvars...)

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
    transport_PL[end] = ((thesebcs[2, 2] # in from upper boundary layer
                          - atmdict[sp][end]*thesebcs[2, 1]) # leaving out the top boundary
                        + (-atmdict[sp][end]*fluxcoefs_bulk_layers[sp][end, 1] # leaving out to layer below
                           +atmdict[sp][end-1]*fluxcoefs_bulk_layers[sp][end-1, 2])) # coming in to top layer from layer below

    # Use these for a sanity check if you like. 
    # println("Activity in the top layer for sp $(sp):")
    # println("In from upper bdy: $(thesebcs[2, 2])")
    # println("Out through top bdy: $(atmdict[sp][end]*thesebcs[2, 1])")
    # println("Down to layer below: $(-atmdict[sp][end]*fluxcoefs_all[sp][end, 1])")
    # println("In from layer below: $(atmdict[sp][end-1]*fluxcoefs_all[sp][end-1, 2])")
    return transport_PL
end

function Keddy(z, nt)
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

function scaleH(z::Vector{Float64}, sp::Symbol, T::Array; globvars...)
    #=
    Input:
        z: Altitudes in cm
        sp: Speciecs to calculate for
        T: temperature array for this species
    Output: 
        species-specific scale height at all altitudess
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
        Mean atmospheric scale height at all altitudes
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :molmass, :n_alt_index])

    mm_vec = meanmass(atmdict; globvars...)#, all_species, mmass) # vector version.

    return @. kB*T/(mm_vec*mH*marsM*bigG)*(((GV.alt+radiusM))^2)
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
    @assert all(x->x in keys(GV), [:Tn, :Tp, :bcdict, :all_species, :neutral_species, :molmass, :alt, :n_alt_index, :polarizability, 
                                   :Tprof_for_diffusion, :q])

    ncur_with_bdys = ncur_with_boundary_layers(atmdict; GV.n_alt_index, GV.all_species)
    
    K = Keddy(GV.alt, n_tot(ncur_with_bdys; GV.all_species, GV.n_alt_index))
    H0_dict = Dict{String, Vector{ftype_ncur}}("neutral"=>scaleH(ncur_with_bdys, GV.Tn; globvars...),
                                               "ion"=>scaleH(ncur_with_bdys, GV.Tp; globvars...))
    
    # Molecular diffusion is only needed for transport species, though.  
    Dcoef_dict = Dict{Symbol, Vector{ftype_ncur}}([s=>deepcopy(Dcoef!(D_coefs, GV.Tprof_for_diffusion[charge_type(s)], s, ncur_with_bdys; globvars...)) for s in species_list])

    return K, H0_dict, Dcoef_dict
end

function update_transport_coefficients(species_list, atmdict::Dict{Symbol, Vector{ftype_ncur}}, activellsp, D_coefs; globvars...) 
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
    @assert all(x->x in keys(GV), [:Tn, :Tp, :Hs_dict, :bcdict, :all_species, :neutral_species, :transport_species, :molmass, :alt, :n_alt_index, :polarizability, 
                                   :num_layers, :dz, :Tprof_for_diffusion, :n_all_layers, :q])
    
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

    # transport coefficients for boundary layers
    bc_dict = boundaryconditions(fluxcoefs_all; globvars...)
    tlower = permutedims(reduce(hcat, [bc_dict[sp][1,:] for sp in GV.transport_species]))
    tupper = permutedims(reduce(hcat, [bc_dict[sp][2,:] for sp in GV.transport_species]))

    return tlower, tup, tdown, tupper
end

# **************************************************************************** #
#                                                                              #
#                             Chemistry Functions                              #
#                                                                              #
# **************************************************************************** #

# Network formatting and loading ====================================================

function format_neutral_network(reactions_spreadsheet, used_species; verbose=false)
    #=
    Inputs:
        reactions_spreadsheet: Properly formatted .xlsx of chemical reactions with 
                               sheet names: Neutral reactions, Ion reactions.
        used_species: List of species in use, really just used once to find unneeded
                      reactions
    Outputs:
        lists of several reaction types, all neutrals, in form [[R1, R2], [P1, P2], :(k)].
        3 reactants or products are also possible.
    =#
    n_table_raw = DataFrame(XLSX.readtable(reactions_spreadsheet, "Neutral reactions")...)
    
    # Replace missing values with proper stuff 
    replace!(n_table_raw."R2", missing=>"none");
    replace!(n_table_raw."R3", missing=>"none");
    replace!(n_table_raw."P2", missing=>"none");
    replace!(n_table_raw."P3", missing=>"none");
    replace!(n_table_raw."M2", missing=>1);
    replace!(n_table_raw."M1", missing=>1);
    replace!(n_table_raw."pow", missing=>0);
    replace!(n_table_raw."BR", missing=>1);

    # Get rid of reactions with a species we don't use
    n_table = filter(row -> all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2), Symbol(row.R3)]), n_table_raw)
    if verbose
        println("Removed reactions: $(filter(row->!all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2), Symbol(row.R3)]), n_table_raw))")
    end

    # Set up the array to store neutral network 
    neutral_network = vec(Array{Any}(undef, nrow(n_table)))

    # Just some quick error handling
    @assert all(x->x>=0, n_table."kA")
    @assert all(x->x>=0, n_table."k0A")
    
    # Get all the counts of each type 
    counts = Dict([i=>count(x->x==i, n_table[:, "type"]) for i in 1:6])

    n = 1 # Keeps track of iteration, since some counts may be 0.

    for i in 1:1+counts[1]-1 # Type 1: Pressure independent unimolecular rxns, high pressure limit
        @assert n_table[i, "type"]==1
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]

        neutral_network[i] = [reactants, products, :($(make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn", 
                                                                   n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])))]
        n += 1
    end

    for i in n:n+counts[2]-1 # Type 2: P independent bimolecular; use k (same as *k inf), units cm^3 s^-1
        @assert n_table[i, "type"]==2
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        neutral_network[i] = [reactants, products, make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn",
                                                               n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])]
        n += 1
    end 

    for i in n:n+counts[3]-1 # Type 3: P dependent bimolecular.
        @assert n_table[i, "type"]==3
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        neutral_network[i] = [reactants, products, make_Troe([n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"]], 
                                                             [n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"]],
                                                              n_table[i, "F"], n_table[i, "M2"], n_table[i, "M1"], 
                                                              n_table[i, "pow"], n_table[i, "BR"])]
        n =+ 1
    end 

    for i in n:n+counts[4]-1 # Type 4: P dependent association rxns
        @assert n_table[i, "type"]==4
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]        
        neutral_network[i] = [reactants, products, make_modified_Troe([n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"]], 
                                                                      [n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"]],
                                                                      [n_table[i, "kradA"], n_table[i, "kradB"], n_table[i, "kradC"]],
                                                                       n_table[i, "F"], n_table[i, "M2"], n_table[i, "M1"], 
                                                                       n_table[i, "pow"], n_table[i, "BR"])]
        n += 1
    end 

    for i in n:n+counts[5]-1 # Three body reactions 
        @assert n_table[i, "type"]==5
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        
        k0 = make_k_expr(n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"], "Tn", 
                         n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])             
        kinf = make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn", 
                           n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])

        neutral_network[i] = [reactants, products, :($k0 .* M ./ (1 .+ $k0 .* M ./ $kinf).*0.6 .^ ((1 .+ (log10.($k0 .* M ./ $kinf)) .^2).^-1.0))]
        n += 1
    end 

    for i in n:n+counts[6]-1 # Three body reactions 
        @assert n_table[i, "type"]==6
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        k0 = make_k_expr(n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"], "Tn", 
                         n_table[i, "M2"], n_table[i, "M2"], n_table[i, "pow"], n_table[i, "BR"])             
        kinf = make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn", 
                           n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])

        neutral_network[i] = [reactants, products, :($k0 ./ (1 .+ $k0 ./ ($kinf ./ M)).*0.6 .^ ((1 .+ (log10.($k0 ./ ($kinf .* M))) .^2).^-1.0))]
        n += 1
    end

    return neutral_network, n_table."Num"
end

function format_ion_network(reactions_spreadsheet, used_species; verbose=false)
    #=
    Inputs:
        reactions_spreadsheet: Properly formatted .xlsx of chemical reactions with 
                               sheet names: Neutral reactions, Ion reactions.
    Outputs:
        ion_network, list of all the ion reactions in format [[R1, R2], [P1, P2], :(k)].
        3 reactants or products are also possible.
    =#
    
    ion_table_raw = DataFrame(XLSX.readtable(reactions_spreadsheet, "Ion reactions")...)
    
    # Replace missing values with proper stuff 
    replace!(ion_table_raw."R2", missing=>"none");
    replace!(ion_table_raw."P2", missing=>"none");
    replace!(ion_table_raw."P3", missing=>"none");
    replace!(ion_table_raw."M2", missing=>1);
    replace!(ion_table_raw."M1", missing=>1);
    replace!(ion_table_raw."pow", missing=>0);
    replace!(ion_table_raw."BR", missing=>1);

    # Filter out reactions with reactants we don't use
    ion_table = filter(row->all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2)]), ion_table_raw)
    if verbose
        println("Removed reactions: $(filter(row->!all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2)]), ion_table_raw))")
    end 

    numrows = size(ion_table)[1]
    ion_network = Array{Any}(undef, numrows, 1)

    # Just some quick error handling
    @assert all(x->x>=0, ion_table."kA")
    @assert all(x->x>=0, ion_table."k0A")
    @assert all(x->x>=0, ion_table."kradA")
    
    # Fill in the network
    for i in range(1, stop=numrows)

        # Set the temperature term
        Tstr = ion_table[i, "type"] == -4 ? "Te" : "Ti"
            
        kA = ion_table[i, "kA"]
        kB = ion_table[i, "kB"]
        kC = ion_table[i, "kC"]
        M2 = ion_table[i, "M2"]
        M1 = ion_table[i, "M1"]
        pow = ion_table[i, "pow"]
        BR = ion_table[i, "BR"]
       
        # Fill out the entry
        ion_network[i] = [[Symbol(ion_table[i, "R1"]), Symbol(ion_table[i, "R2"])],  # REACTANT LIST
                          [Symbol(j) for j in filter!(j->j!="none", [ion_table[i, "P1"], ion_table[i, "P2"], ion_table[i, "P3"]])], # PRODUCT LIST
                          make_k_expr(kA, kB, kC, Tstr, M2, M1, pow, BR)] # RATE COEFFICIENT EXPRESSION
    end
        
    return ion_network, ion_table."Num"
end

function load_reaction_network(spreadsheet; ions_on=true, get_inds=false, globvars...)
    #=
    Inputs:
        spreadsheet: path and filename of the reaction network spreadsheet
        Jrl: Jratelist, defined in PARAMETERS.jl
        absorber: Photochemistry absorbing species for each Jrate, defined in PARAMETERS.jl
        photo_prods: list of products for each Jrate, defined in PARAMETERS.jl
    Output:
        reaction_network, a vector of vectors in the format
        [Symbol[reactants...], Symbol[products...], :(rate coefficient expression)
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :Jratelist, :absorber, :photolysis_products])

    if ions_on == true
        ionnet, i_nums = format_ion_network(spreadsheet, GV.all_species)
    else 
        ionnet = []
        i_nums = []
    end 

    neutral_net, n_nums = format_neutral_network(spreadsheet, GV.all_species)
    Jrxns = [[[GV.absorber[Jr]], GV.photolysis_products[Jr], Jr] for Jr in GV.Jratelist]
    whole_network = [Jrxns..., neutral_net..., ionnet...]
    
    if get_inds
        return whole_network, n_nums, i_nums
    else
        return whole_network
    end
end

function make_k_expr(A, B, C, T::String, M2, M1, pow, BR)
    #=
    Constructs the modified Arrhenius equation, with leading multipliers for
    mass scaling and branching ratios.

    Inputs:
        A, B, C: Arrhenius equation parameters
        T: Temperature string to use; can be Tn, Ti, Te.
        M2, M1: Masses of the heavier and lighter isotope respectively. 1 if no mass scaling
        pow: -0.5. Passed in, in case it needs to change later. 
        BR: Branching ratio.
    Output:
        k = BR * (M2/M1)^pow * A * T^B * exp(C/T)
    =#
    if (B == 0) & (C == 0) # constant rate coefficient
        k = :($(((M2/M1)^pow)*BR*A))
    elseif (B != 0) & (C == 0)
        k = :($(((M2/M1)^pow)*BR*A) .* $(Meta.parse(T)) .^ $(B))
    elseif (B != 0) & (C != 0)
        k = :($(((M2/M1)^pow)*BR*A) .* $(Meta.parse(T)) .^ $(B) .* exp.($(C) ./ $(Meta.parse(T))))
    elseif (B == 0) & (C != 0)
        k = :($(((M2/M1)^pow)*BR*A) .* exp.($(C) ./ $(Meta.parse(T))))
    end
    return k
end

function make_Troe(k0_ABC, kinf_ABC, F, M2, M1, pow, BR)
    #=
    Make expression for type 3 reactions, of which there are currently none in the network (Feb 2022).

    Input:
        k0_ABC, kinf_ABC, kR_ABC: A, B, C arrhenius parameters for k0 (low pressure limit), kinf (high pressure limit), krad.
        F: Troe parameter
        M2, M1: Masses for heavier and lighter isotope for mass scaling
        pow: power (usually -0.5) for mass scaling
        BR: Branching ratio 
    Output: a symbolic expression for the rate coefficient
    =#
    k0 = make_k_expr(k0_ABC..., "Tn", M2, M1, pow, BR)
    kinf = make_k_expr(kinf_ABC..., "Tn", M2, M1, pow, BR)
    
    if F == 0
        return :(($(k0) .* $(kinf) .* M) ./ ($(k0) .* M .+ $(kinf)))  # Updated to match Roger's code 4 Feb 2021
    else        
        FF = troe_expr(k0, kinf, F)
        return :($(FF) .* ($(k0) .* $(kinf)) ./ ($(k0) .* M .+ $(kinf)))  # Updated to match Roger's code 5 Feb 2021
    end
end

function make_modified_Troe(k0_ABC, kinf_ABC, kR_ABC, F, M2, M1, pow, BR)
    #=
    "Type 4" pressure dependent association reactions. 
    Roger has the result for F=0 in his code as k_tot = k2 + k0*kinf*M/(k0*M+kinf), with k2 (kR in the appendix from Vuitton)
    being the general form and using the third set of A, B, C coefficients in his file. Those columns are 0 in the file he gave me, but
    I have coded it like this in case I ever get one where it's not 0.
    
    Input:
        k0_ABC, kinf_ABC, kR_ABC: A, B, C arrhenius parameters for k0 (low pressure limit), kinf (high pressure limit), krad.
        F: Troe parameter
        M2, M1: Masses for heavier and lighter isotope for mass scaling
        pow: power (usually -0.5) for mass scaling
        BR: Branching ratio 
    Output: a symbolic expression for the rate coefficient
    =#
    k0 = make_k_expr(k0_ABC..., "Tn", M2, M1, pow, BR)
    kinf = make_k_expr(kinf_ABC..., "Tn", M2, M1, pow, BR)
    kR = make_k_expr(kR_ABC..., "Tn", M2, M1, pow, BR)
    
    if F == 0
        return :($(kR) .+ ($(k0) .* $(kinf) .* M) ./ ($(k0) .* M .+ $(kinf)))
    else
        FF = troe_expr(k0, kinf, F)
        
        # for this one, we need the minimum of k_inf and a more complicated expression.
        # Format for minimum of A and B is: :(min.($:(A), $:(B)))
        return :(min.($:($kinf), $:($(kR) .+ ($(FF) .* $(k0) .* $(kinf) .* M) ./ ($(k0) .* M .+ $(kinf)))))
    end
end

function troe_expr(k0, kinf, F)
    #=
    Put together a really nasty expression
    =#
    outer_numer = :(log10.($(F)))
    inner_numer = :(log10.(($(k0) .* M) ./ ($(kinf))) .- 0.4 .- 0.67 .* log10.($(F)))
    inner_denom = :(0.75 .- 1.27 .* log10.($(F)) .- 0.14 .* (log10.(($(k0) .* M) ./ ($(kinf))) .- 0.4 .- 0.67 .* log10.($(F))))
    outer_denom = :(1 .+ (($(inner_numer)) ./ ($(inner_denom))) .^ 2)
    FF = :(10 .^ (($(outer_numer)) ./ ($(outer_denom))))
    return FF
end

# end formatting ======================================================================

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

# function replace_E(eqn, insertme)
#     i = findall(x->x==:E, eqn)

#     if length(i)!=0
#         deleteat!(eqn, i)
#         for j in 1:length(insertme)
#             insert!(eqn, i[1], insertme[j])
#         end
#     end
    
#     return eqn 
# end

function chemical_jacobian(specieslist, dspecieslist; diff_wrt_e=true, diff_wrt_m=true, globvars...)
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
    @assert all(x->x in keys(GV), [:chem_species, :transport_species, :chemnet, :transportnet])

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)  # this is the active species. 
    ndspecies = length(dspecieslist)  # this is the species with respect to which we differentiate
    if diff_wrt_e==true
        @assert all(x->x in keys(GV), [:ion_species])
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

function eval_rate_coef(atmdict::Dict{Symbol, Vector{ftype_ncur}}, krate::Expr, E_prof; globvars...)
    #=
    Evaluates a chemical reaction rate coefficient, krate, for all levels of the atmosphere. 

    Input:
        atmdict: the atmospheric state dictionary
        krate: rate coefficient for a single reaction
        tn, _i, _e: temperature profiles for neutrals, ions, electrons
    Output:
        rate_coefficient: evaluated rate coefficient at all atmospheric layers
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :all_species, :ion_species])

    # Set stuff up
    eval_k = mk_function(:((Tn, Ti, Te, M, E) -> $krate))

    return eval_k(GV.Tn, GV.Ti, GV.Te, sum([atmdict[sp] for sp in GV.all_species]), E_prof) # sum([atmdict[sp] for sp in GV.ion_species])
end 

function get_column_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}; which="all", sp2=nothing, role="product", startalt_i=1, globvars...)
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
                    useful, for example, for water, which is fixed at 80 km so we don't care what produces/consumes it.
    Output:
        sorted: Total column rates for all reactions of species sp. Sorted, in order of largest rate to smallest. NOT a dictionary.
                sorted[1] is the top production mechanism, e.g.
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :all_species, :ion_species, :rxnnet, :num_layers, :dz])
    
    rxd, coefs = get_volume_rates(sp, atmdict; species_role=role, which=which, globvars...)# Tn[2:end-1], Ti[2:end-1], Te[2:end-1], all_species, ionsp, rxnnet, numlyrs, 
                                   
    # Make the column rates dictionary for production
    columnrate = Dict()
    for k in keys(rxd)
        columnrate[k] = sum(rxd[k][startalt_i:end] .* GV.dz)
    end
    
    # Optionally one can specify a second species to include in the sorted result, i.e. a species' ion.
    if sp2 != nothing
        rxd2, coefs2 = get_volume_rates(sp2, atmdict; species_role=role, which=which, globvars...) # Tn[2:end-1], Ti[2:end-1], Te[2:end-1], all_species, ionsp, rxnnet, numlyrs, 

        columnrate2 = Dict()

        for k in keys(rxd2)
            columnrate2[k] = sum(rxd2[k][startalt_i:end] .* GV.dz)
        end

        colrate_dict = merge(columnrate, columnrate2)
    else
        colrate_dict = columnrate
    end
    
    sorted = sort(collect(colrate_dict), by=x->x[2], rev=true)
    
    return sorted
end

function get_volume_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, E_prof; species_role="both", which="all", globvars...)
    #=
    Input:
        sp: Species name
        atmdict: Present atmospheric state dictionary
        Tn, Ti, Te: temperature arrays
        species_role: whether to look for the species as a reactant, product, or both.  If it has a value, so must species.
        which: "all", "Jrates", "krates". Whether to fill the dictionary with all reactions, only photochemistry/photoionization 
               (Jrates) or only chemistry (krates).
    Output: 
        rxn_dat: Evaluated rates, i.e. k[A][B], units #/cm^3/s for bimolecular rxns
        rate_coefs: Evaluated rate coefficients for each reaction 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :Ti, :Te, :all_species, :ion_species, :rxnnet, :num_layers])

    # Make sure temperatures are correct format
    @assert length(GV.Tn)==GV.num_layers 
    @assert length(GV.Ti)==GV.num_layers
    @assert length(GV.Te)==GV.num_layers

    # Fill in the rate x density dictionary ------------------------------------------------------------------------------
    rxn_dat =  Dict{String, Array{ftype_ncur, 1}}()
    rate_coefs = Dict{String, Array{ftype_ncur, 1}}()

    # Select either photodissociation or bi-/tri-molecular reactions
    if which=="Jrates"
        selected_rxns = filter(x->occursin("J", string(x[3])), GV.rxnnet)
    elseif which=="krates" 
        selected_rxns = filter(x->!occursin("J", string(x[3])), GV.rxnnet)
    elseif which=="all"
        selected_rxns = deepcopy(GV.rxnnet)
    end

    # now select only the reactions with the species in question (there is probably a less repetative way to write this)
    filtered_rxn_list = Any[]

    # need to make a regular expression containing the species that will ignore cases where
    # the species string appears as a subset of another species string, i.e. when O2pl is 
    # detected as being in the reaction CO2 -> CO2pl.
    species_re = r"\b"*string(sp)*r"\b"

    if species_role=="reactant"
        found_rxns = filter(x->(occursin(species_re, string(x[1]))), selected_rxns)
        filtered_rxn_list = vcat(filtered_rxn_list, found_rxns)
    elseif species_role=="product"
        found_rxns = filter(x->(occursin(species_re, string(x[2]))), selected_rxns)
        filtered_rxn_list = vcat(filtered_rxn_list, found_rxns)
    elseif species_role=="both"
        found_rxns = filter(x->(occursin(species_re, string(x[1])) || occursin(species_re, string(x[2]))), selected_rxns)
        filtered_rxn_list = vcat(filtered_rxn_list, found_rxns)
    end

    filtered_rxn_list = unique(filtered_rxn_list)  # gets rid of duplicates since occursin() is greedy
    
    for rxn in filtered_rxn_list
        # get the reactants and products in string form for use in plot labels
        rxn_str = format_chemistry_string(rxn[1], rxn[2])

        # Fill in rate coefficient * species density for all reactions
        if typeof(rxn[3]) == Symbol # for photodissociation
            rxn_dat[rxn_str] = atmdict[rxn[1][1]] .* atmdict[rxn[3]]
            rate_coefs[rxn_str] = atmdict[rxn[3]]
        else                        # bi- and ter-molecular chemistry
            density_prod = reactant_density_product(atmdict, rxn[1]; globvars...)
            thisrate = typeof(rxn[3]) != Expr ? :($rxn[3] + 0) : rxn[3]
            rate_coef = eval_rate_coef(atmdict, thisrate, E_prof; globvars...)

            rxn_dat[rxn_str] = density_prod .* rate_coef # This is k * [R1] * [R2] where [] is density of a reactant. 
            if typeof(rate_coef) == Float64
                rate_coef = rate_coef * ones(GV.num_layers)
            end
            rate_coefs[rxn_str] = rate_coef
        end
    end

    return rxn_dat, rate_coefs
end

function get_volume_rates(sp::Symbol, source_rxn::Vector{Any}, atmdict::Dict{Symbol, Vector{ftype_ncur}}; globvars...)
    #=
    Override to call for a single reaction. Useful for doing non-thermal flux boundary conditions.
    Input:
        sp: Species name
        source_rxn: chemical reaction for which to get the volume rate 
        atmdict: Present atmospheric state dictionary
        Tn, Ti, Te: temperature arrays
      Output: 
        vol_rates: Evaluated rates, i.e. k[A][B], units #/cm^3/s for bimolecular rxns, for the whole atmosphere.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:Tn, :Ti, :Te, :all_species, :ion_species, :num_layers])

    # Make sure temperatures are correct format
    @assert length(GV.Tn)==GV.num_layers 
    @assert length(GV.Ti)==GV.num_layers
    @assert length(GV.Te)==GV.num_layers

    # Fill in the rate x density dictionary ------------------------------------------------------------------------------
    vol_rates = Array{ftype_ncur}(undef, GV.num_layers, 1)

    # get the reactants and products in string form for use in plot labels
    rxn_str = format_chemistry_string(source_rxn[1], source_rxn[2])

    # Fill in rate coefficient * species density for all reactions
    if typeof(source_rxn[3]) == Symbol # for photodissociation
        vol_rates = atmdict[source_rxn[1][1]] .* atmdict[source_rxn[3]]
    else                        # bi- and ter-molecular chemistry
        thisrate = typeof(source_rxn[3]) != Expr ? :($source_rxn[3] + 0) : source_rxn[3]
        rate_coef = eval_rate_coef(atmdict, thisrate; globvars...)#, Tn, Ti, Te, all_species, ionsp, numlyrs)
        vol_rates = reactant_density_product(atmdict, source_rxn[1]; globvars...)#=, all_species, ionsp, numlyrs)=# .* rate_coef # This is k * [R1] * [R2] where [] is density of a reactant. 
    end
    
    return vol_rates
end

function getrate(sp::Symbol; chemistry_on=true, transport_on=true, sepvecs=false, globvars...)
    #=
    Creates a symbolic expression for the rate at which a given species is
    either produced or lost due to chemical reactions or transport.

    Input:
        sp: species for which to get the rate 
        chemnet: chemistry reaction array
        transportnet: transport network array
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
    @assert all(x->x in keys(GV),  [:chemnet, :transportnet, :chem_species, :transport_species])


    # This block will return the total net change for the species, P - L.
    if sepvecs == false
        rate = :(0.0)
        if issubset([sp], GV.chem_species) && chemistry_on
            rate = :($rate 
                     + $(production_rate(sp, GV.chemnet, sepvecs=sepvecs)) 
                     - $(loss_rate(sp, GV.chemnet, sepvecs=sepvecs)) 
                    )
        end
        if issubset([sp], GV.transport_species) && transport_on
            rate = :($rate 
                     + $(production_rate(sp, GV.transportnet, sepvecs=sepvecs)) 
                     - $(loss_rate(sp, GV.transportnet, sepvecs=sepvecs))
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
            transprod_rate = production_rate(sp, GV.transportnet, sepvecs=sepvecs)
            transloss_rate = loss_rate(sp, GV.transportnet, sepvecs=sepvecs)
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

function reactant_density_product(atmdict::Dict{Symbol, Vector{ftype_ncur}}, reactants; globvars...)
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
    @assert all(x->x in keys(GV),  [:all_species, :ion_species, :num_layers])


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

function rxns_where_species_is_observer(sp, chemnet)
    #=
    Input:
        sp: Species which may be an observer
        chemnet: chemistry reaction network
    Output: either:
        nothing: if the species was not found to be an observer in any reaction
        twosides: List of reactions where the species is found on both LHS and RHS
    =#

    krate_rxns = filter(x->!occursin("J", string(x[3])), chemnet)
    
    twosides = []

    flag = false
    for rxn in krate_rxns
        if in(sp, rxn[1]) && in(sp, rxn[2])
            flag = true
            twosides = push!(twosides, rxn)
        end
    end

    if flag == true 
        return twosides
    else
        return nothing
    end
end

# **************************************************************************** #
#                                                                              #
#                     Photochemical Equilibrium Functions                      #
#                                                                              #
# **************************************************************************** #

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

# **************************************************************************** #
#                                                                              #
#                          Photochemistry Functions                            #
#                                                                              #
# **************************************************************************** #

function binupO2(list)
    #=
    Mike originally wrote this
    =#
    ret = Float64[];
    for i in [176:203;]
        posss = getpos(list[:,1],x->i<x<i+1)
        dl = diff([map(x->list[x[1],1],posss); i])
        x0 = map(x->list[x[1],2], posss)
        x1 = map(x->list[x[1],3], posss)
        x2 = map(x->list[x[1],4], posss)
        ax0 = reduce(+,map(*,x0, dl))/reduce(+, dl)
        ax1 = reduce(+,map(*,x1, dl))/reduce(+, dl)
        ax2 = reduce(+,map(*,x2, dl))/reduce(+, dl)
        append!(ret,[i+0.5, ax0, ax1, ax2])
    end
    return transpose(reshape(ret, 4, 203-176+1))
end

function co2xsect(co2xdata, T)
    #=
    Makes an array of CO2 cross sections at temperature T, of format 
    [wavelength in nm, xsect]. 

    Data are available at 195K and 295K, so crosssections at temperatures between 
    these are created by first calculating the fraction that T is along the 
    line from 195 to 295, and then calculating a sort of weighted average of 
    the crosssections at 195 and 295K based on that information. 
    =#
    clamp(T, 195, 295)
    Tfrac = (T-195)/(295-195)

    arr = [co2xdata[:,1]; (1-Tfrac)*co2xdata[:,2] + Tfrac*co2xdata[:,3]]
    reshape(arr, length(co2xdata[:,1]),2)
end

function h2o2xsect_l(l, T) 
    #=
    from 260-350 the following analytic calculation fitting the
    temperature dependence is recommended by Sander 2011.

    Analytic calculation of H2O2 cross section using temperature dependencies
    l: wavelength in nm
    T: temperature in K
    =#
    l = clamp(l, 260, 350)
    T = clamp(T, 200, 400)

    A = [64761., -921.70972, 4.535649,
         -0.0044589016, -0.00004035101,
         1.6878206e-7, -2.652014e-10, 1.5534675e-13]
    B = [6812.3, -51.351, 0.11522, -0.000030493, -1.0924e-7]

    lpowA = map(n->l^n,[0:7;])
    lpowB = map(n->l^n,[0:4;])

    expfac = 1.0/(1+exp(-1265/T))

    return 1e-21*(expfac*reduce(+, map(*, A, lpowA))+(1-expfac)*reduce(+, map(*, B, lpowB)))
end

function h2o2xsect(h2o2xdata, T) 
    #=
    stitches together H2O2 cross sections, some from Sander 2011 table and some
    from the analytical calculation recommended for 260-350nm recommended by the
    same.
    T: temperature in K
    =#
    retl = h2o2xdata[:,1]
    retx = 1e4*h2o2xdata[:,2] # factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    return reshape([retl; retx], length(retl), 2)
end

function hdo2xsect(hdo2xdata, T) 
    #=
    Currently this is a direct copy of h2o2xsect because no HDO2 crosssections
    exist. In the future this could be expanded if anyone ever makes those
    measurements.
    =#
    retl = hdo2xdata[:,1]
    retx = 1e4*hdo2xdata[:,2] # factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    reshape([retl; retx], length(retl), 2)
end

function ho2xsect_l(l) 
    #= 
    compute HO2 cross-section as a function of wavelength l in nm, as given by 
    Sander 2011 JPL Compilation 
    =#
    a = 4.91
    b = 30612.0
    sigmamed = 1.64e-18
    vmed = 50260.0
    v = 1e7/l;
    if 190<=l<=250
        return HO2absx = sigmamed / ( 1 - b/v ) * exp( -a * log( (v-b)/(vmed-b) )^2 )
    else
        return 0.0
    end
end

function o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, T)
    #=
    For O2, used in quantum yield calculations,
    including temperature-dependent Schumann-Runge bands.

    ...details missing
    =# 
    
    o2x = deepcopy(o2xdata);
    # fill in the schumann-runge bands according to Minschwaner 1992
    T = clamp(T, 130, 500)
    if 130<=T<190
        o2schr = o2schr130K
    elseif 190<=T<280
        o2schr = o2schr190K
    else
        o2schr = o2schr280K
    end

    del = ((T-100)/10)^2

    for i in [176.5:203.5;]
        posO2 = something(findfirst(isequal(i), o2x[:, 1]), 0)
        posschr = something(findfirst(isequal(i), o2schr[:, 1]), 0)
        o2x[posO2, 2] += 1e-20*(o2schr[posschr, 2]*del^2
                                + o2schr[posschr, 3]*del
                                + o2schr[posschr, 4])
    end

    # add in the herzberg continuum (though tiny)
    # measured by yoshino 1992
    for l in [192.5:244.5;]
        posO2 = something(findfirst(isequal(l), o2x[:, 1]), 0)
        o2x[posO2, 2] += 1e-24*(-2.3837947e4
                            +4.1973085e2*l
                            -2.7640139e0*l^2
                            +8.0723193e-3*l^3
                            -8.8255447e-6*l^4)
    end

    return o2x
end

function O3O1Dquantumyield(lambda, temp) 
    #=
    Quantum yield for O(1D) from O3 photodissociation. 
    "The quantum yield of O1D from ozone photolysis is actually well-studied! 
    This adds some complications for processing." - Mike
    =#
    if lambda < 306. || lambda > 328.
        return 0.
    end
    temp=clamp(temp, 200, 320)#expression is only valid in this T range

    X = [304.225, 314.957, 310.737];
    w = [5.576, 6.601, 2.187];
    A = [0.8036, 8.9061, 0.1192];
    v = [0.,825.518];
    c = 0.0765;
    R = 0.695;
    q = exp.(-v/(R*temp))
    qrat = q[1]/(q[1]+q[2])

    (q[1]/sum(q)*A[1]*exp.(-((X[1]-lambda)/w[1])^4.)
     +q[2]/sum(q)*A[2]*(temp/300.)^2 .* exp.(-((X[2]-lambda)/w[2])^2.)
     +A[3]*(temp/300.)^1.5*exp.(-((X[3]-lambda)/w[3])^2.)
     +c)
end

function padtosolar(solarflux, crosssection::Array{Float64, 2}) 
    #=
    a function to take an Nx2 array crosssection and pad it with zeroes until it's the
    same length as the solarflux array. Returns the cross sections only, as
    the wavelengths are shared by solarflux 
    =#
    positions = map(x->something(findfirst(isequal(x), solarflux[:,1]), 0), crosssection[:,1])
    retxsec = fill(0.,length(solarflux[:,1]))
    retxsec[positions] = crosssection[:,2]
    return retxsec
end

function populate_xsect_dict(pd_dataf; ion_xsects=true, globvars...)
    #=
    Creates a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    Jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, xsect_dict has already folded in quantum
    efficiency considerations (branching ratios).

    Input
        pd_dataf: Dictionary of filenames associated with a given species for its photodissociation reactions.
        Tn_array: Neutral temperatures for certain photolysis processes.
        ion_xsects: Whether to fill in crosssection information for ion species
        Jrates: the list of Jrates as defined in the PARAMETER file.
    Output
        crosssections: dictionary of cross sections by wavelength for each species. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:Tn, :n_all_layers, :Jratelist])
    

    # Set up =======================================================================
    xsect_dict = Dict{Symbol, Array{Array{Float64}}}()

    # Loading Data =================================================================
    # CO2 photodissociation --------------------------------------------------------
    # temperature-dependent between 195-295K
    co2xdata = readdlm(xsecfolder*pd_dataf[:CO2]["main"],'\t', Float64, comments=true, comment_char='#')

    # H2O & HDO --------------------------------------------------------------------
    h2oxdata = readdlm(xsecfolder*pd_dataf[:H2O]["main"],'\t', Float64, comments=true, comment_char='#')

    # These xsect_dicts for HDO are for 298K.
    hdoxdata = readdlm(xsecfolder*pd_dataf[:HDO]["main"],'\t', Float64, comments=true, comment_char='#')

    # H2O2 + HDO2 ------------------------------------------------------------------
    # the data in the following table cover the range 190-260nm
    h2o2xdata = readdlm(xsecfolder*pd_dataf[:H2O2]["main"],'\t', Float64, comments=true, comment_char='#')
    hdo2xdata = readdlm(xsecfolder*pd_dataf[:HDO2]["main"],'\t', Float64, comments=true, comment_char='#')

    # O3 ---------------------------------------------------------------------------
    # including IR bands which must be resampled from wavenumber
    o3xdata = readdlm(xsecfolder*pd_dataf[:O3]["main"],'\t', Float64, comments=true, comment_char='#')
    global o3ls = o3xdata[:,1]
    global o3xs = o3xdata[:,2]
    o3chapxdata = readdlm(xsecfolder*pd_dataf[:O3]["chapman"],'\t', Float64, comments=true, comment_char='#')
    o3chapxdata[:,1] = map(p->1e7/p, o3chapxdata[:,1])
    for i in [round(Int, floor(minimum(o3chapxdata[:,1]))):round(Int, ceil(maximum(o3chapxdata))-1);]
        posss = getpos(o3chapxdata, x->i<x<i+1)
        dl = diff([map(x->o3chapxdata[x[1],1], posss); i])
        x = map(x->o3chapxdata[x[1],2],posss)
        ax = reduce(+,map(*,x, dl))/reduce(+,dl)
        global o3ls = [o3ls; i+0.5]
        global o3xs = [o3xs; ax]
    end
    o3xdata = reshape([o3ls; o3xs],length(o3ls),2)
    

    # O2 ---------------------------------------------------------------------------
    o2xdata = readdlm(xsecfolder*pd_dataf[:O2]["main"],'\t', Float64, comments=true, comment_char='#')
    o2schr130K = readdlm(xsecfolder*pd_dataf[:O2]["schr_short"],'\t', Float64, comments=true, comment_char='#')
    o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
    o2schr130K = binupO2(o2schr130K)
    o2schr190K = readdlm(xsecfolder*pd_dataf[:O2]["schr_mid"],'\t', Float64, comments=true, comment_char='#')
    o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
    o2schr190K = binupO2(o2schr190K)
    o2schr280K = readdlm(xsecfolder*pd_dataf[:O2]["schr_long"],'\t', Float64, comments=true, comment_char='#')
    o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
    o2schr280K = binupO2(o2schr280K)

    # HO2 & DO2 --------------------------------------------------------------------
    ho2xsect = [190.5:249.5;]
    ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)
    do2xsect = deepcopy(ho2xsect)

    # H2 & HD ----------------------------------------------------------------------
    h2xdata = readdlm(xsecfolder*pd_dataf[:H2]["main"],',',Float64, comments=true, comment_char='#')
    hdxdata = readdlm(xsecfolder*pd_dataf[:HD]["main"],',',Float64, comments=true, comment_char='#')

    # OH & OD ----------------------------------------------------------------------
    ohxdata = readdlm(xsecfolder*pd_dataf[:OH]["main"],',',Float64, comments=true, comment_char='#')
    ohO1Dxdata = readdlm(xsecfolder*pd_dataf[:OH]["O1D+H"],',',Float64, comments=true, comment_char='#')
    odxdata = readdlm(xsecfolder*pd_dataf[:OD]["main"],',',Float64, comments=true, comment_char='#')


    # Populating the dictionary ======================================================
    #CO2+hv->CO+O
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))),
              map(t->co2xsect(co2xdata, t), GV.Tn)), get_Jrate_symb(["CO2", "CO", "O"], GV.Jratelist))
    #CO2+hv->CO+O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))),
              map(t->co2xsect(co2xdata, t), GV.Tn)), get_Jrate_symb(["CO2", "CO", "O1D"], GV.Jratelist))

    # O2 photodissociation ---------------------------------------------------------
    #O2+hv->O+O
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), GV.Tn)),
              get_Jrate_symb(["O2", "O", "O"], GV.Jratelist))
    #O2+hv->O+O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), GV.Tn)),
              get_Jrate_symb(["O2", "O", "O1D"], GV.Jratelist))

    # O3 photodissociation ---------------------------------------------------------
    # O3+hv->O2+O
    setindex!(xsect_dict,
              map(t->quantumyield(o3xdata,
                                  (
                                   (l->l<193, 1-(1.37e-2*193-2.16)),
                                   (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
                                   (l->225<=l<306, 0.1),
                                   (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
                                   (l->328<=l<340, 0.92),
                                   (l->340<=l, 1.0)
                                  )), GV.Tn), get_Jrate_symb(["O3", "O2", "O"], GV.Jratelist))
    # O3+hv->O2+O1D
    setindex!(xsect_dict,
              map(t->quantumyield(o3xdata,
                                  (
                                   (l->l<193, 1.37e-2*193-2.16),
                                   (l->193<=l<225, l->(1.37e-2*l-2.16)),
                                   (l->225<=l<306, 0.9),
                                   (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                                   (l->328<=l<340, 0.08),
                                   (l->340<=l, 0.0)
                                  )), GV.Tn), get_Jrate_symb(["O3", "O2", "O1D"], GV.Jratelist))
    # O3+hv->O+O+O
    setindex!(xsect_dict,
              fill(quantumyield(o3xdata,((x->true, 0.),)),GV.n_all_layers),
              get_Jrate_symb(["O3", "O", "O", "O"], GV.Jratelist))

    # H2 and HD photodissociation --------------------------------------------------
    # H2+hv->H+H
    setindex!(xsect_dict, fill(h2xdata, GV.n_all_layers), get_Jrate_symb(["H2", "H", "H"], GV.Jratelist))
    # HD+hν -> H+D 
    setindex!(xsect_dict, fill(hdxdata, GV.n_all_layers), get_Jrate_symb(["HD", "H", "D"], GV.Jratelist))

    # OH and OD photodissociation --------------------------------------------------
    # OH+hv->O+H
    setindex!(xsect_dict, fill(ohxdata, GV.n_all_layers), get_Jrate_symb(["OH", "O", "H"], GV.Jratelist))
    # OH + hv -> O(¹D) + H
    setindex!(xsect_dict, fill(ohO1Dxdata, GV.n_all_layers), get_Jrate_symb(["OH", "O1D", "H"], GV.Jratelist))
    # OD + hv -> O+D  
    setindex!(xsect_dict, fill(odxdata, GV.n_all_layers), get_Jrate_symb(["OD", "O", "D"], GV.Jratelist))
    # OD + hν -> O(¹D) + D 
    setindex!(xsect_dict, fill(ohO1Dxdata, GV.n_all_layers), get_Jrate_symb(["OD", "O1D", "D"], GV.Jratelist))

    # HO2 and DO2 photodissociation ------------------------------------------------
    # HO2 + hν -> OH + O
    setindex!(xsect_dict, fill(ho2xsect, GV.n_all_layers), get_Jrate_symb(["HO2", "OH", "O"], GV.Jratelist))
    # DO2 + hν -> OD + O
    setindex!(xsect_dict, fill(do2xsect, GV.n_all_layers), get_Jrate_symb(["DO2", "OD", "O"], GV.Jratelist))

    # H2O and HDO photodissociation ------------------------------------------------
    # H2O+hv->H+OH
    setindex!(xsect_dict,
              fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),GV.n_all_layers),
              get_Jrate_symb(["H2O", "H", "OH"], GV.Jratelist))

    # H2O+hv->H2+O1D
    setindex!(xsect_dict,
              fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),GV.n_all_layers),
              get_Jrate_symb(["H2O", "H2", "O1D"], GV.Jratelist))

    # H2O+hv->H+H+O
    setindex!(xsect_dict,
              fill(quantumyield(h2oxdata,((x->true, 0),)),GV.n_all_layers),
              get_Jrate_symb(["H2O", "H", "H", "O"], GV.Jratelist))

    # HDO + hν -> H + OD
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),GV.n_all_layers),
              get_Jrate_symb(["HDO", "H", "OD"], GV.Jratelist))

    # HDO + hν -> D + OH
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),GV.n_all_layers),
              get_Jrate_symb(["HDO", "D", "OH"], GV.Jratelist))

    # HDO + hν -> HD + O1D
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),GV.n_all_layers),
              get_Jrate_symb(["HDO", "HD", "O1D"], GV.Jratelist))

    # HDO + hν -> H + D + O
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->true, 0),)),GV.n_all_layers),
              get_Jrate_symb(["HDO", "H", "D", "O"], GV.Jratelist))


    # H2O2 and HDO2 photodissociation ----------------------------------------------
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
              map(t->h2o2xsect(h2o2xdata, t), GV.Tn)), get_Jrate_symb(["H2O2", "2OH"], GV.Jratelist))

    # H2O2+hv->HO2+H
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))),
              map(t->h2o2xsect(h2o2xdata, t), GV.Tn)), get_Jrate_symb(["H2O2", "HO2", "H"], GV.Jratelist))

    # H2O2+hv->H2O+O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(h2o2xdata, t),
              GV.Tn)), get_Jrate_symb(["H2O2", "H2O", "O1D"], GV.Jratelist))

    # HDO2 + hν -> OH + OD
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
              map(t->hdo2xsect(hdo2xdata, t), GV.Tn)), get_Jrate_symb(["HDO2", "OH", "OD"], GV.Jratelist))

    # HDO2 + hν-> DO2 + H
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
              map(t->hdo2xsect(hdo2xdata, t), GV.Tn)), get_Jrate_symb(["HDO2", "DO2", "H"], GV.Jratelist))

    # HDO2 + hν-> HO2 + D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
              map(t->hdo2xsect(hdo2xdata, t), GV.Tn)), get_Jrate_symb(["HDO2", "HO2", "D"], GV.Jratelist))

    # HDO2 + hν -> HDO + O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(hdo2xdata, t),
              GV.Tn)), get_Jrate_symb(["HDO2", "HDO", "O1D"], GV.Jratelist))

    if ion_xsects == true
        # NEW: CO2 photodissociation ---------------------------------------------------------
        # Source: Roger Yelle
        # CO₂ + hν -> C + O + O; JCO2toCpOpO
        thisjr = get_Jrate_symb(["CO2", "C", "O", "O"], GV.Jratelist)
        CO2_totaldiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_totaldiss_data, GV.n_all_layers), thisjr)

        # CO2 + hν -> C + O₂; JCO2toCpO2
        thisjr = get_Jrate_symb(["CO2", "C", "O2"], GV.Jratelist)
        CO2_diss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_diss_data, GV.n_all_layers), thisjr)

        # NEW: CO photodissociation ---------------------------------------------------------
        # Source: Roger Yelle

        # CO + hν -> C + O; JCOtoCpO
        thisjr = get_Jrate_symb(["CO", "C", "O"], GV.Jratelist)
        CO_diss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_diss_data, GV.n_all_layers), thisjr)


        # NEW: Nitrogen species photodissociation --------------------------------------------
        # Source: Roger Yelle

        # N₂ + hν -> N₂ + O(¹D); JN2OtoN2pO1D
        thisjr = get_Jrate_symb(["N2O", "N2", "O1D"], GV.Jratelist)
        N2_diss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2_diss_data, GV.n_all_layers), thisjr)

        # NO₂ + hν -> NO + O; JNO2toNOpO
        thisjr = get_Jrate_symb(["NO2", "NO", "O"], GV.Jratelist)
        NO2_diss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO2_diss_data, GV.n_all_layers), thisjr)

        # NO + hν -> N + O; JNOtoNpO
        thisjr = get_Jrate_symb(["NO", "N", "O"], GV.Jratelist)
        NO_diss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO_diss_data, GV.n_all_layers), thisjr)

        # Photoionization or ionizing dissociation reactions ============================================

        # NEW: CO₂ ionization -----------------------------------------------------------------
        # Source: Roger Yelle

        # CO₂ + hν -> CO₂⁺; JCO2toCO2pl
        thisjr = get_Jrate_symb(["CO2", "CO2pl"], GV.Jratelist)
        CO2_ionize_data = readdlm(xsecfolder*"$(thisjr).csv", ',', Float64, comments=true, comment_char='#')  # NOTE: replaced with Mike's file 19-Jan-2021.
        setindex!(xsect_dict, fill(CO2_ionize_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> CO₂²⁺; JCO2toCO2plpl (even though we don't track doubly ionized CO₂)
        thisjr = get_Jrate_symb(["CO2", "CO2plpl"], GV.Jratelist)
        CO2_doubleion_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_doubleion_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> C²⁺ + O₂; JCO2toCplplpO2
        thisjr = get_Jrate_symb(["CO2", "Cplpl", "O2"], GV.Jratelist)
        CO2_ionC2diss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionC2diss_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> C⁺ + O₂; JCO2toCplpO2
        thisjr = get_Jrate_symb(["CO2", "Cpl", "O2"], GV.Jratelist)
        CO2_ionCdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCdiss_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> CO⁺ + O⁺; JCO2toCOplpOpl
        thisjr = get_Jrate_symb(["CO2", "COpl", "Opl"], GV.Jratelist)
        CO2_ionCOandOdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCOandOdiss_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> CO⁺ + O; JCO2toCOplpO
        thisjr = get_Jrate_symb(["CO2", "COpl", "O"], GV.Jratelist)
        CO2_ionCOdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCOdiss_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> CO + O⁺; JCO2toOplpCO
        thisjr = get_Jrate_symb(["CO2", "Opl", "CO"], GV.Jratelist)
        CO2_ionOdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionOdiss_data, GV.n_all_layers), thisjr)

        # CO₂ + hν -> C⁺ + O⁺ + O; JCO2toOplpCplpO
        thisjr = get_Jrate_symb(["CO2", "Opl", "Cpl", "O"], GV.Jratelist)
        CO2_ionCandOdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCandOdiss_data, GV.n_all_layers), thisjr)

        # NEW: H2O ionization --------------------------------------------------------------
        # Source: Roger Yelle

        # TODO: Left off here 

        # H2O + hν -> H2O⁺; JH2OtoH2Opl
        H2Oionize_jr = get_Jrate_symb(["H2O", "H2Opl"], GV.Jratelist)
        h2o_ionize_data = readdlm(xsecfolder*"$(H2Oionize_jr).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionize_data, GV.n_all_layers), H2Oionize_jr)

        # HDO + hν -> HDO⁺; JHDOtoHDOpl # TODO: replace with HDO photoionization xsects when they exist
        thisjr = get_Jrate_symb(["HDO", "HDOpl"], GV.Jratelist)
        hdo_ionize_data = readdlm(xsecfolder*"$(H2Oionize_jr).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(hdo_ionize_data, GV.n_all_layers), thisjr)

        # H2O + hν -> O⁺ + H2; JH2OtoOplpH2
        thisjr = get_Jrate_symb(["H2O", "Opl", "H2"], GV.Jratelist)
        h2o_ionOdiss_data = readdlm(xsecfolder*"$(thisjr).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionOdiss_data, GV.n_all_layers), thisjr)

        # H2O + hν -> H⁺ + OH; JH2OtoHplpOH
        thisjr = get_Jrate_symb(["H2O", "Hpl", "OH"], GV.Jratelist)
        h2o_ionHdiss_data = readdlm(xsecfolder*"$(thisjr).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionHdiss_data, GV.n_all_layers), thisjr)

        # H2O + hν -> OH⁺ + H; JH2OtoOHplpH
        thisjr = get_Jrate_symb(["H2O", "OHpl", "H"], GV.Jratelist)
        h2o_ionOHdiss_data = readdlm(xsecfolder*"$(thisjr).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionOHdiss_data, GV.n_all_layers), thisjr)

        # NEW: CO ionization ----------------------------------------------------------------
        # Source: Roger Yelle

        # CO + hν -> CO⁺; JCOtoCOpl
        thisjr = get_Jrate_symb(["CO", "COpl"], GV.Jratelist)
        CO_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_ionize_data, GV.n_all_layers), thisjr)

        # CO + hν -> C + O⁺; JCOtoCpOpl
        thisjr = get_Jrate_symb(["CO", "C", "Opl"], GV.Jratelist)
        CO_ionOdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_ionOdiss_data, GV.n_all_layers), thisjr)

        # CO + hν -> C⁺ + O; JCOtoOpCpl
        thisjr = get_Jrate_symb(["CO", "O", "Cpl"], GV.Jratelist)
        CO_ionCdiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_ionCdiss_data, GV.n_all_layers), thisjr)

        # NEW: Nitrogen species ionization --------------------------------------------------
        # Source: Roger Yelle

        # N₂ + hν -> N₂⁺; JN2toN2pl
        thisjr = get_Jrate_symb(["N2", "N2pl"], GV.Jratelist)
        N2_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2_ionize_data, GV.n_all_layers), thisjr)

        # N₂ + hν -> N⁺ + N; JN2toNplpN
        thisjr = get_Jrate_symb(["N2", "Npl", "N"], GV.Jratelist)
        N2_iondiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2_iondiss_data, GV.n_all_layers), thisjr)

        # NO₂ + hν -> NO₂⁺; JNO2toNO2pl
        thisjr = get_Jrate_symb(["NO2", "NO2pl"], GV.Jratelist)
        NO2_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO2_ionize_data, GV.n_all_layers), thisjr)

        # NO + hν -> NO⁺; JNOtoNOpl
        thisjr = get_Jrate_symb(["NO", "NOpl"], GV.Jratelist)
        NO_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO_ionize_data, GV.n_all_layers), thisjr)

        # N₂O + hν -> N₂O⁺; JN2OtoN2Opl
        thisjr = get_Jrate_symb(["N2O", "N2Opl"], GV.Jratelist)
        N2O_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2O_ionize_data, GV.n_all_layers), thisjr)

        # NEW: Molecular and atomic hydrogen ionization -------------------------------------
        # Source: Roger Yelle

        # H + hν -> H⁺; JHtoHpl
        thisjr = get_Jrate_symb(["H", "Hpl"], GV.Jratelist)
        H_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H_ionize_data, GV.n_all_layers), thisjr)

        # H₂ + hν -> H₂⁺; JH2toH2pl
        thisjr = get_Jrate_symb(["H2", "H2pl"], GV.Jratelist)
        H2_ion_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H2_ion_data, GV.n_all_layers), thisjr)

        # HD + hν -> HD⁺; JH2toH2pl # TODO: Load HD crosssections when they exist
        thisjr = get_Jrate_symb(["HD", "HDpl"], GV.Jratelist)
        HD_ion_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(HD_ion_data, GV.n_all_layers), thisjr)

        # H₂ + hν -> H⁺ + H; JH2toHplpH
        thisjr = get_Jrate_symb(["H2", "Hpl", "H"], GV.Jratelist)
        H2_iondiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H2_iondiss_data, GV.n_all_layers), thisjr)

        # H₂O₂ + hν -> H₂O₂⁺; JH2O2toH2O2pl
        thisjr = get_Jrate_symb(["H2O2", "H2O2pl"], GV.Jratelist)
        H2O2_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H2O2_ionize_data, GV.n_all_layers), thisjr)

        # NEW: Oxygen and ozone ionization --------------------------------------------------
        # Source: Roger Yelle

        # O + hν -> O⁺; JOtoOpl
        thisjr = get_Jrate_symb(["O", "Opl"], GV.Jratelist)
        O_iondiss_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(O_iondiss_data, GV.n_all_layers), thisjr)

        # O₂ + hν -> O₂⁺; JO2toO2pl
        thisjr = get_Jrate_symb(["O2", "O2pl"], GV.Jratelist)
        O2_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(O2_ionize_data, GV.n_all_layers), thisjr)

        # # O₃ + hν -> O₃⁺; JO3toO3pl
        thisjr = get_Jrate_symb(["O3", "O3pl"], GV.Jratelist)
        O3_ionize_data = readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(O3_ionize_data, GV.n_all_layers), thisjr)
    end
    
    return xsect_dict
end

function quantumyield(xsect, arr)
    #= 
    function to assemble cross-sections for a given pathway. 

    xsect: an Nx2 array, format [wavelength, crosssection], where N is 
           the number of wavelengths available.

    arr: a tuple of tuples with a condition and a quantum yield multiplicative 
         factor, either constant or a function of wavelength in the given regime. 
         Return is an array with all of the matching wavelengths and the scaled cross-sections.
    =#
    lambdas = Float64[];
    rxs = Float64[];
    for (cond, qeff) in arr
        places = findall(cond, xsect[:,1])  # locate wavelengths that match the cond condition
        append!(lambdas, xsect[places, 1])  # put those wavelengths into a new list
        # if we have a number then map to a function
        isa(qeff, Function) ? (qefffn = qeff) : (qefffn = x->qeff)
        append!(rxs, map(*,map(qefffn, xsect[places, 1]),xsect[places, 2]))
    end

    return reshape([lambdas; rxs],length(lambdas),2)
end

# **************************************************************************** #
#                                                                              #
#                   Temperature & Water profile functions                      #
#                                                                              #
# **************************************************************************** #

function T_updated(z, Tsurf, Tmeso, Texo, sptype::String)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    z_meso_bottom = alt[searchsortednearest(alt, (Tmeso-Tsurf)/(lapserate))]
    z_meso_top = 110e5  # height of the tropopause top
    
    # These are the altitudes at which we "stitch" together the profiles 
    # from fitting the tanh profile in Ergun+2015,2021 to Hanley+2021 DD8
    # ion and electron profiles, and the somewhat arbitary profiles defined for
    # the region roughly between z_meso_top and the bottom of the fitted profiles.
    stitch_alt_electrons = 142e5
    stitch_alt_ions = 135e5

    function T_upper_atmo_neutrals(zee)
        return Texo - (Texo - Tmeso)*exp(-((zee - z_meso_top)^2)/(8e10*Texo))
    end
    
    function bobs_profile(z, TH, TL, z0, H0)
        #=
        Functional form from Ergun+ 2015 and 2021, but fit to data for electrons and ions
        in Hanley+ 2021, DD* data.
        =#
        return ((TH+TL)/2) + ((TH-TL)/2) * tanh(((z/1e5)-z0)/H0)
    end
    
    function T_ions_thermalize_region(z)
        #=
        This is a totally arbitary functional form for the region from [z_meso_top, 138]
        =#
        M = 6.06308414
        B = -538.97784139
        return M*(z/1e5) + B
    end
    
    # In the lower atmosphere, neutrals, ions, and electrons all 
    # have the same temperatures. 
    if z < z_meso_bottom
        return Tsurf + lapserate*z
    elseif z_meso_bottom <= z <= z_meso_top 
        return Tmeso
    
    # Near the top of the isothermal mesosphere, profiles diverge.        
    elseif z > z_meso_top
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="electron"
            # This region connects the upper atmosphere with the isothermal mesosphere
            if z_meso_top <= z < stitch_alt_electrons
                return bobs_profile(z, -1289.05806755, 469.31681082, 72.24740123, -50.84113252)
            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            elseif z >= stitch_alt_electrons
                return bobs_profile(z, 1409.23363494, 292.20319103, 191.39012079, 36.64138724)
            end
        elseif sptype=="ion"
            # This is similar to the electron handling, but for the ion profile in Hanley+2021, DD8.
            if z_meso_top < z <= stitch_alt_ions
                return T_ions_thermalize_region(z) < Tmeso ? Tmeso : T_ions_thermalize_region(z)
            elseif z > stitch_alt_ions
                return bobs_profile(z, 4.87796600e+06, 2.15643719e+02, 7.83610155e+02, 1.16129872e+02)
            end
        end
    end
end

function T_all(z, Tsurf, Tmeso, Texo, sptype::String)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in K
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of mesosphere. This code allows it to vary for surface or 
    # mesosphere experiments.
    ztropo_bot = (Tmeso-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    function find_offset(z1, T1, A, B, C)
        #=
        (z1, T1): the new point that satisfies T1 = T_upper_atmo_neutrals(z1), which we
                  also want to solve T_upper_atmo_ions.
        =#
        Q = A*(exp(-B*exp(-C*(z1/1e5))))
        
        D = T1 - Q
        return D
    end

    function T_upper_atmo_neutrals(zee)
        return Texo - (Texo - Tmeso)*exp(-((zee - ztropo)^2)/(8e10*Texo))
    end

    function T_upper_atmo_electrons(z1, T1)
        # for electrons, this is valid above 130 km in height.
        A = 2.24100983e+03
        B = 1.84024165e+01
        C = 1.44238590e-02
        D = find_offset(z1, T1, A, B, C)  
        return A*(exp(-B*exp(-C*z/1e5))) + D
    end

    function T_upper_atmo_ions(z1, T1)
        # valid above 160 km 
        A = 9.96741381e2
        B = 1.66317054e2
        C = 2.49434339e-2
        D = find_offset(z1, T1, A, B, C)  #
        return A*(exp(-B*exp(-C*(z/1e5 + 0)))) + D
    end

    # Used for correctly configuring the ion and electron profiles
    z1_ion = 160e5
    T1_ion = T_upper_atmo_neutrals(z1_ion)
    z1_e = 130e5
    T1_e = T_upper_atmo_neutrals(z1_e)

    # this is the WOOO--OOORRRST but someday I'll make it better MAYBE???
    if z >= 160e5 
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="ion"
            return T_upper_atmo_ions(z1_ion, T1_ion) # only valid up here
        elseif sptype=="electron"
            T_upper_atmo_electrons(z1_e, T1_e)
        end
    elseif 130e5 < z < 160e5
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="ion"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="electron"
            T_upper_atmo_electrons(z1_e, T1_e)
        end
    elseif ztropo < z <= 130e5   # upper atmosphere until neutrals, ions, electrons diverge
        return Texo - (Texo - Tmeso)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif (ztropo - ztropowidth) < z <= ztropo  # tropopause
        return Tmeso
    elseif z <= ztropo-ztropowidth  # lower atmosphere; <= makes it work for isothermal atm
        return Tsurf + lapserate*z
    end
end

function Tpiecewise(z, Tsurf, Tmeso, Texo)
    #= 
    FOR USE WITH FRACTIONATION FACTOR PROJECT ONLY.

    DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 "half-Gaussian" function for temperatures 
    altitudes above the tropopause, with a constant lapse rate (1.4K/km) 
    in the lower atmosphere. The tropopause width is allowed to vary
    in certain cases.

    z: altitude above surface in cm
    Tsurf: Surface temperature in K
    Tmeso: tropopause/mesosphere tempearture
    Texo: exobase temperature
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of mesosphere. This code allows it to vary for surface or 
    # mesosphere experiments.
    ztropo_bot = (Tmeso-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Tmeso)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Tmeso
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

# 1st term is a conversion factor to convert to (#/cm^3) from Pa. Source: Marti & Mauersberger 1993
Psat(T) = (1e-6/(kB_MKS * T))*(10^(-2663.5/T + 12.537))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# However, this function is defined on the offchance someone studies HDO.
Psat_HDO(T) = (1e-6/(kB_MKS * T))*(10^(-2663.5/T + 12.537))

end
