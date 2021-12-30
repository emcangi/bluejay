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


export # Basic utility functions
       charge_type, create_folder, deletefirst, find_nonfinites, fluxsymbol, format_chemistry_string, format_sec_or_min, getpos, input, nans_present, 
       next_in_loop, searchdir, search_subfolders, subtract_difflength,
       # Plotting functions
       get_colors, get_grad_colors, plot_atm, plot_bg, plot_extinction, plot_Jrates, plot_rxns, plot_temp_prof, plot_water_profile,                 
       # Reaction rate functions
       eval_rate_coef, get_column_rates, get_volume_rates, reactant_density_product, 
       # Atmosphere array manipulation                                                       
       atm_dict_to_matrix, atm_matrix_to_dict, column_density, flatten_atm, get_ncurrent, ncur_with_boundary_layers, n_tot, precip_microns, unflatten_atm, write_atmosphere,
       # Boundary condition functions                                                   
       boundaryconditions, effusion_velocity,
       # transport functions                                                                           
       Dcoef, Dcoef!, fluxcoefs, flux_param_arrays, flux_pos_and_neg, get_flux, Keddy, scaleH, update_diffusion_and_scaleH, update_transport_coefficients,                      
       # Chemistry functions
       calculate_stiffness,check_jacobian_eigenvalues, chemical_jacobian, getrate, loss_equations, loss_rate, make_chemjac_key, 
       make_net_change_expr, meanmass, production_equations, production_rate, rxns_where_species_is_observer, 
       # Photochemical equilibrium functions
       choose_solutions, construct_quadratic, group_terms, loss_coef, linear_in_species_density,
       # Photochemistry functions
       binupO2, co2xsect, h2o2xsect_l, h2o2xsect, hdo2xsect, ho2xsect_l, o2xsect, O3O1Dquantumyield, padtosolar, populate_xsect_dict, quantumyield, 
       # Temperature functions
       T_all, Tpiecewise,                                                                                                      
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
    return (10^y for y in range(log10(x1), log10(x2), length=n))
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

function column_density(n)
    #=
    Input
        n: species number density (#/cm³) by altitude
    Output
        Column density (#/cm²)
    =#
    return sum(n .* dz)
end

function flatten_atm(atmdict::Dict{Symbol, Vector{ftype_ncur}}, species_list, numlyrs::Int64) # params to pass
    #=
    Input:
        atmdict: atmospheric densities by altitude
        species_list: Included species which will have profiles flattened
    Output:
        Vector of form [n_sp1(z=0), n_sp2(z=0)...n_sp1(z=zmax)...n_spN(z=zmax)]
    
    This function is the reverse of unflatten_atm. 
    =#

    return deepcopy(ftype_ncur[[atmdict[sp][ialt] for sp in species_list, ialt in 1:numlyrs]...])
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

function meanmass(atmdict::Dict{Symbol, Vector{ftype_ncur}}, z, allsp)   # params to pass
    #= 
    find the mean molecular mass at a given altitude z

    atmdict: species number density by altitude
    z: float; altitude in atmosphere in cm

    return: mean molecular mass in amu
    =#
    thisaltindex = n_alt_index[z]
    c = [atmdict[sp][thisaltindex] for sp in allsp]
    m = [molmass[sp] for sp in allsp]
    return sum(c.*m)/sum(c)
end

function meanmass(atmdict::Dict{Symbol, Vector{ftype_ncur}}, allsp, mmass::Dict{Symbol, Int64}) # params to pass
    #= 
    Override for vector form. Calculates mean molecular mass at all atmospheric layers.

    atmdict: Array; species number density by altitude
    returns: mean molecular mass in amu for all atmospheric layers.
    =#

    # Gets the atmosphere as a matrix with rows = altitudes and cols = species
    # so we can do matrix multiplication.
    n_mat = transpose(atm_dict_to_matrix(atmdict, allsp))

    m = [mmass[sp] for sp in allsp] # this will always be 1D

    weighted_mm = zeros(size(n_mat)[1]) # This will store the result

    # Multiply densities of each species by appropriate molecular mass 
    mul!(weighted_mm, n_mat, m)

    return weighted_mm ./ n_tot(atmdict, allsp)
end

function ncur_with_boundary_layers(atmdict_no_bdys::Dict{Symbol, Vector{ftype_ncur}}, allsp) # params to pass
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
    
    # This gets a sorted list of the clamped indices, so it's [1, 1, 2, 3...end-1, end, end].
    clamped_n_alt_index = sort(collect(values(n_alt_index)))
    
    atmdict_with_bdy_layers = Dict{Symbol, Vector{ftype_ncur}}()
    
    # Fill the dictionary with the profile. This duplicates the lowest and highest altitude values.
    for i in 1:length(allsp)
        atmdict_with_bdy_layers[allsp[i]] = atmdict_no_bdys[allsp[i]][clamped_n_alt_index]
    end
    return atmdict_with_bdy_layers
end

function n_tot(atmdict::Dict{Symbol, Vector{ftype_ncur}}, z, allsp) # params to pass
    #= 
    Calculates total atmospheric density at altitude z.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        z: altitude, in cm
    Output: 
        Density of the atmosphere at altitude z
    =#
    thisaltindex = n_alt_index[z]
    return sum( [atmdict[s][thisaltindex] for s in allsp] )
end

function n_tot(atmdict::Dict{Symbol, Vector{ftype_ncur}}, allsp) # params to pass
    #= 
    Override to calculate total atmospheric density at all altitudes.

    Input: 
        atmdict: dictionary of atmospheric density profiles by altitude
        z: altitude, in cm
    Output: 
        Density of the atmosphere at altitude z

    This function is agnostic as to the size of the altitude grid--
    it collects it directly from atmdict.
    =#

    ndensities = zeros(length(allsp), length(atmdict[collect(keys(atmdict))[1]]))
    for i in 1:length(allsp)
        ndensities[i, :] = atmdict[allsp[i]]
    end

    # returns the sum over all species at each altitude as a vector.
    return vec(sum(ndensities, dims=1)) 
end

function precip_microns(sp, sp_profile, allsp, dz, alts, mmass::Dict{Symbol, Int64})
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
    col_abundance = column_density(sp_profile)
    cc_per_g = mmass[sp] / mmass[:H2O] # Water is 1 g/cm^3. Scale appropriately.

    #pr μm = (#/cm²) * (1 mol/molecules) * (g/1 mol) * (1 cm^3/g) * (10^4 μm/cm)
    pr_microns = col_abundance * (1/6.02e23) * (mmass[sp]/1) * (cc_per_g / 1) * (1e4/1)
    return pr_microns
end

function unflatten_atm(n_vec, species_list, numlyrs) # params to pass
    #=
    Input:
        n_vec: flattened density vector for the species in species_list: [n_sp1(z=0), n_sp2(z=0)...n_sp1(z=250)...n_spN(z=250)] 
    Output:
        dictionary of atmospheric densities by altitude with species as keys 

    This function is the reverse of flatten_atm.
    =#
    n_matrix = reshape(n_vec, (length(species_list), numlyrs))

    return atm_matrix_to_dict(n_matrix, species_list)
end

function write_atmosphere(atmdict::Dict{Symbol, Vector{ftype_ncur}}, filename::String, alts, numlyrs::Int64) # params to pass 
    #=
    Writes out the current atmospheric state to an .h5 file

    Input: 
        atmdict: atmospheric density profile dictionary
        filename: filename to write to
    Output: none
    =# 
    atm_mat = Array{Float64}(undef, numlyrs, length(collect(keys(atmdict))));
    for ispecies in [1:length(collect(keys(atmdict)));]
        for ialt in [1:numlyrs;]
            atm_mat[ialt, ispecies] = convert(Float64, atmdict[collect(keys(atmdict))[ispecies]][ialt])
        end
    end
    delete_old_h5file(filename)
    h5open(filename, "w") do f # this syntax is ok because we never write multiple times to a file.
        write(f, "n_current/n_current_mat", atm_mat)
        write(f, "n_current/alt", alts)
        write(f, "n_current/species", map(string, collect(keys(atmdict))))
    end
end

# **************************************************************************** #
#                                                                              #
#                             Plotting Functions                               #
#                                                                              #
# **************************************************************************** #

function get_colors(L::Int64, cmap)
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
    colors = [Cmap(x) for x in range(0, stop=1, length=L)]
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

function plot_atm(atmdict::Dict{Symbol, Vector{ftype_ncur}}, species_lists::Vector{Vector{Any}}, savepath::String, # params to pass
                  plotalts, scolor::Dict{Symbol, String}, sstyle::Dict{Symbol, String}, zmax, atol; 
                  t="", showonly=false, xlab="", xlim_1=(1e-12, 1e18), xlim_2=(1e-5, 1e5))
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
                                :N2,:N2pl,
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
    if length(species_lists)==2  # neutrals and ions 
        
        # set up the overall plot -------------------------------------------------------------
        atm_fig, atm_ax = subplots(3, 2, sharex=false, sharey=true, figsize=(14, 16))
        subplots_adjust(wspace=0, hspace=0)
        tight_layout()
                
        # only the neutral-col axes
        atm_ax[1, 1].set_title("Neutrals")
        for i in 1:3
            plot_bg(atm_ax[i, 1])
            atm_ax[i, 1].set_xlim(xlim_1[1], xlim_1[2])
            atm_ax[i, 1].fill_betweenx(plotalts, xlim_1[1] .* ones(size(plotalts)), x2=atol, alpha=0.1, color=medgray, zorder=10)
            atm_ax[i, 1].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end

        x_axis_label = xlab == "" ? L"Species concentration (cm$^{-3}$)" : xlab
        atm_ax[3, 1].set_xlabel(x_axis_label)
        
        # only the ion-col axes
        atm_ax[1, 2].set_title("Ions")
        for i in 1:3
            plot_bg(atm_ax[i, 2])
            atm_ax[i, 2].set_xlim(xlim_2[1], xlim_2[2])
            atm_ax[i, 2].fill_betweenx(plotalts, xlim_2[1] .* ones(size(plotalts)), x2=atol, alpha=0.1, color=medgray, zorder=10)
            atm_ax[i, 2].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end
        atm_ax[3, 2].set_xlabel(x_axis_label)
        
        # y axes labels
        for j in 1:3
            atm_ax[j, 1].set_ylabel("Altitude [km]")
        end
        
        # plot the neutrals according to logical groups -------------------------------------------------------
        for sp in species_lists[1]
            atm_ax[axes_by_sp[sp], 1].plot(convert(Array{Float64}, atmdict[sp]), plotalts, color=get(scolor, sp, "black"),
                                           linewidth=2, label=sp, linestyle=get(sstyle, sp, "-"), zorder=10)
        end
        
        # plot the ions according to logical groups ------------------------------------------------------------
        for sp in species_lists[2]
            atm_ax[axes_by_sp[sp], 2].plot(convert(Array{Float64}, atmdict[sp]), plotalts, color=get(scolor, sp, "black"),
                                           linewidth=2, label=sp, linestyle=get(sstyle, sp, "-"), zorder=10)
        end

        # stuff that applies to all axes
        for a in atm_ax
            a.set_ylim(0, zmax/1e5)
            a.set_xscale("log")
            handles, labels = a.get_legend_handles_labels()
            if isempty(handles) == false
                a.legend(handles, labels, fontsize=8)#bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)
            end
        end

    # Plot only neutrals - to support the fractionation factor project ==========================================
    elseif length(species_lists)==1
        atm_fig, atm_ax = subplots(figsize=(16,6))
        tight_layout()
        for sp in species_lists[1]
            atm_ax.plot(convert(Array{Float64}, atmdict[sp]), plotalts, color=get(scolor, sp, "black"),
                        linewidth=2, label=sp, linestyle=get(sstyle, sp, "-"), zorder=1)
            atm_ax.set_xlim(xlim_1[1], xlim_1[2])
            atm_ax.set_ylabel("Altitude [km]")
            # atm_ax.set_title("Neutrals")
        end
        atm_ax.tick_params(which="both", labeltop=true, top=true)
        plot_bg(atm_ax)
        atm_ax.set_ylim(0, zmax/1e5)
        atm_ax.set_xscale("log")
        atm_ax.set_xlabel(x_axis_label)
        atm_ax.legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0, fontsize=16)
    else
        throw("unexpected number of species lists to plot")
    end

    suptitle(t, y=1.05)

    if showonly==false  
        atm_fig.savefig(savepath, bbox_inches="tight", dpi=300)
        # show()
    else
        show()
    end
end

function plot_bg(axob)
    #=
    Make plots not look ugly. 

    intended to immitate seaborn
    =#
    axob.set_facecolor("#ededed")
    axob.grid(zorder=-5, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function plot_extinction(solabs, path::String, dt, iter::Int64, zmax=zmax; # Not currently called. params to pass
                         tauonly=false, xsect_info=nothing, solflux=nothing) 
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
    
    fig, ax = subplots()

    num_x = length(solabs[1])  # extinction is in an awkward format
    num_y = size(solabs)[1]
    X = 1:num_x  # start with 0 because we are sending it into a python plotting library, not julia 
    Y = 0:2:(zmax/1e5 - 2)

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
            # println(maximum(xsect))
            z = z .* xsect
            titlestr = L"e^{-\tau}\sigma_{\lambda}" * " for $(xsect_sp)"
            savestr = "extinction_$(dt)_$(iter)_$(xsect_sp).png"
        end

        if solflux != nothing
            solflux = reshape(solflux, (1, length(solflux)))
            # println(maximum(solflux))
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
    yticks(0:50:(zmax/1e5))

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(25))

    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(25))

    title(titlestr)
    ax.set_xlim(0, 300)  # turn off to see full range to 2000 nm.
    savefig(path*savestr, bbox_inches="tight")
    close(fig)
end

function plot_Jrates(sp, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, bcdict::Dict{Symbol, Matrix{Any}}, # Not currently called. params to pass
                     allsp, ionsp, rxnnet, numlyrs::Int64, plotalts, results_dir::String; 
                     filenameext="")
    #=
    Plots the Jrates for each photodissociation or photoionizaiton reaction. Override for small groups of species.

    species: species for which to plot the Jrates

    TODO: Needs to be worked on especially get_volume_rates works as expected
    =#

    # --------------------------------------------------------------------------------
    # calculate reaction rates x density of the species at each level of the atmosphere.
    
    rxd_prod = get_volume_rates(sp, atmdict, Tn, Ti, Te, bcdict, allsp, ionsp, rxnnet, numlyrs, which="Jrates")
    rxd_loss = get_volume_rates(sp, atmdict, Tn, Ti, Te, bcdict, allsp, ionsp, rxnnet, numlyrs, which="Jrates")

    # ---------------------------------------------------------------------------------
    # Plot reaction rates and transport rates by altitude
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

   
    fig, ax = subplots(figsize=(8,6))
    plot_bg(ax)

    minx = 1e10
    maxx = 0
    
    # Collect chem production equations and total 
    for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
        lbl = "$(kv[1])"
        ax.semilogx(kv[2], plotalts, linestyle="-", linewidth=1, label=lbl)
        
        # set the xlimits
        if minimum(kv[2]) <= minx
            minx = minimum(kv[2])
        end
        if maximum(kv[2]) >= maxx
            maxx = maximum(kv[2])
        end
    end

    suptitle("J rates", fontsize=20)
    ax.set_ylabel("Altitude (km)")
    ax.legend(bbox_to_anchor=(1.01, 1))
    ax.set_xlabel("Chemical reaction rate [check units] ("*L"cm^{-3}s^{-1})")
    # labels and such 
    if minx < 1e-20
        minx = 1e-20
    end
    maxx = 10^(ceil(log10(maxx)))
    ax.set_xlim([minx, maxx])
    savefig(results_dir*"J_rates_$(filenameext).png", bbox_inches="tight", dpi=300)
    show()
end

function plot_rxns(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, Tp, # params to pass
                   bcdict::Dict{Symbol, Matrix{Any}}, rxnnet, 
                   allsp, ionsp, transsp, chemsp, mmass::Dict{Symbol, Int64}, alts, n_alt_index::Dict, dz, 
                   polar::Dict{Symbol, Float64}, numlyrs::Int64, plotalts, 
                   results_dir::String, T_for_Hs::Dict{String, Vector{Any}}, T_for_diff::Dict{String, Vector{Any}}; 
                   shown_rxns=nothing, subfolder="", plotsfolder="", dt=nothing, num="", extra_title="", 
                   plot_timescales=false, plot_total_rate_coefs=false, showonly=false)
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

    # Arrays to store the total reactions per second for this species of interest
    total_prod_rate = zeros(numlyrs)
    total_loss_rate = zeros(numlyrs)

    # Arrays to hold the total chemical production and loss 
    total_chem_prod = zeros(numlyrs)
    total_chem_loss = zeros(numlyrs)
    total_chem_prod_ratecoef = zeros(numlyrs)
    total_chem_loss_ratecoef = zeros(numlyrs)

    if sp in chemsp
        # --------------------------------------------------------------------------------
        # calculate reaction rates x density for all reactions of the species at each level of the atmosphere.
        # Temperature arrays include boundary layers so that diffusion can be calculated later, so we have to 
        # index these in this way to make the evaluation of chemistry reaction rate coefficients work. 
        rxd_prod, rate_coefs_prod = get_volume_rates(sp, atmdict, Tn[2:end-1], Ti[2:end-1], Te[2:end-1], bcdict, allsp, ionsp, rxnnet, numlyrs, species_role="product")
        rxd_loss, rate_coefs_loss = get_volume_rates(sp, atmdict, Tn[2:end-1], Ti[2:end-1], Te[2:end-1], bcdict, allsp, ionsp, rxnnet, numlyrs, species_role="reactant")

        # Water is turned off in the lower atmosphere, so we should represent that.
        if sp in [:H2O, :HDO] 
            for (prod_k, loss_k) in zip(keys(rxd_prod), keys(rxd_loss))
                rxd_prod[prod_k][1:upper_lower_bdy_i] .= NaN
                rxd_loss[loss_k][1:upper_lower_bdy_i] .= NaN
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
                    ax[1].semilogx(kv[2], plotalts, linestyle=ls[ls_i], marker=9, markevery=20, color=cols[col_i], linewidth=1, label=kv[1])
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
                    ax[1].semilogx(kv[2], plotalts, linestyle=ls[ls_i], marker=8, markevery=20, color=cols[col_i], linewidth=1, label=kv[1])
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
        ax[1].semilogx(total_chem_prod, plotalts, color="xkcd:forest green", linestyle=(0, (4,2)), marker=9, markevery=20, linewidth=2, label="Total chemical production", zorder=5)
        ax[1].semilogx(total_chem_loss, plotalts, color="xkcd:shamrock", linestyle=(0, (4,2)), marker=8, markevery=20, linewidth=2, label="Total chemical loss", zorder=5)

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
            ax_1_2.semilogx(total_chem_prod_ratecoef, plotalts, color="xkcd:royal purple", linestyle=":", linewidth=3, label="Total chemical production rate coef", zorder=6)
            ax_1_2.semilogx(total_chem_loss_ratecoef, plotalts, color="xkcd:lavender", linestyle=":", linewidth=3, label="Total chemical loss rate coef", zorder=6)
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
    if sp in transsp
        transportPL = get_transport_PandL_rate(sp, atmdict, Tn, Ti, Te, Tp, bcdict, allsp, setdiff(allsp, ionsp), transsp, 
                                                              mmass, alts, n_alt_index, polar, dz, T_for_Hs, T_for_diff)
        # now separate into two different arrays for ease of addition.
        production_i = transportPL .>= 0  # boolean array for where transport entries > 0 (production),
        loss_i = transportPL .< 0 # and for where transport entries < 0 (loss).
        total_transport_prod = production_i .* transportPL
        total_transport_loss = loss_i .* abs.(transportPL)

        if sp in [:H2O, :HDO] # Water is turned off in the lower atmosphere, so we should represent that.
            total_transport_prod[1:upper_lower_bdy_i] .= NaN
            total_transport_loss[1:upper_lower_bdy_i] .= NaN
        end

        # set the x lims for transport axis. Special because total_transport_prod, and etc are incomplete arrays.
        minx[2] = 10^(floor(log10(minimum(abs.(transportPL)))))
        maxx[2] = 10^(ceil(log10(maximum(abs.(transportPL)))))

        # Plot the transport production and loss without the boundary layers
        ax[2].scatter(total_transport_prod, plotalts, color="red", marker=9, label="Total gain this layer", zorder=4)
        ax[2].scatter(total_transport_loss, plotalts, color="blue", marker=8, label="Total loss this layer", zorder=4)
        ax[2].set_xscale("log")

        ax[2].legend(fontsize=12)

        plottitle_ext = " by chemistry & transport"
    else
        total_transport_prod = zeros(numlyrs)
        total_transport_loss = zeros(numlyrs)
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

    ax[3].semilogx(total_prod_rate, plotalts, color="black", marker=9, markevery=15, linewidth=2, label="Total production", zorder=3) #linestyle=(2, (1,2)), 
    ax[3].semilogx(total_loss_rate, plotalts, color="gray", marker=8, markevery=15, linewidth=2, label="Total loss", zorder=3) #linestyle=(2, (1,2)),

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
        ax[total_ax].plot(timescale_of_change, plotalts, color=c, zorder=5)
    end

    # Final plotting tasks ============================================================
    # check for and correct any ridiculously low limits
    for i in 1:length(minx)
        if minx[i] < 1e-12
            minx[i] = 1e-12
        end
    end

    titles = ["Chemistry", "Transport", "Total chem+trans"]
    notes = ["Chemical production/loss is off for $(sp) below $(Int64(upper_lower_bdy /1e5)) km.", 
             "Vertical transport is off for $(sp) below $(Int64(upper_lower_bdy /1e5)) km.", 
             "Abundance of $(sp) below $(Int64(upper_lower_bdy /1e5)) km is fixed."] 
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

function plot_temp_prof(n_temps, alts; # params to pass
                        i_temps=nothing, e_temps=nothing, savepath=nothing, showonly=false) 
    #=
    Creates a .png image of the tepmeratures plotted by altitude in the atmosphere

    Inputs:
        n_temps: an array of neutral temperature by altitude
        ion_temps: same, but for the ions
        e_temps: same but for the electrons
        savepath: where to save the resulting .png image
        showonly: whether to just show() the figure. If false, must be accompanied by a value for savepath.
    =#

    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    fig, ax = subplots(figsize=(4,6))
    plot_bg(ax)

    plot(n_temps, alts./1e5, label="Neutrals", color=medgray)

    if i_temps != nothing
        ax.plot(i_temps, alts./1e5, label="Ions", color="xkcd:bright orange")
        ax.legend(fontsize=16)
        ax.set_xscale("log")
    end
    if e_temps != nothing
        ax.plot(e_temps, alts./1e5, label="Electrons", color="cornflowerblue")
        ax.legend(fontsize=16, loc="center right")
        ax.set_xscale("log")
    end

    # plot the control temps
    ax.scatter(n_temps[1], 0, marker="o", color=medgray, zorder=10)
    ax.text(n_temps[1], 0, L"\mathrm{T}_{\mathrm{surface}}"*" = $(Int64(round(n_temps[1], digits=0))) K ")

    middle_ind = Int64(length(n_temps)/2)
    ax.scatter(n_temps[middle_ind-1], 75, marker="o", color=medgray, zorder=10)
    ax.text(n_temps[middle_ind]+5, 75, L"\mathrm{T}_{\mathrm{meso}}"*" = $(Int64(round(n_temps[middle_ind-1], digits=0))) K ")

    ax.scatter(n_temps[end], 250, marker="o", color=medgray, zorder=10)
    ax.text(n_temps[end]*1.05, 240, L"\mathrm{T}_{\mathrm{exo}}"*" = $(Int64(round(n_temps[end], digits=0))) K ")
    
    ax.set_ylabel("Altitude [km]")
    ax.set_yticks(collect(0:50:alts[end]/1e5))
    ax.set_yticklabels(collect(0:50:alts[end]/1e5))
    ax.set_xlabel("Temperature [K]")
    ax.set_xlim(100, 2e3)

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

function plot_water_profile(H2Oinitf, HDOinitf, nH2O, nHDO, savepath::String, # params to pass
                            plotalts; 
                            showonly=false, watersat=nothing) 
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

    fig = figure(figsize=(4,6))
    ax = gca()
    plot_bg(ax)
    ax.tick_params(axis="x", which="minor", bottom=true, top=true)

    
    ax1col = "#88527F"
    
    # mixing ratio in PPM axis
    ax.semilogx(convert(Array{Float64}, H2Oinitf)/1e-6, plotalts, color=ax1col, linewidth=2)
    ax.semilogx(convert(Array{Float64}, HDOinitf)/1e-6, plotalts, color=ax1col, linestyle="--", linewidth=2)
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
    ax2.semilogx(convert(Array{Float64}, nH2O), plotalts, color=ax2col, linewidth=2, label=L"H$_2$O")
    ax2.semilogx(convert(Array{Float64}, nHDO), plotalts, color=ax2col, linestyle="--", linewidth=2, label="HDO")
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
        semilogx(convert(Array{Float64}, H2Oinitf), plotalts, color=ax1col, linewidth=3, label=L"H$_2$O initial fraction")
        semilogx(convert(Array{Float64}, watersat[2:end-1]), plotalts, color="black", alpha=0.5, linewidth=3, label=L"H$_2$O saturation")
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
function boundaryconditions(fluxcoef_dict::Dict, bcdict::Dict{Symbol, Matrix{Any}}, allsp, notranssp, dz) # params to pass
    #= 
    Inputs:
        fluxcoef_dict: a dictionary containing the K and D flux coefficients for every species throughout
                       the atmosphere. Format species=>Array(length(all_species), length(alt)).
        bcdict: Boundary conditions dictionary specified in parameters file
    Outputs:
        boundary conditions for species in the format:
        [n_1 -> n_0, n_0 -> n_1;      n_(num_layers) -> n_(num_layers+1), n_(num_layers+1) -> n_(num_layers)]
        where n_0 is the boundary layer from [-1, 1], n_1 is the first atmospheric layer describing [1 km, 3 km],
        n_(num_layers) is the topmost atmospheric layer, and n_(num_layers+1) is the top boundary layer.
    =#
    
    # This is where we will store all the boundary condition numbers in the format that the code
    # understands.
    bc_dict = Dict{Symbol, Array{ftype_ncur}}([s=>[0 0; 0 0] for s in allsp])

    for s in allsp
        # retrieves user-supplied boundary conditions. When unsupplied, falls back to flux = 0 at both boundaries.
        bcs = get(bcdict, s, ["f" 0.; "f" 0.]) 
        # Flux is also 0 at both boundaries for species disallowed from transport.
        if issubset([s], notranssp)
            bcs = ["f" 0.; "f" 0.]
        end

        # first element returned corresponds to lower BC, second to upper
        # BC transport rate. Within each element, the two rates correspond
        # to the two equations
        # n_b  -> NULL (first rate, depends on species concentration)
        # NULL -> n_b  (second rate, independent of species concentration)
        # note these are chemical equations not indicating atmospheric layers!
        bcvec = ftype_ncur[0 0; 0 0]
        # the signs on these are the opposite of what you expect because it is necessary to 
        # encode loss in production equations and vice versa.
        # in first element, negative means you're adding something, and vice versa.
        # in second element, positive means you're adding something, and vice versa.

        # Note: Because of the above, all the comments below talking about the layers are probably worng.
        # LOWER
        if bcs[1, 1] == "n"
            # Density conditions at the surface.
            # n_1 -> n_0 is just the combined transport coefficient (D+K)/Δz² + ⌷/(2Δz), directed downward (the [1]).
            # n_0 -> n_1 condition is the upward transport coefficient, but multiplied by the boundary condition since [-1, 1] is not a real cell.
            bcvec[1,:]=[fluxcoef_dict[s][2, :][1], # indices: s = species, alt = 2 km, 1 = into the boundary layer
                        fluxcoef_dict[s][1, :][2]*bcs[1,2]] # indices: species, alt = 0 km, 2 = out of bdy layer "lower up"
                         # Units: (#/s * #/cm³ = #/cm³/s)
        elseif bcs[1, 1] == "f"
            # No condition for n_1 -> n_0 flux. Specified condition for n_0 -> n_1 flux.
            bcvec[1,:] = [0.0, -bcs[1, 2]/dz]  # division by dz converts to #/cm^3. NOTE: this negative sign added 11 Feb 2021, shouldn't affect anything as of now.
        elseif bcs[1, 1] == "v"
            # there is a condition for n_1 -> n_0 ("depositional velocity"), but no condition for n_0 -> n_1 velocity.
            bcvec[1,:] = [bcs[1, 2]/dz, 0.0] # division by dz converts velocity to #/s
        else
            throw("Improper lower boundary condition!")
        end

        # UPPER
        if bcs[2, 1] == "n"
            # Similar conditions for density here, but for the top of the atmosphere.
            bcvec[2,:] = [fluxcoef_dict[s][end-1, :][2],
                          fluxcoef_dict[s][end, :][1]*bcs[2,2]]
        elseif bcs[2, 1] == "f"  
            # n_(top+1) -> n_top flux is constrained to this value (negative because of the way transport is encoded as chemistry.)
            # allows terms not proportional to density - requires another reaction.
            bcvec[2,:] = [0.0,-bcs[2, 2]/dz] # this represents removal from the atmosphere, while going upwards. sign not directional. 
        elseif bcs[2, 1] == "v"
            # Sets a velocity condition for particles leaving n_top to the exosphere, but no condition for particles coming in from the exosphere.
            bcvec[2,:] = [bcs[2, 2]/dz, 0.0] 
        else
            throw("Improper upper boundary condition!")
        end
        
        bc_dict[s] .= bcvec
    end

    # remember each element is 2x2 array where the first row is for the lower boundary, and 2nd for upper boundary.
    return bc_dict
end

function effusion_velocity(Texo, m, zmax) # params to pass
    #=
    Returns effusion velocity for a species in cm/s

    Inputs:
        Texo: temperature of the exobase (upper boundary) in K
        m: mass of one molecule of species in amu
        zmax: max altitude in cm
    Outputs:
        v: effusion velocity for species of mass m 
    =#
    
    # lambda is the Jeans parameter (Gronoff 2020), basically the ratio of the 
    # escape velocity GmM/z to the thermal energy, kT.
    lambda = (m*mH*bigG*marsM)/(kB*Texo*(radiusM+zmax))
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

function Dcoef!(D_arr, T_arr, sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, bcdict::Dict{Symbol, Matrix{Any}}, # params to pass
                allsp, mmass::Dict{Symbol, Int64}, polarizability::Dict{Symbol, Float64}, neutralsp) 
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
   
    # Calculate as if it was a neutral
    D_arr[:] .= (diffparams(sp)[1] .* 1e17 .* T_arr .^ (diffparams(sp)[2])) ./ n_tot(atmdict, allsp)

    # a place to store the density array
    species_density = zeros(size(T_arr))
    
    # If an ion, overwrite with the ambipolar diffusion
    if charge_type(sp) == "ion"
        sum_nu_in = zeros(size(T_arr))

        mi = mmass[sp] .* mH
        # create the sum of nu_in. Note that this depends on density, but we only have density for the real layers,
        # so we have to assume the density at the boundary layers is the same as at the real layers.
        for n in neutralsp
            species_density = atmdict[n]

            # This sets the species density to a boundary condition if it exists. 
            bccheck = get(bcdict, n, ["f" 0.; "f" 0.])
            if bccheck[1,1] == "n"
                species_density[1] = bccheck[1,2]
            end
            if bccheck[2,1] == "n"  # currently this should never apply.
                species_density[end] = bccheck[2,2]
            end

            mu_in = (1 ./ mi .+ 1 ./ (mmass[n] .* mH)) .^ (-1) # reduced mass in g
            sum_nu_in .+= 2 .* pi .* (((polarizability[n] .* q .^ 2) ./ mu_in) .^ 0.5) .* species_density

        end
        
        D_arr .= (kB .* T_arr) ./ (mi .* sum_nu_in)
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

function fluxcoefs(sp::Symbol, z, dz, Kv, Dv, Tv_n, Tv_p, # Checked - global free
                   Hsv, H0v) 
    #= 
    base function to generate flux coefficients of the transport network. 
    
    For all the arrays, length = num_layers 

    Inputs:
        sp: species symbol 
        z: altitude in cm.
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

    # Initialize arrays for downward (i to i-1) and upward (i to i+1) coefficients
    Dl = zeros(length(z))
    Kl = zeros(length(z))
    Tl_n = zeros(length(z))
    Tl_p = zeros(length(z))
    dTdzl_n = zeros(length(z))
    dTdzl_p = zeros(length(z))
    Hsl = zeros(length(z))
    H0l = zeros(length(z))

    Du = zeros(length(z))
    Ku = zeros(length(z))
    Tu_n = zeros(length(z))
    Tu_p = zeros(length(z))
    dTdzu_n = zeros(length(z))
    dTdzu_p = zeros(length(z))
    Hsu = zeros(length(z))
    H0u = zeros(length(z))

    # Calculate the coefficients between this layer and the lower layer. 
    Dl[2:end] = @. (Dv[sp][1:end-1] + Dv[sp][2:end]) /  2.0
    Kl[2:end] = @. (Kv[1:end-1] + Kv[2:end]) / 2.0
    Tl_n[2:end] = @. (Tv_n[1:end-1] + Tv_n[2:end]) / 2.0
    Tl_p[2:end] = @. (Tv_p[1:end-1] + Tv_p[2:end]) / 2.0
    dTdzl_n[2:end] = @. (Tv_n[2:end] - Tv_n[1:end-1]) / dz
    dTdzl_p[2:end] = @. (Tv_p[2:end] - Tv_p[1:end-1]) / dz
    Hsl[2:end] = @. (Hsv[sp][1:end-1] + Hsv[sp][2:end]) / 2.0
    H0l[2:end] = @. (H0v[charge_type(sp)][1:end-1] + H0v[charge_type(sp)][2:end]) / 2.0

    # Handle the lower boundary layer:
    Dl[1] = @. (1 + Dv[sp][1]) /  2.0
    Kl[1] = @. (1 + Kv[1]) / 2.0
    Tl_n[1] = @. (1 + Tv_n[1]) / 2.0
    Tl_p[1] = @. (1 + Tv_p[1]) / 2.0
    dTdzl_n[1] = @. (Tv_n[1] - 1) / dz
    dTdzl_p[1] = @. (Tv_p[1] - 1) / dz
    Hsl[1] = @. (1 + Hsv[sp][1]) / 2.0
    H0l[1] = @. (1 + H0v[charge_type(sp)][1]) / 2.0

    # Now the coefficients between this layer and upper layer
    Du[1:end-1] = @. (Dv[sp][1:end-1] + Dv[sp][2:end]) /  2.0
    Ku[1:end-1] = @. (Kv[1:end-1] + Kv[2:end]) / 2.0
    Tu_n[1:end-1] = @. (Tv_n[1:end-1] + Tv_n[2:end]) / 2.0
    Tu_p[1:end-1] = @. (Tv_p[1:end-1] + Tv_p[2:end]) / 2.0
    dTdzu_n[1:end-1] = @. (Tv_n[2:end] - Tv_n[1:end-1]) / dz
    dTdzu_p[1:end-1] = @. (Tv_p[2:end] - Tv_p[1:end-1]) / dz
    Hsu[1:end-1] = @. (Hsv[sp][1:end-1] + Hsv[sp][2:end]) / 2.0
    H0u[1:end-1] = @. (H0v[charge_type(sp)][1:end-1] + H0v[charge_type(sp)][2:end]) / 2.0

    # Handle upper boundary layer:
    Du[end] = @. (Dv[sp][end] + 1) /  2.0
    Ku[end] = @. (Kv[end] + 1) / 2.0
    Tu_n[end] = @. (Tv_n[end] + 1) / 2.0
    Tu_p[end] = @. (Tv_p[end] + 1) / 2.0
    dTdzu_n[end] = @. (1 - Tv_n[end]) / dz
    dTdzu_p[end] = @. (1 - Tv_p[end]) / dz
    Hsu[end] = @. (Hsv[sp][end] + 1) / 2.0
    H0u[end] = @. (H0v[charge_type(sp)][end] + 1) / 2.0


    # two flux terms: eddy diffusion and gravity/thermal diffusion.
    # these are found in line 5 of Mike's transport_as_chemistry.pdf:
    # sumeddy = (D+K)/(Δz²), gravthermal = ☐/(2Δz), where ☐ = {D(1/H + 1+(α/T)(dT/dz)) + K(1/H_H + (1/T)(dT/dz))}
    sumeddyl = @. (Dl+Kl)/dz/dz
    if charge_type(sp) == "neutral"
        gravthermall = @. (Dl*((1/Hsl) + ((1+thermaldiff(sp))/Tl_n)*dTdzl_n) +
                        Kl*((1/H0l) + (1/Tl_n)*dTdzl_n))/(2*dz)
    elseif charge_type(sp) == "ion"
        gravthermall = @. (Dl*((1/Hsl) + ((1+thermaldiff(sp))/Tl_p)*dTdzl_p) +
                        Kl*((1/H0l) + (1/Tl_n)*dTdzl_n))/(2*dz)
    elseif charge_type(sp) == "electron"
        throw("Electrons not handled as individual species")
    end

    sumeddyu = @. (Du+Ku)/dz/dz  # this is the line where we divide by cm^2
    if charge_type(sp) == "neutral"
        gravthermalu = @. (Du*((1/Hsu) + ((1 + thermaldiff(sp))/Tu_n)*dTdzu_n) +
                        Ku*((1/H0u) + (1/Tu_n)*dTdzu_n))/(2*dz)
    elseif charge_type(sp) == "ion"
        gravthermalu = @. (Du*((1/Hsu) + ((1 + thermaldiff(sp))/Tu_p)*dTdzu_p) +
                        Ku*((1/H0u) + (1/Tu_n)*dTdzu_n))/(2*dz)
    elseif charge_type(sp) == "electron"
        throw("Electrons not handled as individual species")
    end
    
    # this results in the following coupling coefficients; sumeddy + gravthermal = (D+K)/(Δz²) + ☐/(2Δz), units 1/s <-----_!!!!! important
    # first return is this term between layer i and i-1 for whole atmosphere.
    # second return is between layer i and i+1
    return sumeddyl .+ gravthermall,  # down
            sumeddyu .- gravthermalu # up; negative because gravity points down. I think that's why.
end

function fluxcoefs(species_list::Vector{Any}, T_neutral, T_plasma, K, D, # params to pass
                   H0, Hs, alts, dz) 
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
                transport species during th emain simulation run, and for all species when trying to plot 
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
    
    # the return dictionary: Each species has 2 entries for every layer of the atmosphere.
    fluxcoef_dict = Dict{Symbol, Array{ftype_ncur}}([s=>fill(0., length(alts), 2) for s in species_list])

    for s in species_list
        layer_below_coefs, layer_above_coefs = fluxcoefs(s, alts, dz, K, D, T_neutral, T_plasma, Hs, H0) 
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

function get_flux(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, Tp,  # params to pass
                  bcdict::Dict{Symbol, Matrix{Any}}, 
                  allsp, neutralsp, transsp, mmass::Dict{Symbol,Int64}, alts, n_alt_index::Dict, 
                  polar::Dict{Symbol,Float64}, numlyrs::Int64, dz, T_for_Hs::Dict{String,Vector{Any}}, T_for_diff::Dict{String,Vector{Any}}) 
    #=

    TODO: Having problems calling this? You moved sp to the first argument, that's probably why. 

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
    
    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(allsp, atmdict, Tn, Tp, D_arr, bcdict, allsp, neutralsp, 
                                                                 mmass, alts, n_alt_index, polar, T_for_diff)
    Hs_dict = Dict{Symbol, Vector{ftype_ncur}}([sp=>scaleH(alt, sp, T_for_Hs[charge_type(sp)], molmass) for sp in allsp])  

    fluxcoefs_all = fluxcoefs(allsp, Tn, Tp, Keddy_arr, Dcoef_dict, H0_dict, Hs_dict, alts, dz)
    bc_dict = boundaryconditions(fluxcoefs_all, bcdict, allsp, setdiff(allsp, transsp), dz)

    # each element in bulk_layer_coefs has the format [downward flow (i to i-1), upward flow (i to i+1)].  units 1/s
    bulk_layer_coefs = fluxcoefs_all[sp][2:end-1]

    bcs = bc_dict[sp]
    
    net_bulk_flow = fill(convert(ftype_ncur, NaN), length(alt)-1)  # units #/cm^3/s; tracks the cell boundaries, of which there are length(alt)-1

    # We will calculate the net flux across each boundary, with sign indicating direction of travel.
    # First boundary is the boundary between the surface layer and the first atmospheric layer (alt = 1 km)
    net_bulk_flow[1] = (bcs[1, 2]  # increase of the lowest atmospheric layer's density. Will always be 0 unless the species has a density or flux condition
                       - atmdict[sp][1]*bcs[1, 1]) # lowest atmospheric layer --> surface ("depositional" term)
                        
    for ialt in 2:numlyrs  # now iterate through every cell boundary within the atmosphere. boundaries at 3 km, 5...247. 123 elements.
        net_bulk_flow[ialt] = (atmdict[sp][ialt-1]*bulk_layer_coefs[ialt-1][2]   # coming up from below: cell i-1 to cell i. Should be positive * positive
                              - atmdict[sp][ialt]*bulk_layer_coefs[ialt][1])     # leaving to the layer below: downwards: cell i to cell i-1
    end

    # now the top boundary - between 124th atmospheric cell (alt = 249 km)
    net_bulk_flow[end] = (atmdict[sp][end]*bcs[2, 1] # into exosphere from the cell
                         - bcs[2, 2]) # into top layer from exosphere. do not question this
                
    return net_bulk_flow .* dz # now it is a flux. hurrah.
end

function get_transport_PandL_rate(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, Tp, bcdict::Dict{Symbol, Matrix{Any}}, # params to pass
                                  allsp, neutralsp, transsp, mmass::Dict{Symbol,Int64}, alts, n_alt_index::Dict, polar::Dict{Symbol,Float64}, dz, 
                                  T_for_Hs, T_for_diff::Dict{String,Vector{Any}}) 
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
    
    # Generate the fluxcoefs dictionary and boundary conditions dictionary
    D_arr = zeros(size(Tn))
    Keddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(allsp, atmdict, Tn, Tp, D_arr, bcdict, allsp, neutralsp, 
                                                                 mmass, alts, n_alt_index, polar, T_for_diff)

    Hs_dict = Dict{Symbol, Vector{ftype_ncur}}([sp=>scaleH(alt, sp, T_for_Hs[charge_type(sp)], molmass) for sp in allsp]) 


    fluxcoefs_all = fluxcoefs(allsp, Tn, Tp, Keddy_arr, Dcoef_dict, H0_dict, Hs_dict, alts, dz)

    # For the bulk layers only to make the loops below more comprehendable: 
    fluxcoefs_bulk_layers = Dict([s=>fluxcoefs_all[s][2:end-1, :] for s in keys(fluxcoefs_all)])

    bc_dict = boundaryconditions(fluxcoefs_all, bcdict, allsp, setdiff(allsp, transsp), dz)  # here we are passing in no_transport_species

    # each element in thesebcs has the format [downward, upward]
    thesebcs = bc_dict[sp]

    # Fill array 
    transport_PL = fill(convert(ftype_ncur, NaN), length(alts)-2)

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

# Species-specific scale height
function scaleH(z::Float64, sp::Symbol, controltemps::Array, mmass) # params to pass
    #=
    Input:
        z: ONE altitude in cm
        sp: species to calculate for
        controltemps: temperatures at surface, mesosphere, exobase
    Output:
        species-specific scale height 
    =#  

    T = T_all(z, controltemps[1], controltemps[2], controltemps[3], charge_type(sp))
    mm = mmass[sp]
    return kB*T/(mm*mH*marsM*bigG)*(((z+radiusM))^2)
end

# VECTORIZED scale heights
function scaleH(z::Vector{Float64}, sp::Symbol, T::Array, mmass) # params to pass
    #=
    Input:
        z: Altitudes in cm
        sp: Speciecs to calculate for
        T: temperature array for this species
    Output: 
        species-specific scale height at all altitudess
    =#  

    mm = mmass[sp]
    return @. kB*T/(mm*mH*marsM*bigG)*(((z+radiusM))^2)
end

function scaleH(z::Vector{Float64}, atmdict::Dict{Symbol, Vector{ftype_ncur}}, T::Vector{Any}, allsp::Array, mmass) # params to pass
    #= 
    Input:
        z: altitudes in cm
        atmdict: Present atmospheric state dictionary
        T: temperature array for the neutral atmosphere
    Output:
        Mean atmospheric scale height at all altitudes
    =#
    mm_vec = meanmass(atmdict, allsp, mmass) # vector version.

    return @. kB*T/(mm_vec*mH*marsM*bigG)*(((z+radiusM))^2)
end

# thermal diffusion factors
thermaldiff(sp) = get(Dict(:H=>-0.25, :H2=>-0.25, :D=>-0.25, :HD=>-0.25,
                                :He=>-0.25, 
                                :Hpl=>-0.25, :H2pl=>-0.25, :Dpl=>-0.25, :HDpl=>-0.25,
                                :Hepl=>-0.25), sp, 0)

function update_diffusion_and_scaleH(species_list, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Tp, # params to pass
                                     D_coefs, bcdict, 
                                     allsp, neutralsp, mmass, alts, n_alt_index, polar, T_for_diff) 
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

    ncur_with_bdys = ncur_with_boundary_layers(atmdict, allsp)
    
    K = Keddy(alts, n_tot(ncur_with_bdys, allsp))
    H0_dict = Dict{String, Vector{ftype_ncur}}("neutral"=>scaleH(alts, ncur_with_bdys, Tn, allsp, mmass),
                                            "ion"=>scaleH(alts, ncur_with_bdys, Tp, allsp, mmass))
    
    # Molecular diffusion is only needed for transport species, though.  
    Dcoef_dict = Dict{Symbol, Vector{ftype_ncur}}([s=>deepcopy(Dcoef!(D_coefs, T_for_diff[charge_type(s)], s, ncur_with_bdys, bcdict, 
                                                                   allsp, mmass, polar, neutralsp)) for s in species_list])

    return K, H0_dict, Dcoef_dict
end

function update_transport_coefficients(species_list, atmdict::Dict{Symbol, Vector{ftype_ncur}}, activellsp, Tn, Tp, # params to pass
                                       D_coefs, Hs, bcdict, 
                                       allsp, neutralsp, transsp, mmass, n_alt_index, polar, alts, numlyrs, dz, T_for_diff) 
    #=
    Input:
        species_list: Species which will have transport coefficients updated
        atmdict: Atmospheric state dictionary for bulk layers
        activellsp: active long-lived species symbol list
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
    
    # Update the diffusion coefficients and scale heights
    K_eddy_arr, H0_dict, Dcoef_dict = update_diffusion_and_scaleH(species_list, atmdict, Tn, Tp, D_coefs, bcdict, allsp, neutralsp, 
                                                                  mmass, alts, n_alt_index, polar, T_for_diff)
    
    # Get flux coefficients
    fluxcoefs_all = fluxcoefs(species_list, Tn, Tp, K_eddy_arr, Dcoef_dict, H0_dict, Hs, alts, dz)
    
    # Transport coefficients, non-boundary layers
    tup = fill(-999., length(activellsp), numlyrs)
    tdown = fill(-999., length(activellsp), numlyrs)
    for (i, s) in enumerate(activellsp)
        tup[i, :] .= fluxcoefs_all[s][2:end-1, 2]
        tdown[i, :] .= fluxcoefs_all[s][2:end-1, 1]
    end

    # transport coefficients for boundary layers
    bc_dict = boundaryconditions(fluxcoefs_all, bcdict, allsp, setdiff(allsp, transsp), dz)
    tlower = permutedims(reduce(hcat, [bc_dict[sp][1,:] for sp in activellsp]))
    tupper = permutedims(reduce(hcat, [bc_dict[sp][2,:] for sp in activellsp]))

    return tlower, tup, tdown, tupper
end

# **************************************************************************** #
#                                                                              #
#                             Chemistry Functions                              #
#                                                                              #
# **************************************************************************** #

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

function chemical_jacobian(chemnet, transportnet, specieslist, dspecieslist, chemsp, transportsp; chem_on=true, trans_on=true)
    #= 
    Compute the symbolic chemical jacobian of a supplied chemnet and transportnet
    for the specified specieslist. 

    Input:
        chemnet: chemical reactions
        transportnet: transport equations
        specieslist: list of species to calculate for; i.e. the species whose equation is differentiated
        dspecieslist: species with respect to which the equations are differentiated
        chemsp: list of chemistry species
        transportsp: list of transport species
    Output:
        three arrays suitable for constructing a sparse matrix: 
            i: row indices
            j: column indices
            v: values to place at (i, j)
    =#

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)  # this is the active species. 
    ndspecies = length(dspecieslist)  # this is the species with respect to which we differentiate

    for i in 1:nspecies # for each species
        ispecies = specieslist[i]
        # get the production and loss equations
        peqn = []
        leqn = []
        if issubset([ispecies], chemsp)
            peqn = [peqn; production_equations(ispecies, chemnet)] 
            leqn = [leqn; loss_equations(ispecies, chemnet)]
        end
        if issubset([ispecies],transportsp)
            peqn = [peqn; production_equations(ispecies, transportnet)]
            leqn = [leqn; loss_equations(ispecies, transportnet)]
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
    return (ivec, jvec, tvec) # Expr(:vcat, tvec...))  # 
end

function eval_rate_coef(atmdict::Dict{Symbol, Vector{ftype_ncur}}, krate::Expr, tn, ti, te, # params to pass 
                        allsp, ionsp, numlyrs::Int64) 
    #=
    Evaluates a chemical reaction rate coefficient, krate, for all levels of the atmosphere. 

    Input:
        atmdict: the atmospheric state dictionary
        krate: rate coefficient for a single reaction
        tn, _i, _e: temperature profiles for neutrals, ions, electrons
    Output:
        rate_coefficient: evaluated rate coefficient at all atmospheric layers
    =#

    # Set stuff up
    rate_coefficient = zeros(numlyrs)
    M_by_alt = sum([atmdict[sp] for sp in allsp]) 
    E_by_alt = sum([atmdict[sp] for sp in ionsp])
    @eval ratefunc(Tn, Ti, Te, M, E) = $krate
    rate_coefficient .= Base.invokelatest(ratefunc, tn, ti, te, M_by_alt, E_by_alt)

    return rate_coefficient
end 

function get_column_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, bcdict::Dict{Symbol, Matrix{Any}}, # params to pass
                          allsp, ionsp, dz, rxnnet, numlyrs::Int64; 
                          which="all", sp2=nothing, role="product", startalt_i=1) 
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
    
    rxd, coefs = get_volume_rates(sp, atmdict, Tn[2:end-1], Ti[2:end-1], Te[2:end-1], bcdict, allsp, ionsp, rxnnet, numlyrs, 
                                   species_role=role, which=which);
    
    # Make the column rates dictionary for production
    columnrate = Dict()
    for k in keys(rxd)
        columnrate[k] = sum(rxd[k][startalt_i:end] .* dz)
    end
    
    # Optionally one can specify a second species to include in the sorted result, i.e. a species' ion.
    if sp2 != nothing
        rxd2, coefs2 = get_volume_rates(sp2, atmdict, Tn[2:end-1], Ti[2:end-1], Te[2:end-1], bcdict, allsp, ionsp, rxnnet, numlyrs, 
                                         species_role=role, which=which);

        columnrate2 = Dict()

        for k in keys(rxd2)
            columnrate2[k] = sum(rxd2[k][startalt_i:end] .* dz)
        end

        colrate_dict = merge(columnrate, columnrate2)
    else
        colrate_dict = columnrate
    end
    
    sorted = sort(collect(colrate_dict), by=x->x[2], rev=true)
    
    return sorted
end

function get_volume_rates(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, Tn, Ti, Te, bcdict::Dict{Symbol, Matrix{Any}}, 
                           allsp, ionsp, rxnnet, numlyrs::Int64; 
                           species_role="both", which="all") 
    #=
    Input:
        sp: Species name
        atmdict: Present atmospheric state dictionary
        Tn, Ti, Te: temperature arrays
        bcdict: Boundary conditions dictionary specified in parameters file
        species_role: whether to look for the species as a reactant, product, or both.  If it has a value, so must species.
        which: "all", "Jrates", "krates". Whether to fill the dictionary with all reactions, only photochemistry/photoionization 
               (Jrates) or only chemistry (krates).
    Output: 
        rxn_dat: Evaluated rates, i.e. k[A][B], units #/cm^3/s for bimolecular rxns
        rate_coefs: Evaluated rate coefficients for each reaction 
    =#

    # Fill in the rate x density dictionary ------------------------------------------------------------------------------
    rxn_dat =  Dict{String, Array{ftype_ncur, 1}}()
    rate_coefs = Dict{String, Array{ftype_ncur, 1}}()

    # Select either photodissociation or bi-/tri-molecular reactions
    if which=="Jrates"
        selected_rxns = filter(x->occursin("J", string(x[3])), rxnnet)
    elseif which=="krates" 
        selected_rxns = filter(x->!occursin("J", string(x[3])), rxnnet)
    elseif which=="all"
        selected_rxns = deepcopy(rxnnet)
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
        reactants = rxn[1]
        products = rxn[2]  # vector of the product symbols

        # get the reactants and products in string form for use in plot labels
        rxn_str = format_chemistry_string(rxn[1], rxn[2])

        # Fill in rate coefficient * species density for all reactions
        if typeof(rxn[3]) == Symbol # for photodissociation
            rxn_dat[rxn_str] = atmdict[rxn[1][1]] .* atmdict[rxn[3]]
            rate_coefs[rxn_str] = atmdict[rxn[3]]
        else                        # bi- and ter-molecular chemistry
            density_prod = reactant_density_product(atmdict, rxn[1], allsp, ionsp, numlyrs)
            rate_coef = eval_rate_coef(atmdict, rxn[3], Tn, Ti, Te, allsp, ionsp, numlyrs)

            rxn_dat[rxn_str] = density_prod .* rate_coef
            rate_coefs[rxn_str] = rate_coef
        end
    end

    return rxn_dat, rate_coefs
end

function getrate(sp::Symbol, chemnet, transportnet, chemsp, transportsp; chemistry_on=true, transport_on=true, sepvecs=false) 
    #=
    Creates a symbolic expression for the rate at which a given species is
    either produced or lost due to chemical reactions or transport.

    Input:
        sp: species for which to get the rate 
        chemnet: chemistry reaction array
        transportnet: transport network array
        chemsp: species with active chemistry
        transportsp: species which transport
        chemistry_on: set to false to disallow chemical changes to species
        transport_on: set to false to disallow transport of a species
        sepvecs: Allows this function to return a vector of expressions for chemical production and loss 
                 and transport production and loss
    Output: either
        rate: the final value of dn/dt for sp from all processes, or
        chemprod_rate, chemloss_rate, transprod_rate, transloss_rate: dn/dt for sp due to these processes, calculated separately.
    =#

    # This block will return the total net change for the species, P - L.
    if sepvecs == false
        rate = :(0.0)
        if issubset([sp], chemsp) && chemistry_on
            rate = :($rate 
                     + $(production_rate(sp, chemnet, sepvecs=sepvecs)) 
                     - $(loss_rate(sp, chemnet, sepvecs=sepvecs)) 
                    )
        end
        if issubset([sp],transportsp) && transport_on
            rate = :($rate 
                     + $(production_rate(sp, transportnet, sepvecs=sepvecs)) 
                     - $(loss_rate(sp, transportnet, sepvecs=sepvecs))
                    )
        end
        return rate
    else  # if we want a vector of expressions for each production and loss (4 terms, 2 each for chemistry and transport)
        if issubset([sp],chemsp) && chemistry_on
            chemprod_rate = production_rate(sp, chemnet, sepvecs=sepvecs)
            chemloss_rate = loss_rate(sp, chemnet, sepvecs=sepvecs)
        else
            chemprod_rate = [:(0.0 + 0.0)]  # Doing it this way because it's the easiest way to make a vector of one expression that's just 0
            chemloss_rate = [:(0.0 + 0.0)]
        end
        
        if issubset([sp],transportsp) && transport_on
            transprod_rate = production_rate(sp, transportnet, sepvecs=sepvecs)
            transloss_rate = loss_rate(sp, transportnet, sepvecs=sepvecs)
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

function reactant_density_product(atmdict::Dict{Symbol, Vector{ftype_ncur}}, reactants, 
                                  allsp, ionsp, numlyrs::Int64)
    #=
    Calculates the product of all reactant densities for a chemical reaction for the whole atmosphere, 
    i.e. for A + B --> C + D, return n_A * n_B.

    Input:
        atmdict: the atmospheric state dictionary
        reactants: a list of reactant symbols.
    Output: 
        density_product: returns n_A * n_B for all altitudes for the reaction A + B --> ...
    =#

    density_product = ones(numlyrs)
    for r in reactants
        if r != :M && r != :E
            # species densities by altitude
            density_product .*= atmdict[r]  # multiply by each reactant density
        elseif r == :M
            density_product .*= sum([atmdict[sp] for sp in allsp]) 
        elseif r == :E
            density_product .*= sum([atmdict[sp] for sp in ionsp])
        else
            throw("Got an unknown symbol in a reaction rate: $(r)")
        end
    end

    return density_product 
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

function loss_coef(leqn_vec, sp::Symbol; calc_tau_chem=false)  # Checked - global-free
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

function populate_xsect_dict(controltemps::Array, alts; ion_xsects=true) # TODO: Put xsect file names into a large list argument. 
    #=
    Creates a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    Jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, xsect_dict has already folded in quantum
    efficiency considerations (branching ratios).

    ion_xsects: Whether to fill in crosssection information for ion species
    =#

    # Set up =======================================================================
    Temp_n(z) = T_all(z, controltemps[1], controltemps[2], controltemps[3], "neutral")
    xsect_dict = Dict{Symbol, Array{Array{Float64}}}()

    # Loading Data =================================================================
    # CO2 photodissociation --------------------------------------------------------
    # temperature-dependent between 195-295K
    co2xdata = readdlm(xsecfolder*co2file,'\t', Float64, comments=true, comment_char='#')

    # CO2 photoionization (used to screen high energy sunlight)
    # co2exdata = readdlm(xsecfolder*co2exfile,',',Float64, comments=true, comment_char='#')

    # H2O & HDO --------------------------------------------------------------------
    h2oxdata = readdlm(xsecfolder*h2ofile,'\t', Float64, comments=true, comment_char='#')

    # These xsect_dicts for HDO are for 298K.
    hdoxdata = readdlm(xsecfolder*hdofile,'\t', Float64, comments=true, comment_char='#')

    # H2O2 + HDO2 ------------------------------------------------------------------
    # the data in the following table cover the range 190-260nm
    h2o2xdata = readdlm(xsecfolder*h2o2file,'\t', Float64, comments=true, comment_char='#')
    hdo2xdata = readdlm(xsecfolder*hdo2file,'\t', Float64, comments=true, comment_char='#')

    # O3 ---------------------------------------------------------------------------
    # including IR bands which must be resampled from wavenumber
    o3xdata = readdlm(xsecfolder*o3file,'\t', Float64, comments=true, comment_char='#')
    global o3ls = o3xdata[:,1]
    global o3xs = o3xdata[:,2]
    o3chapxdata = readdlm(xsecfolder*o3chapfile,'\t', Float64, comments=true, comment_char='#')
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
    o2xdata = readdlm(xsecfolder*o2file,'\t', Float64, comments=true, comment_char='#')
    o2schr130K = readdlm(xsecfolder*o2_130_190,'\t', Float64, comments=true, comment_char='#')
    o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
    o2schr130K = binupO2(o2schr130K)
    o2schr190K = readdlm(xsecfolder*o2_190_280,'\t', Float64, comments=true, comment_char='#')
    o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
    o2schr190K = binupO2(o2schr190K)
    o2schr280K = readdlm(xsecfolder*o2_280_500,'\t', Float64, comments=true, comment_char='#')
    o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
    o2schr280K = binupO2(o2schr280K)

    # HO2 & DO2 --------------------------------------------------------------------
    ho2xsect = [190.5:249.5;]
    ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)
    do2xsect = deepcopy(ho2xsect)

    # H2 & HD ----------------------------------------------------------------------
    h2xdata = readdlm(xsecfolder*h2file,',',Float64, comments=true, comment_char='#')
    hdxdata = readdlm(xsecfolder*hdfile,',',Float64, comments=true, comment_char='#')

    # OH & OD ----------------------------------------------------------------------
    ohxdata = readdlm(xsecfolder*ohfile,',',Float64, comments=true, comment_char='#')
    ohO1Dxdata = readdlm(xsecfolder*oho1dfile,',',Float64, comments=true, comment_char='#')
    odxdata = readdlm(xsecfolder*odfile,',',Float64, comments=true, comment_char='#')


    # Populating the dictionary ======================================================
    #CO2+hv->CO+O
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))),
              map(t->co2xsect(co2xdata, t), map(Temp_n, alts))), :JCO2toCOpO)
    #CO2+hv->CO+O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))),
              map(t->co2xsect(co2xdata, t), map(Temp_n, alts))), :JCO2toCOpO1D)

    # O2 photodissociation ---------------------------------------------------------
    #O2+hv->O+O
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), map(Temp_n, alts))),
              :JO2toOpO)
    #O2+hv->O+O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), map(Temp_n, alts))),
              :JO2toOpO1D)

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
                                  )), map(Temp_n, alts)), :JO3toO2pO)
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
                                  )), map(Temp_n, alts)), :JO3toO2pO1D)
    # O3+hv->O+O+O
    setindex!(xsect_dict,
              fill(quantumyield(o3xdata,((x->true, 0.),)),length(alts)),
              :JO3toOpOpO)

    # H2 and HD photodissociation --------------------------------------------------
    # H2+hv->H+H
    setindex!(xsect_dict, fill(h2xdata, length(alts)), :JH2toHpH)
    # HD+hν -> H+D 
    setindex!(xsect_dict, fill(hdxdata, length(alts)), :JHDtoHpD)

    # OH and OD photodissociation --------------------------------------------------
    # OH+hv->O+H
    setindex!(xsect_dict, fill(ohxdata, length(alts)), :JOHtoOpH)
    # OH+hv->O1D+H
    setindex!(xsect_dict, fill(ohO1Dxdata, length(alts)), :JOHtoO1DpH)
    # OD + hv -> O+D  
    setindex!(xsect_dict, fill(odxdata, length(alts)), :JODtoOpD)
    # OD + hν -> O(¹D) + D 
    setindex!(xsect_dict, fill(ohO1Dxdata, length(alts)), :JODtoO1DpD)

    # HO2 and DO2 photodissociation ------------------------------------------------
    # HO2 + hν -> OH + O
    setindex!(xsect_dict, fill(ho2xsect, length(alts)), :JHO2toOHpO)
    # DO2 + hν -> OD + O
    setindex!(xsect_dict, fill(do2xsect, length(alts)), :JDO2toODpO)

    # H2O and HDO photodissociation ------------------------------------------------
    # H2O+hv->H+OH
    setindex!(xsect_dict,
              fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),length(alts)),
              :JH2OtoHpOH)

    # H2O+hv->H2+O1D
    setindex!(xsect_dict,
              fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alts)),
              :JH2OtoH2pO1D)

    # H2O+hv->H+H+O
    setindex!(xsect_dict,
              fill(quantumyield(h2oxdata,((x->true, 0),)),length(alts)),
              :JH2OtoHpHpO)

    # HDO + hν -> H + OD
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alts)),
              :JHDOtoHpOD)

    # HDO + hν -> D + OH
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alts)),
              :JHDOtoDpOH)

    # HDO + hν -> HD + O1D
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alts)),
              :JHDOtoHDpO1D)

    # HDO + hν -> H + D + O
    setindex!(xsect_dict,
              fill(quantumyield(hdoxdata,((x->true, 0),)),length(alts)),
              :JHDOtoHpDpO)


    # H2O2 and HDO2 photodissociation ----------------------------------------------
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
              map(t->h2o2xsect(h2o2xdata, t), map(Temp_n, alts))), :JH2O2to2OH)

    # H2O2+hv->HO2+H
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))),
              map(t->h2o2xsect(h2o2xdata, t), map(Temp_n, alts))), :JH2O2toHO2pH)

    # H2O2+hv->H2O+O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(h2o2xdata, t),
              map(Temp_n, alts))), :JH2O2toH2OpO1D)

    # HDO2 + hν -> OH + OD
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
              map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alts))), :JHDO2toOHpOD)

    # HDO2 + hν-> DO2 + H
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
              map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alts))), :JHDO2toDO2pH)

    # HDO2 + hν-> HO2 + D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
              map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alts))), :JHDO2toHO2pD)

    # HDO2 + hν -> HDO + O1D
    setindex!(xsect_dict,
              map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(hdo2xdata, t),
              map(Temp_n, alts))), :JHDO2toHDOpO1D)

    if ion_xsects == true
        # NEW: CO2 photodissociation ---------------------------------------------------------
        # Source: Roger Yelle
        # CO₂ + hν -> C + O + O
        CO2_totaldiss_data = readdlm(xsecfolder*"$(:JCO2toCpOpO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_totaldiss_data, length(alts)), :JCO2toCpOpO)

        # CO2 + hν -> C + O₂
        CO2_diss_data = readdlm(xsecfolder*"$(:JCO2toCpO2).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_diss_data, length(alts)), :JCO2toCpO2)

        # NEW: CO photodissociation ---------------------------------------------------------
        # Source: Roger Yelle

        # CO + hν -> C + O
        CO_diss_data = readdlm(xsecfolder*"$(:JCOtoCpO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_diss_data, length(alts)), :JCOtoCpO)


        # NEW: Nitrogen species photodissociation --------------------------------------------
        # Source: Roger Yelle

        # N₂ + hν -> N₂ + O(¹D)
        N2_diss_data = readdlm(xsecfolder*"$(:JN2OtoN2pO1D).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2_diss_data, length(alts)), :JN2OtoN2pO1D)

        # NO₂ + hν -> NO + O
        NO2_diss_data = readdlm(xsecfolder*"$(:JNO2toNOpO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO2_diss_data, length(alts)), :JNO2toNOpO)

        # NO + hν -> N + O
        NO_diss_data = readdlm(xsecfolder*"$(:JNOtoNpO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO_diss_data, length(alts)), :JNOtoNpO)

        # Photoionization or ionizing dissociation reactions ============================================

        # NEW: CO₂ ionization -----------------------------------------------------------------
        # Source: Roger Yelle

        # CO₂ + hν -> CO₂⁺
        CO2_ionize_data = readdlm(xsecfolder*"$(:JCO2toCO2pl).csv", ',', Float64, comments=true, comment_char='#')  # NOTE: replaced with Mike's file 19-Jan-2021.
        setindex!(xsect_dict, fill(CO2_ionize_data, length(alts)), :JCO2toCO2pl)

        # CO₂ + hν -> CO₂²⁺  (even though we don't track doubly ionized CO₂)
        CO2_doubleion_data = readdlm(xsecfolder*"$(:JCO2toCO2plpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_doubleion_data, length(alts)), :JCO2toCO2plpl)

        # CO₂ + hν -> C²⁺ + O₂
        CO2_ionC2diss_data = readdlm(xsecfolder*"$(:JCO2toCplplpO2).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionC2diss_data, length(alts)), :JCO2toCplplpO2)

        # CO₂ + hν -> C⁺ + O₂
        CO2_ionCdiss_data = readdlm(xsecfolder*"$(:JCO2toCplpO2).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCdiss_data, length(alts)), :JCO2toCplpO2)

        # CO₂ + hν -> CO⁺ + O⁺
        CO2_ionCOandOdiss_data = readdlm(xsecfolder*"$(:JCO2toCOplpOpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCOandOdiss_data, length(alts)), :JCO2toCOplpOpl)

        # CO₂ + hν -> CO⁺ + O
        CO2_ionCOdiss_data = readdlm(xsecfolder*"$(:JCO2toCOplpO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCOdiss_data, length(alts)), :JCO2toCOplpO)

        # CO₂ + hν -> CO + O⁺
        CO2_ionOdiss_data = readdlm(xsecfolder*"$(:JCO2toOplpCO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionOdiss_data, length(alts)), :JCO2toOplpCO)

        # CO₂ + hν -> C⁺ + O⁺ + O
        CO2_ionCandOdiss_data = readdlm(xsecfolder*"$(:JCO2toOplpCplpO).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO2_ionCandOdiss_data, length(alts)), :JCO2toOplpCplpO)

        # NEW: H2O ionization --------------------------------------------------------------
        # Source: Roger Yelle

        # H2O + hν -> H2O⁺
        h2o_ionize_data = readdlm(xsecfolder*"$(:JH2OtoH2Opl).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionize_data, length(alts)), :JH2OtoH2Opl)

        # H2DO + hν -> HDO⁺ # TODO: replace with HDO photoionization xsects when they exist
        hdo_ionize_data = readdlm(xsecfolder*"$(:JH2OtoH2Opl).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(hdo_ionize_data, length(alts)), :JHDOtoHDOpl)

        # H2O + hν -> O⁺ + H2
        h2o_ionOdiss_data = readdlm(xsecfolder*"$(:JH2OtoOplpH2).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionOdiss_data, length(alts)), :JH2OtoOplpH2)

        # H2O + hν -> H⁺ + OH
        h2o_ionHdiss_data = readdlm(xsecfolder*"$(:JH2OtoHplpOH).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionHdiss_data, length(alts)), :JH2OtoHplpOH)

        # H2O + hν -> OH⁺ + H
        h2o_ionOHdiss_data = readdlm(xsecfolder*"$(:JH2OtoOHplpH).csv", ',', Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(h2o_ionOHdiss_data, length(alts)), :JH2OtoOHplpH)

        # NEW: CO ionization ----------------------------------------------------------------
        # Source: Roger Yelle

        # CO + hν -> CO⁺
        CO_ionize_data = readdlm(xsecfolder*"$(:JCOtoCOpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_ionize_data, length(alts)), :JCOtoCOpl)

        # CO + hν -> C + O⁺
        CO_ionOdiss_data = readdlm(xsecfolder*"$(:JCOtoCpOpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_ionOdiss_data, length(alts)), :JCOtoCpOpl)

        # CO + hν -> C⁺ + O
        CO_ionCdiss_data = readdlm(xsecfolder*"$(:JCOtoOpCpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(CO_ionCdiss_data, length(alts)), :JCOtoOpCpl)

        # NEW: Nitrogen species ionization --------------------------------------------------
        # Source: Roger Yelle

        # N₂ + hν -> N₂⁺
        N2_ionize_data = readdlm(xsecfolder*"$(:JN2toN2pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2_ionize_data, length(alts)), :JN2toN2pl)

        # N₂ + hν -> N⁺ + N
        N2_iondiss_data = readdlm(xsecfolder*"$(:JN2toNplpN).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2_iondiss_data, length(alts)), :JN2toNplpN)

        # NO₂ + hν -> NO₂⁺
        NO2_ionize_data = readdlm(xsecfolder*"$(:JNO2toNO2pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO2_ionize_data, length(alts)), :JNO2toNO2pl)

        # NO + hν -> NO⁺
        NO_ionize_data = readdlm(xsecfolder*"$(:JNOtoNOpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(NO_ionize_data, length(alts)), :JNOtoNOpl)

        # N₂O + hν -> N₂O⁺
        N2O_ionize_data = readdlm(xsecfolder*"$(:JN2OtoN2Opl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(N2O_ionize_data, length(alts)), :JN2OtoN2Opl)

        # NEW: Molecular and atomic hydrogen ionization -------------------------------------
        # Source: Roger Yelle

        # H + hν -> H⁺
        H_ionize_data = readdlm(xsecfolder*"$(:JHtoHpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H_ionize_data, length(alts)), :JHtoHpl)

        # H₂ + hν -> H₂⁺
        H2_ion_data = readdlm(xsecfolder*"$(:JH2toH2pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H2_ion_data, length(alts)), :JH2toH2pl)

        # HD + hν -> HD⁺ # TODO: Load HD crosssections when they exist
        HD_ion_data = readdlm(xsecfolder*"$(:JH2toH2pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(HD_ion_data, length(alts)), :JHDtoHDpl)

        # H₂ + hν -> H⁺ + H
        H2_iondiss_data = readdlm(xsecfolder*"$(:JH2toHplpH).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H2_iondiss_data, length(alts)), :JH2toHplpH)

        # H₂O₂ + hν -> H₂O₂⁺
        H2O2_ionize_data = readdlm(xsecfolder*"$(:JH2O2toH2O2pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(H2O2_ionize_data, length(alts)), :JH2O2toH2O2pl)

        # NEW: Oxygen and ozone ionization --------------------------------------------------
        # Source: Roger Yelle

        # O + hν -> O⁺
        O_iondiss_data = readdlm(xsecfolder*"$(:JOtoOpl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(O_iondiss_data, length(alts)), :JOtoOpl)

        # O₂ + hν -> O₂⁺
        O2_ionize_data = readdlm(xsecfolder*"$(:JO2toO2pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(O2_ionize_data, length(alts)), :JO2toO2pl)

        # # O₃ + hν -> O₃⁺
        O3_ionize_data = readdlm(xsecfolder*"$(:JO3toO3pl).csv",',',Float64, comments=true, comment_char='#')
        setindex!(xsect_dict, fill(O3_ionize_data, length(alts)), :JO3toO3pl)
    end
    
    return xsect_dict
end

function quantumyield(xsect, arr) # Checked - global free
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
        D = find_offset(z1, T1, A, B, C)  #3.66763690e+00
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

function Tpiecewise(z, Tsurf, Tmeso, Texo) # Checked - global-free
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
