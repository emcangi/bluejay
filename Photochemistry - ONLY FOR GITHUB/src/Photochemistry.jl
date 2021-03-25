__precompile__()

module Photochemistry

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils

export charge_type, create_folder, deletefirst, fluxsymbol, getpos, input, next_in_loop, searchdir, search_subfolders,         # Basic utility functions
       get_colors, get_grad_colors, plot_atm, plot_bg, plot_extinction, plot_Jrates, plot_rxns, plot_temp_prof,                # plotting functions
            plot_water_profile,     
       get_column_rates, make_ratexdensity, rxn_chem, rxn_photo,                                                               # Reaction rate functions
       get_ncurrent, write_ncurrent, n_tot,                                                                                    # Atmosphere array manipulation
       boundaryconditions, effusion_velocity,                                                                                  # Boundary condition functions
       Dcoef, fluxcoefs, flux_pos_and_neg, get_flux, Keddy, lower_up, scaleH, upper_down, #total_H_or_D_flux,                  # transport functions
       chemical_jacobian, getrate, lossequations, loss_rate, meanmass, production_equations, production_rate,                  # Chemistry functions
       binupO2, co2xsect, h2o2xsect_l, h2o2xsect, hdo2xsect, ho2xsect_l, o2xsect, O3O1Dquantumyield, padtosolar, quantumyield, # Photochemistry functions
       T_all, Tpiecewise,                                                                                                      # Temperature functions
       Psat, Psat_HDO                                                                                                          # Water profile functions                                                                       

# Load the parameter file ==========================================================
# TODO: Make this smarter somehow
include("/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/PARAMETERS.jl")     
println("NOTICE: Parameter file in use is /home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/PARAMETERS.jl")

# include("/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/PARAMETERS_neutralsonly.jl")
# println("NOTICE: Parameter file in use is /home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/PARAMETERS_neutralsonly.jl")

# include("/home/emc/GDrive-CU/Research-Modeling/FractionationFactor/Code/PARAMETERS.jl")
# println("NOTICE: Parameter file in use is /home/emc/GDrive-CU/Research-Modeling/FractionationFactor/Code/PARAMETERS.jl")

# basic utility functions ==========================================================

function charge_type(spsymb::Symbol)
    if occursin("pl", String(spsymb))
        return "ion"
    elseif spsymb==:E
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
    println()
end

function delete_old_h5file(filename)
    #= 
    the HDF5 module doesn't allow implicit deletion of existing names in the function
    h5write, so I have to create this function to clean up before writing.
    =#
    current_folder = pwd()
    dircontents = readdir(current_folder)
    if filename in dircontents
        println("File $(filename) exists, deleting")
        rm(filename)
    end
end

function deletefirst(A::Array, v)
    #=
    returns: list A with its first element equal to v removed.
    =#
    index = something(findfirst(isequal(v), A), 0)  # this horrible syntax introduced by Julia.
    keep = setdiff([1:length(A);],index)
    return A[keep]
end

function fluxsymbol(x)
    #= 
    Converts string x to a symbol. f for flux. 
    =#
    return Symbol(string("f",string(x)))
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

function next_in_loop(i::Int64, n::Int64)
    #=
    returns i+1, restricted to values within [1, n]. 
    ret_i(n) returns 1.
    =#
    return i % n + 1
end

# searches path for key
searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function search_subfolders(path, key; type="folders")
    #=
    path: a folder containing subfolders and files.
    key: the text pattern in folder names that you wish to find.

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

# plot functions ===============================================================

function get_colors(L::Int64, cmap)
    #=
    Generates some colors based on a non-gradient color map for use in plotting a 
    bunch of lines all at once.
    L: number of colors to generate.
    cmap: color map name

    NOTE: This function refuses to work unless "using PlotUtils" is in the master
          module file. It's not enough to have "using PlotUtils" appear in the 
          same submodule file in which this function is defined for some reason.
          Thus I have this function here in the main file.
    =#

    Cmap = get_cmap(Symbol(cmap))
    colors_dumb = [Cmap(x) for x in range(0, stop=1, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = colors_dumb[i][1]
        c[i, 2] = colors_dumb[i][2]
        c[i, 3] = colors_dumb[i][3]
    end
    return c
end

function get_grad_colors(L::Int64, cmap; strt=0, stp=1)
    #=
    Generates some colors based on a GRADIENT color map for use in plotting a 
    bunch of lines all at once.
    L: number of colors to generate.
    cmap: color map name

    AVAILABLE MAPS: blues, viridis, pu_or, magma, plasma, inferno
    =#

    colors_dumb = [cgrad(Symbol(cmap))[x] for x in range(strt, stop=stp, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors_dumb))
        c[i, 1] = red(colors_dumb[i])
        c[i, 2] = green(colors_dumb[i])
        c[i, 3] = blue(colors_dumb[i])
    end
    return c
end

function plot_atm(ncur, spclists, savepath::String; t=nothing, iter=nothing)
    #=
    Makes a "spaghetti plot" of the species concentrations by altitude in the
    atmosphere. 

    ncur: dictionary of species densities by altitude
    spclists: a list of lists [neutrals, ions]
    savepath: path and name for saving resulting .png file
    t: timestep, for plotting the atmosphere during convergence
    iter: iteration, for plotting the atmosphere during convergence
    =#

    # Establish logical groups for axes =======================================================

    # This dictionary was created by hand. Species may be missing. If so, an error will automatically be thrown by the code
    logical_groups = Dict(2=>[:H, :H2, :H2O, :H2O2, :HO2, :HOCO, :HCO, :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OH, :OD, 
                              :Hpl, :H2pl, :H2Opl, :H3pl, :H3Opl, :HO2pl, :HCO2pl, :HCOpl, :Dpl, :DCOpl, :HDpl, :HD2pl, :H2Dpl, :H2DOpl, :OHpl, :ArHpl, :ArDpl],
                          1=>[:C,   :CH,          :CO,   :CO2,    :Ar,   :O, :O1D, :O2, :O3, :N2, 
                              :Cpl, :CHpl, :CDpl, :COpl, :CO2pl,  :Arpl, :Opl,     :O2pl, :N2pl, :NDpl],
                          3=>[:CN, :HCN, :HNO, :N, :N2O, :NH, :NH2, :NO, :NO2, 
                              :CNpl, :HCNpl, :HNOpl, :HCNHpl, :HN2Opl, :HOCpl, :Npl,  :NHpl, :NH2pl, :N2Hpl, :NOpl, :NO2pl, :N2Opl, :NH3pl, :N2Dpl, ])

    axes_by_sp = Dict()

    for k in keys(logical_groups)
        for sp in logical_groups[k]
            axes_by_sp[sp] = k
        end
    end

    # Plot neutrals and ions together =========================================================
    if length(spclists)==2  # neutrals and ions 
        
        # set up the overall plot -------------------------------------------------------------
        atm_fig, atm_ax = subplots(3, 2, sharex=false, sharey=true, figsize=(14, 16))
        subplots_adjust(wspace=0, hspace=0)
        tight_layout()
                
        # only the neutral-col axes
        atm_ax[1, 1].set_title("Neutrals")
        for i in 1:3
            plot_bg(atm_ax[i, 1])
            atm_ax[i, 1].set_xlim(1e-12, 1e18)
            atm_ax[i, 1].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end
        atm_ax[3, 1].set_xlabel(L"Species concentration (cm$^{-3}$)")
        
        # only the ion-col axes
        atm_ax[1, 2].set_title("Ions")
        for i in 1:3
            plot_bg(atm_ax[i, 2])
            atm_ax[i, 2].set_xlim(1e-2, 1e5)
            atm_ax[i, 2].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end
        atm_ax[3, 2].set_xlabel(L"Species concentration (cm$^{-3}$)")
        
        # y axes labels
        for j in 1:3
            atm_ax[j, 1].set_ylabel("Altitude [km]")
        end
        
        # plot the neutrals according to logical groups -------------------------------------------------------
        for sp in spclists[1]
            atm_ax[axes_by_sp[sp], 1].plot(ncur[sp], plot_grid, color=get(speciescolor, sp, "black"),
                                           linewidth=2, label=sp, linestyle=get(speciesstyle, sp, "-"), zorder=10)
        end
        
        # plot the ions according to logical groups ------------------------------------------------------------
        for sp in spclists[2]
            atm_ax[axes_by_sp[sp], 2].plot(ncur[sp], plot_grid, color=get(speciescolor, sp, "black"),
                                           linewidth=2, label=sp, linestyle=get(speciesstyle, sp, "-"), zorder=10)
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
    elseif length(spclists)==1
        atm_fig, atm_ax = subplots(figsize=(16,6))
        tight_layout()
        for sp in spclists[1]
            atm_ax.plot(ncur[sp], plot_grid, color=get(speciescolor, sp, "black"),
                        linewidth=2, label=sp, linestyle=get(speciesstyle, sp, "-"), zorder=1)
            atm_ax.set_xlim(1e-15, 1e18)
            atm_ax.set_ylabel("Altitude [km]")
            # atm_ax.set_title("Neutrals")
        end
        atm_ax.tick_params(which="both", labeltop=true, top=true)
        plot_bg(atm_ax)
        atm_ax.set_ylim(0, zmax/1e5)
        atm_ax.set_xscale("log")
        atm_ax.set_xlabel(L"Species concentration (cm$^{-3}$)")
        atm_ax.legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0, fontsize=16)
    else
        throw("unexpected number of species lists to plot")
    end

    if t!=nothing # there will always be a timestep
        titlestr = "time=$(t)"
        if iter!=nothing  # and there may be an iteration associated.
            titlestr = titlestr * ", iter=$(iter)"
        end
        suptitle(titlestr, y=1.05)
    else
        suptitle("Initial conditions", y=1.05)
    end

    atm_fig.savefig(savepath, bbox_inches="tight")
    close(atm_fig)
end

function plot_bg(axob)
    #=
    Takes a PyPlot axis object axob and
    1) sets the background color to a light gray
    2) creates a grid for the major tick marks
    3) turns off all the axis borders

    intended to immitate seaborn
    =#
    axob.set_facecolor("#ededed")
    axob.grid(zorder=-5, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function plot_extinction(solabs, path::String, dt::Float64, iter::Int64; tauonly=false, xsect_info=nothing, solflux=nothing)
    #=
    Used to plot the extinction in the atmosphere by wavelength and altitude

    solabs: a ROW VECTOR (sadly) of solar absorption in the atmosphere. rows = altitudes, "columns" = wavelengths,
            but being a row vector it supposedly has only one column, but each element is a row.
    path: a path at which to save the plot
    dt: timestep, used for making a unique image file name
    iter: iteration number, same purpose as dt
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

function plot_Jrates(ncur, controltemps, speciesbclist; filenameext="")
    #=
    Plots the Jrates for each photodissociation or photoionizaiton reaction.

    ncur: the atmospheric state to plot 
    t: the non-mean value of the temperature for T_exptype, i.e. might be 150 for T_surf instead of the mean 216.
    exptype: whether the "surf", "tropo" or "exo" temperature is the parameter being examined
    filenameext: something to append to the image filename
    =#

    # --------------------------------------------------------------------------------
    # calculate reaction rates x density of the species at each level of the atmosphere.
    rxd_prod = make_ratexdensity(ncur, controltemps, which="Jrates", speciesbclist)
    rxd_loss = make_ratexdensity(ncur, controltemps, which="Jrates", speciesbclist)

    # ---------------------------------------------------------------------------------
    # Plot reaction rates and transport rates by altitude
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

    whichax = Dict("HO2"=>1, "H2"=>1, "H2O"=>1, "H2O2"=>1, "H"=>1, "OH"=>1, 
                   "HDO"=>2, "HDO2"=>2, "HD"=>2, "DO2"=>2, "OD"=>2, 
                   "CO"=>3, "CO2"=>3, "O"=>3, "O3"=>3, "O2"=>3, 
                   "N2"=>4, "NO2"=>4, "N2O"=>4, "NO"=>4)
    
    fig, ax = subplots(1, 4, sharey=true, figsize=(36,6))
    subplots_adjust(wspace=0.7, bottom=0.25)
    for a in ax
        plot_bg(a)
    end
    
    # x lim starting values
    minx = Array{Float64, 1}([1e10, 1e10, 1e10, 1e10])
    maxx = Array{Float64, 1}([0, 0, 0, 0])

    # color options - distinct colors generated by colorgorical
    colopts = ["#68affc", "#99ea40", "#ea3ffc", "#1c9820", "#8a0458", "#32e195", "#eb36a0", "#20f53d", "#f7393a", "#0b5313", "#fbacf6", 
               "#7b9b47", "#442edf", "#dcda5e", "#66457a", "#f7931e", "#11677e", "#d47767", "#20d8fd", "#842411", "#b8deb7", "#0064e1", 
               "#fec9af", "#604020", "#a88c65"]
    
    i_col = [1, 1, 1, 1]  # iterators to pick from the colors for ions
    n_col = [1, 1, 1, 1]  # and neutrals, separate for each axis 1-4.
    # Collect chem production equations and total 
    for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
        axnum = whichax[match(r".+(?=\ -->)", kv[1]).match] # different species go onto different axes to keep things tidy
        lbl = "$(kv[1])"
        if occursin("pl", kv[1])
            ls = "--"
            c = colopts[i_col[axnum]]
            i_col[axnum] += 1
        else
            ls = "-"
            c = colopts[n_col[axnum]]
            n_col[axnum] += 1
        end
        ax[axnum].semilogx(kv[2], plot_grid, linestyle=ls, linewidth=1, label=lbl, color=c)
        
        # set the xlimits
        if minimum(kv[2]) <= minx[axnum]
            minx[axnum] = minimum(kv[2])
        end
        if maximum(kv[2]) >= maxx[axnum]
            maxx[axnum] = maximum(kv[2])
        end
    end
    
    suptitle("J rates", fontsize=20)
    ax[1].set_ylabel("Altitude (km)")
    for i in range(1, stop=4)
        ax[i].legend(bbox_to_anchor=(1.01, 1))
        ax[i].set_xlabel("Chemical reaction rate [check units] ("*L"cm^{-3}s^{-1})")
        # labels and such 
        if minx[i] < 1e-20
            minx[i] = 1e-20
        end
        maxx[i] = 10^(ceil(log10(maxx[i])))
        ax[i].set_xlim([minx[i], maxx[i]])
        end
    ax[1].set_title("H")
    ax[2].set_title("D")
    ax[3].set_title("C and O species")
    ax[4].set_title("N species")
    savefig(results_dir*"J_rates_"*filenameext*".png", bbox_inches="tight", dpi=300)
    # show()
    close(fig)
end

function plot_Jrates(ncur, controltemps, spc_list, speciesbclist; filenameext="")
    #=
    Plots the Jrates for each photodissociation or photoionizaiton reaction. Override for small groups of species.

    species: species for which to plot the Jrates
    =#

    # --------------------------------------------------------------------------------
    # calculate reaction rates x density of the species at each level of the atmosphere.
    rxd_prod = make_ratexdensity(ncur, controltemps, species=spc_list, which="Jrates", speciesbclist)
    rxd_loss = make_ratexdensity(ncur, controltemps, species=spc_list, which="Jrates", speciesbclist)

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
        ax.semilogx(kv[2], plot_grid, linestyle="-", linewidth=1, label=lbl)
        
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

function plot_rxns(sp::Symbol, ncur, controltemps::Array, speciesbclist; plot_indiv_rxns=false, thresh=1e-8, subfolder="", plotsfolder="", dt=nothing, num="", trans=true, extra_title=nothing, plot_timescales=false)
    #=
    Plots the production and loss rates by altitude for a given species, sp, at a given 
    snapshot in time of the atmosphere, ncur. 

    sp: species to make plot for
    ncur: the atmospheric state to plot, a dictionary of format :Sp => [densities by alt]
    controltemps: an array containing the surface, tropopause, and exobase temperatures.
    plot_indiv_rxns: whether to plot lines for the individual chemical reactions. if false, only the total will be plotted.
    thresh: a threshhold of reaction rate. will only plot the reactions that have values above this threshhold
    subfolder: an optional subfolder within results_dir
    plotsfolder: specific folder within which the plots will go.
    dt: timestep, optional, to be included in the plot title
    num: optional string to add to filename
    trans: whether to include physical transport in the total rates
    extra_title: extra text for the title, e.g. "converged"
    =#

    # --------------------------------------------------------------------------------
    # calculate reaction rates x density of the species at each level of the atmosphere.
    rxd_prod = make_ratexdensity(ncur, controltemps, species=sp, species_role="product", speciesbclist)
    rxd_loss = make_ratexdensity(ncur, controltemps, species=sp, species_role="reactant", speciesbclist)

    # ---------------------------------------------------------------------------------
    # Set up the plot and totals arrays
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

    # Calculate the total reactions per second for this species of interest
    total_prod_rate = Array{Float64}(undef, num_layers)
    total_prod_rate .= 0
    total_loss_rate = Array{Float64}(undef, num_layers)
    total_loss_rate .= 0

    # Arrays to hold the total chemical production and loss 
    total_chem_prod = deepcopy(total_prod_rate)
    total_chem_loss = deepcopy(total_loss_rate)

    # Arrays for transport don't have to be pre-declared because we only enter values into them once.

    # this will be used to determine the x limits. the three entries are the chemistry panel, transport panel, total panel
    # these are just starting values. 
    minx = [1e-8, 1e-8, 1e-8]  
    maxx = [1e2, 1e2, 1e2]

    # set up a selection of colors and linestyles so that individual reaction rate lines are distinguishable
    cols = ["#cc156d", "#7e5fe7", "#b26145", "#257950", "#265582", "#fc5931"]  # colorgorical
    ls = ["-", "--", (0, (3, 1, 1, 1)), ":"]

    # ---------------------------------------------------------------------------------
    col_i = 1
    ls_i = 1
    
    # Chemical production - add up total, and plot individual reactions if needed
    for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
        lbl = "$(kv[1])"
        if plot_indiv_rxns == true
            if any(x->x>thresh, kv[2])   # looks for values above the threshhold in the lower atmosphere (1:40 is approx surface to 80 km)
                ax[1].semilogx(kv[2], plot_grid, linestyle=ls[ls_i], marker=9, markevery=20, color=cols[col_i], linewidth=1, label=lbl)
                col_i = next_in_loop(col_i, length(cols))
                ls_i = next_in_loop(ls_i, length(ls))
            end
        end
        total_chem_prod += kv[2]
    end

    # Chemical loss 
    for kv in rxd_loss
        lbl = "$(kv[1])"
        if plot_indiv_rxns == true
            if any(x->x>thresh, kv[2])
                ax[1].semilogx(kv[2], plot_grid, linestyle=ls[ls_i], marker=8, markevery=20, color=cols[col_i], linewidth=1, label=lbl)
                col_i = next_in_loop(col_i, length(cols))
                ls_i = next_in_loop(ls_i, length(ls))
            end
        end
        total_chem_loss += kv[2]
    end

    # ---------------------------------------------------------------------------------
    # Calculate the fluxes for the species
    plottitle_ext = "" # no extra info in the plot title if flux==false
    if trans==true
        transportPL = get_transport_production_and_loss_rate(ncur, sp, controltemps, speciesbclist)
        # now separate into two different arrays for ease of addition.
        production_i = transportPL .>= 0
        loss_i = transportPL .< 0
        total_transport_prod = production_i .* transportPL
        total_transport_loss = loss_i .* abs.(transportPL)

        total_prod_rate = total_transport_prod .+ total_chem_prod
        total_loss_rate = total_transport_loss .+ total_chem_loss

        plottitle_ext = " by chemistry & transport"
    end

    # ---------------------------------------------------------------------------------
    # calculate the timescale of change
    change_timescale = ncur[sp] ./ abs.(total_prod_rate .- total_loss_rate)

    # ---------------------------------------------------------------------------------
    # plot the totals 
    # chemistry
    ax[1].semilogx(total_chem_prod, plot_grid, color="xkcd:forest green", linestyle=(0, (4,2)), marker=9, markevery=20, linewidth=2, label="Total chemical production", zorder=5)
    ax[1].semilogx(total_chem_loss, plot_grid, color="xkcd:shamrock", linestyle=(0, (4,2)), marker=8, markevery=20, linewidth=2, label="Total chemical loss", zorder=5)

    # transport - must scatter because not every altitude has production and not every altitude has loss.
    ax[2].scatter(total_transport_prod, plot_grid, color="red", marker=9, label="Total gain this layer", zorder=4)#linestyle=(1, (3,1)), markevery=20, linewidth=2)
    ax[2].scatter(total_transport_loss, plot_grid, color="blue", marker=8, label="Total loss this layer", zorder=4)#linestyle=(1, (3,1)), markevery=20, linewidth=2)
    ax[2].set_xscale("log")
    # all together
    # c_all = "black"
    ax[3].semilogx(total_prod_rate, plot_grid, color="black", marker=9, markevery=15, linewidth=2, label="Total production", zorder=3) #linestyle=(2, (1,2)), 
    ax[3].semilogx(total_loss_rate, plot_grid, color="gray", marker=8, markevery=15, linewidth=2, label="Total loss", zorder=3) #linestyle=(2, (1,2)),

    # ---------------------------------------------------------------------------------
    # Plot the timescale of dn/dt
    if plot_timescales==true
        c = "teal"
        ax[total_ax].set_xscale("log")
        ax[total_ax].plot(change_timescale, plot_grid, color=c, zorder=5)
    end

    # ---------------------------------------------------------------------------------
    # Handle labels and plot limits for the primary axis
    # set the x lims for chem axis
    minx[1] = minimum([minimum(total_chem_prod), minimum(total_chem_loss)])
    maxx[1] = maximum([maximum(total_chem_prod), maximum(total_chem_loss)])
    minx[1] = 10^(floor(log10(minx[1])))
    maxx[1] = 10^(ceil(log10(maxx[1])))

    # set the x lims for transport axis
    minx[2] = minimum(abs.(transportPL))  # this one is special because total_transport_prod, and etc are incomplete arrays.
    maxx[2] = maximum(abs.(transportPL))
    minx[2] = 10^(floor(log10(minx[2])))
    maxx[2] = 10^(ceil(log10(maxx[2])))

    # set the x lims for total prod/loss axis 
    minx[3] = minimum([minimum(total_prod_rate), minimum(total_loss_rate)])
    maxx[3] = maximum([maximum(total_prod_rate), maximum(total_loss_rate)])
    minx[3] = 10^(floor(log10(minx[3])))
    maxx[3] = 10^(ceil(log10(maxx[3])))

    # check for and correct any ridiculously low limits
    for i in 1:length(minx)
        if minx[i] < 1e-12
            minx[i] = 1e-12
        end
    end

    # labels and such 
    titles = ["Chemistry", "Transport", "Total chem+trans"]
    if plot_timescales==true
        titles = ["Chemistry", "Transport", "Total chem+trans", "Timescale of dn/dt (s)"]
    end
    dtstr = dt == nothing ? "" : ", $(dt)"
    titlestr = extra_title == nothing ? "" : ", "*extra_title
    for i in 1:length(ax)
        ax[i].legend(fontsize=12)
        ax[i].set_title(titles[i])
        ax[i].set_xlim(minx[i], maxx[i])
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

    savefig(savepathname, bbox_inches="tight", dpi=300)
    close(fig)
end

function plot_temp_prof(n_temps::Array{Float64,1}, savepath::String; i_temps=nothing, e_temps=nothing)
    #=
    Creates a .png image of the tepmeratures plotted by altitude in the atmosphere

    n_temps: an array of neutral temperature by altitude
    ion_temps: same, but for the ions
    e_temps: same but for the electrons
    savepath: where to save the resulting .png image
    =#

    fig, ax = subplots(figsize=(4,6))
    plot_bg(ax)

    plot(n_temps, alt./1e5, label="Neutrals", color=medgray)

    if i_temps != nothing
        ax.plot(i_temps, alt./1e5, label="Ions", color="xkcd:bright orange")
        ax.legend(fontsize=16)
        ax.set_xscale("log")
    end
    if e_temps != nothing
        ax.plot(e_temps, alt./1e5, label="Electrons", color="cornflowerblue")
        ax.legend(fontsize=16)
        ax.set_xscale("log")
    end

    # plot the control temps
    # ax.text(n_temps[end]*0.9, 185, L"T_{exo}")
    # ax.text(n_temps[Int64(length(n_temps)/2)]+5, 75, L"T_{tropo}")
    # ax.text(n_temps[1], 10, L"T_{surface}")

    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(collect(0:50:alt[end]/1e5))
    ax.set_yticklabels(collect(0:50:alt[end]/1e5))
    ax.set_xlabel("Temperature (K)")

    savefig(savepath*"/temp_profile.png", bbox_inches="tight") 
end

function plot_water_profile(H2Oinitf::Array, HDOinitf::Array, nH2O::Array, nHDO::Array, savepath::String; watersat=nothing)
    #=
    Plots the water profile in mixing ratio and number densities, in two panels.

    H2Oinitf: initial fraction of H2O in the atmosphere
    HDOinitf: same for HDO
    nH2O: number density of H2O
    nHDO: same for HDO
    watersat: optional. must be a list of the saturation fractions with HDO second, 
              i.e. [H2Osatfrac, HDOsatfrac]
    =#

    fig, ax = subplots(1, 2, sharey=true, sharex=false, figsize=(12,8))
    for a in ax
        plot_bg(a)
        a.tick_params(axis="x", which="minor", bottom=true, top=true)
    end
    
    # mixing ratio in PPM axis
    ax[1].semilogx(H2Oinitf/1e-6, plot_grid, color="cornflowerblue", linewidth=2)
    ax[1].semilogx(HDOinitf/1e-6, plot_grid, color="cornflowerblue", linestyle="--", linewidth=2)
    ax[1].set_xlabel("Volume Mixing Ratio [ppm]")
    ax[1].set_ylabel("Altitude [km]")
    ax[1].set_title("Water mixing ratios")
    ax[1].legend()

    # number density axis
    ax[2].semilogx(nH2O, plot_grid, color="cornflowerblue", linewidth=2, label=L"H$_2$O")
    ax[2].semilogx(nHDO, plot_grid, color="cornflowerblue", linestyle="--", linewidth=2, label="HDO")
    ax[2].set_xlabel(L"Species density [cm$^{-2}$")
    ax[2].set_title("Water number density")
    ax[2].legend()

    suptitle(L"H$_2$O and HDO model profiles")
    # save it
    savefig(savepath*"/water_profiles.png")
    close(fig)

    if watersat != nothing
        fig, ax = subplots(figsize=(6,9))
        plot_bg(ax)
        semilogx(H2Oinitf, plot_grid, color="cornflowerblue", linewidth=3, label=L"H$_2$O initial fraction")
        semilogx(watersat[2:end-1], plot_grid, color="black", alpha=0.5, linewidth=3, label=L"H$_2$O saturation")
        xlabel("Mixing ratio", fontsize=18)
        ylabel("Altitude [km]", fontsize=18)
        title(L"H$_2$O saturation fraction", fontsize=20)
        ax.tick_params("both",labelsize=16)
        legend()
        savefig(savepath*"/water_MR_and_saturation.png")
    end
end

# Reaction rate functions ======================================================
function get_column_rates(sp, ncur, controltemps, bcdict)
    #=
    I can't believe I haven't written this function earlier, but it calculates total column rates 
    for a given reaction. Bonus: the returned arrays are sorted, so you can just look at
    sorted_prod[1] to get the top production mechanism, etc. Note that the returns are NOT dictionaries!!!
    
    sp: species for which to search for reactions
    ncur: the present atmospheric state to calculate on
    controltemps: [T_surf, T_tropo, T_exo]
    bcdict: dictionary of boundary conditions
    =#
    
    rxd_prod = make_ratexdensity(ncur, controltemps, species=sp, species_role="product", bcdict);
    rxd_loss = make_ratexdensity(ncur, controltemps, species=sp, species_role="reactant", bcdict)
    
    # Make the column rates dictionary for production
    columnrate_prod = Dict()
    for k in keys(rxd_prod)
        percm2 = rxd_prod[k] .* dz
        columnrate_prod[k] = sum(percm2)
    end
    sorted_prod = sort(collect(columnrate_prod), by=x->x[2], rev=true)
    
    # Make the column rates dictionary for loss
    columnrate_loss = Dict()
    for k in keys(rxd_loss)
        percm2 = rxd_loss[k] .* dz
        columnrate_loss[k] = sum(percm2)
    end
    sorted_loss = sort(collect(columnrate_loss), by=x->x[2], rev=true)
    
    return sorted_prod, sorted_loss
end

function make_ratexdensity(ncur, controltemps::Array, speciesbclist; species=Nothing, species_role="both", which="all")
    #=
    ncur: a given result file for a converged atmosphere
    t: a specified temperature parameter, either T_surf, T_tropo, or T_tropoexo; which is identified by exptype.
    exptype: "surf", "tropo", "exo", just allows for specifiying the temperature profile.
    species: Symbol; only reactions including this species will be plotted. If it has a value, so must species_role.
    species_role: whether to look for the species as a reactant, product, or both.  If it has a value, so must species.
    which: "all", "Jrates", "krates". Whether to fill the dictionary with all reactions, only photochemistry/photoionization 
           (Jrates) or only chemistry (krates).
    =#

    # TODO: Put in some logic to catch problems where species or species_role is specified but not the other.

    # Construct temperature profile (by altitude) ------------------------------------------------------------------------------
    Temp_n(z::Float64) = T_all(z, controltemps[1], controltemps[2], controltemps[3], "neutral")
    Temp_i(z::Float64) = T_all(z, controltemps[1], controltemps[2], controltemps[3], "ion")
    Temp_e(z::Float64) = T_all(z, controltemps[1], controltemps[2], controltemps[3], "electron")

    temps_neutrals = map(Temp_n, non_bdy_layers)
    temps_ions = map(Temp_i, non_bdy_layers)
    temps_electrons = map(Temp_e, non_bdy_layers)

    # Fill in the rate x density dictionary ------------------------------------------------------------------------------
    rxn_dat =  Dict{String,Array{Float64, 1}}()

    # TODO: select here only the rates we want
    if which=="Jrates"
        selected_rxns = filter(x->occursin("J", string(x[3])), reactionnet)
    elseif which=="krates" 
        selected_rxns = filter(x->!occursin("J", string(x[3])), reactionnet)
    elseif which=="all"
        selected_rxns = deepcopy(reactionnet)
    end

    # now select only the reactions with the species in question (there is probably a less repetative way to write this)
    if species==Nothing
        filtered_rxn_list = selected_rxns
    else
        filtered_rxn_list = Any[]

        # need to make a regular expression containing the species that will ignore cases where
        # the species string appears as a subset of another species string, i.e. when O2pl is 
        # detected as being in the reaction CO2 -> CO2pl.
        species_re = r"\b"*string(species)*r"\b"

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
    end

    for rxn in filtered_rxn_list
        reactants = rxn[1]
        products = rxn[2]  # vector of the product symbols

        # get the reactants and products in string form for use in plot labels
        rxn_str = string(join(rxn[1], " + ")) * " --> " * string(join(rxn[2], " + "))

        # calculate the reaction strength = rate coefficient * species density.
        if typeof(rxn[3]) == Symbol # for photodissociation
            rxn_dat[rxn_str] = rxn_photo(ncur, rxn[1][1], rxn[3])
        else                        # bi- and ter-molecular chemistry
            rxn_dat[rxn_str] = rxn_chem(ncur, rxn[1], rxn[3], temps_neutrals, temps_ions, temps_electrons)
        end
    end

    return rxn_dat
end

function rxn_chem(ncur, reactants::Array, krate, temps_n::Array, temps_i::Array, temps_e::Array)
    #=
    Calculates total rate = (r)ate coefficient (x) species de(n)sity at each level of the atmosphere
    for bi- or tri-molecular chemical reactions. rate = k[A][B]...
    Thus, the return value is an array by altitude of a reaction rate, i.e. in #/cm^3/s for bimolecular
    and #/cm^6/s for trimolecular.

    ncur: current atmospheric state dictionary
    reactants: a list of species symbols that are reactants
    krate: the chemical reaction rate
    temps_n, _i, _e: neutral, ion, electron temperature profiles by altitude
    =#

    if in(:Opl, keys(ncur)) # this just checks if the simulation has ions - because this function is sometimes called for the neutral-only model.
        modeltype = "ions"
    else
        modeltype = "neutral"
    end
    # Calculate the sum of all bodies in a layer (M) for third body reactions. 
    # This does it in an array so we can easily plot.
    M_by_alt = sum([ncur[sp] for sp in fullspecieslist]) 
    if modeltype == "ions"
        E_by_alt = sum([ncur[sp] for sp in ionlist])
    end

    # multiply the result array by all reactant densities
    rate_x_density = ones(num_layers)
    for r in reactants
        if r != :M && r != :E
            # species densities by altitude
            rate_x_density .*= ncur[r]  # multiply by each reactant density
        elseif r == :M
            rate_x_density .*= M_by_alt
        elseif r == :E
            rate_x_density .*= E_by_alt
        else
            throw("Got an unknown symbol in a reaction rate: $(r)")
        end
    end

    # WARNING: invokelatest is key to making this work. I don't really know how. At some point I did. WITCHCRAFT
    if modeltype == "ions"
        @eval ratefunc(Tn, Ti, Te, M, E) = $krate
        rate_arr = Base.invokelatest(ratefunc, temps_n, temps_i, temps_e, M_by_alt, E_by_alt)
    else  # neutrals-only model
        @eval ratefunc(T, M) = $krate
        rate_arr = Base.invokelatest(ratefunc, temps_n, M_by_alt)
    end
    
    rate_x_density .*= rate_arr  # this is where we multiply the product of species densities by the reaction rate 
    return rate_x_density
end

function rxn_photo(ncur, reactant::Symbol, Jrate) 
    #=
    calculates J(r)ate (x) species de(n)sity at each level of the atmosphere
    for photochemical or photoionization reactions.
    Thus, the return value is an array by altitude of reaction rates in #/s.

    ncur: current atmospheric state dictionary
    reactant: a species symbol
    Jrate: Jrate symbol
    returns: an array of rate x n_reactant by altitude.
    =#
    n_by_alt = ncur[reactant]
    rate_by_alt = ncur[Jrate]
    return rate_x_density = n_by_alt .* rate_by_alt
end

# Atmosphere array manipulation ================================================
function get_ncurrent(readfile::String)
    #=
    Retrieves the matrix of species concentrations by altitude from an HDF5
    file, readfile, containing a converged atmosphere.
    =#
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()

    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies], length(n_current_mat[:, ispecies]))
    end
    return n_current
end

function n_tot(ncur, z, n_alt_index)
    #= 
    Calculates the total number density in #/cm^3 at a given altitude.

    ncur: the array of the current atmospheric state in species density by altitude
    z: altitude in question
    n_alt_index: index array, shouldn't really need to be passed in because it's global
    =#
    thisaltindex = n_alt_index[z]
    return sum( [ncur[s][thisaltindex] for s in specieslist] )
end

function n_tot(ncur, z)
    #= 
    Override for when n_alt_index is global
    =#
    thisaltindex = n_alt_index[z]
    return sum( [ncur[s][thisaltindex] for s in specieslist] )
end

function write_ncurrent(ncur, filename::String)
    #=
    Writes out the current state of species density by altitude to an .h5 file
    (for converged atmosphere). 

    ncur: the array of the current atmospheric state in species density by altitude
    filename: whatever you want to call the filename when it's saved as .h5
    =# 
    ncur_mat = Array{Float64}(undef, num_layers, length(collect(keys(ncur))));
    for ispecies in [1:length(collect(keys(ncur)));]
        for ialt in [1:num_layers;]
            ncur_mat[ialt, ispecies] = ncur[collect(keys(ncur))[ispecies]][ialt]
        end
    end
    delete_old_h5file(filename)
    h5open(filename, "w") do f # this syntax is ok because we never write multiple times to a file.
        write(f, "n_current/n_current_mat", ncur_mat)
        write(f, "n_current/alt", alt)
        write(f, "n_current/species", map(string, collect(keys(ncur))))
    end
end

# boundary condition functions =================================================
function boundaryconditions(species::Symbol, dz::Float64, ncur, controltemps::Array, speciesbclist)
    #= 
    Returns boundary conditions for species in the format:
    [n_1 -> n_0, n_0 -> n_1;      n_(num_layers) -> n_(num_layers+1), n_(num_layers+1) -> n_(num_layers)]
    where n_0 is the boundary layer from [-1, 1], n_1 is the first atmospheric layer from [1, 3],
    n_(num_layers) is the topmost atmospheric layer, and n_(num_layers+1) is the top boundary layer.

    dz: atmospheric layer thickness in cm
    ncur: the atmospheric state dictionary
    controltemps: [T_surf, T_tropo, T_exo], required to retrieve transport coefficients
    speciesbclist: dictionay of boundary conditions.
    =#
    bcs = get(speciesbclist, species, ["f" 0.; "f" 0.])
    if issubset([species],notransportspecies)
        bcs = ["f" 0.; "f" 0.]
    end

    # first element returned corresponds to lower BC, second to upper
    # BC transport rate. Within each element, the two rates correspond
    # to the two equations
    # n_b  -> NULL (first rate, depends on species concentration)
    # NULL -> n_b  (second rate, independent of species concentration)
    # note these are chemical equations not indicating atmospheric layers!
    bcvec = Float64[0 0; 0 0]
    # the signs on these are the opposite of what you expect because it is necessary to 
    # encode loss in production equations and vice versa.
    # in first element, negative means you're adding something, and vice versa.
    # in second element, positive means you're adding something, and vice versa.

    # LOWER
    if bcs[1, 1] == "n"
        # Density conditions at the surface.
        # n_1 -> n_0 is just the combined transport coefficient (D+K)/Δz² + ⌷/(2Δz), directed downward (the [1]).
        # n_0 -> n_1 condition is the upward transport coefficient, but multiplied by the boundary condition since [-1, 1] is not a real cell.
        bcvec[1,:]=[fluxcoefs(alt[2], dz, species, ncur, controltemps)[1], # 
                    lower_up(alt[1], dz, species, ncur, controltemps)*bcs[1, 2]]  # (#/s * #/cm³ = #/cm³/s)
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
        bcvec[2,:] = [fluxcoefs(alt[end-1], dz, species, ncur, controltemps)[2],
                    upper_down(alt[end], dz, species, ncur, controltemps)*bcs[1, 2]]
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

    return bcvec
end

function effusion_velocity(Texo::Float64, m::Float64, zmax)
    #=
    Returns effusion velocity for a species in cm/s

    Texo: temperature of the exobase (upper boundary) in K
    m: mass of one molecule of species in amu
    zmax: max altitude in cm
    =#
    
    # lambda is the Jeans parameter (Gronoff 2020), basically the ratio of the 
    # escape velocity GmM/z to the thermal energy, kT.
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    vth = sqrt(2*boltzmannK*Texo/(m*mH))  # this one is in m/s
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)  # this is in cm/s
    return v
end

# Transport functions ==========================================================

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

function Dcoef(T, n::Real, species::Symbol)
    #=
    Calculates molecular diffusion coefficient for a particular slice of the
    atmosphere using D = AT^s/n, from Banks and Kockarts Aeronomy, part B,
    pg 41, eqn 15.30 and table 15.2 footnote

    T: temperature (K)
    n: number of molecules (all species) this altitude, usually calculated by n_tot
    species: whichever species we are calculating for
    =#
    dparms = diffparams(species)
    return dparms[1]*1e17*T^(dparms[2])/n
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
diffparams(species) = get(Dict(:H=>[8.4, 0.597], :H2=>[2.23, 0.75],
                               :D=>[5.98, 0.597], :HD=>[1.84, 0.75],
                               :Hpl=>[8.4, 0.597], :H2pl=>[2.23, 0.75],
                               :Dpl=>[5.98, 0.597], :HDpl=>[1.84, 0.75]),
                               species,[1.0, 0.75])

function fluxcoefs(z, dz, Kv, Dv, Tv, Hsv, H0v, species::Symbol)
    #= 
    base function to generate coefficients of the transport network. 

    z: Float64; altitude in cm.
    dz: Float64; altitude layer thickness ("resolution")
    Kv: Array; 3 elements (lower, this, and upper layer). eddy diffusion coefficient
    Dv: Array; 3 elements, same as K. molecular diffusion coefficient
    Tv: Array; 3 elements, same as K. temperature
    Hsv: Array; 3 elements, same as K. scale height by species
    H0v: Array; 3 elements, same as K. mean atmospheric scale height
    species: species symbol 

    v refers to "vector"
    u refers to "upper" (the upper parcel)
    l refers to "lower" (the lower parcel)

    Returns: a list of the coefficients for [downward, upward] flux for a given atmospheric layer.
    Note that even though it's defined as being between a layer and the one above or below, the value is 
    evaluated at the center of the layer 
    units are in 1/s.
    =#

    # Calculate the coefficients between this layer and the lower layer. 
    Dl = (Dv[1] + Dv[2])/2.0
    Kl = (Kv[1] + Kv[2])/2.0
    Tl = (Tv[1] + Tv[2])/2.0
    dTdzl = (Tv[2] - Tv[1])/dz
    Hsl = (Hsv[1] + Hsv[2])/2.0
    H0l = (H0v[1] + H0v[2])/2.0

    # two flux terms: eddy diffusion and gravity/thermal diffusion.
    # these are found in line 5 of Mike's transport_as_chemistry.pdf:
    # sumeddy = (D+K)/(Δz²), gravthermal = ☐/(2Δz), where ☐ = {D(1/H + 1+(α/T)(dT/dz)) + K(1/H_H + (1/T)(dT/dz))}
    sumeddyl = (Dl+Kl)/dz/dz
    gravthermall = (Dl*(1/Hsl + (1+thermaldiff(species))/Tl*dTdzl) +
                    Kl*(1/H0l + 1/Tl*dTdzl))/(2*dz)

    # Now the coefficients between this layer and upper layer.
    Du = (Dv[2] + Dv[3])/2.0
    Ku = (Kv[2] + Kv[3])/2.0
    Tu = (Tv[2] + Tv[3])/2.0
    dTdzu = (Tv[3] - Tv[2])/dz
    Hsu = (Hsv[2] + Hsv[3])/2.0
    H0u = (H0v[2] + H0v[3])/2.0

    sumeddyu = (Du+Ku)/dz/dz  # this is the line where we divide by cm^2
    gravthermalu = (Du*(1/Hsu + (1 + thermaldiff(species))/Tu*dTdzu) +
                    Ku*(1/H0u + 1/Tu*dTdzu))/(2*dz)

    # this results in the following coupling coefficients; sumeddy + gravthermal = (D+K)/(Δz²) + ☐/(2Δz), units 1/s <-----_!!!!! important
    # first row is this term between layer i and i-1, second row between layer i and i+1
    return [sumeddyl+gravthermall, # down
            sumeddyu-gravthermalu] # up; negative because gravity points down. I think that's why.
end

function fluxcoefs(z, dz, species::Symbol, ncur, controltemps::Array)
    #=
    NEW VERSION 

    generates the coefficients K, D, T, Hs if they are not supplied (most common)

    z: a specific altitude in cm
    dz: thickness of an altitude later (2 km, but in cm)
    species: the species for which to calculate the coefficient. Symbol
    ncur: array of species densities by altitude, the current state of the atmosphere
    controltemps: T_surf, T_tropo, T_exo 

    p: upper layer ("plus")
    0: this layer
    m: lower layer ("minus")
    =#

    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    Tp = T_all(z+dz, controltemps[1], controltemps[2], controltemps[3], species_type)
    T0 = T_all(z, controltemps[1], controltemps[2], controltemps[3], species_type)
    Tm = T_all(z-dz, controltemps[1], controltemps[2], controltemps[3], species_type)

    ntp = n_tot(ncur, z+dz)
    nt0 = n_tot(ncur, z)
    ntm = n_tot(ncur, z-dz)
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)

    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = scaleH(z+dz, species, controltemps)
    Hs0 = scaleH(z, species, controltemps)
    Hsm = scaleH(z-dz, species, controltemps)
    H0p = scaleH(z+dz, Tp, ncur)
    H00 = scaleH(z, T0, ncur)
    H0m = scaleH(z-dz, Tm, ncur)

    # return the coefficients
    return fluxcoefs(z, dz, [Km , K0, Kp], [Dm , D0, Dp], [Tm , T0, Tp],
                     [Hsm, Hs0, Hsp], [H0m, H00, H0p], species)
end

function flux_pos_and_neg(fluxarr)
    #=
    Receives the output of get_flux and generates two arrays, one with the positive flux
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

function get_flux(ncur, species::Symbol, controltemps::Array, speciesbclist)
    #=
    NEW VERSION : THIS IS THE BETTER VERSION NOW! But only for fluxes.
    
    Returns a 1D array of flux (#/cm²/s) for a given species at each boundary between layers of the
    atmosphere, including the extreme boundaries (at 1 km and 249 km) and the between-air-layers boundaries
    (3, 5...247 km).
    
    ncur: Array; species number density by altitude
    species: Symbol
    controltemps: [T_surf, T_tropo, T_exo] in use, which are needed for getting the transport coefficients.
    speciesbclist: the boundary condition dictionary.

    returns: Array of bulk flow values (#/cm³/s) at each altitude layer boundary.  
             i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
             and the end of the array is the boundary at 249 km.
    =#
    
    # each element in thesecoefs has the format [downward flow (i to i-1), upward flow (i to i+1)]. 
    # units 1/s, 
    thesecoefs = [fluxcoefs(a, dz, species, ncur, controltemps) for a in non_bdy_layers] 

    bcs = boundaryconditions(species, dz, ncur, controltemps, speciesbclist)
    
    net_bulk_flow = fill(convert(Float64, NaN), length(alt)-1)  # units #/cm^3/s; tracks the cell boundaries, of which there are length(alt)-1

    # We will calculate the net flux across each boundary, with sign indicating direction of travel.
    # First boundary is the boundary between the surface layer and the first atmospheric layer (alt = 1 km)
    net_bulk_flow[1] = (bcs[1, 2]  # increase of the lowest atmospheric layer's density. Will always be 0 unless the species has a density or flux condition
                       - ncur[species][1]*bcs[1, 1]) # lowest atmospheric layer --> surface ("depositional" term)
                        
    for ialt in 2:num_layers  # now iterate through every cell boundary within the atmosphere. boundaries at 3 km, 5...247. 123 elements.
        net_bulk_flow[ialt] = (ncur[species][ialt-1]*thesecoefs[ialt-1][2]   # coming up from below: cell i-1 to cell i
                              - ncur[species][ialt]*thesecoefs[ialt][1])     # leaving to the layer below: downwards: cell i to cell i-1
    end

    # now the top boundary - between 124th atmospheric cell (alt = 249 km)
    net_bulk_flow[end] = (ncur[species][end]*bcs[2, 1] # into exosphere from the cell
                         - bcs[2, 2]) # into top layer from exosphere. do not question this
                
    return net_bulk_flow .* dz # now it is a flux. hurrah.
end

function get_transport_production_and_loss_rate(ncur, species::Symbol, controltemps::Array, speciesbclist)
    #=
    Returns a 1D array of production and loss by transport for a given species at each boundary between layers of the
    atmosphere, including the extreme boundaries (at 1 km and 249 km) and the between-air-layers boundaries
    (3, 5...247 km).
    
    ncur: Array; species number density by altitude
    species: Symbol
    controltemps: [T_surf, T_tropo, T_exo] in use, which are needed for getting the transport coefficients.
    speciesbclist: the boundary condition dictionary.

    returns: Array of production and loss (#/cm³/s) at each altitude layer boundary.  
             i = 1 in the net_bulk_flow array corresponds to the boundary at 1 km,
             and the end of the array is the boundary at 249 km.
    =#
    
    # each element in thesecoefs has the format [downward, upward]
    thesecoefs = [fluxcoefs(a, dz, species, ncur, controltemps) for a in alt[2:end-1]]

    # thesebcs has the format [lower bc; upper bc], where each row contains a 
    # character showing the type of boundary condition, and a number giving its value
    thesebcs = boundaryconditions(species, dz, ncur, controltemps, speciesbclist)

    transport_PL = fill(convert(Float64, NaN),length(intaltgrid))

    # These are the derivatives, which should be what we want (check math)
    transport_PL[1] = ((ncur[species][2]*thesecoefs[2][1]  # in from layer above
                        -ncur[species][1]*thesecoefs[1][2]) # out to layer above
                    +(-ncur[species][1]*thesebcs[1, 1] # out to boundary layer
                      +thesebcs[1, 2])) # in from the boundary layer
    for ialt in 2:length(intaltgrid)-1
        transport_PL[ialt] = ((ncur[species][ialt+1]*thesecoefs[ialt+1][1]  # coming in from above
                               -ncur[species][ialt]*thesecoefs[ialt][2])    # leaving out to above layer
                             +(-ncur[species][ialt]*thesecoefs[ialt][1]     # leaving to the layer below
                               +ncur[species][ialt-1]*thesecoefs[ialt-1][2]))  # coming in from below
    end
    transport_PL[end] = ((thesebcs[2, 2]
                          - ncur[species][end]*thesebcs[2, 1])
                        + (-ncur[species][end]*thesecoefs[end][1]
                           +ncur[species][end-1]*thesecoefs[end-1][2]))
    return transport_PL
end

function Keddy(z::Real, nt::Real)
    #=
    eddy diffusion coefficient, stolen from Krasnopolsky (1993).
    Scales as the inverse sqrt of atmospheric number density

    z: some altitude in cm.
    nt: number total of species at this altitude (I think)
    =#
    z <= 60.e5 ? 10^6 : 2e13/sqrt(nt)
end

function Keddy(ncur, z)
    #=
    Override for when nt is not provided, and only the current atmospheric state
    (ncur) is
    =#
    z <= 60.e5 ? 10^6 : 2e13/sqrt(n_tot(ncur, z))
end

function lower_up(z, dz, species::Symbol, ncur, controltemps::Array)
    #= 
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one above. 
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus") 

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    ncur: Array; species number density by altitude

    returns: array of fluxcoefs
    =#

    species_type = charge_type(species)

    Tp = T_all(z+dz, controltemps[1], controltemps[2], controltemps[3], species_type)
    T0 = T_all(z, controltemps[1], controltemps[2], controltemps[3], species_type)
    Tm = 1

    ntp = n_tot(ncur, z+dz)
    nt0 = n_tot(ncur, z)
    ntm = 1
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z,nt0)
    Km = 1

    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = 1
    Hsp = scaleH(z+dz, species, controltemps)
    Hs0 = scaleH(z, species, controltemps)
    Hsm = 1
    H0p = scaleH(z+dz, Tp, ncur)
    H00 = scaleH(z, T0, ncur)
    H0m = 1

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[2]
end

function scaleH(z, T::Float64, mm::Real)
    #= 
    Computes the species-specific scale height (cm) of the atmosphere for the species with mass multiplier mm (* hydrogen mass) 
    at altitude z

    z: Float or Int; unit: cm. altitude in atmosphere at which to calculate scale height
    T: temperature in Kelvin
    mm: species mass multiplier in amu 
    =#
    return boltzmannK*T/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
    # constants are in MKS. Convert to m and back to cm.
end

function scaleH(z, T::Float64, species::Symbol)
    #=
    Overload: Scale height for a particular species
    =#
    mm = speciesmolmasslist[species]
    scaleH(z, T, mm)
end

function scaleH(z, T::Float64, ncur)
    #= 
    Overload: scale height when mean mass is not provided but general atmospheric
              state (ncur) is
    =#
    mm = meanmass(ncur, z)
    scaleH(z, T, mm)
end

function scaleH(z, species::Symbol, controltemps::Array)
    #=
    NEW VERSION for make_ratexdensity

    =#  

    # since we need to pass in temperatures.
    #Th is normally T, but I changed the temperature function to T, so now it's Th, the h being for scale height.
    Th = T_all(z, controltemps[1], controltemps[2], controltemps[3], charge_type(species))
    mm = speciesmolmasslist[species]
    return boltzmannK*Th/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
end

# thermal diffusion factors (all verified with Krasnopolsky 2002)
thermaldiff(species) = get(Dict(:H=>-0.25, :H2=>-0.25, :D=>-0.25, :HD=>-0.25,
                                :He=>-0.25, 
                                :Hpl=>-0.25, :H2pl=>-0.25, :Dpl=>-0.25, :HDpl=>-0.25,
                                :Hepl=>-0.25), species, 0)

function upper_down(z, dz, species::Symbol, ncur, controltemps::Array)
    #= 
    NEW VERSION 

    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one below. 
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus") 

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    ncur: Array; species number density by altitude

    returns: return of fluxcoefs
    =#

    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    Tp = 1
    T0 = T_all(z, controltemps[1], controltemps[2], controltemps[3], species_type)
    Tm = T_all(z-dz, controltemps[1], controltemps[2], controltemps[3], species_type)

    ntp = 1
    nt0 = n_tot(ncur, z)
    ntm = n_tot(ncur, z-dz)
    Kp = 1
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)

    Dp = 1
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = 1
    Hs0 = scaleH(z, species, controltemps)
    Hsm = scaleH(z-dz, species, controltemps)
    H0p = 1
    H00 = scaleH(z, T0, ncur)
    H0m = scaleH(z-dz, Tm, ncur)

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[1]
end

# function total_H_or_D_flux(species, readfile, controltemps, speciesbclist; oflux=1.2e8, repro=false, therm_only=false)
#     #=
#     Retrieves the flux for either H or D at the top of the equilibrated 
#     atmosphere. This function is redundant with get_flux so it is commented out.
    
#     species: species in question, :H or :D. no error control right now
#     readfile: the file with simulation results
#     oflux: flux of O in /cm^2s. 
#     temps: array of [Ts, Tt, Te]
#     repro: whether flux is being calculated for reproduction of a past study
#     therm_only: whether to return flux_t only, false by default
#     =#
#     n_current = get_ncurrent(readfile)

#     # the species which can have ither D or H in them for loss.
#     bearer = Dict(:D=>[:D, :HD], :H=>[:H, :HD, :H2])
#     num_D_or_H = Dict(:D=>[1, 1], :H=>[1, 1, 2])

#     # this dict keeps track of loss due to each species. order: H, D, H2, HD
#     contrib_t = Dict(:H=>0., :D=>0., :H2=>0., :HD=>0.)
#     contrib_nt = Dict(:H=>0., :D=>0., :H2=>0., :HD=>0.)

#     # total thermal and non-thermal fluxes, separated.
#     flux_t = 0
#     flux_nt = 0

#     # Calculate the thermal escape
#     for (s, m) in zip(bearer[species], num_D_or_H[species])
#         println("the boundary condition being multiplied for $(s): $(get(speciesbclist, s, ["f" 0.; "f" 0.])[2,2])")
#         this_species_t_flux = m * n_current[s][end]*get(speciesbclist, s, ["f" 0.; "f" 0.])[2,2]
#         flux_t += this_species_t_flux
#         contrib_t[s] += this_species_t_flux
#     end
    
#     # If inclusion of non-thermal escape is requested (therm_only==false), this section
#     # will add the effect of non-thermal escape. 
#     if therm_only==false
#         # Nonthermal ecsape velocities for temperatures: T_exo = 150K, 205K, 250K,
#         # using ratios of thermal/nonthermal from Kras 2002, in cm/s.
#         inds = Dict(150=>1, 205=>2, 250=>3) # indices for different exobase temps
#         i = inds[Int(controltemps[3])]             # convert exobase temp to an index 
#         v_nt = Dict(:H => [26.7, 34.4, 45.1], :H2 => [8.05, 10.2, 13.26], # assuming 2nd order polynomial
#                     :D => [10.65, 13.47, 17.47], :HD => [5.65, 6.00, 7.42])  # in cm/s.

#         for (s, m) in zip(bearer[species], num_D_or_H[species])
#             this_species_nt_flux = m * n_current[s][end] * v_nt[s][i]
#             flux_nt += this_species_nt_flux
#             contrib_nt[s] += this_species_nt_flux
#         end

#         return flux_t, flux_nt, contrib_t, contrib_nt
#     else
#         return flux_t, contrib_t
#     end
# end

# chemistry functions ==========================================================

function chemical_jacobian(chemnetwork, transportnetwork, specieslist, dspecieslist)
    #= 
    Compute the symbolic chemical jacobian of a supplied chemnetwork and transportnetwork
    for the specified specieslist. Returns three arrays suitable for
    constructing a sparse matrix: lists of the first and second indices
    and the symbolic value to place at that index.
    =#

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)  # this is generally the active species. 
    ndspecies = length(dspecieslist)  # this is generally the active species at either the same layer, above, or below.

    for i in 1:nspecies # for each species
        ispecies = specieslist[i]
        # get the production and loss equations
        peqn = []
        leqn = []
        if issubset([ispecies],chemspecies)
            peqn = [peqn; production_equations(chemnetwork, ispecies)] 
            leqn = [leqn; loss_equations(chemnetwork, ispecies)]
        end
        if issubset([ispecies],transportspecies)
            peqn = [peqn; production_equations(transportnetwork, ispecies)]
            leqn = [leqn; loss_equations(transportnetwork, ispecies)]
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

function getrate(chemnet, transportnet, species::Symbol)
    #=
    Creates a symbolic expression for the rate at which a given species is
    either produced or lost. Production is from chemical reaction yields or
    entry from other atmospheric layers. Loss is due to consumption in reactions
    or migration to other layers.
    =#
    rate = :(0.0)
    if issubset([species],chemspecies)
        rate = :($rate
               + $(production_rate(chemnet, species))
               - $(      loss_rate(chemnet, species)))
    end
    if issubset([species],transportspecies)
        rate = :($rate
               + $(production_rate(transportnet, species))
               - $(      loss_rate(transportnet, species)))
    end

    return rate
end

function loss_equations(network, species::Symbol)
    #=  
    given a network of equations in the form of reactionnet, this
    function returns the loss equations and rate coefficient for all
    reactions where the supplied species is consumed, in the form of an array
    where each entry is of the form [reactants, rate] 
    =#

    # get list of all chemical reactions species participates in:
    speciespos = getpos(network, species)
    # find pos where species is on LHS but not RHS:
    lhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(lhspos, rhspos)
        lhspos = deletefirst(lhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    losseqns=map(x->vcat(Any[network[x][1]...,network[x][3]]), lhspos)
    # automatically finds a species where it occurs twice on the LHS
end

function loss_rate(network, species::Symbol)
    #= 
    return a symbolic expression for the loss rate of species in the
    supplied reaction network. Format is a symbolic expression containing a sum
    of reactants * rate. 
    =#
    leqn=loss_equations(network, species) # get the equations
    lval=:(+($( # and add the products together
               map(x->:(*($(x...))) # take the product of the
                                    # concentrations and coefficients
                                    # for each reaction
                   ,leqn)...)))
end

function meanmass(ncur, z)
    #= 
    find the mean molecular mass at a given altitude 

    ncur: Array; species number density by altitude
    z: Float64; altitude in atmosphere in cm

    return: mean molecular mass in amu
    =#
    thisaltindex = n_alt_index[z]
    c = [ncur[sp][thisaltindex] for sp in specieslist]
    m = [speciesmolmasslist[sp] for sp in specieslist]
    return sum(c.*m)/sum(c)
end

function production_equations(network, species::Symbol)
    #= 
    given a network of equations in the form of reactionnet, this
    function returns the production equations and rate coefficient for all
    reactions where the supplied species is produced, in the form of an array
    where each entry is of the form [reactants, rate] 
    =#

    speciespos = getpos(network, species) # list of all reactions where species is produced
    # find pos where species is on RHS but not LHS
    lhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], #select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(rhspos, lhspos)
        rhspos = deletefirst(rhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    prodeqns = map(x->vcat(Any[network[x][1]...,network[x][3]]),
                 # automatically finds and counts duplicate
                 # production for each molecule produced
                 rhspos)

    return prodeqns
end

function production_rate(network, species::Symbol)
    #= 
    return a symbolic expression for the loss rate of species in the
    supplied reaction network.
    =#

    # get the reactants and rate coefficients
    peqn = production_equations(network, species)

    # add up and take the product of each set of reactants and coeffecient
    pval = :(+ ( $(map(x -> :(*($(x...))), peqn) ...) ))
end

# photochemistry functions ========================================================
function binupO2(list)
    #=
    For O2

    ...details on this missing
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

function co2xsect(co2xdata, T::Float64)
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

function h2o2xsect_l(l::Float64, T::Float64)
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

function h2o2xsect(h2o2xdata, T::Float64)
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

function hdo2xsect(hdo2xdata, T::Float64)
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

function ho2xsect_l(l::Float64)
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

function o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, T::Float64)
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

function quantumyield(xsect::Array, arr)
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

# Temperature functions ========================================================

function T_all(z::Float64, Tsurf::Float64, Ttropo::Float64, Texo::Float64, sptype::String)
    #= 
    The temperature function to END ALL TEMPERATURE FUNCTIONS!!
    Handles Neutrals!
    Ions!
    ELECTRONS!!
    WHAT WILL THEY THINK OF NEXT?

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 "half-Gaussian" function for temperatures 
    altitudes above the tropopause, with a constant lapse rate (1.4K/km) 
    in the lower atmosphere. The tropopause width is allowed to vary
    in certain cases.

    z: altitude above surface in cm
    Tsurf: Surface temperature in K
    Tropo: tropopause tempearture
    Texo: exobase temperature
    sptype: "neutral", "ion" or "electron". NECESSARY!
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of tropopause. This code allows it to vary for surface or 
    # tropopause experiments.
    ztropo_bot = (Ttropo-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    function T_upper_atmo_electrons()
        # for electrons, this is valid above 130 km in height.
        A = 2.24100983e+03
        B = 1.84024165e+01
        C = 1.44238590e-02
        D = 3.66763690e+00
        return A*(exp(-B*exp(-C*z/1e5))) + D
    end

    function T_upper_atmo_neutrals()
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    end

    function T_upper_atmo_ions()
        # valid above 160 km 
        A = 9.96741381e2
        B = 1.66317054e2
        C = 2.49434339e-2
        D = 1.29787053e2
        return A*(exp(-B*exp(-C*(z/1e5 + 0)))) + D
    end


    # this is the WOOO--OOORRRST but someday I'll make it better MAYBE???
    if z > 160e5 
        if sptype=="neutral"
            return T_upper_atmo_neutrals()
        elseif sptype=="ion"
            return T_upper_atmo_ions() # only valid up here
        elseif sptype=="electron"
            T_upper_atmo_electrons()
        end
    elseif 130e5 < z <= 160e5
        if sptype=="neutral"
            return T_upper_atmo_neutrals()
        elseif sptype=="ion"
            return T_upper_atmo_neutrals()
        elseif sptype=="electron"
            T_upper_atmo_electrons() # only valid up here
        end
    elseif ztropo < z <= 130e5   # upper atmosphere until neutrals, ions, electrons diverge
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo - ztropowidth < z <= ztropo  # tropopause
        return Ttropo
    elseif z < ztropo-ztropowidth  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

function Tpiecewise(z::Float64, Tsurf::Float64, Ttropo::Float64, Texo::Float64)
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
    Tropo: tropopause tempearture
    Texo: exobase temperature
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of tropopause. This code allows it to vary for surface or 
    # tropopause experiments.
    ztropo_bot = (Ttropo-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Ttropo
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

# Water functions ==============================================================
# 1st term is a conversion factor to convert to (#/cm^3) from Pa. Source: Marti & Mauersberger 1993
Psat(T::Float64) = (1e-6/(boltzmannK * T))*(10^(-2663.5/T + 12.537))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# However, this function is defined on the offchance someone studies HDO.
Psat_HDO(T::Float64) = (1e-6/(boltzmannK * T))*(10^(-2663.5/T + 12.537))

end