module Analysis

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils

# using Photochemistry

include("/home/emc/GDrive-CU/Research-FF/PARAMETERS.jl")

export get_ncurrent, write_ncurrent, n_tot, 
       effusion_velocity, speciesbcs, 
       areadensity_to_micron_atm, molec_to_GEL, GEL_to_molecule, 
       get_flux, calculate_f, 
       search_subfolders, create_folder, searchdir, input,
       get_colors, get_grad_colors, plot_bg,
       Tpiecewise, Psat, Psat_HDO

# Array manipulation ===========================================================
function get_ncurrent(readfile)
    #=
    Retrieves the matrix of species concentrations by altitude from an HDF5
    file containing a converged atmosphere.
    =#
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()

    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    return n_current
end

function write_ncurrent(n_current, filename, alt)
    #=
    Writes out the current state of species concentrations by altitude to a file
    (for converged atmosphere). 
    =# 
    n_current_mat = Array{Float64}(undef, length(alt)-2, 
                                   length(collect(keys(n_current))));
    for ispecies in [1:length(collect(keys(n_current)));]
        for ialt in [1:length(alt)-2;]
            n_current_mat[ialt, ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end
    h5write(filename,"n_current/n_current_mat",n_current_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string, collect(keys(n_current))))
end

function write_ncurrent(n_current, filename)
    #=
    Writes out the current state of species concentrations by altitude to a file
    (for converged atmosphere). 
    =# 
    n_current_mat = Array{Float64}(undef, length(alt)-2, 
                                   length(collect(keys(n_current))));
    for ispecies in [1:length(collect(keys(n_current)));]
        for ialt in [1:length(alt)-2;]
            n_current_mat[ialt, ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end
    h5write(filename,"n_current/n_current_mat",n_current_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string, collect(keys(n_current))))
end

function n_tot(n_current, z, n_alt_index)
    #= get the total number , density at a given altitude =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function n_tot(n_current, z)
    #= get the total number density at a given altitude =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

# boundary condition functions =================================================

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

function speciesbcs(species, surface_watersat, v_eff, oflux)
    H2Osat = surface_watersat["H2O"]
    HDOsat = surface_watersat["HDO"]

    speciesbclist=Dict(
                        :CO2=>["n" 2.1e17; "f" 0.],
                        :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                        :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                        :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                        :HDO=>["n" HDOsat[1]; "f" 0.],
                        :O=>["f" 0.; "f" oflux],
                        :H2=>["f" 0.; "v" v_eff["H2"]],
                        :HD=>["f" 0.; "v" v_eff["HD"]],
                        :H=>["f" 0.; "v" v_eff["H"]],
                        :D=>["f" 0.; "v" v_eff["D"]],
                      );
    get(speciesbclist, species, ["f" 0.; "f" 0.])
end

# Conversions ==================================================================
function areadensity_to_micron_atm(numpercm2)
    # #/cm^2 * (cm^2/m^2) * (10μm-am/2.687e20 #/m^2)
    return numpercm2 * (1e4/1) * (10/2.687e20)
end

function molec_to_GEL(molecules, HorH2O)
    #=
    Converts molecules of H2O per cm^2 to a global equivalent layer in meters. 
    molecules: number of molecules of species HorH2O
    HorH2O: string, either "H" or "H2O" to specify which type the molecule is
    =#
    if HorH2O == "H"
        molecules = molecules / 2
    end
        
    return molecules * 1e4 * 18 * 1.67e-27 / 999.89
end

function GEL_to_molecule(GEL, HorH2O)
    #=
    Converts a global equivalent layer in meters to molecules per cm^2. 
    GEL: layer of water in meters
    HorH2O: string, either "H" or "H2O" to specify which molecule to return
    =#
    
    molec_H2O = (GEL * 999.89) / (1e4 * 18 * 1.67e-27)
    
    if HorH2O == "H"
        return 2 * molec_H2O
    elseif HorH2O == "H2O"
        return molec_H2O
    end
end

# Flux and f ===================================================================
function get_flux(species, readfile, oflux, temps; repro=false, therm_only=false)
    #=
    Retrieves the flux for either H or D at the top of the equilibrated 
    atmosphere.

    species: species in question, :H or :D. no error control right now
    readfile: the file with simulation results
    oflux: flux of O in /cm^2s. 
    temps: array of [Ts, Tt, Te]
    repro: whether flux is being calculated for reproduction of a past study
    therm_only: whether to return flux_t only, false by default
    =#
    n_current = get_ncurrent(readfile)

    # find the experiment type for the purpose of setting the boundary 
    # conditions correcftly 
    exptype = match(r"[a-z]{0,5}(?=_.+)",readfile).match

    # the species which can have ither D or H in them for loss.
    bearer = Dict(:D=>[:D, :HD], :H=>[:H, :HD, :H2])
    num_D_or_H = Dict(:D=>[1, 1], :H=>[1, 1, 2])

    # this dict keeps track of loss due to each species. order: H, D, H2, HD
    contrib_t = Dict(:H=>0., :D=>0., :H2=>0., :HD=>0.)
    contrib_nt = Dict(:H=>0., :D=>0., :H2=>0., :HD=>0.)

    flux_t = 0
    flux_nt = 0

    # set things up for accurate boundary conditions
    Temp(z::Float64) = Tpiecewise(z, temps[1], temps[2], temps[3])
    Temp_keepSVP(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe)
    if exptype=="temp"
        H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt))
        HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))
    else
        H2Osat = map(x->Psat(x), map(Temp, alt))
        HDOsat = map(x->Psat_HDO(x), map(Temp, alt))
    end
    surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
    H_effusion_velocity = effusion_velocity(Temp(zmax), 1.0, zmax)
    H2_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
    D_effusion_velocity = effusion_velocity(Temp(zmax), 2.0, zmax)
    HD_effusion_velocity = effusion_velocity(Temp(zmax), 3.0, zmax)

    # Used for passing a variable speciesbcs function
    v_eff = Dict("H"=>H_effusion_velocity, "D"=>D_effusion_velocity, 
                 "H2"=>H2_effusion_velocity, "HD"=>HD_effusion_velocity)

    # Calculate the thermal escape
    for (s, m) in zip(bearer[species], num_D_or_H[species])
        this_species_t_flux = m * n_current[s][end]*speciesbcs(s, surface_watersat, v_eff, oflux)[2,2]
        flux_t += this_species_t_flux
        contrib_t[s] += this_species_t_flux
    end
    
    # Nonthermal ecsape velocities for temperatures: T_exo = 158K, 205K, 264K. 
    # using ratios of thermal/nonthermal from Kras 2010, in cm/s.
    if therm_only==false
        if repro==false
            inds = Dict(150=>1, 205=>2, 250=>3) # indices for different exobase temps
            i = inds[Int(temps[3])]             # convert exobase temp to an index 
            v_nt = Dict(:H => [3.8, 49.3, 106.5], :H2 => [1, 5.04, 11.5], 
                        :D => [8.6, 15.5, 24.2], :HD => [0.19, 3.9, 8.5])  #in cm/s. TODO: are these values suspicious?????
        else
            # If reproducing past studies, need the nonthermal escape from Kras 2002.
            inds = Dict(200 => 1, 270 => 2, 350 => 3)
            i = inds[temps[3]]
            # Nonthermal: cm/s. Each species has a value recorded at T = 200K, 
            # 270K, and 350K.
            v_nt = Dict(:H => [38, 56, 89], :H2 => [12.9, 18.2, 28], :D => [17, 24, 37],
                        :HD => [8.2, 11.5, 17.7])  
        end

        for (s, m) in zip(bearer[species], num_D_or_H[species])
            this_species_nt_flux = m * n_current[s][end] * v_nt[s][i]
            flux_nt += this_species_nt_flux
            contrib_nt[s] += this_species_nt_flux
        end
    end

    if therm_only==true
        return flux_t, contrib_t
    else
        return flux_t, flux_nt, contrib_t, contrib_nt
    end
end

function calculate_f(thefile, flux_type, temps, oflux; reprod=false)
    #=
    A function to calculate f or a single simulation.

    thefile: an equilibrated atmosphere simulation for which to calculate f.

    =#
    ncur = get_ncurrent(thefile)

    # contrib dictionaries (of how each bearer species contributes to escape)
    # are not used in this function.

    t_flux_H, nt_flux_H, contrib_t_H, contrib_nt_H = get_flux(:H, thefile, oflux, temps, repro=reprod)
    t_flux_D, nt_flux_D, contrib_t_D, contrib_nt_D = get_flux(:D, thefile, oflux, temps, repro=reprod)

    if flux_type=="thermal"
        Hf = t_flux_H
        Df = t_flux_D
    elseif flux_type=="both"
        Hf = t_flux_H + nt_flux_H
        Df = t_flux_D + nt_flux_D
    elseif flux_type=="nonthermal"
        Hf = nt_flux_H
        Df = nt_flux_D
    else
        println("Invalid escape type: $(flux_type)")
    end
  
    return 2*(Df/Hf) / (ncur[:HDO][1]/ncur[:H2O][1])
end

# Plot Hacks ===================================================================
function get_colors(L, cmap)
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

function plot_bg(axob)
    axob.set_facecolor("#ededed")
    axob.grid(zorder=0, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
    # return axob
end

function get_grad_colors(L, cmap; strt=0, stp=1)
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

# TEMPERATURE ==================================================================

function Tpiecewise(z::Float64, Tsurf, Ttropo, Texo, E="")
    #= DO NOT MODIFY! If you want to change the temperature, define a
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
    E: type of experiment, used for determining if mesopause width will vary 
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of tropopause. It varies unless we're only varying the 
    # exobase temperature.
    if (E=="tropo") || (E=="surf")
        ztropo_bot = (Ttropo-Tsurf)/(lapserate)
        ztropowidth = ztropo - ztropo_bot
    else
        ztropo_bot = (Ttropo-Tsurf)/(lapserate)
        ztropowidth = ztropo - ztropo_bot
    end

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Ttropo
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

# WATER ========================================================================
# 1st term is a conversion factor to convert to (#/cm^3) from Pa. Source: Marti & Mauersberger 1993
Psat(T::Float64) = (1e-6/(boltzmannK * T))*(10^(-2663.5/T + 12.537))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# However, this function is defined on the offchance someone studies HDO.
Psat_HDO(T::Float64) = (1e-6/(boltzmannK * T))*(10^(-2663.5/T + 12.537))


# Utility ======================================================================
searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function search_subfolders(path, key)
    #=
    path: a folder containing subfolders and files.
    key: the pattern you wish to find.

    Searches the top level subfolders within path for all folders matching a 
    certain regex given by key. Does not search files or sub-subfolders.
    =#
    folders = []
    for (root, dirs, files) in walkdir(path)
        if root==path
            for dir in dirs
                push!(folders, joinpath(root, dir)) # path to directories
            end
        end
    end

    wfolders = filter(x->occursin(key, x), folders)
    return wfolders
end

function create_folder(foldername, parentdir)
    #=
    Creates a folder, or if it already exists, notifies the user.
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

function input(prompt::String="")::String
   print(prompt)
   return chomp(readline())
end

end