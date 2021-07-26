# progress bar witchcraft
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

t1 = time()
using Revise
using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
# using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Photochemistry  # custom module for this project
using OrdinaryDiffEq
using Sundials

t2 = time()

println("Time to load modules: $(round(t2-t1, digits=1)) seconds")
include("PARAMETERS-CO2Only.jl")

################################################################################
#                                   FUNCTIONS                                  #
################################################################################

#=
These functions are required to be in this file for one of two reasons:
1) Because they call a function that is couched in an @eval statement, which 
   cannot be relocated,
2) They are the main routine functions and will not be shared by any other scripts.

For additional functions, see the Photochemistry module.
=#

# chemistry functions ==========================================================

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] .= ratefn_local([nthismat[:,1]; 
                                     nthismat[:,2];
                                     fill(1.0, length(activespecies));
                                     inactivemat[:,1]; 
                                     Jrates[:,1]; 
                                     Tn[1]; Ti[1]; Te[1];
                                     tup[:,1];
                                     tlower[:,1];
                                     tdown[:,2];
                                     tlower[:,2]]...)

    # iterate through other altitudes in the lower atmosphere
    for ialt in 2:(num_layers-1)
        returnrates[:,ialt] .= ratefn_local([nthismat[:,ialt];
                                            nthismat[:,ialt+1];
                                            nthismat[:,ialt-1];
                                            inactivemat[:,ialt];
                                            Jrates[:,ialt];
                                            Tn[ialt]; Ti[ialt]; Te[ialt];
                                            tup[:,ialt]; 
                                            tdown[:,ialt];
                                            tdown[:,ialt+1]; 
                                            tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude
    returnrates[:,end] .= ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       Tn[end]; Ti[end]; Te[end];
                                       tupper[:,1]; 
                                       tdown[:,end];
                                       tupper[:,2]; 
                                       tup[:,end-1]]...)

    return [returnrates...;]
end


function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    #=
    Collects coordinate tuples of (I, J, V) [row index, column index, value] for a sparse matrix
    representing the chemical jacobian of the atmospheric system. 

    nthis: The atmospheric densities array, but flattened, in the form [n_CO(z=0), n_CO2(z=0)...n_N2Dpl(z=0), n_CO(z=2)...n_N2Dpl(z=250)]
    inactive: A flattened array of the atmospheric densities of any inactive species, same format as nthis. Functionally constant.
    activespecies: List of active species (const)
    inactivespecies: List of inactivespecies (const)
    Jrates: Flattened array of Jrates, same format as nthis.
    Tn, Ti, Te: Temperature-vs-altitude arrays for neutrals, ions, electrons. (const)
    tup, tdown: Transport coefficients
    tlower, tupper: Transport coefficients

    dt: Can be specified to manually control the timestep by which the derivatives are multiplied.
        If not supplied, then dt=1 so that the external solver can manage the multiplication by time.
    =#

    nthismat = reshape(nthis, (length(activespecies), num_layers))
    inactivemat = reshape(inactive, (length(inactivespecies), num_layers))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    # tc___ are the coordinate tuples containing (I, J, V) to be used to fill a sparse matrix.
    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1]; nthismat[:,2]; fill(1.0, length(activespecies));
                                                  inactivemat[:,1]; Jrates[:,1];
                                                  Tn[1]; Ti[1]; Te[1]; 
                                                  tup[:,1]; tlower[:,1];
                                                  tdown[:,2]; tlower[:,2]]...)


    # add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])

    # and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(num_layers-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt]; Tn[ialt]; Ti[ialt]; Te[ialt];
                                                      tup[:,ialt];
                                                      tdown[:,ialt];
                                                      tdown[:,ialt+1];
                                                      tup[:,ialt-1]]...) 
        # add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclocal[2].+(ialt-1)*length(activespecies))
        append!(chemJval, tclocal[3])
        # and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tcupper[2].+(ialt  )*length(activespecies))
        append!(chemJval, tcupper[3])
        # and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclower[2].+(ialt-2)*length(activespecies))
        append!(chemJval, tclower[3])
    end

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,end];
                                              fill(1.0, length(activespecies));
                                              nthismat[:,end-1];
                                              inactivemat[:,end];
                                              Jrates[:,end];
                                              Tn[end]; Ti[end]; Te[end]; 
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1]]...)

    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(num_layers-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(num_layers-1)*length(activespecies))
    append!(chemJval, tclocal[3])

    # and the lower densities
    append!(chemJi, tclower[1].+(num_layers-1)*length(activespecies))
    append!(chemJj, tclower[2].+(num_layers-2)*length(activespecies))
    append!(chemJval, tclower[3])

    return sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +);
end

function update_Jrates!(n_cur_densities::Dict{Symbol, Array{Float64, 1}})
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#

    # Initialize an array, length=number of active layers
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, num_layers)
    for i in range(1, length=num_layers)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_cur_densities[species][ialt]*dz

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
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1,
                       solarabs[ialt],1)
        end
    end

    # solarabs now records the total optical depth of the atmosphere at
    # each wavelength and altitude

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2].*exp.(-solarabs[ialt])
    end

    # each species absorbs according to its cross section at each
    # altitude times the actinic flux.
    # BLAS.dot includes an integration (sum) across wavelengths, i.e:
    # (a·b) = aa + ab + ab + bb etc that kind of thing
    for j in Jratelist
        n_cur_densities[j] = similar(n_cur_densities[:CO2])
        n_cur_densities[j] .= 0
        for ialt in [1:nalt;]
            n_cur_densities[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1,
                                                 crosssection[j][ialt+1], 1)
        end
    end
end

# main routine functions =======================================================
function get_transport_and_J_rates(n, inactive, activesp, inactivesp)
    #=
    This takes the current densities of all active species (n), transforms
    back into an atmospheric state dictionary, calculates the current Jrates, 
    and then returns the transport coefficients, boundary layer transport
    coefficients, and Jrates needed to run ratefn and chemical_jacobian.
    =#
    
    # transform n vector back into n_current so we can operate on it --------------------
    n_cur_active = unflatten_atm(n, activesp)
    n_cur_inactive = unflatten_atm(inactive, inactivesp)
    
    n_cur_all = Dict(vcat([k=>n_cur_active[k] for k in keys(n_cur_active)],
                           [k=>n_cur_inactive[k] for k in keys(n_cur_inactive)]))
    
    # Calculate things which will be passed into ratefn() -------------------------------
    #=
    Calculates things which are needed for both the production and loss equation and
    the chemical jacobian.
    =#
    Jrates = deepcopy(Float64[n_cur_all[sp][ialt] for sp in Jratelist, ialt in 1:length(non_bdy_layers)])
    
    # transport coefficients:
    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    tup = Float64[issubset([sp], notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_cur_all, [T_surf, T_tropo, T_exo])[2] for sp in specieslist, a in non_bdy_layers]
    tdown = Float64[issubset([sp], notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_cur_all, [T_surf, T_tropo, T_exo])[1] for sp in specieslist, a in non_bdy_layers]
    tlower = permutedims(reduce(hcat, [boundaryconditions(sp, dz, n_cur_all, [T_surf, T_tropo, T_exo], speciesbclist)[1,:] for sp in specieslist]))
    tupper = permutedims(reduce(hcat, [boundaryconditions(sp, dz, n_cur_all, [T_surf, T_tropo, T_exo], speciesbclist)[2,:] for sp in specieslist]))
    
    return Jrates, tup, tdown, tlower, tupper
end

function PnL_eqn(dn, n, p, t)
    #=
    This is the production and loss equation, non-discretized because the
    solver will do that. It just returns the rates of change for each species,
    i.e. P-L = dn/dt.
    =#

    n[n .< 0] .= 0

    # Unpack the parameters needed to call ratefn 
    inactive, inactivesp, activesp, Tn, Ti, Te = p 
    
    Jrates, tup, tdown, tlower, tupper = get_transport_and_J_rates(n, inactive, activesp, inactivesp)

    dn .= ratefn(n, inactive, inactivesp, activesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
end

function make_jacobian_wrapper(J, n, p, t)
    #=
    dear god
    =#
    J .= make_jacobian(n, p, t)
    nothing
end

function make_jacobian(n, p, t)
    #=
    Constructs the chemical jacobian in the normal way, including stuff to calculate parameters for chemJmat.
    =#
    
    n[n .< 0] .= 0
   
    # Unpack the parameters ---------------------------------------------------------------
    inactive, inactivesp, activesp, Tn, Ti, Te = p

    Jrates, tup, tdown, tlower, tupper = get_transport_and_J_rates(n, inactive, activesp, inactivesp)
    
    return chemJmat(n, inactive, activesp, inactivesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
end

################################################################################
#                                  MAIN SETUP                                  #
################################################################################

# Set up simulation files and experiment type ==================================

# Note: directory paths are in PARAMETERS.jl
# get command line arguments for experiment type. format:
# <experiment type> <parameters> <solar cycle type>
# examples: 
# temp Tsurf Ttropo Texo mean; water <mixing ratio> mean; dh <multiplier> mean; 
# Oflux <cm^-2s^-1> mean
# examples: temp 190 110 200 mean; water 1e-3 mean; dh 8 mean; Oflux 1.2e8 mean.
# last argument is solar cycle: min, mean, or max.
args = Any["temp", "200", "200", "200", "mean"]#Any[ARGS[i] for i in 1:1:length(ARGS)]

# Establish a pattern for filenames. FNext = filename extension
# TODO: We probably won't be running all these experiments this time, so this section can probably go away.
# we may not even have to run from command line at all...
if args[1]=="temp" # TODO: this needs to handle the ion and electron temps.
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    FNext = "water_$(args[2])"
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
    DH = parse(Float64, args[2]) * 1.6e-4
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
else
    throw("Error! Bad experiment type")
end

user_input_folder_name = input("Enter a name for the results folder or press enter to use default: ")
sim_folder_name = user_input_folder_name == "" ? FNext : user_input_folder_name
# Set up the folder if it doesn't exist
create_folder(sim_folder_name, results_dir)

n_current = Dict{Symbol,Array{Float64,1}}()
n_current[:CO2] = fill(2e17, num_layers)

# This flattened array is used elsewhere, but never changes.
const inactive = flatten_atm(n_current, inactivespecies)

# Whether to initialize new species as zeros or not
use_nonzero_initial_profiles = true

# Set whether neutrals or ions are being converged. The other group must be entered in nochemspecies, notransportspecies.
converge_which = input("Converging ions, neutrals or both?: ")
while converge_which != "neutrals" && converge_which != "ions" && converge_which != "both"
    println("Bad entry! Please enter ions or neturals or both.")
    global converge_which = input("Converging ions, neutrals or both?: ")
end
println()


# turn chemistry on or off
do_chem = input("Allow chemistry? (on/off): ")
while do_chem != "on" && do_chem != "off" 
    println("Bad entry! Please enter on or off")
    global do_chem = input("Allow chemistry? (on/off): ")
end
do_chem = do_chem=="on" ? true : false
println()

do_trans = input("Allow transport? (on/off): ")
while do_trans != "on" && do_trans != "off" 
    println("Bad entry! Please enter on or off")
    global do_chem = input("Allow transport? (on/off): ")
end
do_trans = do_trans=="on" ? true : false
println()

# Zero out the new Jrates
for nj in newJrates
    n_current[nj] = zeros(num_layers)
end

# Set up the timesteps
mindt = dt_min_and_max[converge_which][1]
maxdt = dt_min_and_max[converge_which][2]
total_timesteps = abs(mindt) + maxdt + 1 # the 1 is to account for 10^0. Used in some plotting codes. TODO: Check if not used anymore.

# Set solar cycle file and alert user 
cycle = args[end]
solar_data_file = Dict("max"=>"marssolarphotonflux_solarmax.dat", 
                       "mean"=>"marssolarphotonflux_solarmean.dat", 
                       "min"=>"marssolarphotonflux_solarmin.dat")
solarfile = solar_data_file[cycle]

# Convert the arguments to numbers so we can use them to do maths
for i in 2:1:length(args)-1
    args[i] = parse(Float64, args[i])
end

# Let the user know what is being done 
if cycle != "mean"
    println("ALERT: Solar $(cycle) data being used")
end

# Plot styles ==================================================================
rcParams = PyCall.PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Louis George Caf?"]
rcParams["font.monospace"] = ["FreeMono"]
rcParams["font.size"] = 22
rcParams["axes.labelsize"]= 24
rcParams["xtick.labelsize"] = 22
rcParams["ytick.labelsize"] = 22

################################################################################
#                       TEMPERATURE/PRESSURE PROFILES                          #
################################################################################
println("Setting up temperature functions...")
# If changes to the temperature are needed, they should be made here 
if args[1]=="temp"
    global T_surf = args[2]
    global T_tropo = args[3]
    global T_exo = args[4]
    Temp_n(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "neutral")
    Temp_i(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "ion")
    Temp_e(z::Float64) = T_all(z, T_surf, T_tropo, T_exo, "electron")
    Temp_keepSVP(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral") # for testing temp without changing SVP. #TODO: adjust if needed for ions?
else
    global T_surf = meanTs
    global T_tropo = meanTt
    global T_exo = meanTe
    Temp_n(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "neutral")
    Temp_i(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "ion")
    Temp_e(z::Float64) = T_all(z, meanTs, meanTt, meanTe, "electron")
end

# set temperature arrays to be used everywhere
Tn_arr = fill(200., num_layers)
Ti_arr = Float64[Temp_i(a) for a in non_bdy_layers]
Te_arr = Float64[Temp_e(a) for a in non_bdy_layers]

# plot all 3 profiles on top of each other
# plot_temp_prof([Temp_n(a) for a in alt], savepath=results_dir*sim_folder_name, i_temps=[Temp_i(a) for a in alt], e_temps=[Temp_e(a) for a in alt])


################################################################################
#                               WATER PROFILES                                 #
################################################################################
println("Setting up the water profile...")
# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
else
    fix_SVP = false
    H2Osat = map(x->Psat(x), map(Temp_n, alt)) # array in #/cm^3 by altitude
    HDOsat = map(x->Psat_HDO(x), map(Temp_n, alt))
end

# H2O Water Profile ============================================================
surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
# H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo
# # set H2O SVP fraction to minimum for all alts above first time min is reached
# H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
# H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
#                fill(minimum(H2Osatfrac), num_layers-length(H2Oinitfrac))]

# make profile constant in the lower atmosphere (well-mixed).
# when doing water experiments, the temperature profile is such that the minimum
# in the SVP curve occurs at about 55 km alt. This was found by manual tweaking.
# thus we need to only set the lower atmo mixing ratio below that point, or 
# there will be a little spike in the water profile.
# if args[1] == "water"
#     H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= args[2]
#     MR = args[2] # mixing ratio 
# else
#     MR = MR_mean_water
#     H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR # 10 pr μm
# end

# for i in [1:length(H2Oinitfrac);]
#     H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
# end

MR = MR_mean_water

# set the water profiles =======================================================
# n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers)
# n_current[:HDO] = 2 * DH * n_current[:H2O] # This should be the correct way to set the HDO profile.
# n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), non_bdy_layers) # OLD WAY that is wrong.


# We still have to calculate the HDO initial fraction in order to calculate the pr um 
# and make water plots.
# HDOinitfrac = n_current[:HDO] ./ map(z->n_tot(n_current, z), non_bdy_layers)  

# ADD EXCESS WATER AS FOR DUST STORMS. TODO: revert
# H2Oppm = 1e-6*map(x->250 .* exp(-((x-42)/12.5)^2), non_bdy_layers/1e5) + H2Oinitfrac  # 250 ppm at 42 km (peak)
# HDOppm = 1e-6*map(x->0.350 .* exp(-((x-38)/12.5)^2), non_bdy_layers/1e5) + HDOinitfrac  # 350 ppb at 38 km (peak)
# n_current[:H2O] = H2Oppm .* map(z->n_tot(n_current, z), non_bdy_layers)
# n_current[:HDO] = HDOppm .* map(z->n_tot(n_current, z), non_bdy_layers)

# Compute total water column for logging and checking that we did things right =
# H2O #/cm^3 (whole atmosphere) = sum(MR * n_tot) for each alt
# H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
# HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))

# pr μm = (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23, 
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
# H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
# HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)


################################################################################
#                             BOUNDARY CONDITIONS                              #
################################################################################
println("Setting boundary conditions and creating metaprogramming...")
H_veff = effusion_velocity(Temp_n(zmax), 1.0, zmax)
H2_veff = effusion_velocity(Temp_n(zmax), 2.0, zmax)
D_veff = effusion_velocity(Temp_n(zmax), 2.0, zmax)
HD_veff = effusion_velocity(Temp_n(zmax), 3.0, zmax)

#=
    boundary conditions for each species (mostly from Nair 1994, Yung 1988). For 
    most species, default boundary condition is zero (f)lux at top and bottom. 
    Atomic/molecular hydrogen and deuterated analogues have a nonzero effusion 
    (v)elocity at the upper layer of the atmosphere. Some species have a (n)umber 
    density boundary condition.
=#

if args[1] == "Oflux"
    global const Of = args[2]
else
    global const Of = 1.2e8
end

# This has to be defined here, because it uses as a boundary condition the H2O 
# and HDO saturation at the surface.
global const speciesbclist=Dict(
                :CO2=>["n" 2.1e17; "f" 0.],
                :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                :HDO=>["n" HDOsat[1]; "f" 0.],
                :O=>["f" 0.; "f" Of],
                :H2=>["f" 0.; "v" H2_veff],
                :HD=>["f" 0.; "v" HD_veff],
                :H=>["f" 0.; "v" H_veff],
                :D=>["f" 0.; "v" D_veff],
                # TODO: Ion boundary conditions?
               );

################################################################################
#                       COMBINED CHEMISTRY AND TRANSPORT                       #
################################################################################

#=
    We now have objects that return the list of indices and coefficients
    for transport, assuming no other species in the atmosphere
    (transportmat), and for chemistry, assuming no other altitudes
    (chemical_jacobian). We need to perform a kind of outer product on
    these operators, to determine a fully coupled set of equations for
    all species at all altitudes.

    the rates at each altitude can be computed using the reaction network
    already in place, plus additional equations describing the transport
    to and from the cells above and below:
=#
upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
          Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
          for s in specieslist]

downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
            Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
            for s in specieslist]

local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in specieslist]
                          [Symbol("t"*string(s)*"_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_above_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_below_up") for s in specieslist]]...;]

transportnet = [[upeqns...;]; [downeqns...;]]

# define names for all the species active in the coupled rates:
# DO NOT CHANGE!!
const active_above = [Symbol(string(s)*"_above") for s in activespecies]
const active_below = [Symbol(string(s)*"_below") for s in activespecies]

# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x, chem_on=do_chem, trans_on=do_trans), activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

# TODO: Turn these notes into a test that looks for "+()" in rates_local instead of debugging lines. 
# -----------------------------------------------------------------------------------------------
# These lines are useful for troubleshooting if you're getting weird errors with ratefn.
# "no method matching +()" means the chemical system is unbalanced and you've got a
# species that has production but no consumption, or vice versa.
# other errors may occur. Uncomment these 3 lines to inspect what goes into rates_local.
# println("The contents of rates_local: ")
# println(rates_local)
# println("That was rates_local")


arglist_local = [activespecies; active_above; active_below; inactivespecies;
                 Jratelist; :Tn; :Ti; :Te; local_transport_rates]

arglist_local_typed = [:($s::Float64) for s in arglist_local]

# These expressions which are evaluated below enable a more accurate assessment of M and E values
# by calculating at the time of being called rather than only at each timestep.
Mexpr = Expr(:call, :+, fullspecieslist...)


# NOTE: These functions within @eval cannot be moved. Do not move them.
@eval begin
    function ratefn_local($(arglist_local_typed...))

        # M and E are calculated here to ensure that the right number of ions/electrons
        # is used. It is for only the altitude at which this function was called 
        # (i.e. all the arguments to the function, when it's called, are for only one altitude)
        M = $Mexpr
        $rates_local # evaluates the rates_local expression
    end
end

@eval begin
    function chemJmat_local($(arglist_local_typed...))

        M = $Mexpr

        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])

        # section to populate localchemJval element by element,
        # which is having problems with compiler errors when ----->
        # localchemJval = similar($chemJ_local[3], Float64)
        #     # unroll the loop and evaluate each expression 1 by 1 thank you 张实唯 of SE
        #     $((quote
        #         localchemJval[$i] = $(chemJ_local[3][i])
        #     end for i in 1:length(chemJ_local[3]))...)
        # localchemJval = -dt*localchemJval  # the old code that throws compiler error:*$(chemJ_local[3])  
        # <----- end section
        localchemJval = $(Expr(:vcat, chemJ_local[3]...)) 

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = $(Expr(:vcat, chemJ_above[3]...))

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = $(Expr(:vcat, chemJ_below[3]...))

        # return the actual values of I, J, and V (indices and numerical value):
        return ((localchemJi, localchemJj, localchemJval),
                (abovechemJi, abovechemJj, abovechemJval),
                (belowchemJi, belowchemJj, belowchemJval))
    end
end


################################################################################
#                         PHOTOCHEMICAL CROSS SECTIONS                         #
################################################################################
println("Populating cross section dictionary...")

crosssection = populate_xsect_dict([T_surf, T_tropo, T_exo])#, o3xdata)

# Solar Input ==================================================================

const solarflux=readdlm(research_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]
solarflux[:,2] = solarflux[:,2]/2  # To roughly put everything at an SZA=60° (from a Kras comment)

lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

# if !(setdiff(solarflux[:,1],lambdas)==[])
#     throw("Need a broader range of solar flux values!")
# end

# pad all cross-sections to solar
for j in Jratelist, ialt in 1:length(alt)
    crosssection[j][ialt] = padtosolar(solarflux, crosssection[j][ialt])
end

# this is the unitialized array for storing values
solarabs = fill(fill(0.,size(solarflux, 1)), num_layers);

################################################################################
#                                 LOGGING                                      #
################################################################################
println("Creating the simulation log file...")
# crosssection dict for logging purposes =======================================
xsect_dict = Dict("CO2"=>[co2file],
                  "H2O, HDO"=>[h2ofile, hdofile],
                  "H2O2, HDO2"=>[h2o2file, hdo2file],
                  "O3"=>[o3file, o3chapfile],
                  "O2"=>[o2file, o2_130_190, o2_190_280, o2_280_500],
                  "H2, HD"=>[h2file, hdfile],
                  "OH, OD"=>[ohfile, oho1dfile, odfile])

# Log temperature and water parameters =========================================
if args[1]=="temp"
    input_string = "T_s=$(args[2]), T_tropo=$(args[3]), T_exo=$(args[4])" *
                   "\nwater init=$(MR)\nDH=5.5 \nOflux=1.2e8"
elseif args[1]=="water"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\n" *
                   "water init=$(args[2])\nDH=5.5\nOflux=1.2e8\n"
elseif args[1]=="dh"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)\n" *
                   "DH=$(args[2]) \nOflux=1.2e8\n"
elseif args[1]=="Oflux"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)" *
                   "\nDH=5.5\nOflux=$(Of)\n"
end

# Write the log ================================================================
f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "w")
write(f, "$(args[1]) experiment: \n")
write(f, input_string)
write(f, "\n")

# Mean temperatures
write(f, "Mean temperatures used:\n")
write(f, "Surface: $(meanTs) K, Tropopause: $(meanTt) K, Exobase: $(meanTe) K\n\n")

# which species are turned on
write(f, "All species: $(join([string(i) for i in fullspecieslist], ", "))\n")
write(f, "No-chem species: $(join([string(i) for i in nochemspecies], ", "))\n")
write(f, "No-transport species: $(join([string(i) for i in notransportspecies], ", "))\n")
write(f, "Active species: $(join(sort([string(i) for i in activespecies]), ", "))\n\n")

# Which ion reactions are turned on 
write(f, "New J rates that are turned on: \n")
for rxn in [string(i) for i in newJrates]
    write(f, rxn*"\n")
end
write(f, "\n")

# gets just chemistry reactions with ions as reactants or products
ion_chem_rxns = filter(x->(occursin("pl", string(x[1])) || occursin("pl", string(x[2]))), filter(x->!occursin("J", string(x[3])), reactionnet)) 
write(f, "Ion chemistry reactions: \n")
for rxn in [string(i) for i in ion_chem_rxns]
    write(f, rxn*"\n")
end
write(f, "\n")

# Whether nonzero initial profiles were used
write(f, "Non-zero initial profiles used: $(use_nonzero_initial_profiles)\n\n")

# Solar cycle and SVP stuff
write(f, "\nSVP fixed: $(fix_SVP)\n\n")
write(f, "Solar cycle status: solar $(cycle)\n\n")

# Water profile information
# write(f, "Water profile information: \n")
# write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
# write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
# write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
# write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
# write(f, "Total H2O (pr μm): $(H2Oprum)\n")
# write(f, "Total HDO (pr μm): $(HDOprum)\n")
# write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)\n\n")


# cross sections
write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
write(f, "\n")

# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="),
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end
write(f, "\n")

# timestep etc

close(f)


function evolve_atmosphere(atm_init::Dict{Symbol, Array{Float64, 1}}, t_to_save) 
    #=
    Sets up the initial conditions for the simulation and calls the ODE solver. 
    t_to_save: timesteps at which to save a snapshot of the atmosphere.
    
    plotJratesflag: if set to true, it will plot the Jrates.
    =#

    # Set up the initial state and important simulation parameters
    nstart = flatten_atm(atm_init, activespecies)
    find_nonfinites(nstart, collec_name="nstart")
    tspan = (10.0^(mindt), 10.0^(maxdt))
    params = [inactive, inactivespecies, activespecies, Tn_arr, Ti_arr, Te_arr]

    dummy_jac = make_jacobian(nstart, params, tspan[1])
    find_nonfinites(dummy_jac, collec_name="dummy_jac")
    println("Sparsity of the chemical jacobian: $(round(length(dummy_jac.nzval)*100/(size(dummy_jac)[1] * size(dummy_jac)[2]), digits=2))%")
    # figure()
    # scatter(findnz(dummy_jac)[1], findnz(dummy_jac)[2], findnz(dummy_jac)[3])
    # title("Sparsity pattern")ratefn(n, inactive, inactivesp, activesp, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper)
    # show()
    
    f = ODEFunction(PnL_eqn, jac=make_jacobian_wrapper, jac_prototype=dummy_jac)
    prob = ODEProblem(f, nstart, tspan, params)
    println("Starting the solver...")
    sol = solve(prob, CVODE_BDF(linear_solver=:KLU), saveat=t_to_save, progress=false, dt=(10.0^(mindt-3)), save_start=true, abstol=1e-3, reltol=1e-3)
    return sol
end

################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# do the convergence ===========================================================
println("Beginning Convergence")
ti = time()
times_to_save = [10.0^t for t in mindt:maxdt]
atm_soln = evolve_atmosphere(n_current, times_to_save)
tf = time()
println("Finished convergence in $((tf-ti)/60) minutes")


# Save the results =============================================================

L = length(times_to_save)
println(L)
i = 1
# Write all states to individual files.
for (timestep, atm_state) in zip(atm_soln.t, atm_soln.u)
    # This is the contents of unflatten_atm.

    nc = unflatten_atm(atm_state, fullspecieslist)

    if i == L
        plot_atm(nc, [neutrallist, ionlist], results_dir*sim_folder_name*"/final_atmosphere.png")
    end
    global i += 1
    println(i)

    filepath = results_dir*sim_folder_name*"/atm_state_t_$(timestep).h5"
    write_ncurrent(nc, filepath)   
end

println("Wrote all states to files")

println()

f = open(results_dir*sim_folder_name*"/simulation_params_"*FNext*".txt", "a")

println("Finished")
close(f)