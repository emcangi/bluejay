###############################################################################
# converge_atmo_instantiate_ions.jl
# TYPE: (1) Model files - required
# DESCRIPTION: does initial convergence for a photochemistry experiment of the
# Martian atmosphere; this file adds in initial abundances for ions and converges
# the atmosphere so that they will take on their equilibrium profiles
#
# Eryn Cangi
# Created 2018
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
###############################################################################

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Photochemistry  # custom module for this project
# The following functions are imported directly in order to overload them, making
# multiple methods available.
import Photochemistry.fluxcoefs, Photochemistry.scaleH


include("PARAMETERS.jl")

################################################################################
#                                   FUNCTIONS                                  #
################################################################################

#=
These functions are required to be in this file for one of three reasons:
1) because they, or one of their overrides, call on the dynamically-defined
   Temp(z) function,
2) Because they call a function that is couched in an @eval statement, which
   cannot be relocated,
3) They are the main routine functions and will not be shared by any other scripts.

The Temp(z) function is defined during each run as a shortcut to Tpiecewise(),
which takes as arguments the temperatures at the surface, tropopause, and exobase.
Temp(z) is defined so that those arguments need not be passed constantly. This
may be fixed in the future.

For additional functions, see the Photochemistry module.
=#


# transport/scale height =======================================================

function scaleH(z, species::Symbol)
    #=
    Same as first scaleH, but for a particular atomic/molecular species.
    =#

    sptype = charge_type(species)
    if sptype == "ion"
        T = Temp_i(z)
    elseif sptype == "electron"
        T = Temp_e(z)
    elseif sptype == "neutral"
        T = Temp_n(z)
    end

    mm = speciesmolmasslist[species]
    return scaleH(z, T, mm)
end

# transport functions ==========================================================

#=
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

function fluxcoefs(z, dz, species, n_current)
    #=
    overload to generate the coefficients K, D, T, Hs if they are not supplied (most common)

    z: altitude in cm
    dz: altitudinal layer thickness in cm
    species: species symbol for which flux coefficients are calculated
    n_current: current atmospheric state array

    p: upper layer ("plus")
    0: this layer
    m: lower layer ("minus")
    =#

    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    if species_type == "ion"
        Tp = Temp_i(z+dz)
        T0 = Temp_i(z)
        Tm = Temp_i(z-dz)
    elseif species_type == "electron"
        Tp = Temp_e(z+dz)
        T0 = Temp_e(z)
        Tm = Temp_e(z-dz)
    elseif species_type == "neutral"
        Tp = Temp_n(z+dz)
        T0 = Temp_n(z)
        Tm = Temp_n(z-dz)
    end

    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)

    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz, [Km , K0, Kp], [Dm , D0, Dp], [Tm , T0, Tp],
                     [Hsm, Hs0, Hsp], [H0m, H00, H0p], species)
end

function lower_up(z, dz, species, n_current)
    #=
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one above.
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus")

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    n_current: Array; species number density by altitude

    returns: return of fluxcoefs
    =#

    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    if species_type == "ion"
        Tp = Temp_i(z+dz)
        T0 = Temp_i(z)
        Tm = 1
    elseif species_type == "electron"
        Tp = Temp_e(z+dz)
        T0 = Temp_e(z)
        Tm = 1
    elseif species_type == "neutral"
        Tp = Temp_n(z+dz)
        T0 = Temp_n(z)
        Tm = 1
    end

    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = 1
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z,nt0)
    Km = 1

    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = 1
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z, species)
    Hsm = 1
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
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

function upper_down(z, dz, species, n_current)
    #=
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one below.
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus")

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    n_current: Array; species number density by altitude

    returns: return of fluxcoefs
    =#

    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    if species_type == "ion"
        Tp = 1
        T0 = Temp_i(z)
        Tm = Temp_i(z-dz)
    elseif species_type == "electron"
        Tp = 1
        T0 = Temp_e(z)
        Tm = Temp_e(z-dz)
    elseif species_type == "neutral"
        Tp = 1
        T0 = Temp_n(z)
        Tm = Temp_n(z-dz)
    end

    ntp = 1
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = 1
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    # Tp = 1
    # T0 = Temp(z)
    # Tm = Temp(z-dz)
    Dp = 1
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = 1
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = 1
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[1]
end

function boundaryconditions(species, speciesbclist, dz, n_current)
    #=
    returns the symbolic transport coefficients that encode the
    boundary conditions for the null-pointing equations

    n_1->NULL t_lower_bc
    n_f->NULL t_upper_bc

    this defines two additional symbols for each species that need
    to be resolved in the function call macro:
                tspecies_lower_up and
                tspecies_upper_down
    these are found by passing the appropriate values to fluxcoefs
    and selecting the correct output.

    species: Symbol
    speciesbclist: Has to be passed in to pass on to speciesbcs. Because water
                   saturation at the surface can vary, but is a bc.
    dz: Float64; layer thickness in cm
    n_current: Array; species number density by altitude

    returns: 2x2 boundary condition array where the first row is for the surface
             layer and second row is for the top of the atmosphere.
    =#

    bcs = speciesbcs(species, speciesbclist)
    if issubset([species],notransportspecies)
        bcs = ["f" 0.; "f" 0.]
    end

    # first element returned corresponds to lower BC, second to upper
    # BC transport rate. Within each element, the two rates correspond
    # to the two equations
    # n_b  -> NULL (first rate, depends on species concentration)
    # NULL -> n_b  (second rate, independent of species concentration
    bcvec = Float64[0 0;0 0]

    # LOWER
    if bcs[1, 1] == "n"
        bcvec[1,:]=[fluxcoefs(alt[2], dz, species, n_current)[1],
                    lower_up(alt[1], dz, species, n_current)*bcs[1, 2]]
    elseif bcs[1, 1] == "f"
        bcvec[1,:] = [0.0, bcs[1, 2]/dz]
    elseif bcs[1, 1] == "v"
        bcvec[1,:] = [bcs[1, 2]/dz, 0.0]
    else
        throw("Improper lower boundary condition!")
    end

    # UPPER
    if bcs[2, 1] == "n"
        bcvec[2,:] = [fluxcoefs(alt[end-1],dz, species, n_current)[2],
                    upper_down(alt[end],dz, species, n_current)*bcs[1, 2]]
    elseif bcs[2, 1] == "f"
            bcvec[2,:] = [0.0,-bcs[2, 2]/dz]
    elseif bcs[2, 1] == "v"
        bcvec[2,:] = [bcs[2, 2]/dz, 0.0]
    else
        throw("Improper upper boundary condition!")
    end

    return bcvec
end

# chemistry functions ==========================================================

function reactionrates(n_current)
    #=
    Creates an array of size length(intaltgrid) x (number of reactions).
    Populated with chemical reaction rates for each reaction based on species
    populations.
    =#
    
    theserates = fill(convert(Float64, NaN), (length(intaltgrid), length(reactionnet)))
    for ialt in 1:length(intaltgrid)
        theserates[ialt,:] = reactionrates_local([n_current[sp][ialt] for sp in specieslist], 
                                                 [n_current[J][ialt] for J in Jratelist], 
                                                 Temp_n(alt[ialt+1]), Temp_i(alt[ialt+1]), Temp_e(alt[ialt+1]))  # check that this works

                                                #([[n_current[sp][ialt] for sp in specieslist];
                                                #[n_current[J][ialt] for J in Jratelist];
                                                #Temp_n(alt[ialt+1]); Temp_i(alt[ialt+1]); Temp_e(alt[ialt+1])]...)
                                                #n_tot(n_current, alt[ialt+1])]...)   # this is the old argument to reactionrates_local 
    end
    return theserates
end

function ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper) # DEBUG: add M after Te
    #=
    at each altitude, get the appropriate group of concentrations,
    coefficients, and rates to pass to ratefn_local
    =#
    println("ratefn has been called")
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive,(length(inactivespecies),length(intaltgrid)))
    returnrates = zero(nthismat)

    println("The argument to ratefn_local:")
    thing1 = [nthismat[:,1]; nthismat[:,2]; fill(1.0, length(activespecies)); inactivemat[:,1]; Jrates[:,1]; Tn[1]; Ti[1]; Te[1]; # M[1]; E[1]; # DEBUG: replace if needed
              tup[:,1]; tlower[:,1]; tdown[:,2]; tlower[:,2]]

    println(thing1)
    println("Length of argument to ratefn_local: $(length(thing1))")
    println()
    # fill the first altitude entry with information for all species
    returnrates[:,1] = ratefn_local([nthismat[:,1]; nthismat[:,2];
                                    fill(1.0, length(activespecies));
                                    inactivemat[:,1]; Jrates[:,1]; Tn[1]; Ti[1]; Te[1]; # M[1]; E[1]; # DEBUG: replace if needed
                                    tup[:,1]; tlower[:,1]; tdown[:,2];
                                    tlower[:,2]])#...) # DEBUG: replace splat after last ] if need be.


    # iterate through other altitudes except the last level, filling the info in
    for ialt in 2:(length(intaltgrid)-1)
        returnrates[:,ialt] = ratefn_local([nthismat[:,ialt];
                                          nthismat[:,ialt+1];
                                          nthismat[:,ialt-1];
                                          inactivemat[:,ialt];
                                          Jrates[:,ialt];
                                          Tn[ialt]; Ti[ialt]; Te[ialt]; # M[ialt]; E[ialt];# DEBUG: replace if needed
                                          tup[:,ialt]; tdown[:,ialt];
                                          tdown[:,ialt+1]; tup[:,ialt-1]])#...)# DEBUG: replace splat after last ] if need be.
    end

    # fill in the last level of altitude (200 km)
    returnrates[:,end] = ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       Tn[end]; Ti[end]; Te[end]; #M[end]; E[end]; # DEBUG: replace if needed
                                       tupper[:,1]; tdown[:,end];
                                       tupper[:,2]; tup[:,end-1]])#...) # DEBUG: replace splat if problem
    return [returnrates...;]
end

function chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, tup, tdown, tlower, tupper, dt) # DEBUG: M and E after Te if need be
    #=
    TODO: docstring
    =#

    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive, (length(inactivespecies), length(intaltgrid)))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]

    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1]; nthismat[:,2];
                                                  fill(1.0, length(activespecies));
                                                  inactivemat[:,1]; Jrates[:,1];
                                                  Tn[1]; Ti[1]; Te[1]; #M[1]; E[1]; 
                                                  tup[:,1]; tlower[:,1];
                                                  tdown[:,2]; tlower[:,2];dt])#...) # DEBUG: splat
    #add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])
    #and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2] .+ length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(length(intaltgrid)-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt]; Tn[ialt]; Ti[ialt]; Te[ialt];
                                                      #M[ialt]; E[ialt]; 
                                                      tup[:,ialt];
                                                      tdown[:,ialt];
                                                      tdown[:,ialt+1];
                                                      tup[:,ialt-1]; dt])#...) # DEBUG: splat
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
                                              Tn[end]; Ti[end]; Te[end]; #M[end]; E[end];
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1]; dt])#...) # DEBUG: splat
    # add the influence of the local densities
    append!(chemJi, tclocal[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJval, tclocal[3])
    # and the lower densities
    append!(chemJi, tclower[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclower[2].+(length(intaltgrid)-2)*length(activespecies))
    append!(chemJval, tclower[3])

    # make sure to add 1's along the diagonal
    append!(chemJi,[1:length(nthis);])
    append!(chemJj,[1:length(nthis);])
    append!(chemJval, fill(1.0, length(nthis)))

    sparse(chemJi, chemJj, chemJval, length(nthis), length(nthis), +);
end

# main routine functions =======================================================
function update_Jrates!(n_current::Dict{Symbol, Array{Float64, 1}})
    #=
    this function updates the photolysis rates stored in n_current to
    reflect the altitude distribution of absorbing species
    =#

    # Initialize an array, length=num of altitude levels - 2.
    # Each sub-array is an array of length 2000, corresponding to 2000 wavelengths.
    solarabs = Array{Array{Float64}}(undef, length(alt)-2)
    for i in range(1, length=length(alt)-2)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_current[species][ialt]*dz
            # if jspecies==:JO2toOpO
            #     println(string("At alt = ",alt[ialt+1],
            #                    ", n_",species," = ",n_current[species][ialt],
            #                    ", jcolumn = ",jcolumn))
            #     println("and solarabs[ialt] is $(solarabs[ialt]) before we do axpy")
            #     readline(STDIN)
            # end

            # add the total extinction to solarabs:
            # multiplies air column density at all wavelengths by crosssection
            # to get optical depth. This is an override of axpy! to use the
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

    #solarabs now records the total optical depth of the atmosphere at
    #each wavelength and altitude

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
        for ialt in [1:nalt;]
            n_current[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1,
                                          crosssection[j][ialt+1], 1)
        end
       # this section for testing sensitivity to J rates
        # if contains(string(j), "H2O") | contains(string(j), "HDO")
        # if contains(string(j), "CO2toCOpO")
           # n_current[j] = n_current[j] ./ 10
        # end
    end
end

function timeupdate(mytime, imgpath)
    println("Entering timeupdate")
    if mytime > 1000
        numiters = 25
    else
        numiters = 15
    end

    printE = Dict(15=>true, 25=>true)

    for i = 1:numiters
        plotatm(n_current, [neutrallist, ionlist], imgpath, t=mytime, iter=i)
        # write_ncurrent(n_current, results_dir*"ncurrent_$(mytime)_$(i).h5") #debug
        # println("dt: $(mytime)")

        printEflag = get(printE, i, false)
        update!(n_current, mytime, printEflag) # TODO: remove this weird flag when done checking electrons
    end
    # show()
    # yield() # this is to make the plot appear in the window and refresh
end

function next_timestep(nstart::Array{Float64, 1}, nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1}, activespecies, inactivespecies,
                       Jrates::Array{Float64, 2},
                       Tn::Array{Float64, 1}, Ti::Array{Float64, 1}, Te::Array{Float64, 1},
                       #M::Array{Float64, 1}, E::Array{Float64, 1},
                       tup::Array{Float64, 2}, tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2}, tupper::Array{Float64, 2},
                       dt::Float64)
    #=
    moves to the next timestep using Newton's method on the linearized
    coupled transport and chemical reaction network.
    =#
    println("in next_timestep")
    eps = 1.0 # ensure at least one iteration
    iter = 0
    while eps>1e-8 
        nold = deepcopy(nthis)

        # stuff concentrations into update function and jacobian.
        println("attempting fval calculation")
        fval = nthis - nstart - dt*ratefn(nthis, inactive, inactivespecies, activespecies, Jrates, Tn, Ti, Te, #M, E, 
                                          tup, tdown, tlower, tupper)
        
        println("Now we will break because updatemat isn't defined")
        error("Updatemat currently commented out to test code")
        # updatemat = chemJmat(nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, #M, E, 
        #                      tup, tdown, tlower, tupper, dt)

        # update
        nthis = nthis - (updatemat \ fval) # backslash - matrix solution operator
        # check relative size of update
        eps = maximum(abs.(nthis-nold)./nold)
        iter += 1
        if iter>1e3
            throw("too many iterations in next_timestep! eps: $(eps)")
        end
    end

    return nthis
end

function update!(n_current::Dict{Symbol, Array{Float64, 1}}, dt, printEflag)
    #=
    update n_current using the coupled reaction network, moving to
    the next timestep
    =#

    println("entered update!")
    #set auxiliary (not solved for in chemistry) species values, photolysis rates
    inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(intaltgrid)]...])
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(intaltgrid)])

    # extract concentrations
    nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(intaltgrid)]...])
    # M = sum([n_current[sp] for sp in fullspecieslist])
    # E = sum([n_current[sp] for sp in ionlist])  # Electron density calculated as the density of total ions to simplify things.
    if printEflag
        # write(Edensityfile, "electron density at time $(dt) on last iter: $(E)\n")
        plot_Jrates(n_current, 216.0, "surf", "time=$(dt)")
    end

    # set temperatures
    Tn = Float64[Temp_n(a) for a in alt[2:end-1]]
    Ti = Float64[Temp_i(a) for a in alt[2:end-1]]
    Te = Float64[Temp_e(a) for a in alt[2:end-1]]

    # take initial guess
    nthis = deepcopy(nstart)

    # these are the sum of the transport flux coefficients D+K, divided by Δz², units 1/s
    tup = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[2] for sp in specieslist, a in alt[2:end-1]]
    tdown = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[1] for sp in specieslist, a in alt[2:end-1]]

    # put the lower layer and upper layer boundary conditions in separate arrays; but they are not the
    # right shape!
    tlower_temporary = [boundaryconditions(sp, speciesbclist, dz, n_current)[1,:] for sp in specieslist]
    tupper_temporary = [boundaryconditions(sp, speciesbclist, dz, n_current)[2,:] for sp in specieslist]

    # reshape tlower and tupper into 2x2 arrays
    tlower = zeros(Float64, length(tlower_temporary), 2)
    tupper = zeros(Float64, length(tupper_temporary), 2)

    # tlower_temporary & tupper_temporary have same length; OK to use lower for the range
    for r in range(1, length=length(tlower_temporary))
        tlower[r, :] = tlower_temporary[r]
        tupper[r, :] = tupper_temporary[r]
    end

    println("calling next_timestep")
    nthis = next_timestep(nstart, nthis, inactive, activespecies, inactivespecies, Jrates, Tn, Ti, Te, #M, E, 
                          tup, tdown, tlower, tupper, dt)
    nthismat = reshape(nthis,(length(activespecies),length(intaltgrid)))

    # write found values out to n_current
    for s in 1:length(activespecies)
        for ia in 1:length(intaltgrid)
            tn = nthismat[s, ia]
            n_current[activespecies[s]][ia] = tn > 0. ? tn : 0.
        end
    end

    update_Jrates!(n_current)
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
args = Any[ARGS[i] for i in 1:1:length(ARGS)]

# Establish a pattern for filenames. FNext = filename extension
# TODO: We probably won't be running all these experiments this time, so this section can probably go away.
# we may not even have to run from command line at all...
if args[1]=="temp"
    FNext = "temp_$(args[2])_$(args[3])_$(args[4])"
elseif args[1]=="water"
    FNext = "water_$(args[2])"
elseif args[1]=="dh"
    FNext = "dh_$(args[2])"
elseif args[1]=="Oflux"
    FNext = "Oflux_$(args[2])"
else
    throw("Error! Bad experiment type")
end

# Set up the folder if it doesn't exist
create_folder(FNext, results_dir)

# Case where this file will be used to converge an atmosphere of a different extent
# make_new_alt_grid = input("Would you like to use the script to converge a new atmosphere of a different extent? (y/n): ")
make_new_alt_grid = "n" # TODO: revert
if make_new_alt_grid=="y"
    readfile = research_dir*"converged_200km_atmosphere.h5"
    const alt = convert(Array, (0:2e5:200e5))
    n_current = get_ncurrent(readfile)

    new_zmax = parse(Int64, input("Enter the new top of the atmosphere in km: "))
    extra_entries = Int64((new_zmax - 200)/(dz/1e5))

    # Extend the grid
    for (k,v) in zip(keys(n_current), values(n_current))
       append!(v, fill(v[end], extra_entries))  # repeats the last value in the array for the upper atmo as an initial value.
    end

    const alt = convert(Array, (0:2e5:new_zmax*10^5))

    if new_zmax != 250
        println("Warning: You entered $(new_zmax) for the new max altitude but
                 250 is hard-coded in the PARAMETERS.jl file. I haven't made this
                 general yet. So the code is probably about to break")
    end
elseif make_new_alt_grid=="n"
    # Set up the converged file to read from and load the simulation state at init.
    file_to_use = input("Enter the name of a file containing a converged, 250 km atmosphere to use (press enter to use default): ")
    readfile = file_to_use == "" ? "converged_250km_atmosphere.h5" : file_to_use
    n_current = get_ncurrent(readfile)
else
    error("Didn't understand response")
end

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
println("ALERT: running sim for $(FNext)")
println("ALERT: Using file: ", readfile)
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
#                             ADD DEUTERATED SPECIES                           #
################################################################################

if args[1] == "dh"
    DH = parse(Float64, args[2]) * 1.6e-4
end

# Set all existing neutrals to zero to let simulation build up  # DEBUG
# for k in keys(n_current)
#     if k!=:CO2 && k!=:N2 && k!=:Ar
#         n_current[k] .= 0
#     end
# end


# modify n_current with deuterated species profiles ============================
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
# Commented out 9 December after successfully converging a file with these new neutrals but not ions.
# (probably, anyway, I didn't notice I hadn't entered a comment of when I commented these out)
# n_current[:HDO] = n_current[:H2O] * DH
# n_current[:OD] = n_current[:OH] * DH
# n_current[:HDO2] = n_current[:H2O2] * DH
# n_current[:D] = n_current[:H] * DH
# n_current[:DO2] = n_current[:HO2] * DH
# n_current[:HD] = n_current[:H2] * DH
# n_current[:DOCO] = n_current[:HOCO] * DH

# new neutrals from Roger's model (this screwy code loads the initial profile from a
# file and then puts it in the right shape because Julia is stupid about arrays)
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
# Commented out 9 December after successfully converging a file with these new neutrals but not ions.
# n_current[:C] = reshape(readdlm("../Resources/initial_profiles/C_initial_profile.txt", '\r'), (124,))
# n_current[:CH] = reshape(readdlm("../Resources/initial_profiles/CH_initial_profile.txt", '\r'), (124,))
# n_current[:CN] = reshape(readdlm("../Resources/initial_profiles/CN_initial_profile.txt", '\r'), (124,))
# n_current[:HCN] = reshape(readdlm("../Resources/initial_profiles/HCN_initial_profile.txt", '\r'), (124,))
# n_current[:HCO] = reshape(readdlm("../Resources/initial_profiles/HCO_initial_profile.txt", '\r'), (124,))
# n_current[:HNO] = reshape(readdlm("../Resources/initial_profiles/HNO_initial_profile.txt", '\r'), (124,))
# n_current[:N] = reshape(readdlm("../Resources/initial_profiles/N_initial_profile.txt", '\r'), (124,))
# n_current[:NH] = reshape(readdlm("../Resources/initial_profiles/NH_initial_profile.txt", '\r'), (124,))
# n_current[:NH2] = reshape(readdlm("../Resources/initial_profiles/NH2_initial_profile.txt", '\r'), (124,))
# n_current[:N2O] = reshape(readdlm("../Resources/initial_profiles/N2O_initial_profile.txt", '\r'), (124,))
# n_current[:NO] = reshape(readdlm("../Resources/initial_profiles/NO_initial_profile.txt", '\r'), (124,))
# n_current[:NO2] = reshape(readdlm("../Resources/initial_profiles/NO2_initial_profile.txt", '\r'), (124,))

# new ions from Roger's model
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
n_current[:Arpl] = reshape(readdlm("../Resources/initial_profiles/Arpl_initial_profile.txt", '\r'), (124,))
n_current[:ArHpl] = reshape(readdlm("../Resources/initial_profiles/ArHpl_initial_profile.txt", '\r'), (124,))
n_current[:Cpl] = reshape(readdlm("../Resources/initial_profiles/Cpl_initial_profile.txt", '\r'), (124,))
n_current[:CHpl] = reshape(readdlm("../Resources/initial_profiles/CHpl_initial_profile.txt", '\r'), (124,))
n_current[:CNpl] = reshape(readdlm("../Resources/initial_profiles/CNpl_initial_profile.txt", '\r'), (124,))
n_current[:COpl] = reshape(readdlm("../Resources/initial_profiles/COpl_initial_profile.txt", '\r'), (124,))
n_current[:CO2pl] = reshape(readdlm("../Resources/initial_profiles/CO2pl_initial_profile.txt", '\r'), (124,)) # Nair minimal ionosphere
n_current[:Hpl] = reshape(readdlm("../Resources/initial_profiles/Hpl_initial_profile.txt", '\r'), (124,))
n_current[:H2pl] = reshape(readdlm("../Resources/initial_profiles/H2pl_initial_profile.txt", '\r'), (124,))
n_current[:H3pl] = reshape(readdlm("../Resources/initial_profiles/H3pl_initial_profile.txt", '\r'), (124,))
n_current[:H2Opl] = reshape(readdlm("../Resources/initial_profiles/H2Opl_initial_profile.txt", '\r'), (124,))
n_current[:H3Opl] = reshape(readdlm("../Resources/initial_profiles/H3Opl_initial_profile.txt", '\r'), (124,))
n_current[:HCNpl] = reshape(readdlm("../Resources/initial_profiles/HCNpl_initial_profile.txt", '\r'), (124,))
n_current[:HCNHpl] = reshape(readdlm("../Resources/initial_profiles/HCNHpl_initial_profile.txt", '\r'), (124,))
n_current[:HCOpl] = reshape(readdlm("../Resources/initial_profiles/HCOpl_initial_profile.txt", '\r'), (124,))
n_current[:HCO2pl] = reshape(readdlm("../Resources/initial_profiles/HCO2pl_initial_profile.txt", '\r'), (124,)) # Nair minimal ionosphere
n_current[:HNOpl] = reshape(readdlm("../Resources/initial_profiles/HNOpl_initial_profile.txt", '\r'), (124,))
n_current[:HN2Opl] = reshape(readdlm("../Resources/initial_profiles/HN2Opl_initial_profile.txt", '\r'), (124,))
n_current[:HOCpl] = reshape(readdlm("../Resources/initial_profiles/HOCpl_initial_profile.txt", '\r'), (124,))
n_current[:HO2pl] = reshape(readdlm("../Resources/initial_profiles/HO2pl_initial_profile.txt", '\r'), (124,))
n_current[:Npl] = reshape(readdlm("../Resources/initial_profiles/Npl_initial_profile.txt", '\r'), (124,))
n_current[:N2pl] = reshape(readdlm("../Resources/initial_profiles/N2pl_initial_profile.txt", '\r'), (124,))
n_current[:NHpl] = reshape(readdlm("../Resources/initial_profiles/NHpl_initial_profile.txt", '\r'), (124,))
n_current[:NH2pl] = reshape(readdlm("../Resources/initial_profiles/NH2pl_initial_profile.txt", '\r'), (124,))
n_current[:NH3pl] = reshape(readdlm("../Resources/initial_profiles/NH3pl_initial_profile.txt", '\r'), (124,))
n_current[:N2Hpl] = reshape(readdlm("../Resources/initial_profiles/N2Hpl_initial_profile.txt", '\r'), (124,))
n_current[:N2Opl] = reshape(readdlm("../Resources/initial_profiles/N2Opl_initial_profile.txt", '\r'), (124,))
n_current[:NOpl] = reshape(readdlm("../Resources/initial_profiles/NOpl_initial_profile.txt", '\r'), (124,))
n_current[:NO2pl] = reshape(readdlm("../Resources/initial_profiles/NO2pl_initial_profile.txt", '\r'), (124,))
n_current[:Opl] = reshape(readdlm("../Resources/initial_profiles/Opl_initial_profile.txt", '\r'), (124,)) # Nair minimal ionosphere
n_current[:O2pl] = reshape(readdlm("../Resources/initial_profiles/O2pl_initial_profile.txt", '\r'), (124,)) # Nair minimal ionosphere
n_current[:OHpl] = reshape(readdlm("../Resources/initial_profiles/OHpl_initial_profile.txt", '\r'), (124,))

# new D neutrals: TODO: Find?

# new D ions:    # TODO: turn this on after testing that Roger's network works.
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
# n_current[:ArDpl] = n_current[:ArHpl] * DH
# n_current[:CDpl] = n_current[:CHpl] * DH
# n_current[:HDpl] = n_current[:H2pl] * DH
# n_current[:H2Dpl] = n_current[:H3pl] * DH
# n_current[:HD2pl] = n_current[:H3pl] * DH
# n_current[:DCOpl] = n_current[:HCOpl] * DH
# n_current[:Dpl] = n_current[:Hpl] * DH
# n_current[:N_2Dpl] = n_current[:N2Hpl] * DH
# n_current[:NDpl] = n_current[:NHpl] * DH
# n_current[:ODpl] = n_current[:OHpl] * DH
# n_current[:NHDpl] = zeros(length(alt))
# n_current[:NHD2pl] = zeros(length(alt))
# n_current[:ND3pl] = zeros(length(alt))
# n_current[:DOCOpl] = zeros(length(alt))


# add the new Jrates --the values will get calculated ================================
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
# Commented out 22 December 2020 after realizing they were still turned on and probably shouldn't be. 
# n_current[:JHDOtoHpOD] = zeros(length(alt)-2)
# n_current[:JHDOtoDpOH] = zeros(length(alt)-2)
# n_current[:JHDO2toOHpOD] = zeros(length(alt)-2)
# n_current[:JHDOtoHDpO1D] = zeros(length(alt)-2)
# n_current[:JHDOtoHpDpO] = zeros(length(alt)-2)
# n_current[:JODtoOpD] = zeros(length(alt)-2)
# n_current[:JHDtoHpD] = zeros(length(alt)-2)
# n_current[:JDO2toODpO] = zeros(length(alt)-2)
# n_current[:JHDO2toDO2pH] = zeros(length(alt)-2)
# n_current[:JHDO2toHO2pD] = zeros(length(alt)-2)
# n_current[:JHDO2toHDOpO1D] = zeros(length(alt)-2)
# n_current[:JODtoO1DpD] = zeros(length(alt)-2)

# NEW: photodissociation and photoionization reactions from Roger to be initialized.
# INITIAL CONVERGENCE ONLY: remove after converging the first equilibrated atmosphere.
n_current[:JCO2toCO2pl] = zeros(length(alt)-2) # Nair minimal ionosphere
n_current[:JCO2toOplpCO] = zeros(length(alt)-2) # Nair minimal ionosphere
n_current[:JCO2toCO2plpl] = zeros(length(alt)-2)
n_current[:JCO2toCpO2] = zeros(length(alt)-2)
n_current[:JCO2toCplplpO2] = zeros(length(alt)-2)
n_current[:JCO2toCOplpOpl] = zeros(length(alt)-2)
n_current[:JCO2toOplpCplpO] = zeros(length(alt)-2)
n_current[:JCO2toCOplpO] = zeros(length(alt)-2)
n_current[:JCO2toCpOpO] = zeros(length(alt)-2)
n_current[:JCO2toCplpO2] = zeros(length(alt)-2)

n_current[:JCOtoCpO] = zeros(length(alt)-2)
n_current[:JCOtoCOpl] = zeros(length(alt)-2)
n_current[:JCOtoCpOpl] = zeros(length(alt)-2)
n_current[:JCOtoOpCpl] = zeros(length(alt)-2)

n_current[:JHtoHpl] = zeros(length(alt)-2)
n_current[:JH2toH2pl] = zeros(length(alt)-2)
n_current[:JH2toHplpH] = zeros(length(alt)-2)
n_current[:JH2OtoH2Opl] = zeros(length(alt)-2) 
n_current[:JH2OtoOplpH2] = zeros(length(alt)-2)
n_current[:JH2OtoHplpOH] = zeros(length(alt)-2)
n_current[:JH2OtoOHplpH] = zeros(length(alt)-2)
n_current[:JH2O2toH2O2pl] = zeros(length(alt)-2)

n_current[:JN2toN2pl] = zeros(length(alt)-2)
n_current[:JN2toNplpN] = zeros(length(alt)-2)
n_current[:JN2OtoN2pO1D] = zeros(length(alt)-2)
n_current[:JN2OtoN2Opl] = zeros(length(alt)-2)
n_current[:JNOtoNOpl] = zeros(length(alt)-2)
n_current[:JNOtoNpO] = zeros(length(alt)-2)
n_current[:JNO2toNOpO] = zeros(length(alt)-2)
n_current[:JNO2toNO2pl] = zeros(length(alt)-2)

n_current[:JOtoOpl] = zeros(length(alt)-2) # Nair minimal ionosphere
n_current[:JO2toO2pl] = zeros(length(alt)-2) # Nair minimal ionosphere
n_current[:JO3toO3pl] = zeros(length(alt)-2) # Nair minimal ionosphere

################################################################################
#                       TEMPERATURE/PRESSURE PROFILES                          #
################################################################################

println("As a safety check, here are the mean temperatures as entered in code:")
println("Surface: $(meanTs), Tropopause: $(meanTt), Exobase: $(meanTe)")

#If changes to the temperature are needed, they should be made here
if args[1]=="temp"
    Temp_n(z::Float64) = T(z, args[2], args[3], args[4], "neutral")
    Temp_i(z::Float64) = T(z, args[2], args[3], args[4], "ion")
    Temp_e(z::Float64) = T(z, args[2], args[3], args[4], "electron")
    Temp_keepSVP(z::Float64) = T(z, meanTs, meanTt, meanTe, "neutral") # for testing temp without changing SVP. #TODO: adjust if needed for ions?
else
    Temp_n(z::Float64) = T(z, meanTs, meanTt, meanTe, "neutral")
    Temp_i(z::Float64) = T(z, meanTs, meanTt, meanTe, "ion")
    Temp_e(z::Float64) = T(z, meanTs, meanTt, meanTe, "electron")
end

# plot all 3 profiles on top of each other
plot_nie_temp_profs([Temp_n(a) for a in alt], [Temp_i(a) for a in alt], [Temp_e(a) for a in alt], results_dir*FNext, alt)

################################################################################
#                               WATER PROFILES                                 #
################################################################################

# set SVP to be fixed or variable with temperature
if args[1] == "temp"
    fix_SVP = true
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed
    println("ALERT: SVP will be fixed to that of the mean temperature profile ")
else
    fix_SVP = false
    H2Osat = map(x->Psat(x), map(Temp_n, alt)) # array in #/cm^3 by altitude
    HDOsat = map(x->Psat_HDO(x), map(Temp_n, alt))
end

# H2O Water Profile ============================================================
surface_watersat = Dict("H2O"=>H2Osat[1], "HDO"=>HDOsat[1])
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z), alt)  # get SVP as fraction of total atmo
# set H2O SVP fraction to minimum for all alts above first time min is reached
H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
               fill(minimum(H2Osatfrac), length(alt)-2-length(H2Oinitfrac))]

# make profile constant in the lower atmosphere (well-mixed).
# when doing water experiments, the temperature profile is such that the minimum
# in the SVP curve occurs at about 55 km alt. This was found by manual tweaking.
# thus we need to only set the lower atmo mixing ratio below that point, or
# there will be a little spike in the water profile.
if args[1] == "water"
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= args[2]
    MR = args[2] # mixing ratio
else  # i believe 30 km is supposed to be approximately the hygropause.
    MR = MR_mean_water
    H2Oinitfrac[findall(x->x<hygropause_alt, alt)] .= MR # 10 pr μm
end

for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
end

# HDO water profile ============================================================
HDOsatfrac = HDOsat./map(z->n_tot(n_current, z), alt)
# use D/H ratio to set population of HDO
HDOinitfrac = H2Oinitfrac * DH  # initial profile for HDO

# Compute total water column in pr μm starting with mixing ratio array =========
# H2O #/cm^3 = sum(MR * n_tot) for each alt
H2O_per_cc = sum([MR; H2Oinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))
HDO_per_cc = sum([MR*DH; HDOinitfrac] .* map(z->n_tot(n_current, z), alt[1:end-1]))

# pr μm = (H2O #/cm^3) * cm * (mol/#) * (H2O g/mol) * (1 cm^3/g) * (10^4 μm/cm)
# where the lone cm is a slice of the atmosphere of thickness dz, #/mol=6.023e23,
# H2O or HDO g/mol = 18 or 19, cm^3/g = 1 or 19/18 for H2O or HDO.
# written as conversion factors for clarity.
H2Oprum = (H2O_per_cc * dz) * (18/1) * (1/6.02e23) * (1/1) * (1e4/1)
HDOprum = (HDO_per_cc * dz) * (19/1) * (1/6.02e23) * (19/18) * (1e4/1)

# Write out water content to a file ============================================
f = open(results_dir*FNext*"/water_column_"*FNext*".txt", "w")
write(f, "Total H2O col: $(H2O_per_cc*2e5)\n")
write(f, "Total HDO col: $(HDO_per_cc*2e5)\n")
write(f, "Total water col: $((H2O_per_cc + HDO_per_cc)*2e5)\n")
write(f, "H2O+HDO at surface: $((H2O_per_cc[1] + HDO_per_cc[1])*2e5)\n")
write(f, "Total H2O (pr μm): $(H2Oprum)\n")
write(f, "Total HDO (pr μm): $(HDOprum)\n")
write(f, "Total H2O+HDO, no enhancement: $(H2Oprum + HDOprum)")
close(f)

# Plot the water profiles (as raw mixing ratio) ================================
# fig, ax = subplots(figsize=(6,9))
# plot_bg(ax)
# semilogx(H2Oinitfrac, alt[2:end-1]/1e5, color="cornflowerblue", linewidth=3,
#          label=L"H$_2$O")
# semilogx(HDOinitfrac, alt[2:end-1]/1e5, color="cornflowerblue", linestyle="--",
#          linewidth=3, label="HDO")
# semilogx(H2Osatfrac, alt[1:end]/1e5, color="black", alpha=0.5, linewidth=3,
#          label=L"H$_2$O saturation")
# semilogx(HDOsatfrac, alt[1:end]/1e5, color="black", alpha=0.5, linestyle="--",
#          linewidth=3, label="HDO saturation")
# xlabel("Mixing ratio", fontsize=18)
# ylabel("Altitude [km]", fontsize=18)
# title(L"H$_2$O and HDO model profiles", fontsize=20)
# ax.tick_params("both",labelsize=16)
# legend()
# savefig(results_dir*FNext*"/water_profile_rawMR.png") # DEBUG: uncomment

################################################################################
#                             BOUNDARY CONDITIONS                              #
################################################################################

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
    We now have plot_bgobjects that return the list of indices and coefficients
    for transport, assuming no other species in the atmosphere
    (transportmat), and for chemistry, assuming no other altitudes
    (chemical_jacobian). We need to perform a kind of outer product on
    these operators, to determine a fully coupled set of equations for
    all species at all altitudes.
=#

# need to get a list of all species at all altitudes to iterate over

#=
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
active_above = [Symbol(string(s)*"_above") for s in activespecies]
active_below = [Symbol(string(s)*"_below") for s in activespecies]


# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x), activespecies)...);

# These lines are useful for troubleshooting if you're getting weird errors with ratefn.
# "no method matching +()" means the chemical system is unbalanced and you've got a
# species that has production but no consumption, or vice versa.
# other errors may occur. Uncomment these 3 lines to inspect what goes into rates_local.
# println("The contents of rates_local: ")
# println(rates_local)
# println("That was rates_local")

const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

arglist_local = [activespecies; active_above; active_below; inactivespecies;
                 Jratelist; :Tn; :Ti; :Te; #:M; :E;  DEBUG: add :M and :E back if needed
                 local_transport_rates; :dt]
                 #

arglist_local_typed = [:($s::Float64) for s in arglist_local]

@eval begin
    function ratefn_local(arglist_vector_typed::Array{Float64, 1})#($(arglist_local_typed[1:end-1]...))
        #=
        Produces the symbolic expression for production and loss of every species
        via both transport and chemistry
        =#

        # stuff mike suggested for debugging compilter errors:
        # use this for the argument: (arglist_vector_typed::Array{Float64, 1})
        # unpack vector 
        $(arglist_local_typed[1:end-1]...) = arglist_vector_typed
        # $(rates_local) 
        # end mike's suggestions

        # M and E are calculated here to ensure that the right number of ions/electrons
        # is used.
        M = eval(Expr(:call, :+, $(fullspecieslist...)))
        E = eval(Expr(:call, :+, $(ionlist...)))
        $rates_local # evaluates the rates_local expression
    end
end

@eval begin
    function chemJmat_local(arglist_vector_typed::Array{Float64, 1})#($(arglist_local_typed...))

        # mike's suggestions:
        #(arglist_vector_typed::Array{Float64, 1})
        # unpack vector
        $(arglist_local_typed...) = arglist_vector_typed

        M = eval(Expr(:call, :+, $(fullspecieslist...)))
        E = eval(Expr(:call, :+, $(ionlist...)))

        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = -dt*$(chemJ_local[3])

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = -dt*$(chemJ_above[3])

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(chemJ_below[3])

        ((localchemJi, localchemJj, localchemJval),
         (abovechemJi, abovechemJj, abovechemJval),
         (belowchemJi, belowchemJj, belowchemJval))
    end
end

# TODO: this function also probably needs edits with argument number
@eval begin
    function reactionrates_local(spclist_vector, Jratelist_vector, Tn, Ti, Te)#($(specieslist...), $(Jratelist...), Tn, Ti, Te)#, M, E) # DEBUG: revert as necessary
        #= 
        a function to return chemical reaction rates by multiplying the reaction rate coefficient
        (represented by x[3]) by the concentrations of all reactants (represented by x[1]...).
        the concentrations and Jrates are contained in specieslist and Jratelist, which are flattened over altitude
        using the splats (I think).
        =#

        # unpack the vectors
        $(specieslist...) = spclist_vector
        $(Jratelist...) = Jratelist_vector

        M = eval(Expr(:call, :+, $(fullspecieslist...)))
        E = eval(Expr(:call, :+, $(ionlist...)))
        $(Expr(:vcat, map(x->Expr(:call,:*,x[1]..., x[3]), reactionnet)...))
    end
end

################################################################################
#                         PHOTOCHEMICAL CROSS SECTIONS                         #
################################################################################

xsecfolder = research_dir * "uvxsect/";

# Crosssection Files ===========================================================
co2file = "CO2.dat"
co2exfile = "binnedCO2e.csv" # added to shield short λ of sunlight in upper atmo
h2ofile = "h2oavgtbl.dat"
hdofile = "HDO.dat"#"HDO_250K.dat"#
h2o2file = "H2O2.dat"
hdo2file = "H2O2.dat" #TODO: do HDO2 xsects exist?
o3file = "O3.dat"
o3chapfile = "O3Chap.dat"
o2file = "O2.dat"
o2_130_190 = "130-190.cf4"
o2_190_280 = "190-280.cf4"
o2_280_500 = "280-500.cf4"
h2file = "binnedH2.csv"
hdfile = "binnedH2.csv" # TODO: change this to HD file if xsects ever exist
ohfile = "binnedOH.csv"
oho1dfile = "binnedOHo1D.csv"
odfile = "OD.csv"

# NEW: files for crosssections of specific photodissociation/photoionization
# Source: Roger Yelle
# NOTE: Order has changed compared to the order in the Jrates list and the absorber dict,
# for ease of reading.
H2O_ionize_file = "JH2OtoH2Opl.dat"    # DEBUG
H2O_ionOdiss_file = "JH2OtoOplpH2.dat"
H2O_ionHdiss_file = "JH2OtoHplpOH.dat"
H2O_ionOHdiss_file = "JH2OtoOHplpH.dat"

CO_diss_file = "JCOtoCpO.dat"
CO_ionize_file = "JCOtoCOpl.dat"
CO_ionOdiss_file = "JCOtoCpOpl.dat"
CO_ionCdiss_file = "JCOtoOpCpl.dat"

N2_diss_file = "JN2OtoN2pO1D.dat"
N2_ionize_file = "JN2toN2pl.dat"
N2_iondiss_file = "JN2toNplpN.dat"
NO2_diss_file = "JNO2toNOpO.dat"
NO2_ionize_file = "JNO2toNO2pl.dat"
NO_diss_file = "JNOtoNpO.dat"
NO_ionize_file = "JNOtoNOpl.dat"
N2O_ionize_file = "JN2OtoN2Opl.dat"

H_ionize_file = "JHtoHpl.dat"
H2_ion_file = "JH2toH2pl.dat"
H2_iondiss_file = "JH2toHplpH.dat"

CO2_totaldiss_file = "JCO2toCpOpO.dat"
CO2_diss_file = "JCO2toCpO2.dat"
CO2_ionize_file = "JCO2toCO2pl.dat"
CO2_doubleion_file = "JCO2toCO2plpl.dat"
CO2_ionC2diss_file = "JCO2toCplplpO2.dat"
CO2_ionCdiss_file = "JCO2toCplpO2.dat"
CO2_ionCOandOdiss_file = "JCO2toCOplpOpl.dat"
CO2_ionCOdiss_file = "JCO2toCOplpO.dat"
CO2_ionOdiss_file = "JCO2toOplpCO.dat"
CO2_ionCandOdiss_file = "JCO2toOplpCplpO.dat"

H2O2_ionize_file = "JH2O2toH2O2pl.dat"

O_iondiss_file = "JOtoOpl.dat"
O2_ionize_file = "JO2toO2pl.dat"
O3_ionize_file = "JO3toO3pl.dat"

# Loading Data =================================================================
# CO2 photodissociation --------------------------------------------------------
# temperature-dependent between 195-295K
co2xdata = readdlm(xsecfolder*co2file,'\t', Float64, comments=true, comment_char='#')

# CO2 photoionization (used to screen high energy sunlight)
co2exdata = readdlm(xsecfolder*co2exfile,',',Float64, comments=true, comment_char='#')

# H2O & HDO --------------------------------------------------------------------
h2oxdata = readdlm(xsecfolder*h2ofile,'\t', Float64, comments=true, comment_char='#')

# These crosssections for HDO are for 298K.
hdoxdata = readdlm(xsecfolder*hdofile,'\t', Float64, comments=true, comment_char='#')

# H2O2 + HDO2 ------------------------------------------------------------------
# the data in the following table cover the range 190-260nm
h2o2xdata = readdlm(xsecfolder*h2o2file,'\t', Float64, comments=true, comment_char='#')
hdo2xdata = readdlm(xsecfolder*hdo2file,'\t', Float64, comments=true, comment_char='#')

# O3 ---------------------------------------------------------------------------
# including IR bands which must be resampled from wavenumber
o3xdata = readdlm(xsecfolder*o3file,'\t', Float64, comments=true, comment_char='#')
o3ls = o3xdata[:,1]
o3xs = o3xdata[:,2]
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

# PHOTODISSOCIATION ============================================================

const solarflux=readdlm(research_dir*solarfile,'\t', Float64, comments=true, comment_char='#')[1:2000,:]  # cut off at 2000 nm
solarflux[:,2] = solarflux[:,2]/2  # To roughly put everything at an SZA=60° (from a Kras comment)

absorber = Dict(:JCO2ion =>:CO2,
                :JCO2toCOpO=>:CO2,
                :JCO2toCOpO1D=>:CO2,
                :JO2toOpO=>:O2,
                :JO2toOpO1D=>:O2,
                :JO3toO2pO=>:O3,
                :JO3toO2pO1D=>:O3,
                :JO3toOpOpO=>:O3,
                :JH2toHpH=>:H2,
                :JHDtoHpD=>:HD,
                :JOHtoOpH=>:OH,
                :JOHtoO1DpH=>:OH,
                :JODtoOpD=>:OD,
                :JODtoO1DpD=>:OD,
                :JHO2toOHpO=>:HO2,
                :JDO2toODpO=>:DO2,
                :JH2OtoHpOH=>:H2O,
                :JH2OtoH2pO1D=>:H2O,
                :JH2OtoHpHpO=>:H2O,
                :JH2O2to2OH=>:H2O2,
                :JH2O2toHO2pH=>:H2O2,
                :JH2O2toH2OpO1D=>:H2O2,
                :JHDO2toHDOpO1D=>:HDO2,
                :JHDOtoHpOD=>:HDO,
                :JHDO2toOHpOD=>:HDO2,
                :JHDO2toDO2pH=>:HDO2,
                :JHDO2toHO2pD=>:HDO2,
                :JHDOtoDpOH=>:HDO,
                :JHDOtoHpDpO=>:HDO,
                :JHDOtoHDpO1D=>:HDO,
                # NEW: reactions from Roger's model. 
                :JH2OtoH2Opl=>:H2O,   # DEBUG
                :JH2OtoOplpH2=>:H2O,
                :JCOtoCpO=>:CO,
                :JCOtoCOpl=>:CO,
                :JN2OtoN2pO1D =>:N2O,
                :JH2toH2pl=>:H2,
                :JCOtoCpOpl=>:CO,
                :JNO2toNOpO=>:NO2,
                :JCO2toCpO2=>:CO2,
                :JCO2toCplplpO2=>:CO2,
                :JNOtoNOpl=>:NO,
                :JH2toHplpH=>:H2,
                :JH2OtoHplpOH=>:H2O,
                :JH2O2toH2O2pl=>:H2O2,
                :JN2toN2pl=>:N2,
                :JCO2toCOplpOpl=>:CO2,
                :JCOtoOpCpl=>:CO,
                :JCO2toOplpCplpO=>:CO2,
                :JNOtoNpO=>:NO,
                :JCO2toCplpO2=>:CO2,
                :JCO2toCO2pl=>:CO2,
                :JOtoOpl=>:O,
                :JH2OtoOHplpH=>:H2O,
                :JNO2toNO2pl=>:NO2,
                :JCO2toCOplpO=>:CO2,
                :JN2toNplpN=>:N2,
                :JCO2toCpOpO=>:CO2,
                :JCO2toCO2plpl=>:CO2,
                :JCO2toOplpCO=>:CO2,
                :JO2toO2pl=>:O2,
                :JHtoHpl=>:H,
                :JN2OtoN2Opl=>:N2O,
                :JO3toO3pl=>:O3
                );

# Enter crosssections into master dictionary ==================================
#=
    this is a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, crosssection has already folded in quantum
    efficiency considerations.
=#
crosssection = Dict{Symbol, Array{Array{Float64}}}()

# CO2 photodissociation --------------------------------------------------------
setindex!(crosssection, fill(co2exdata, length(alt)), :JCO2ion)
#CO2+hv->CO+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))),
          map(t->co2xsect(co2xdata, t),map(Temp_n, alt))), :JCO2toCOpO)
#CO2+hv->CO+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))),
          map(t->co2xsect(co2xdata, t),map(Temp_n, alt))), :JCO2toCOpO1D)

# O2 photodissociation ---------------------------------------------------------
#O2+hv->O+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), map(Temp_n, alt))),
          :JO2toOpO)
#O2+hv->O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), map(Temp_n, alt))),
          :JO2toOpO1D)

# O3 photodissociation ---------------------------------------------------------
# O3+hv->O2+O
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1-(1.37e-2*193-2.16)),
                               (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
                               (l->225<=l<306, 0.1),
                               (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
                               (l->328<=l<340, 0.92),
                               (l->340<=l, 1.0)
                              )), map(Temp_n, alt)), :JO3toO2pO)
# O3+hv->O2+O1D
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1.37e-2*193-2.16),
                               (l->193<=l<225, l->(1.37e-2*l-2.16)),
                               (l->225<=l<306, 0.9),
                               (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                               (l->328<=l<340, 0.08),
                               (l->340<=l, 0.0)
                              )), map(Temp_n, alt)), :JO3toO2pO1D)
# O3+hv->O+O+O
setindex!(crosssection,
          fill(quantumyield(o3xdata,((x->true, 0.),)),length(alt)),
          :JO3toOpOpO)

# H2 and HD photodissociation --------------------------------------------------
# H2+hv->H+H
setindex!(crosssection, fill(h2xdata, length(alt)), :JH2toHpH)
# HD+hν -> H+D
setindex!(crosssection, fill(hdxdata, length(alt)), :JHDtoHpD)

# OH and OD photodissociation --------------------------------------------------
# OH+hv->O+H
setindex!(crosssection, fill(ohxdata, length(alt)), :JOHtoOpH)
# OH+hv->O1D+H
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JOHtoO1DpH)
# OD + hv -> O+D
setindex!(crosssection, fill(odxdata, length(alt)), :JODtoOpD)
# OD + hν -> O(¹D) + D
setindex!(crosssection, fill(ohO1Dxdata, length(alt)), :JODtoO1DpD)

# HO2 and DO2 photodissociation ------------------------------------------------
# HO2 + hν -> OH + O
setindex!(crosssection, fill(ho2xsect, length(alt)), :JHO2toOHpO)
# DO2 + hν -> OD + O
setindex!(crosssection, fill(do2xsect, length(alt)), :JDO2toODpO)

# H2O and HDO photodissociation ------------------------------------------------
# H2O+hv->H+OH
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),length(alt)),
          :JH2OtoHpOH)

# H2O+hv->H2+O1D
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JH2OtoH2pO1D)

# H2O+hv->H+H+O
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->true, 0),)),length(alt)),
          :JH2OtoHpHpO)

# HDO + hν -> H + OD
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoHpOD)

# HDO + hν -> D + OH
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),length(alt)),
          :JHDOtoDpOH)

# HDO + hν -> HD + O1D
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JHDOtoHDpO1D)

# HDO + hν -> H + D + O
setindex!(crosssection,
          fill(quantumyield(hdoxdata,((x->true, 0),)),length(alt)),
          :JHDOtoHpDpO)


# H2O2 and HDO2 photodissociation ----------------------------------------------
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->h2o2xsect(h2o2xdata, t), map(Temp_n, alt))), :JH2O2to2OH)

# H2O2+hv->HO2+H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))),
          map(t->h2o2xsect(h2o2xdata, t), map(Temp_n, alt))), :JH2O2toHO2pH)

# H2O2+hv->H2O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(h2o2xdata, t),
          map(Temp_n, alt))), :JH2O2toH2OpO1D)

# HDO2 + hν -> OH + OD
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alt))), :JHDO2toOHpOD)

# HDO2 + hν-> DO2 + H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alt))), :JHDO2toDO2pH)

# HDO2 + hν-> HO2 + D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))),
          map(t->hdo2xsect(hdo2xdata, t), map(Temp_n, alt))), :JHDO2toHO2pD)

# HDO2 + hν -> HDO + O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(hdo2xdata, t),
          map(Temp_n, alt))), :JHDO2toHDOpO1D)

# NEW: CO2 photodissociation ---------------------------------------------------------
# Source: Roger Yelle
# DEBUG: uncomment
# CO₂ + hν -> C + O + O
CO2_totaldiss_data = readdlm(xsecfolder*CO2_totaldiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_totaldiss_data, length(alt)), :JCO2toCpOpO)

# CO2 + hν -> C + O₂
CO2_diss_data = readdlm(xsecfolder*CO2_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_diss_data, length(alt)), :JCO2toCpO2)

# # NEW: CO photodissociation ---------------------------------------------------------
# # Source: Roger Yelle

# CO + hν -> C + O
CO_diss_data = readdlm(xsecfolder*CO_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_diss_data, length(alt)), :JCOtoCpO)

# # NEW: Nitrogen species photodissociation --------------------------------------------
# # Source: Roger Yelle

# N₂ + hν -> N₂ + O(¹D)
N2_diss_data = readdlm(xsecfolder*N2_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2_diss_data, length(alt)), :JN2OtoN2pO1D)

# NO₂ + hν -> NO + O
NO2_diss_data = readdlm(xsecfolder*NO2_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO2_diss_data, length(alt)), :JNO2toNOpO)

# NO + hν -> N + O
NO_diss_data = readdlm(xsecfolder*NO_diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO_diss_data, length(alt)), :JNOtoNpO)

# # Photoionization or ionizing dissociation reactions ============================================

# # NEW: CO₂ ionization -----------------------------------------------------------------
# # Source: Roger Yelle

# CO₂ + hν -> CO₂⁺
CO2_ionize_data = readdlm(xsecfolder*CO2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionize_data, length(alt)), :JCO2toCO2pl)

# # CO₂ + hν -> CO₂²⁺  (even though we don't track doubly ionized CO₂)
CO2_doubleion_data = readdlm(xsecfolder*CO2_doubleion_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_doubleion_data, length(alt)), :JCO2toCO2plpl)

# # CO₂ + hν -> C²⁺ + O₂
CO2_ionC2diss_data = readdlm(xsecfolder*CO2_ionC2diss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionC2diss_data, length(alt)), :JCO2toCplplpO2)

# # CO₂ + hν -> C⁺ + O₂
CO2_ionCdiss_data = readdlm(xsecfolder*CO2_ionCdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCdiss_data, length(alt)), :JCO2toCplpO2)

# # CO₂ + hν -> CO⁺ + O⁺
CO2_ionCOandOdiss_data = readdlm(xsecfolder*CO2_ionCOandOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCOandOdiss_data, length(alt)), :JCO2toCOplpOpl)

# # CO₂ + hν -> CO⁺ + O
CO2_ionCOdiss_data = readdlm(xsecfolder*CO2_ionCOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCOdiss_data, length(alt)), :JCO2toCOplpO)

# # CO₂ + hν -> CO + O⁺
CO2_ionOdiss_data = readdlm(xsecfolder*CO2_ionOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionOdiss_data, length(alt)), :JCO2toOplpCO)

# # CO₂ + hν -> C⁺ + O⁺ + O
CO2_ionCandOdiss_data = readdlm(xsecfolder*CO2_ionCandOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO2_ionCandOdiss_data, length(alt)), :JCO2toOplpCplpO)

# # NEW: H2O ionization --------------------------------------------------------------
# # Source: Roger Yelle

# H2O + hν -> H2O⁺
h2o_ionize_data = readdlm(xsecfolder*H2O_ionize_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionize_data, length(alt)), :JH2OtoH2Opl)

# H2O + hν -> O⁺ + H2
h2o_ionOdiss_data = readdlm(xsecfolder*H2O_ionOdiss_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionOdiss_data, length(alt)), :JH2OtoOplpH2)

# # H2O + hν -> H⁺ + OH
h2o_ionHdiss_data = readdlm(xsecfolder*H2O_ionHdiss_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionHdiss_data, length(alt)), :JH2OtoHplpOH)

# # H2O + hν -> OH⁺ + H
h2o_ionOHdiss_data = readdlm(xsecfolder*H2O_ionOHdiss_file, ',', Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(h2o_ionOHdiss_data, length(alt)), :JH2OtoOHplpH)

# # NEW: CO ionization ----------------------------------------------------------------
# # Source: Roger Yelle

# # CO + hν -> CO⁺
CO_ionize_data = readdlm(xsecfolder*CO_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_ionize_data, length(alt)), :JCOtoCOpl)

# # CO + hν -> C + O⁺
CO_ionOdiss_data = readdlm(xsecfolder*CO_ionOdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_ionOdiss_data, length(alt)), :JCOtoCpOpl)

# # CO + hν -> C⁺ + O
CO_ionCdiss_data = readdlm(xsecfolder*CO_ionCdiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(CO_ionCdiss_data, length(alt)), :JCOtoOpCpl)

# # NEW: Nitrogen species ionization --------------------------------------------------
# # Source: Roger Yelle

# # N₂ + hν -> N₂⁺
N2_ionize_data = readdlm(xsecfolder*N2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2_ionize_data, length(alt)), :JN2toN2pl)

# # N₂ + hν -> N⁺ + N
N2_iondiss_data = readdlm(xsecfolder*N2_iondiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2_iondiss_data, length(alt)), :JN2toNplpN)

# # NO₂ + hν -> NO₂⁺
NO2_ionize_data = readdlm(xsecfolder*NO2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO2_ionize_data, length(alt)), :JNO2toNO2pl)

# # NO + hν -> NO⁺
NO_ionize_data = readdlm(xsecfolder*NO_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(NO_ionize_data, length(alt)), :JNOtoNOpl)

# # N₂O + hν -> N₂O⁺
N2O_ionize_data = readdlm(xsecfolder*N2O_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(N2O_ionize_data, length(alt)), :JN2OtoN2Opl)

# # NEW: Molecular and atomic hydrogen ionization -------------------------------------
# # Source: Roger Yelle

# # H + hν -> H⁺
H_ionize_data = readdlm(xsecfolder*H_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H_ionize_data, length(alt)), :JHtoHpl)

# # H₂ + hν -> H₂⁺
H2_ion_data = readdlm(xsecfolder*H2_ion_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H2_ion_data, length(alt)), :JH2toH2pl)

# # H₂ + hν -> H⁺ + H
H2_iondiss_data = readdlm(xsecfolder*H2_iondiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H2_iondiss_data, length(alt)), :JH2toHplpH)

# # H₂O₂ + hν -> H₂O₂⁺
H2O2_ionize_data = readdlm(xsecfolder*H2O2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(H2O2_ionize_data, length(alt)), :JH2O2toH2O2pl)

# # NEW: Oxygen and ozone ionization --------------------------------------------------
# # Source: Roger Yelle

# O + hν -> O⁺
O_iondiss_data = readdlm(xsecfolder*O_iondiss_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(O_iondiss_data, length(alt)), :JOtoOpl)

# O₂ + hν -> O₂⁺
O2_ionize_data = readdlm(xsecfolder*O2_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(O2_ionize_data, length(alt)), :JO2toO2pl)

# # O₃ + hν -> O₃⁺
O3_ionize_data = readdlm(xsecfolder*O3_ionize_file,',',Float64, comments=true, comment_char='#')
setindex!(crosssection, fill(O3_ionize_data, length(alt)), :JO3toO3pl)


# Solar Input ==================================================================
lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

if !(setdiff(solarflux[:,1],lambdas)==[])
    throw("Need a broader range of solar flux values!")
end

# pad all cross-sections to solar
for j in Jratelist, ialt in 1:length(alt)
    crosssection[j][ialt] = padtosolar(solarflux, crosssection[j][ialt])
end

################################################################################
#                             CONVERGENCE CODE                                 #
################################################################################

# set the water profiles =======================================================
n_current[:H2O] = H2Oinitfrac.*map(z->n_tot(n_current, z), alt[2:end-1])
n_current[:HDO] = HDOinitfrac.*map(z->n_tot(n_current, z), alt[2:end-1])

# Plot initial water profile ===================================================
# ...in ppm --------------------------------------------------------------------
fig, ax = subplots(figsize=(6,8))
plot_bg(ax)
ax.tick_params(axis="x", which="minor", bottom=true, top=true)
semilogx(H2Oinitfrac/1e-6, alt[2:end-1]/1e5, color="cornflowerblue", linewidth=2,
         label=L"H$_2$O")
semilogx(HDOinitfrac/1e-6, alt[2:end-1]/1e5, color="cornflowerblue",
         linestyle="--", linewidth=2, label="HDO")
xlabel("Volume Mixing Ratio [ppm]")
ylabel("Altitude [km]")
title(L"H$_2$O and HDO model profiles")
legend()
savefig(results_dir*FNext*"/initial_water_MR.png")
close(fig)

# ...and raw abundance ---------------------------------------------------------
fig, ax = subplots(figsize=(6,8))
plot_bg(ax)
ax.tick_params(axis="x", which="minor", bottom=true, top=true)
semilogx(n_current[:H2O], alt[2:end-1]/1e5, color="cornflowerblue", linewidth=2,
         label=L"H$_2$O")
semilogx(n_current[:HDO], alt[2:end-1]/1e5, color="cornflowerblue",
         linestyle="--", linewidth=2, label="HDO")
xlabel(L"Species density [cm$^{-2}$")
ylabel("Altitude [km]")
title(L"H$_2$O and HDO model profiles")
legend()
savefig(results_dir*FNext*"/initial_water_number.png")
close(fig)

# do the convergence ===========================================================

# convfig.show()
# convfig.canvas.draw()

# plot the initial guess for the atmosphere
println("Plotting the initial condition")
plotatm(n_current, [neutrallist, ionlist], results_dir*FNext*"/initial_atmosphere.png")

# Edensityfile = open(results_dir*"/Edensiteis.txt", "w") # TODO: Remove

println("Beginning Convergence")
img_savepath = results_dir*FNext*"/converged_atm_"*FNext*".png"
initdt = -10
@showprogress 0.1 "Converging over 10 My..." [timeupdate(t, img_savepath) for t in [10.0^(1.0*i) for i in initdt:14]] # TODO: change range to -3:14 (usual) if needed

@showprogress 0.1 "Last convergence steps..." for i in 1:100  # TODO: uncomment after debugging plots
    plotatm(n_current, [neutrallist, ionlist], img_savepath, t="1e14", iter=i)
    # println("dt: 1e14 iter $(i)")
    update!(n_current, 1e14, false)
end

# close(Edensityfile) # TODO: remove

# write out the new converged file to matching folder.
towrite = results_dir*FNext*"/converged_"*FNext*".h5"
write_ncurrent(n_current, towrite)
println("Wrote $(towrite)")

# save the figure
# savefig(img_savepath, bbox_inches="tight")
println("Saved figure to same folder")

################################################################################
#                                 LOGGING                                      #
################################################################################

# crosssection dict for logging purposes =======================================
xsect_dict = Dict("CO2"=>[co2file, co2exfile],
                  "H2O, HDO"=>[h2ofile, hdofile],
                  "H2O2, HDO2"=>[h2o2file, hdo2file],
                  "O3"=>[o3file, o3chapfile],
                  "O2"=>[o2file, o2_130_190, o2_190_280, o2_280_500],
                  "H2, HD"=>[h2file, hdfile],
                  "OH, OD"=>[ohfile, oho1dfile, odfile])

# Log temperature and water parameters =========================================
if args[1]=="temp"
    input_string = "T_0=$(args[2]), T_tropo=$(args[3]), T_exo=$(args[4])" *
                   "\nwater init=$(MR)\nDH=5.5 \nOflux=1.2e8" #*
                   #"\nlapse rate=$(lapserate_logme)\n"
elseif args[1]=="water"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\n" *
                   "water init=$(args[2])\nDH=5.5\nOflux=1.2e8\n" #*
                   #"lapse rate=$(lapserate_logme)\n"
elseif args[1]=="dh"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)\n" *
                   "DH=$(args[2]) \nOflux=1.2e8\n"#lapse rate=$(lapserate_logme)\n"
elseif args[1]=="Oflux"
    input_string = "T_s=$(meanTs), T_tropo=$(meanTt), T_exo=$(meanTe)\nwater=(MR)" *
                   "\nDH=5.5\nOflux=$(Of)\n"#lapse rate=$(lapserate_logme)\n"
end

# Write the log ================================================================
f = open(results_dir*FNext*"/convergence_data_"*FNext*".txt", "w")
write(f, "Finished convergence for $(args[1]) experiment with control parameters: \n")
write(f, input_string)
write(f, "\nSVP fixed: $(fix_SVP)\n")
write(f, "Solar cycle status: solar $(cycle)\n")

# cross sections
write(f, "\nCROSS SECTIONS: \n")
for k in keys(xsect_dict)  # cross sections
    write(f, k*": "*join(xsect_dict[k], ", ")*"\n")
end
# boundary conditions
write(f, "\nBOUNDARY CONDITIONS: \n")
write(f, "n: number density at surface, f: flux at top, v: velocity at top\n")
for k2 in keys(speciesbclist)
    bcstring = join([join(speciesbclist[k2][1, :], "="),
                     join(speciesbclist[k2][2, :], "=")], ", ")
    write(f, string(k2)*": "*bcstring*"\n")
end

# timestep etc
write(f, "Initial atmosphere state: $(readfile)\n")
write(f, "Initial dt size: $(initdt)\n")


close(f)

println("ALERT: Finished")
println()
