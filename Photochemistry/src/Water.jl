# **************************************************************************** #
#                                                                              #
#                       Functions related to water profile                     #
#                                                                              #
# **************************************************************************** #

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
    @assert all(x->x in keys(GV),  [:molmass])

    col_abundance = column_density(sp_profile)
    cc_per_g = GV.molmass[sp] / GV.molmass[:H2O] # Water is 1 g/cm^3. Scale appropriately.

    #pr μm = (#/cm²) * (1 mol/molecules) * (g/1 mol) * (1 cm^3/g) * (10^4 μm/cm)
    pr_microns = col_abundance * (1/6.02e23) * (GV.molmass[sp]/1) * (cc_per_g / 1) * (1e4/1)
    return pr_microns
end

# 1st term is a conversion factor to convert to (#/cm^3) from Pa. Source: Marti & Mauersberger 1993
Psat(T) = (1e-6 ./ (kB_MKS .* T)) .* (10 .^ (-2663.5 ./ T .+ 12.537))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# However, this function is defined on the offchance someone studies HDO.
Psat_HDO(T) = (1e-6/(kB_MKS * T))*(10^(-2663.5/T + 12.537))

function setup_water_profile!(atmdict; dust_storm_on=false, make_sat_curve=false, tanh_prof="standard", fixed_bdy_i=n_alt_index[72e5], 
                                       showonly=false, hygropause_alt=40e5, globvars...)
    #=
    Sets up the water profile as a fraction of the initial atmosphere. 
    Input:
        atmdict: dictionary of atmospheric density profiles by altitude
        Optional:
            dust_storm_on: whether to add an extra parcel of water at a certain altitude.
            tanh_prof: "low", "standard", or "high" to choose 1/10, mean, or 10x as much water in the atmosphere.
            hygropause_alt: altitude at which the water will switch from well-mixed to following the saturation vapor pressure curve.
    Output: 
        atmdict: Modified in place with the new water profile. 
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :num_layers, :DH, :alt, :plot_grid, :n_alt_index,
                                   :non_bdy_layers, :H2Osat, :water_mixing_ratio,
                                   :results_dir, :sim_folder_name]) # not actually required: :H2O_excess, :HDO_excess, :ealt, 

    # H2O Water Profile ================================================================================================================
    H2Osatfrac = GV.H2Osat ./ map(z->n_tot(atmdict, z; GV.all_species, GV.n_alt_index), GV.alt)  # get SVP as fraction of total atmo
    # set H2O SVP fraction to minimum for all alts above first time min is reached
    H2Oinitfrac = H2Osatfrac[1:something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)]
    H2Oinitfrac = [H2Oinitfrac;   # ensures no supersaturation
                   fill(minimum(H2Osatfrac), GV.num_layers-length(H2Oinitfrac))]

    # Set lower atmospheric water to be well-mixed (constant with altitude) below the hygropause
    H2Oinitfrac[findall(x->x<hygropause_alt, GV.alt)] .= GV.water_mixing_ratio

    for i in [1:length(H2Oinitfrac);]
        H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1]
    end

    # For doing highand low water cases ================================================================================================
    if tanh_prof == "high"
        H2Oinitfrac[1:fixed_bdy_i] = H2Oinitfrac[1:fixed_bdy_i] .* water_tanh_prof(GV.non_bdy_layers./1e5; z0=GV.ealt, f=50)[1:fixed_bdy_i]
    elseif tanh_prof=="low"
        H2Oinitfrac[1:fixed_bdy_i] = H2Oinitfrac[1:fixed_bdy_i] .* water_tanh_prof(GV.non_bdy_layers./1e5; z0=GV.ealt, f=0.05)[1:fixed_bdy_i]
    elseif tanh_prof=="standard"
        H2Oinitfrac=H2Oinitfrac
    end

    # set the water profiles ===========================================================================================================
    atmdict[:H2O] = H2Oinitfrac.*n_tot(atmdict; GV.n_alt_index, GV.all_species)
    atmdict[:HDO] = 2 * GV.DH * atmdict[:H2O] 
    HDOinitfrac = atmdict[:HDO] ./ n_tot(atmdict; GV.n_alt_index, GV.all_species)  # Needed to make water plots.

    
    # ADD EXCESS WATER AS FOR DUST STORMS.
    if dust_storm_on
        sigma = 12.5
        H2Oppm = 1e-6*map(z->GV.H2O_excess .* exp(-((z-GV.ealt)/sigma)^2), GV.non_bdy_layers/1e5) + H2Oinitfrac 
        HDOppm = 1e-6*map(z->GV.HDO_excess .* exp(-((z-GV.ealt)/sigma)^2), GV.non_bdy_layers/1e5) + HDOinitfrac
        atmdict[:H2O][1:fixed_bdy_i] = (H2Oppm .* n_tot(atmdict; GV.n_alt_index, GV.all_species))[1:fixed_bdy_i]
        atmdict[:HDO][1:fixed_bdy_i] = (HDOppm .* n_tot(atmdict; GV.all_species))[1:fixed_bdy_i]
    end

    # Plot the water profile ===========================================================================================================
    if make_sat_curve
        satarray = H2Osatfrac
    else
        satarray = nothing 
    end
    plot_water_profile(H2Oinitfrac, HDOinitfrac, atmdict[:H2O], atmdict[:HDO], GV.results_dir*GV.sim_folder_name, watersat=satarray, plot_grid=GV.plot_grid, showonly=showonly)
end 

function water_tanh_prof(z; f=10, z0=62, dz=11)
    #=
    Apply a tanh_prof to the water init fraction to add or subtract water from the atmosphere.
    =#

    # subtract_me = f > 1 ? 1 : f/2

    return ((f .- 1)/2) * (tanh.((z .- z0) ./ dz) .+ 1) .+ 1
end
