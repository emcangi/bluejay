################################################################################
# MODEL_SETUP.jl
# DESCRIPTION: Sets up variables, constants, etc. that typically don't require
# checking by the user before each run but change depending on the type of 
# model being run.
# 
# Eryn Cangi
# Created April 2024
# Last edited: August 2024
# Currently tested for Julia: 1.8.5
################################################################################

using DataFrames
using DoubleFloats
using CSV

# First do some error checking
if (special_seasonal_case!=nothing) & (exp_type=="all")
    throw("Only use exp_type='all' when special_seasonal_case != nothing")
end

# ***************************************************************************************************** #
#                                                                                                       #
#                                       Planet-dependent constants                                      #
#                                                                                                       #
# ***************************************************************************************************** #
# Yeah ok they're not really constants but this is the easiest way to manage it because they
# are used in the photochemistry module, which gets loaded before the parameter/constant files, and which
# cannot load something based on a variable defined in the enclosing scope. 
const M_P = Dict( # Planetary mass in g 
                 "Mars"=>0.1075*5.972e27, 
                 "Venus"=>4.867e27
                )[planet]
const R_P = Dict( # Planetary radius in cm
                 "Mars"=>3396e5, 
                 "Venus"=>6050e5
                )[planet] 
const DH = Dict( # Atmospheric D/H ratio 
                "Mars"=>5.5 * SMOW, # Yung 1988
                "Venus"=>240 * SMOW, # Fedorova 2008
               )[planet]
const sol_in_sec = Dict(
                        "Mars"=>88775.2438,   # One Mars sol in seconds
                        "Venus"=>2.09968e7
                       )[planet]
const season_in_sec = Dict(
                           "Mars"=>1.4838759e7,
                           "Venus"=>4.8535373e6
                          )[planet]
const g = bigG * M_P / (R_P^2);
const SA = 4*pi*(R_P)^2 # cm^2
 

# ***************************************************************************************************** #
#                                                                                                       #
#              Model species, Jrate lists, and lists of chemistry/transport species                     #
#                                                                                                       #
# ***************************************************************************************************** #

# Minor species that have reactions available in the network files, but aren't used. These are just for
# reference: [:CNpl,:HCNpl,:HCNHpl,:HN2Opl,:NH2pl,:NH3pl,:CH,:CN,:HCN,:HNO,:NH,:NH2,:HD2pl]
 
#                                     Full species list
# =======================================================================================================
const neutral_species = [conv_neutrals[planet]..., new_neutrals...];
const ion_species = [conv_ions[planet]..., new_ions...]
const new_species = [new_neutrals..., new_ions...]  # Needed later to be excluded from n_tot() as called 
                                                    # in the water saturation calculation, in the case
                                                    # where new species are being added.
const all_species = [neutral_species..., ion_species...];

#                        Photolysis and Photoionization rate symbol lists 
# =======================================================================================================
const nontherm = ions_included==true ? true : false   # whether to do non-thermal escape. Must be here, sued in call to format Jrates
const conv_Jrates, newJrates = format_Jrates(reaction_network_spreadsheet, all_species, "Jratelist"; ions_on=ions_included, hot_atoms=nontherm)
const Jratelist = [conv_Jrates..., newJrates...];

# These dictionaries specify the species absorbing a photon for each J rate, and the products of the reaction.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])


#                               Miscellaneous logical groupings
# =======================================================================================================
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
const D_bearing_species = get_deuterated(all_species)
const D_ions = get_deuterated(ion_species)
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];


#                             Define short- and long-lived species
# =======================================================================================================

# Short lived species, whose chemical lifetime is << diffusion timescale ------- #
const short_lived_species = [];# technically shortlived but count as longlived: :CH, :HCO, :HO2, :O3, :OH, :O1D, :DO2, :OD...
if assume_photochem_eq
    append!(short_lived_species, [:NO2, :CN, :HNO, :NH, :NH2, :C, :CH])
    append!(short_lived_species, ion_species)
end

# Long lived species 
const long_lived_species = setdiff(all_species, short_lived_species)


#                         Species participating in chemistry and transport
# =======================================================================================================

no_chem_species = []; 
no_transport_species = [];

# Fixed species (Densities don't update)
# -------------------------------------------------------------------
for s in dont_compute_either_chem_or_transport
    push!(no_chem_species, s)
    push!(no_transport_species, s)
end

for s in dont_compute_chemistry
    if ~(s in no_chem_species)
        push!(no_chem_species, s)
    end
end 
for s in dont_compute_transport
    if ~(s in no_transport_species)
        push!(no_transport_species, s)
    end
end


# Chemistry and transport participants
# -------------------------------------------------------------------
if converge_which == "neutrals"
    println("Note: Still removing nitrogen neutrals from the converged species. May want to change this.")
    append!(no_chem_species, union(conv_ions[planet], N_neutrals)) # This is because the N chemistry is highly coupled to the ions.
    append!(no_transport_species, union(conv_ions[planet], N_neutrals, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, conv_neutrals[planet])
    append!(no_transport_species, conv_neutrals[planet])
elseif converge_which == "both"
    append!(no_transport_species, short_lived_species)
elseif converge_which == "ions+nitrogen"
    append!(no_chem_species, setdiff(conv_neutrals[planet], N_neutrals))
    append!(no_transport_species, setdiff(conv_neutrals[planet], N_neutrals))
end

# Disallow transport and/or chemistry if the appropriate setting is toggled
# TODO: Currently not working :(
if do_trans==false
    append!(no_transport_species, all_species)
end

if do_chem==false
    append!(no_chem_species, all_species)
end

const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species 
# -------------------------------------------------------------------
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)

# Sort name lists created here
# -------------------------------------------------------------------
sort!(all_species)
sort!(neutral_species)
sort!(ion_species)
sort!(active_species)
sort!(inactive_species)
sort!(short_lived_species)
sort!(long_lived_species)
sort!(active_longlived)
sort!(active_shortlived)
sort!(chem_species)
sort!(transport_species)
sort!(no_chem_species)
sort!(no_transport_species)
sort!(D_bearing_species)
sort!(D_ions)
sort!(N_neutrals)


# ***************************************************************************************************** #
#                                                                                                       #
#                                        Atmospheric setup                                              #
#                                                                                                       #
# ***************************************************************************************************** #
data_alt = []
data_T = []
water_T = []

if ingest_data
    # load file
    thedata =  CSV.read(data_loc*ingest_scenario*".csv", DataFrame)
    data_alt = thedata[!, "alt"]
    data_T = thedata[!, "T"]
    data_water = thedata[!, "water"]
end

#                             Ion chemistry and non-thermal escape
# =======================================================================================================

const e_profile_type = ions_included==true ? "quasineutral" : "none" 
    # OPTIONS: 
    # "quasineutral" - n_e = sum of all the ion densities; PREFERRED
    # "O2+" - n_e = sum(n_O2pl)
    # "constant" - n_e set to some constant value which I believe is 1e5.
    # "none" - no electrons, for use in neutral-only models

#                                       Altitude grid                  
# =======================================================================================================
const zmin = Dict("Venus"=>90e5, "Mars"=>0.)[planet]
const dz = 2e5  # Discretized layer thickness
const zmax = 250e5  # Top altitude (cm)
const alt = convert(Array, (zmin:dz:zmax)) # These are the layer centers.
const n_all_layers = length(alt)
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]; # the altitude grid CELLS but in integers.
const non_bdy_layers = alt[2:end-1]  # all layers, centered on 2 km, 4...248. Excludes the boundary layers which are [-1, 1] and [249, 251].
const num_layers = length(non_bdy_layers) # there are 124 non-boundary layers.
const plot_grid = non_bdy_layers ./ 1e5;  # for plotting. Points located at atmospheric layer cell centers and in units of km.
const n_alt_index=Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])
const hygropause_alt = 40e5  # Location of the hygropause

#                              Temperature profile construction                      
# =======================================================================================================

# These are used to smooth out the profiles later. Has to be here cause it depends on the variable temp_scenario.
const fudge_factors = Dict("orbit12807_MD"=>15, "perihelion_MD"=>15., "baseline_perihelion_MD"=>15., "aphelion_MD"=>-5., "baseline_MD"=>0., )

const Tsurf = Dict("Mars"=>Dict("mean"=>230., 
                                # Options for data ingestion from MAVEN ('MD' for maven data) in collaboration with Mike Stevens. 
                                "orbit12807_MD"=>285, "perihelion_MD"=>273., "baseline_perihelion_MD"=>273., "aphelion_MD"=>255., "baseline_MD"=>230.,
                                # Special case where we do a seasonal model but vary 3 thigns at once.
                                "inclusive-ap"=>230., "inclusive-mean"=>230., "inclusive-peri"=>230.),

                   "Venus"=>Dict("mean"=>735.), 
                  )

const Tmeso = Dict("Mars"=>Dict("mean"=>130., 
                                "orbit12807_MD"=>169.6, "perihelion_MD"=>169.6, "baseline_perihelion_MD"=>169.6, "aphelion_MD"=>221.4, "baseline_MD"=>130.,
                                # Special case where we do a seasonal model but vary 3 thigns at once.
                               "inclusive-ap"=>130., "inclusive-mean"=>130., "inclusive-peri"=>130.), 
                   "Venus"=>Dict("mean"=>170.)
                   )

# Establish the options for controltemps[3]
const Texo = Dict("Mars"=>Dict("min-P2"=>190., "mean-P2"=>210., "max-P2"=>280.,   # These are based on solar min, mean, max.
                               "min"=>175., "mean"=>225., "max"=>275.,   # These are based on solar min, mean, max.
                               "meansundist"=>225., "aphelion"=>225., "perihelion"=>225., 
                               # Options for data ingestion from MAVEN ('MD' for maven data) in collaboration with Mike Stevens. 
                               "orbit12807_MD"=>250, "perihelion_MD"=>250., "baseline_perihelion_MD"=>250., "aphelion_MD"=>176., "baseline_MD"=>250.,
                               # Special case where we do a seasonal model but vary 3 thigns at once.
                               "inclusive-ap"=>175., "inclusive-mean"=>225., "inclusive-peri"=>275.),
                       "Venus"=>Dict("min"=>260., "mean"=>290., "max"=>320.)
                )


# Modify the settings if doing a special isothermal atmosphere.
if temp_scenario=="isothermal"
    const controltemps = [225., 225., 225.]
    const meantemps = [225., 225., 225.] # Used for saturation vapor pressure. DON'T CHANGE!
else # Set the exobase temp according to the temp scenario.
    if planet=="Venus"
        const meantemps = [Tsurf[planet]["mean"], Tmeso[planet]["mean"], Texo[planet]["min"]] # Used for saturation vapor pressure. DON'T CHANGE!
    elseif planet=="Mars"
        const meantemps = [Tsurf[planet]["mean"], Tmeso[planet]["mean"], Texo[planet]["mean"]] # Used for saturation vapor pressure. DON'T CHANGE!
    end

    const controltemps = [Tsurf[planet][temp_scenario], Tmeso[planet][temp_scenario], Texo[planet][temp_scenario]]
end

# Modify the array for the special case where multiple parameters are changed for the seasonal model
# if special_seasonal_case!=nothing 
#     const controltemps = [Tsurf[planet], Tmeso[planet], Texo_inclusive_opts[special_seasonal_case]]
# end

# Now create the actual temperature profiles
if planet=="Mars"
    if occursin("MD", temp_scenario)  # when ingesting maven data....
        T_array_dict = T_MARS_NEW(controltemps[1], controltemps[2], controltemps[3], data_T, data_alt; 
                                  lapserate=-1.4e-5, weird_Tn_param=8, fudge_factor=get(fudge_factors, temp_scenario, 0.), alt, dz)
    else  # normal model runs.......
        const T_array_dict = T_Mars(controltemps[1], controltemps[2], controltemps[3]; alt)
    end
    const Tn_meanSVP = T_Mars(meantemps...; alt)["neutrals"]; # Needed for boundary conditions.
elseif planet=="Venus"
    const T_array_dict = T_Venus(controltemps[1], controltemps[2], controltemps[3], "Venus-Inputs/FoxandSung2001_temps_mike.txt"; alt);
    const Tn_meanSVP = T_Venus(meantemps..., "Venus-Inputs/FoxandSung2001_temps_mike.txt"; alt)["neutrals"]; # Needed for boundary conditions.
end

const Tn_arr = T_array_dict["neutrals"]
const Ti_arr = T_array_dict["ions"]
const Te_arr = T_array_dict["electrons"]

const Tplasma_arr = Ti_arr .+ Te_arr;
# A comment on the plasma temperature: It's more rightly defined as (Te + Ti)/2, and comes into play in the diffusion
# equation for ions. Per Schunk & Nagy 2009, equations 5.55 and 5.56, the ambipolar diffusion coefficient is 
# Da = 2kT_p / (m_i * ν_in). This is the same as our formulation here (see Core.jl, Dcoef!()), which is 
# Da = k(T_e + T_i) / (m_i * ν_in). The diffusion equation in our model includes a term like 1/T * dT/dz (see also
# Schunk and Nagy 2009, equation 5.54, third term in the bracket). Note that the factor of 2 in the official definition of
# the plasma temperature cancels out here. If we defined plasma temperature as the average of T_e and T_i, we would have an 
# extra factor of 1/2 in the scale height (Schunk & Nagy 2009 equation 5.59) as calculated for plasma, unless we wrote a 
# special version of scaleH() (see Core.jl) to handle the plasma temperature. To avoid reformulating our scale height 
# function, we define the plasma temperature for the purposes of the model as simply T_e + T_i.
# It may be warranted in the future to change this so that the definition matches with the literature appropriately.
# --Eryn, 5 September 2024
const Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
const Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr)


#                                      Water profile settings
# =======================================================================================================

if dust_storm_on==true
    const excess_peak_alt = 60
end

# Water mixing ratios to use
if planet=="Mars"
    const water_MRs = Dict("loweratmo"=>Dict("standard"=>1.3e-4, "low"=>0.65e-4, "high"=>2.6e-4), 
                           "mesosphere"=>Dict("standard"=>1.3e-4, "high"=>1.3e-4, "low"=>1.3e-4), 
                           "everywhere"=>Dict("standard"=>1.3e-4, "high"=>1.3e-4, "low"=>1.3e-4))
    const water_mixing_ratio = water_MRs[water_loc][water_case]
elseif planet=="Venus"
    const water_mixing_ratio = Dict("standard"=>1e-6)[water_case]

    # SPECIAL: Crazy water Mahieux & Viscardy 2024
    if venus_special_water
        const h2o_vmr_low = 10^0.3 * 1e-6
        const h2o_vmr_high = 10^0.7 * 1e-6
        const hdo_vmr_low = 10^(-0.5)  * 1e-6
        const hdo_vmr_high = 10^(0.5) * 1e-6
    else
        const h2o_vmr_low = water_mixing_ratio
        const h2o_vmr_high = nothing
        const hdo_vmr_low = 2*DH*water_mixing_ratio
        const hdo_vmr_high = nothing
    end
end

# Whether to install a whole new water profile or just use the initial guess with modifications (for seasonal model)
if planet=="Venus"
    const reinitialize_water_profile = venus_special_water==true ? true : false
elseif planet=="Mars"
    const reinitialize_water_profile = seasonal_cycle==true ? false : true # should be off if trying to run simulations for seasons
end

const update_water_profile = seasonal_cycle==true ? true : false # this is for modifying the profile during cycling, MAY be fixed?
const modified_water_alts = "below fixed point"

# altitude at which to add the extra water -- applies to both dust storm parcels and the tanh profile
const add_water_alt_opts = Dict("low"=>45, "standard"=>60, "high"=>65)
const f_fac_opts = Dict("low"=>0.005, "standard"=>10, "high"=>100) # a parameter named f which helps manipulate the water profile. Not related to any other f.

# Set the saturation vapor pressure curves, used in boundary conditions
const H2Osat = map(x->Psat(x), Tn_meanSVP) # Using this function keeps SVP fixed 
const HDOsat = map(x->Psat_HDO(x), Tn_meanSVP)

# To allow water to be active in the upper atmosphere but not the lower atmosphere, we need 
# its position within the active species vector - these are used later in chemJ_mat.
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

#                              Species-specific scale heights
# =======================================================================================================
const Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)]; molmass, M_P, R_P) for sp in all_species])


#                                     Boundary conditions
# =======================================================================================================
# "n": density boundary condition; "f": flux bc; "v": velocity bc; 
# "see boundaryconditions()" -- nonthermal escape depends on the dynamic density of the
# atmosphere, so it can't be imposed as a constant here and is calculated on the fly.
if planet=="Mars"
    const speciesbclist=Dict(:CO2=>Dict("n"=>[2.1e17, NaN], "f"=>[NaN, 0.]),
                        :Ar=>Dict("n"=>[2.0e-2*2.1e17, NaN], "f"=>[NaN, 0.]),
                        :N2=>Dict("n"=>[1.9e-2*2.1e17, NaN], "f"=>[NaN, 0.]),
                        #:C=>Dict("f"=>[NaN, 4e5]), # NEW: Based on Lo 2021
                        :H2O=>Dict("n"=>[H2Osat[1], NaN], "f"=>[NaN, 0.]), # bc doesnt matter if H2O fixed
                        :HDO=>Dict("n"=>[HDOsat[1], NaN], "f"=>[NaN, 0.]),
                        :O=> Dict("f"=>[0., 1.2e8]),
                        :H2=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; M_P, R_P, zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),  # velocities are in cm/s
                        :HD=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 3.0; M_P, R_P, zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                        :H=> Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 1.0; M_P, R_P, zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                        :D=> Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; M_P, R_P, zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                       );
elseif planet=="Venus"
    const ntot_at_lowerbdy = 9.5e15 # Based on Fox & Sung 2001

    H2O_lowerbdy = h2o_vmr_low * ntot_at_lowerbdy
    HDO_lowerbdy = hdo_vmr_low * ntot_at_lowerbdy
    
    # END SPECIAL
    # MRs
    SO2mr = 1e-7
    # Simon's notes: Belyaev 2012: this was 0.1 ppmv at 165–170 K to 0.5–1 ppmv at 190–192 K; 
    # It said 0.1ppm was related to the most common temperature reading so I went with that 
    # (this is either 1E-7 or 6.79E-8 depending on the calculation)
    N2mr = 0.032
    H2SO4mr = 3e-9
    O2mr = 3e-3
    COmr = 4.5e-6
    CO2mr = 0.965
    HClmr = 3.66e-7
    # Simon's notes: Krasnopolsky, 2010a: this was 400ppb at 74km in altitude, and the actual number is likely lower 
    # (is either 4.0E-7, or 4.8E-7 depending on the calculation); and according to Zhang 2012 it is 3.66e-7
    SOmr = 1e-7

    
    const KoverH_lowerbdy = Keddy([zmin], [ntot_at_lowerbdy]; planet)[1]/scaleH_lowerboundary(zmin, Tn_arr[1]; molmass, M_P, R_P, zmin)
    const manual_speciesbclist=Dict(# major species neutrals at lower boundary (estimated from Fox&Sung 2001, Hedin+1985, agrees pretty well with VIRA)
                                    :CO2=>Dict("n"=>[CO2mr*ntot_at_lowerbdy, NaN], "f"=>[NaN, 0.]),
                                    :Ar=>Dict("n"=>[5e11, NaN], "f"=>[NaN, 0.]),
                                    :CO=>Dict("n"=>[COmr*ntot_at_lowerbdy, NaN], "f"=>[NaN, 0.]),
                                    :O2=>Dict("n"=>[O2mr*ntot_at_lowerbdy, NaN], "f"=>[NaN, 0.]),
                                    :N2=>Dict("n"=>[N2mr*ntot_at_lowerbdy, NaN]),

                                    
                                    :HCl=>Dict("n"=>[HClmr * ntot_at_lowerbdy, NaN]),
                                    :DCl=>Dict("n"=>[HClmr * DH * ntot_at_lowerbdy, NaN]),


                                    :H2SO4=>Dict("n"=>[H2SO4mr * ntot_at_lowerbdy, NaN]), 

                                    :SO2=>Dict("n"=>[SO2mr * ntot_at_lowerbdy, NaN]),
                                    :SO=>Dict("n"=>[SOmr * ntot_at_lowerbdy, NaN]),

                                    # water mixing ratio is fixed at lower boundary
                                    :H2O=>Dict("n"=>[H2O_lowerbdy, NaN], "f"=>[NaN, 0.]),
                                    # we assume HDO has the bulk atmosphere ratio with H2O at the lower boundary, ~consistent with Bertaux+2007 observations
                                    :HDO=>Dict("n"=>[HDO_lowerbdy, NaN], "f"=>[NaN, 0.]),

                                    # atomic H and D escape solely by photochemical loss to space, can also be mixed downward
                                    :H=> Dict("v"=>[-KoverH_lowerbdy, effusion_velocity(Tn_arr[end], 1.0; zmax, M_P, R_P)],
                                                    #                 ^^^ other options here:
                                                    #                 effusion_velocity(Tn_arr[end], 1.0; zmax) # thermal escape, negligible
                                                    #                 100 # representing D transport to nightside, NOT escape
                                                    #                 NaN # No thermal escape to space, appropriate for global average model
                                              "ntf"=>[NaN, "see boundaryconditions()"]),
                                    :D=> Dict("v"=>[-KoverH_lowerbdy, effusion_velocity(Tn_arr[end], 2.0; zmax, M_P, R_P)], # 
                                                    #                 ^^^ other options here:
                                                    #                  effusion_velocity(Tn_arr[end], 2.0; zmax) # thermal escape, negligible
                                                    #                 100 # representing D transport to nightside, NOT escape
                                                    #                 NaN # No thermal escape to space, appropriate for global average model
                                              "ntf"=>[NaN, "see boundaryconditions()"]),

                                    # # H2 mixing ratio at lower boundary adopted from Yung&DeMore1982 as in Fox&Sung2001
                                    # :H2=>Dict("n"=>[1e-7*ntot_at_lowerbdy, NaN],
                                    #           "v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; zmax)],
                                    #           "ntf"=>[NaN, "see boundaryconditions()"]),
                                    # :HD=>Dict("n"=>[DH*1e-7*ntot_at_lowerbdy, NaN],
                                    #           "v"=>[NaN, effusion_velocity(Tn_arr[end], 3.0; zmax)],
                                    #           "ntf"=>[NaN, "see boundaryconditions()"]),

                                    # unusued neutral boundary conditions
                                    #:O=> Dict("v"=>[-KoverH_lowerbdy, NaN], "f"=>[NaN, 0.#=1.2e6=#]), # no effect on O profile
                                    #:N=>Dict("v"=>[-KoverH_lowerbdy, NaN], "f"=>[NaN, 0.]),
                                    #:NO=>Dict("v"=>[-KoverH_lowerbdy, NaN], #="n"=>[3e8, NaN],=# #="n"=>[5.5e-9*ntot_at_lowerbdy, NaN], =# "f"=>[NaN, 0.]),

                                    # assume no ion loss, appropriate for global average and small observed rates
                                    #:Hpl=>Dict("v"=>[-KoverH_lowerbdy, 0.0 #=effusion_velocity(Ti_arr[end], 1.0; zmax)=#]),#, "f"=>[NaN, 1.6e7]),
                                    #:H2pl=>Dict("v"=>[-KoverH_lowerbdy, 0.0 #=effusion_velocity(Ti_arr[end], 2.0; zmax)=#]),#, "f"=>[NaN, 2e5]),
                                    #:Opl=>Dict("v"=>[-KoverH_lowerbdy, 2e5], ), # "f"=>[NaN, 2.1e8] # tends to cause hollowing out of atmosphere
                                    );

    # add in downward mixing velocity boundary condition for all other species
    auto_speciesbclist = Dict()
    for sp in all_species
        if sp in keys(manual_speciesbclist)
            auto_speciesbclist[sp] = manual_speciesbclist[sp]
        else
            auto_speciesbclist[sp] = Dict("v"=>[-KoverH_lowerbdy, 0.0])
        end
    end

    const speciesbclist = deepcopy(auto_speciesbclist)
end

# ***************************************************************************************************** #
#                                                                                                       #
#                         Set up simulation filenames and define input files                            #
#                                                                                                       #
# ***************************************************************************************************** #

# Crosssection file sources
# -------------------------------------------------------------------
const photochem_data_files = Dict(:CO2=>Dict("main"=>"CO2.dat"), 
                                  :H2O=>Dict("main"=>"h2oavgtbl.dat"), 
                                  :HDO=>Dict("main"=>"HDO.dat"), 
                                  :H2O2=>Dict("main"=>"H2O2.dat"), 
                                  :HDO2=>Dict("main"=>"H2O2.dat"), 
                                  :O3=>Dict("main"=>"O3.dat", "chapman"=>"O3Chap.dat"), 
                                  :O2=>Dict("main"=>"O2.dat", "schr_short"=>"130-190.cf4", "schr_mid"=>"190-280.cf4", "schr_long"=>"280-500.cf4"), 
                                  :H2=>Dict("main"=>"binnedH2.csv"), 
                                  :HD=>Dict("main"=>"binnedH2.csv"), 
                                  :OH=>Dict("main"=>"binnedOH.csv", "O1D+H"=>"binnedOHo1D.csv"), 
                                  :OD=>Dict("main"=>"OD.csv"))

# Filename tags and codes
# -------------------------------------------------------------------

extra_str = seasonal_cycle==true ? "cycle" : "eq"

if seasonal_cycle == true
    # folder naming scheme
    filetag = Dict("temperature"=>"seasons_temp$(extra_str)_Texo=$(Int64(controltemps[3]))_$(results_version)",
                   "water"=>"seasons_water$(extra_str)_$(water_case)_$(water_loc)_$(results_version)",
                   "insolation"=>"seasons_insolation_$(solar_scenario)_$(results_version)",)

    if special_seasonal_case != nothing
        const tag = "seasons_$(special_seasonal_case)_$(results_version)"
    else 
        const tag = filetag[exp_type]
    end
else # not a seasonal cycle experiment
    const tag = "eqrun_$(exp_type)_$(results_version)"
end

# Tags, shortcodes, and filenames
# -------------------------------------------------------------------
# The shortcodes provide unique identifiers for a simulation. Necessary because you end up running the model many times...
const hrshortcode, rshortcode = generate_code(ions_included, round(controltemps[1]), round(controltemps[2]), round(controltemps[3]), water_case, solar_scenario)
const sim_folder_name = "$(hrshortcode)_$(rshortcode)_$(tag)"
const used_rxns_spreadsheet_name = "active_rxns.xlsx"


# ***************************************************************************************************** #
#                                                                                                       #
#                           Algorithm, solver, and float type settings                                  #
#                                                                                                       #
# ***************************************************************************************************** #

# Simulation run time and timestep size  
const season_length_in_sec = seasonal_cycle==true ? season_in_sec : 1e16
const maxlogdt = seasonal_cycle==true ? 5 : 16 # simulation will run until dt = 10^maxlogdt seconds
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "ions+nitrogen"=>[-4, 6], "both"=>[-3, maxlogdt])
const timestep_type = seasonal_cycle==true ? "log-linear" : "dynamic-log" 
    # OPTIONS: "static-log": Logarithmically spaced timesteps that are imposed and don't adjust.
    #                        Should basically never be used, but can be used for testing.
    # "dynamic-log": Logarithmically spaced timesteps which can be dynamically adjusted
    #                as the model runs to improve stability.
    # "log-linear": A hybrid timestep that uses logarithmically-spaced timesteps 
    #               at small elapsed t to get the model going, but then switches
    #               to linearly spaced timesteps of ~1 week. This is used exclusively
    #               with the seasonal model, so that the output can be saved 
    #               at every week and end specifically at 1 season.
    

#                                        Solver algorithm type 
# =======================================================================================================
const problem_type = "Gear" 
    # OPTIONS: 
    # "SS": Julia solver SteadyState.
    # "ODE": Julia ODE solver.
    # "Gear": Preferred solver; Mike's hand-coded Gear method.
# for static timesteps:
const n_steps = 800 # Used with timestep_type="static-log"
# for dynamic timesteps:
const dt_incr_factor = 1.5
const dt_decr_factor = 10

# Sets whether photochemical equilibrium is assumed. Aids in converging ions and neutrals
# together. Generally leave it as is so the code determines it, but you can change it
# if need be
if problem_type == "Gear"
    const assume_photochem_eq = false
else # In using the Julia-provided solvers, it was necessary to assume photochemical equilibrium for short-lived species.
    const assume_photochem_eq = converge_which == "both" ? true : false
end

#                Nitty gritty additions that were included to improve model stability
# =======================================================================================================
# whether to include differentiation terms in Jacobian with respect to electron density 
# or generic thirdbody M. After much testing, these were determined to not be necessary.
const ediff = false # true 
const mdiff = false # true 
const error_checking_scheme = "new"
    # OPTIONS: "new", "old" 
    # We had to work on the method of checking error a lot when dealing with stability issues.
    # The old version has been left in for testing purposes just in case problems are ever 
    # encountered again.

#                                            Float types
# =======================================================================================================
# See CONSTANTS.jl for the setting. no, it's not ideal, but it's the best option 
# right now.
# this means this file must be loaded after CONSTANTS.

# Logic to require Double64 when using the Gear solver. Currently off as Doubles 
# are not being used.
# if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
#     throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in MODEL_SETUP.jl")
# elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
#     println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
# end

# ***************************************************************************************************** #
#                                                                                                       #
#                          Create a parameter dataframe for logging ease                                #
#                                                                                                       #
# ***************************************************************************************************** #

PARAMETERS_GEN = DataFrame(Field=[], Value=[])

push!(PARAMETERS_GEN, ("PLANET", planet));
push!(PARAMETERS_GEN, ("M_P", M_P));
push!(PARAMETERS_GEN, ("R_P", R_P));
push!(PARAMETERS_GEN, ("HRSHORTCODE", hrshortcode));
push!(PARAMETERS_GEN, ("RSHORTCODE", rshortcode));
push!(PARAMETERS_GEN, ("VARIED_PARAM", exp_type))
push!(PARAMETERS_GEN, ("INITIAL_ATM", initial_atm_file));
push!(PARAMETERS_GEN, ("RXN_SOURCE", results_dir*sim_folder_name*"/$(used_rxns_spreadsheet_name)"));
push!(PARAMETERS_GEN, ("IONS", ions_included ));
push!(PARAMETERS_GEN, ("CONVERGE", converge_which));
push!(PARAMETERS_GEN, ("NONTHERMAL_ESC", nontherm));
push!(PARAMETERS_GEN, ("SOLAR_SCENARIO", solar_scenario));
push!(PARAMETERS_GEN, ("SOLARFILE", solarfile));
push!(PARAMETERS_GEN, ("ELECTRON_PROF", e_profile_type));
push!(PARAMETERS_GEN, ("EDIFF", ediff));
push!(PARAMETERS_GEN, ("MDIFF", mdiff));
push!(PARAMETERS_GEN, ("DH", DH));
push!(PARAMETERS_GEN, ("AMBIPOLAR_DIFFUSION_ON", use_ambipolar));
push!(PARAMETERS_GEN, ("MOLEC_DIFFUSION_ON", use_molec_diff));

# Log altitude grid so we can avoid loading this very file when doing later analysis.
PARAMETERS_ALTGRID = DataFrame(Alt=alt, # grid
                               # The following syntax is ugly, and required because the XLSX package won't write columns of different lengths,
                               # so shorter lists must be padded with blank lines.
                               # It appears again below where species lists are written out.
                               non_bdy_layers=[[string(a) for a in non_bdy_layers]..., ["" for i in 1:length(alt)-length(non_bdy_layers)]...],
                               ) 

# Various descriptive things about the altitude grid that are single values, not vectors.
PARAMETERS_ALT_INFO = DataFrame(Field=[], Value=[], Unit=[], Desc=[]);
push!(PARAMETERS_ALT_INFO, ("zmin", zmin, "cm", "Min altitude (altitude of lower boundary)"));
push!(PARAMETERS_ALT_INFO, ("dz", dz, "cm", "Height of a discretized altitude layer"));
push!(PARAMETERS_ALT_INFO, ("zmax", zmax, "cm", "Max altitude (altitude of top boundary)"));
push!(PARAMETERS_ALT_INFO, ("n_all_layers", n_all_layers, "", "Number of discretized altitude layers, including boundary layers (i.e. length(alt))"));
push!(PARAMETERS_ALT_INFO, ("num_layers", num_layers, "", "Number of discretized altitude layers, excluding boundary layers (i.e. length(alt)-2)"));
# push!(PARAMETERS_ALT_INFO, ("upper_lower_bdy", upper_lower_bdy, "cm", "Altitude at which water goes from being fixed to calculated"));
# push!(PARAMETERS_ALT_INFO, ("upper_lower_bdy_i", upper_lower_bdy_i, "", "Index of the line above within the alt grid"));

# Atmospheric conditions.
PARAMETERS_CONDITIONS = DataFrame(Field=[], Value=[], Unit=[]);
push!(PARAMETERS_CONDITIONS, ("SZA", SZA, "deg"));
push!(PARAMETERS_CONDITIONS, ("TSURF", controltemps[1], "K"));
push!(PARAMETERS_CONDITIONS, ("TMESO", controltemps[2], "K"));
push!(PARAMETERS_CONDITIONS, ("TEXO", controltemps[3], "K"));
push!(PARAMETERS_CONDITIONS, ("MEAN_TEMPS", join(meantemps, " "), "K"));
push!(PARAMETERS_CONDITIONS, ("WATER_MR", water_mixing_ratio, "mixing ratio"));
push!(PARAMETERS_CONDITIONS, ("WATER_CASE", water_case, "whether running with 10x, 1/10th, or standard water in middle/upper atmo"));

L = max(length(all_species), length(neutral_species), length(ion_species), length(no_chem_species), length(no_transport_species), length(Jratelist))
PARAMETERS_SPLISTS = DataFrame(AllSpecies=[[string(a) for a in all_species]..., ["" for i in 1:L-length(all_species)]...], 
                               Neutrals=[[string(n) for n in neutral_species]..., ["" for i in 1:L-length(neutral_species)]...], 
                               Ions=[[string(i) for i in ion_species]..., ["" for i in 1:L-length(ion_species)]...],
                               NoChem=[[string(nc) for nc in no_chem_species]..., ["" for i in 1:L-length(no_chem_species)]...],
                               NoTransport=[[string(nt) for nt in no_transport_species]..., ["" for i in 1:L-length(no_transport_species)]...],
                               Jratelist=[[string(j) for j in Jratelist]..., ["" for i in 1:L-length(Jratelist)]...]);
PARAMETERS_SOLVER = DataFrame(Field=[], Value=[]);
PARAMETERS_XSECTS = DataFrame(Species=[], Description=[], Filename=[]);
PARAMETERS_BCS = DataFrame(Species=[], Type=[], Lower=[], Upper=[]);

# LOG THE TEMPERATURES
PARAMETERS_TEMPERATURE_ARRAYS = DataFrame(Neutrals=Tn_arr, Ions=Ti_arr, Electrons=Te_arr); 