################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Global constants, simulation parameters.
# 
# Eryn Cangi
# Created December 2019
# Last edited: 2023
# Currently tested for Julia: 1.85
################################################################################

using DataFrames

# **************************************************************************** #
#                                                                              #
#                         Main simulation parameters                           #
#                                                                              #
# **************************************************************************** #

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                      !!!!! SUPER IMPORTANT !!!!!                       !! #
# !!     !!! Modify the following items each time you run the model !!!     !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# Basic model parameters
const optional_logging_note = "test" # Brief summary of simulation goal
const seasonal_cycle = true # for paper 3
const results_version = "v99"  # Helps keep track of attempts if you need to keep changing things

# Input files ---------------------------------------------------- #
const initial_atm_file = "INITIAL_GUESS.h5" 
const reaction_network_spreadsheet = code_dir*"REACTION_NETWORK.xlsx" # "REACTION_NETWORK_MIN_IONOSPHERE.xlsx" #

# DEFINE THE MAIN INDEPENDENT VARIABLE  ---------------------------- #
# This is the one you'd like to change and see how the atmosphere responds.
# If you just want to run the model and get some basic output, 
# just select "temperature".
const season_exp = "temperature" # "water" #  "insolation"#

# SOLAR CASE ---------------------------------------------------- #
const SZA = 60  # Puts the model at dayside mean
const solarcyc = "mean"# 
# AU OPTIONS: "perihelion" # "aphelion" #  "meandist" # Defined for solar mean. TODO: Program the solar spectrum scaling in Julia and set AU as a parameter
# SOLAR CYCLE OPTIONS: "mean" # "max" # "min" # Defined for mean AU. (solar spectra varies; hand collected by Eryn)

# TEMPERATURE CASES --------------------------------------------- #
const tempcyc = "max" # "mean" # "min"#  
# Also an option which will set all temperatures to the same value: "isothermal" 

# WATER CASES --------------------------------------------------- #
# Amount of water in the atmosphere
const water_case = "standard" #"low" #"high" #  
const water_loc = "mesosphere" # Location to alter water if selecting "low" or "high" water case "loweratmo" #  "everywhere" #   

# OPTIONAL: Extra parcel of water to simulate effect of dust storms
const dust_storm_on = false
const H2O_excess = 250 # excess H2O in ppm
const HDO_excess = 0.350 # excess HDO in ppm (divide by 1000 to get ppb)
if dust_storm_on==true
    const excess_peak_alt = 60
end

# Ion chemistry and non-thermal escape ---------------------------- #
const ions_included = true
const nontherm = ions_included==true ? true : false   # whether to do non-thermal escape
const converge_which = "both" # "ions" "neutrals"
const e_profile_type = ions_included==true ? "quasineutral" : "none" # "O2+" # "constant"# Electrons calculated as sum of ions (quasineutral), sum of O2+, constant, or none

# Plotting option ------------------------------------------------- #
const make_P_and_L_plots = false # Turn off to save several minutes of runtime if you don't need to check for equilibrium.

# Lesser-used options --------------------------------------------- #
const do_chem = true   # Turning this or next one of will toggle chemistry or transport.
const do_trans = true  # Often useful for troubleshooting or converging new atmospheres.
const adding_new_species = false # true#  set to true if introducing a new species.
const make_new_alt_grid = false
const use_nonzero_initial_profiles = true   # Can be turned off to initialize species from zero. 


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                      !!!!!    END CHECK    !!!!!                       !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# **************************************************************************** #
#                                                                              #
#                    Species name lists and J rate lists                       #
#                                                                              #
# **************************************************************************** #

const fixed_species = [:Ar] # here you may enter any species that you want to be completely fixed (no updates to densities from chemistry or transport)

remove_unimportant = true # Whether to use a slightly smaller list of species and reactions (removing minor species that Roger had in his model)
unimportant = [:CNpl,:HCNpl,:HCNHpl,:HN2Opl,:NH2pl,:NH3pl,:N2Opl,:NO2pl,:CH,:CN,:HCN,:HNO,:NH,:NH2,:HD2pl]#:N2O,:NO2

# Neutrals --------------------------------------------------------------------
const orig_neutrals = [:Ar, :CO, :CO2, :H, :H2, :H2O, :H2O2, 
                       :HO2, :HOCO, :N2, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD,

                       # Turn these off for minimal ionosphere:
                       :C, :DCO, :HCN, :HCO, :N, :NO, :Nup2D, 
                       ]; 
const conv_neutrals = remove_unimportant==true ? setdiff(orig_neutrals, unimportant) : orig_neutrals
const new_neutrals = [];
const neutral_species = [conv_neutrals..., new_neutrals...];

# Ions ------------------------------------------------------------------------- #
const orig_ions = [:CO2pl, :HCO2pl, :Opl, :O2pl, # Nair minimal ionosphere 
                   :Arpl, :ArHpl, :ArDpl, 
                   :Cpl, :CHpl,  :COpl, 
                   :Hpl, :Dpl, :H2pl, :HDpl, :H3pl, :H2Dpl, :HD2pl, 
                   :H2Opl,  :HDOpl, :H3Opl, :H2DOpl, 
                   :HO2pl, :HCOpl, :DCOpl, :HOCpl, :DOCpl, :DCO2pl, 
                   :HNOpl,   
                   :Npl, :NHpl, :N2pl, :N2Hpl, :N2Dpl, :NOpl,
                   :OHpl, :ODpl];
const new_ions = [];
const ion_species = remove_unimportant==true ? setdiff([orig_ions..., new_ions...], unimportant) : [orig_ions..., new_ions...]

# Full species list ------------------------------------------------------------- #
const all_species = [neutral_species..., ion_species...];

# Sorted name lists created here ------------------------------------------------- #
sort!(all_species)
sort!(neutral_species)
sort!(ion_species)

# Photolysis and Photoionization rate symbol lists ----------------------------#

const conv_Jrates, newJrates = format_Jrates(reaction_network_spreadsheet, all_species, "Jratelist"; ions_on=ions_included, hot_atoms=nontherm)
const Jratelist = [conv_Jrates..., newJrates...];

# These dictionaries specify the species absorbing a photon for each J rate, and the products of the reaction.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])

# **************************************************************************** #
#                                                                              #
#                             Water settings                                   #
#                                                                              #
# **************************************************************************** #

# Water mixing ratios to use for paper3 runs.
const water_MRs = Dict("loweratmo"=>Dict("standard"=>1.3e-4, "low"=>0.65e-4, "high"=>2.6e-4), 
                       "mesosphere"=>Dict("standard"=>1.3e-4, "high"=>1.3e-4, "low"=>1.3e-4), 
                       "everywhere"=>Dict("standard"=>1.3e-4, "high"=>1.3e-4, "low"=>1.3e-4))
const water_mixing_ratio = water_MRs[water_loc][water_case]
const reinitialize_water_profile = seasonal_cycle==true ? false : true # should be off if trying to run simulations for Mars seasons
const update_water_profile = seasonal_cycle==true ? true : false # this is for modifying the profile during cycling, MAY be fixed?
const modified_water_alts = "below fixed point"

# altitude at which to add the extra water -- applies to both dust storm parcels and the tanh profile
const add_water_alt_opts = Dict("low"=>45, "standard"=>60, "high"=>65)
const f_fac_opts = Dict("low"=>0.005, "standard"=>10, "high"=>100) # a parameter named f which helps manipulate the water profile. Not related to any other f.


# **************************************************************************** #
#                                                                              #
#                        Temperature profile construction                      #
#                                                                              #
# **************************************************************************** #

# Other simulation controls ------------------------------------------------------------------------- #
const Texo_opts = Dict("min"=>175., "mean"=>225., "max"=>275., 
                       "equinox"=>225., "aphelion"=>225., "peihelion"=>225.) # Since these solar inputs are defined for solar mean conditions, use solar mean temp for all.

if seasonal_cycle == true
    println("Simulating an annual cycle")
    const controltemps = [230., 130., Texo_opts["mean"]]

    if season_exp=="temperature"
        const controltemps[3] =  Texo_opts[tempcyc]
    end
else
    const controltemps = [230., 130., Texo_opts[solarcyc]]
end

if tempcyc=="isothermal"
    const controltemps = [225., 225., 225.]
    const meantemps = [225., 225., 225.] # Used for saturation vapor pressure. DON'T CHANGE!
else
    const meantemps = [230., 130., Texo_opts["mean"]] # Used for saturation vapor pressure. DON'T CHANGE!
end

const T_surf = controltemps[1]
const T_meso = controltemps[2]
const T_exo = controltemps[3]

T_array_dict = T(T_surf, T_meso, T_exo; alt);
const Tn_arr = T_array_dict["neutrals"]
const Ti_arr = T_array_dict["ions"]
const Te_arr = T_array_dict["electrons"]

const Tplasma_arr = Ti_arr .+ Te_arr;
const Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
const Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr)
const Tn_meanSVP = T(meantemps...; alt)["neutrals"]; # Needed for boundary conditions.


# **************************************************************************** #
#                                                                              #
#                             Boundary conditions                              #
#                                                                              #
# **************************************************************************** #
const H2Osat = map(x->Psat(x), Tn_meanSVP) # Using this function keeps SVP fixed 
const HDOsat = map(x->Psat_HDO(x), Tn_meanSVP)

const speciesbclist=Dict(:CO2=>Dict("n"=>[2.1e17, NaN], "f"=>[NaN, 0.]),
                        :Ar=>Dict("n"=>[2.0e-2*2.1e17, NaN], "f"=>[NaN, 0.]),
                        :N2=>Dict("n"=>[1.9e-2*2.1e17, NaN], "f"=>[NaN, 0.]),
                        :H2O=>Dict("n"=>[H2Osat[1], NaN], "f"=>[NaN, 0.]), # bc doesnt matter if H2O fixed
                        :HDO=>Dict("n"=>[HDOsat[1], NaN], "f"=>[NaN, 0.]),
                        :O=> Dict("f"=>[0., 1.2e8]),
                        :H2=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),  # velocities are in cm/s
                        :HD=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 3.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                        :H=> Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 1.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                        :D=> Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(Tn_arr[end], 2.0; zmax)], "ntf"=>[NaN, "see boundaryconditions()"]),
                       );

# **************************************************************************** #
#                                                                              #
#              Set up simulation filenames and define input files              #
#                                                                              #
# **************************************************************************** #

if seasonal_cycle == true
    if solarcyc in ["equinox", "aphelion", "perihelion"]
        const solarfile = "marssolarphotonflux_$(solarcyc).dat"
    else 
        const solarfile = "marssolarphotonflux_solar$(solarcyc).dat"
    end
    extra_str = seasonal_cycle==true ? "cycle" : "eq"

    if season_exp=="temperature"
        const tag = "paper3_temp$(extra_str)_Texo=$(Int64(controltemps[3]))_$(results_version)"

    elseif season_exp=="insolation"
        # const solarfile = "marssolarphotonflux_$(solarcyc).dat"
        const tag = "paper3_insolation_$(solarcyc)_$(results_version)"

    elseif season_exp=="water"
        const tag = "paper3_water$(extra_str)_$(water_case)_$(water_loc)_$(results_version)"

    else
        throw("Simulation type is seasonal cycle but no valid testing parameter specified")
    end
else
    if solarcyc in ["equinox", "aphelion", "perihelion"]
        const solarfile = "marssolarphotonflux_$(solarcyc).dat"
    else 
        const solarfile = "marssolarphotonflux_solar$(solarcyc).dat"
    end
    const tag = "s$(solarcyc)_$(results_version)"
end

# Tags, shortcodes, and filenames
# The shortcodes provide unique identifiers for a simulation. Necessary because you end up running the model many times...
const hrshortcode, rshortcode = generate_code(ions_included, controltemps[1], controltemps[2], controltemps[3], water_case, solarcyc)
const sim_folder_name = "$(hrshortcode)_$(rshortcode)_$(tag)"
const final_atm_file = "final_atmosphere.h5"

# **************************************************************************** #
#                                                                              #
#                             Algorithm settings                               #
#                                                                              #
# **************************************************************************** #

# Tolerances
const rel_tol = 1e-6
const abs_tol = 1e-12 

# Simulation run time and timestep size  
const season_length_in_sec = seasonal_cycle==true ? 1.4838759e7 : 1e16
const maxlogdt = seasonal_cycle==true ? 5 : 16 # simulation will run until dt = 10^maxlogdt seconds
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-3, maxlogdt])
const timestep_type = seasonal_cycle==true ? "log-linear" : "dynamic-log" # "static-log"  # Static-log should basically never be used, but can be used for testing.

# Solver algorithm type 
const problem_type = "Gear" #"SS" #"ODE" #  
# for static timesteps:
const n_steps = 800 # for static case
# for dynamic timesteps:
const dt_incr_factor = 1.5
const dt_decr_factor = 10

# whether to include differentiation terms in Jacobian with respect to electron density or generic thirdbody M. 
# After much testing, these were determined to not be necessary.
const ediff = false # true 
const mdiff = false # true 
const error_checking_scheme = "new" #"old" 

# Sets whether photochemical equilibrium is assumed. Aids in converging ions and neutrals
# together. Generally leave it as is so the code determines it, but you can change it
# if need be
if problem_type == "Gear"
    const assume_photochem_eq = false
else # In using the Julia-provided solvers, it was necessary to assume photochemical equilibrium for short-lived species.
    const assume_photochem_eq = converge_which == "both" ? true : false
end

# Check that float type is properly set
# if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
#     throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in CUSTOMIZATIONS.jl")
# elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
#     println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
# end

# **************************************************************************** #
#                                                                              #
#                      Miscellaneous logical groupings                         #
#                                                                              #
# **************************************************************************** #
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
const D_bearing_species = get_deuterated(all_species)
const D_ions = get_deuterated(ion_species) #[s for s in ion_species if occursin('D', string(s))];
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];

# Sort name lists created here ------------------------------------------------- #
sort!(D_bearing_species)
sort!(D_ions)
sort!(N_neutrals)

# **************************************************************************** #
#                                                                              #
#                    Define short- and long-lived species                      #
#                                                                              #
# **************************************************************************** #

# Short lived species, whose chemical lifetime is << diffusion timescale ------- #
const short_lived_species = [];# technically shortlived but count as longlived: :CH, :HCO, :HO2, :O3, :OH, :O1D, :DO2, :OD...
if assume_photochem_eq
    append!(short_lived_species, [:NO2, :CN, :HNO, :NH, :NH2, :C, :CH])
    append!(short_lived_species, ion_species)
end

# Long lived species ------------------------------------------------------------ #
const long_lived_species = setdiff(all_species, short_lived_species)

# Sort name lists created here ------------------------------------------------- #
sort!(short_lived_species)
sort!(long_lived_species)

# **************************************************************************** #
#                                                                              #
#               Species participating in chemistry and transport               #
#                                                                              #
# **************************************************************************** #

# Non-participants ------------------------------------------------------------- # 
const no_chem_species = []; 
const no_transport_species = [];

# This will append any species that you have picked to be completely fixed. 
for fs in fixed_species
    push!(no_chem_species, fs)
    push!(no_transport_species, fs)
end

if converge_which == "neutrals"
    append!(no_chem_species, union(conv_ions, N_neutrals)) # This is because the N chemistry is intimiately tied up with the ions.
    append!(no_transport_species, union(conv_ions, N_neutrals, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, setdiff(conv_neutrals, N_neutrals))
    append!(no_transport_species, setdiff(conv_neutrals, N_neutrals))
elseif converge_which == "both"
    append!(no_transport_species, short_lived_species)
end

# Participants ------------------------------------------------------------------ #
const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species --------------------------------------------------- #
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)

# Sort name lists created here ------------------------------------------------- #
sort!(active_species)
sort!(inactive_species)
sort!(active_longlived)
sort!(active_shortlived)
sort!(chem_species)
sort!(transport_species)
sort!(no_chem_species)
sort!(no_transport_species)

# **************************************************************************** #
#                                                                              #
#              Misc. things that depend on things defined above                #
#                                                                              #
# **************************************************************************** #

# Annoyingly, these have to be here because they depend on other things defined above.
# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict(vcat([s=>"--" for s in setdiff(D_bearing_species, [:HD2pl])], [:HD2pl=>":", :Nup2D=>"-."]) )

# Species-specific scale heights - has to be done here instead of in the param file
const Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)]; molmass) for sp in all_species])

# To allow water to be active in the upper atmosphere but not the lower atmosphere, we need 
# its position within the active species vector 
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

# Altitude at which water transitions from fixed to freely solved for
# First we have to calculate a few intermediaries.
H2Osatfrac = H2Osat ./ map(z->n_tot(get_ncurrent(initial_atm_file), z; all_species, n_alt_index), alt)  # get SVP as fraction of total atmo
const upper_lower_bdy = alt[something(findfirst(isequal(minimum(H2Osatfrac)), H2Osatfrac), 0)] # in cm
const upper_lower_bdy_i = n_alt_index[upper_lower_bdy]  # the uppermost layer at which water will be fixed, in cm

# **************************************************************************** #
#                                                                              #
#               CREATE A PARAMETER DATAFRAME FOR LOGGING EASE                  #
#                                                                              #
# **************************************************************************** # 

PARAMETERS_GEN = DataFrame(Field=[], Value=[])

push!(PARAMETERS_GEN, ("HRSHORTCODE", hrshortcode));
push!(PARAMETERS_GEN, ("RSHORTCODE", rshortcode));
push!(PARAMETERS_GEN, ("INITIAL_ATM", initial_atm_file));
push!(PARAMETERS_GEN, ("RXN_SOURCE", reaction_network_spreadsheet));
push!(PARAMETERS_GEN, ("IONS", ions_included ));
push!(PARAMETERS_GEN, ("CONVERGE", converge_which));
push!(PARAMETERS_GEN, ("NONTHERMAL_ESC", nontherm));
push!(PARAMETERS_GEN, ("SOLARCYCLE", solarcyc));
push!(PARAMETERS_GEN, ("SOLARFILE", solarfile));
push!(PARAMETERS_GEN, ("ELECTRON_PROF", e_profile_type));
push!(PARAMETERS_GEN, ("EDIFF", ediff));
push!(PARAMETERS_GEN, ("MDIFF", mdiff));
push!(PARAMETERS_GEN, ("DH", DH));

PARAMETERS_CONDITIONS = DataFrame(Field=[], Value=[], Unit=[]);

push!(PARAMETERS_CONDITIONS, ("SZA", SZA, "deg"));
push!(PARAMETERS_CONDITIONS, ("TSURF", T_surf, "K"));
push!(PARAMETERS_CONDITIONS, ("TMESO", T_meso, "K"));
push!(PARAMETERS_CONDITIONS, ("TEXO", T_exo, "K"));
push!(PARAMETERS_CONDITIONS, ("MEAN_TEMPS", join(meantemps, " "), "K"));
push!(PARAMETERS_CONDITIONS, ("WATER_MR", water_mixing_ratio, "mixing ratio"));
push!(PARAMETERS_CONDITIONS, ("WATER_CASE", water_case, "whether running with 10x, 1/10th, or standard water in middle/upper atmo"));

waterbdy = :H2O in fixed_species ? 250 : upper_lower_bdy/1e5
push!(PARAMETERS_CONDITIONS, ("WATER_BDY", waterbdy, "km"))

# This is so ugly because the XLSX package won't write columns of different lengths, so I have to pad all the shorter lists
# with blanks up to the length of the longest list and also transform all the symbols into strings. 
L = max(length(all_species), length(neutral_species), length(ion_species), length(no_chem_species), length(no_transport_species))
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
