################################################################################
# INPUT_PARAMETERS.jl
# DESCRIPTION: Controls key variables that set up the model, usually modified at 
# each run by the user. Planet to model, logging notes, filenames, whether it's 
# a seasonal cycle or equilibrium run, which atmospheric parameter to vary, 
# solar input, water parameters, species modeled, algorithm tolerances, and stuff
# to control which species densities are updated due to chemistry or transport.
# 
# Eryn Cangi
# Created April 2024
# Last edited: December 2025
# Currently tested for Julia: 1.12.3
################################################################################

# Set the planet 
# =======================================================================================================
const planet = "Venus"
    # OPTIONS: "Mars", "Venus"

# Input and output files, directory
# =======================================================================================================
const results_dir = code_dir*"../Results_$(planet)/"
const initial_atm_file = "$(planet)-Inputs/INITIAL_GUESS_VENUS_oUT0ZbGN.h5"  # File to use to initialize the atmosphere.
    # OPTIONS: 
    # INITIAL_GUESS_MARS.h5 --> Basic Mars starting file.
    # INITIAL_GUESS_MARS_bxz4YnHk.h5 --> A Mars atmosphere that includes N2O, NO2, and their ions;
    #                                    not particularly motivated by any present-day data.                                         
    # INITIAL_GUESS_VENUS_vGFd5b0a.h5 --> Best Venus initial atmosphere  without sulfur/chlorine
    # INITIAL_GUESS_VENUS_oUT0ZbGN.h5 --> Venus initial atmosphere with basic chlorine and sulfur species included.
const final_atm_file = "final_atmosphere.h5"
const reaction_network_spreadsheet = code_dir*"$(planet)-Inputs/REACTION_NETWORK_$(uppercase(planet)).xlsx"
    # OPTIONS: "REACTION_NETWORK_MIN_IONOSPHERE.xlsx", code_dir*"REACTION_NETWORK_$(uppercase(planet)).xlsx"

# Descriptive attributes of this model run
# =======================================================================================================
const short_summary = "co+od_ratecoeff_eq_co+oh" 
      # a short string that will be added to the results folder, to jog your memory of what you did.
      # Recommended not to include spaces. May be blank.
const logged_long_description = "Test a rate coefficient of +0% (*1) of H equivalent for CO + OD -> CO2 + D" 
      # Brief summary of simulation goal, will be in the log file but not the results folder name.
const results_version = "v0"  
      # Helps keep track of attempts if you need to keep changing things. Will be appended to results
      # folder name.

# Set the modifiable atmospheric parameters
# =======================================================================================================
const seasonal_cycle = false  # Whether to simulate one season or run to equilibrium.
    # OPTIONS: 
    # true (short, one season)
    # false (run to equilibrium)
const exp_type = "temperature" # Goes into the filename and log to show which parameter was varied.
    # OPTIONS: 
    # temperature: varying the exobase temperature
    # water: varying the amount of water in the mesosphere 
    # insolation: varying insolation due to either solar cycle or orbital distance
    # all: all three are modified; use only with special_seasonal_case != nothing
const special_seasonal_case = nothing
    # Special case for running the seasonal model with muliple parameters being modified.
    # At the request of a reviewer for Cangi+ 2024. Note: You still have to change the individual parameters below.
    # OPTIONS: 
    # nothing - if not doing a seasonal run (note, this is not a string, it's the Julia nothing, similar to Python None)
    # "inclusive-peri" - Perihelion conditions, modifying all of exobase temperature, water in mesosphere, insolation.
    # "inclusive-mean" - Same as above but for the time halfway between perihelion and aphelion, either side.
    # "inclusive-ap" - Same as above, but for aphelion.
    
# Temperature
# -------------------------------------------------------------------
# This sets just the exobase temperature. Typically, the exobase temperature tracks the
# solar cycle. You don't HAVE to make them match, but if you don't, you need to be able to 
# interpret that choice and justify your choices when you write up the paper.
const temp_scenario = "mean"  # Temperature selection for the seasonal model run.
    # OPTIONS: 
    # min-P2, mean-P2, max-P2: uses exobase temps from Cangi+2023 (190, 210, 280) K.
    # mean, min, max: temps as in Cangi+ 2024 (175, 225, 275) K - goal for this one was evenly spaced.
    # isothermal (this will set a constant temperature at all alts - 225 K)

# Solar case
# -------------------------------------------------------------------
const SZA = 60  # Puts the model at dayside mean. Enter in degrees please.
const solar_scenario = "solarmean" 
    # Solar scenario definition. You can choose from different Mars-sun distances or parts of the solar cycle.
    # ORBITAL DISTANCE OPTIONS: "perihelion" #  "meansundist" # "aphelion"; these are defined at solar mean. 
    #     NOTE: these options are only available for Mars at present.
    # SOLAR CYCLE OPTIONS: "solarmean" # "solarmax" # "solarmin"; these are defined at the mean AU. 
    # NOTE: Not currently possible to mix orbital distance and solar cycle as I would have to generate more spectra.

    # TODO: Program the solar spectrum scaling in Julia and set AU as a parameter
const solarfile = "$(planet)-Inputs/$(lowercase(planet))solarphotonflux_$(solar_scenario).dat"

# Water
# -------------------------------------------------------------------
const water_case = "standard" # Amount of water in the atmosphere
    # OPTIONS: low, standard, high
    # for Venus: standard
const water_loc = "mesosphere" # Location to modify water abundance, if selecting "low" or "high" water case 
    # OPTIONS: "loweratmo", "everywhere", "mesosphere"

# VENUS OPTIONS:
const venus_special_water = false# true
# MARS DUST OPTIONS:
const dust_storm_on = false  # This adds a parcel of water at some altitude, sort of like if there was a dust storm. Haven't published on this.
    # OPTIONS: True, false
const H2O_excess = 250  # excess H2O in ppm
    # OPTIONS: Typical amounts are 100-300 ppm.
const HDO_excess = 0.350 # excess HDO in ppm (divide by 1000 to get ppb)
    # OPTIONS: 0.350 is reasonable


# Control which species (atoms, molecules) are modeled 
# =======================================================================================================

const ions_included = true
const converge_which =  "both" 
    # OPTIONS: "ions", "neutrals", "both", "ions+nitrogen"

# Species lists
# -------------------------------------------------------------------
# Convention: Alphabetized, except D-bearing species should be mixed in after their H-bearing isotopologue. 
# this makes it easier to see which species have isotopologues.
# 
# WHAT DO THESE LISTS MEAN? 
# 
# new_[neutrals, ions]: These are species you want to add to the modeled atmosphere on THIS RUN.
#                       They should not have a density vector in the initial guess file.
#                       Their initial density guess will be either all 0, or loaded from a file 
#                       depending on what you specify for the use_nonzero_initial_profiles variable.
#                       Note that after you complete a successful run with new species,
#                       you must modify this file manually and enter those species into the 
#                       conv_[neutrals,ions] variable. No way around it!
#
# conv_[neutrals,ions]: These are species which have already been incorporated into the model atmosphere.
#                       They should already have a density vector in your initial guess file. 
#                       Typically you won't need to change these at each run, UNLESS you added new species 
#                       on the previous run.
const new_neutrals = [];
const new_ions = [];

const conv_neutrals = Dict("Mars"=>[:Ar, :C, :CO, :CO2, # Argon and carbon species
                                    :H, :D, :H2, :HD, :H2O, :HDO,  # H and D species
                                    :HCO, :DCO, :HO2, :DO2,        
                                    :H2O2, :HDO2, :HOCO, :DOCO, 
                                    :N, :N2, :NO, :Nup2D, :N2O, :NO2, # Nitrogen species
                                    :O, :O1D, :O2, :O3, :OH, :OD], # Oxygen species
                           "Venus"=>[:Ar, :C, :CO, :CO2, 
                                     :Cl, :ClO, :ClCO, :HCl, :DCl,  # Chlorine species
                                     :H, :D, :H2, :HD, :H2O, :HDO,  # H and D species
                                     :HCO, :DCO, :HO2, :DO2,        
                                     :H2O2, :HDO2, :HOCO, :DOCO, 
                                     :N, :N2, :NO, :Nup2D, :N2O, :NO2,
                                     :O, :O1D, :O2, :O3, :OH, :OD,
                                     :S, :SO, :SO2, :SO3, :H2SO4, :HDSO4] # Sulfur species
                           ); 

const conv_ions = Dict("Mars"=>[:Arpl, :ArHpl, :ArDpl, 
                                :Cpl, :CHpl, :COpl, :CO2pl, 
                                :Dpl, :DCOpl, :DOCpl, :DCO2pl, 
                                :Hpl,  :H2pl, :HDpl, :H3pl, :H2Dpl, 
                                :H2Opl, :HDOpl, :H3Opl, :H2DOpl, 
                                :HO2pl, :HCOpl, :HCO2pl, :HOCpl, :HNOpl,   
                                :Npl, :NHpl, :N2pl, :N2Hpl, :N2Dpl, :NOpl, :N2Opl, :NO2pl,
                                :Opl, :O2pl, :OHpl, :ODpl],
                       "Venus"=>[:Arpl, :ArHpl, :ArDpl, 
                                :Cpl, :CHpl, :COpl, :CO2pl, 
                                :Dpl, :DCOpl, :DOCpl, :DCO2pl, 
                                :Hpl,  :H2pl, :HDpl, :H3pl, :H2Dpl, 
                                :H2Opl, :HDOpl, :H3Opl, :H2DOpl, 
                                :HO2pl, :HCOpl, :HCO2pl, :HOCpl, :HNOpl,   
                                :Npl, :NHpl, :N2pl, :N2Hpl, :N2Dpl, :NOpl, :N2Opl, :NO2pl,
                                :Opl, :O2pl, :OHpl, :ODpl]
                      );

# More specific settings for controling the modeling of species
# -------------------------------------------------------------------
# Chemical species which should never update their densities, but may be chemical reactants.
const dont_compute_chemistry = []
const dont_compute_transport = []
const dont_compute_either_chem_or_transport = []  
    # OPTIONS: Any species included in the model. 
if planet=="Mars" # To avoid convergence problems
    append!(dont_compute_either_chem_or_transport,[:Ar])
elseif planet=="Venus"
    append!(dont_compute_chemistry,[:Ar])
end
    
const assume_photochem_eq = false # whether to turn on photochemical equilibrium for short-lived species

# Turn plots off and on
# =======================================================================================================
const make_P_and_L_plots = true  # Makes a 3-panel plot showing production and loss due to 1) chemistry, 
                                 # 2) transport, 3) sum of both. Turn off to save several minutes of 
                                 # runtime if you don't need to check for equilibrium.
    # OPTIONS: True, false

# Algorithm tolerances
# =======================================================================================================
const rel_tol = planet == "Venus" ? 1e-6 : 1e-3  # Venus: original values, Mars: relaxed for multicolumn stability
const abs_tol = planet == "Venus" ? 1e-12 : 1e-9  # Venus: original values, Mars: relaxed for multicolumn stability

# Helpful options for adding new things to the model 
# =======================================================================================================
const do_chem = true   # Turning this or next one of will toggle chemistry or transport.
const do_trans = true  # Often useful for troubleshooting or converging new atmospheres.
const adding_new_species = false
const make_new_alt_grid = false  # Set to true if extending the altitude grid. TODO: Need to re-write that code.
const use_nonzero_initial_profiles = true
    # OPTIONS: 
    # true -- uses initial guess densities for species based on previous model output.
    # false -- sets species to zero density and lets the chemistry and transport build them up.
const use_ambipolar = true # Toggle ambipolar diffusion for ions.
const use_molec_diff = true # Toggle molecular diffusion. If turned off, eddy diffusion remains active.

# Number of vertical columns in the simulation. Set this to 1 for a single-column run or >1 for a multicolumn model.
const n_horiz = 2

# Cross-terminator (day-to-night) thermospheric transport timescales at Venus are around 23 to 44 hours,
# corresponding to wind speeds of 230 to 120 m/s (2.3e4 to 1.2e4 cm/s), with 30 hours being typical.
# This assumes semi-circumference of Venus ~19,000 km as the characteristic width for day-to-night transport.

# Horizontal column width in cm. This determines the physical scale of horizontal transport.
# For day-night transport setups (n_horiz=2), use larger values for physically realistic transport rates.
# Constant-mode width is surface half-circumference πR using the same R values
# as in MODEL_SETUP.jl (Mars=3396e5 cm, Venus=6050e5 cm).
const horiz_column_width = π * Dict("Mars"=>3396e5, "Venus"=>6050e5)[planet]

# If true, horizontal transport uses altitude-dependent width dx(z) = π(R + z).
# If false, it uses constant `horiz_column_width`.
const use_altitude_dependent_horiz_dx = true

# Horizontal transport timescale in hours. This determines the wind speed via: wind_speed = horiz_column_width / (timescale * 3600)
# For Venus: 23-44 hours corresponds to 230-120 m/s wind speeds
# For Mars: Set to 0 for no horizontal transport
const horiz_transport_timescale = planet == "Venus" ? 30.0 : 0.0  # Baseline shared value (hours)
# Split neutral/ion timescales to reflect faster ion coupling on Venus (literature: ~10–20 h ions, ~20–40 h neutrals)
const horiz_transport_timescale_neutral = planet == "Venus" ? 30.0 : horiz_transport_timescale
const horiz_transport_timescale_ion = planet == "Venus" ? 15.0 : horiz_transport_timescale

# Horizontal wind speed in cm/s calculated from timescale: wind_speed = width / (timescale * 3600)
# This is used to initialize wind profiles in `MODEL_SETUP.jl`
const horiz_wind_speed = horiz_transport_timescale > 0 ? horiz_column_width / (horiz_transport_timescale * 3600) : 0.0
const horiz_wind_speed_neutral = horiz_transport_timescale_neutral > 0 ? horiz_column_width / (horiz_transport_timescale_neutral * 3600) : 0.0
const horiz_wind_speed_ion = horiz_transport_timescale_ion > 0 ? horiz_column_width / (horiz_transport_timescale_ion * 3600) : 0.0

# Whether to allow horizontal transport between columns. When set to `false`
# the model does not compute horizontal advection or diffusion, matching the
# behaviour of the single-column set-up even when multiple columns are present.
const enable_horiz_transport = true

# Optional toggle for horizontal diffusion terms only. Keep off by default so
# horizontal exchange is advective unless explicitly enabled.
const enable_horiz_diffusion = false

# If true, horizontal neighbors wrap periodically (column 1 is behind column n_horiz).
# If false, the model uses explicit back/front edge boundary conditions.
const horiz_transport_cyclic = true

# =======================================================================================================
# Horizontal Column Scenario Configuration
# =======================================================================================================
# Define the characteristics of each horizontal column.
# Use a named mapping first (e.g., "day", "night"), then derive the ordered
# per-column list used by the solver loops.
# Each entry specifies: temperature scenario, solar zenith angle, and optional
# boundary condition modifiers.
# Optional keys:
#   Texo_override (numeric K) overrides Texo_key for that column
#   CO2_lowerbdy_vmr sets a per-column lower boundary CO2 mixing ratio
# Available Texo_key options (from Texo_opts in MODEL_SETUP.jl):
#   Mars: "min-P2", "mean-P2", "max-P2", "min-P3", "mean-P3", "max-P3"
#   Venus: "min", "mean", "max"
#
# SZA notes: When SZA > 90°, the sun is below the horizon and solar flux will be set to zero

const horiz_column_scenario_by_name = if n_horiz == 1
    # Single column: standard uniform setup
    OrderedDict(
        "uniform" => Dict(
            "Texo_key" => planet == "Venus" ? "mean" : "mean-P2",  # Options: "min", "mean", "max", "min-P2", "mean-P2", "max-P2", etc.
            "SZA" => SZA                                           # Solar zenith angle (degrees)
        )
    )
elseif planet == "Venus" && n_horiz == 2
    # Venus day-night test: 2 columns
    OrderedDict(
        "day" => Dict(
            "Texo_key" => "max",       # Keeps the key for reference; Texo_override below is applied when present
            "Texo_override" => 340.0,  # Dayside Texo (K); typical VTGCM/VTS3 dayside 320-360 K
            "CO2_lowerbdy_vmr" => 0.965,
            "SZA" => 0.0               # Subsolar: SZA = 0° (maximum solar flux)
        ),
        "night" => Dict(
            "Texo_key" => "min",
            "Texo_override" => 210.0,  # Nightside Texo (K); typical nightside 180-230 K
            "CO2_lowerbdy_vmr" => 0.960,
            "SZA" => 180.0             # Anti-solar: SZA = 180° (no solar flux)
        )
    )
elseif planet == "Mars" && n_horiz == 2
    # Mars day-night test: 2 columns
    OrderedDict(
        "day" => Dict(
            "Texo_key" => "mean-P2",
            "SZA" => 45.0
        ),
        "night" => Dict(
            "Texo_key" => "min-P2",
            "SZA" => 110.0              # Nightside: SZA > 90° means zero solar flux
        )
    )
else
    # Default: replicate uniform conditions across all columns
    OrderedDict(
        ["col_$i" => Dict(
             "Texo_key" => planet == "Venus" ? "mean" : "mean-P2",
             "SZA" => SZA
         ) for i in 1:n_horiz]
    )
end

# Preserve an ordered list for existing solver loops, while keeping named access available.
const horiz_column_names = collect(keys(horiz_column_scenario_by_name))
const horiz_column_scenario = [merge(Dict("name" => col_name), horiz_column_scenario_by_name[col_name]) for col_name in horiz_column_names]

# Validate that the scenario matches n_horiz
@assert length(horiz_column_scenario_by_name) == n_horiz "horiz_column_scenario_by_name must have exactly n_horiz=$n_horiz entries, but got $(length(horiz_column_scenario_by_name))"
@assert length(horiz_column_scenario) == n_horiz "horiz_column_scenario must have exactly n_horiz=$n_horiz entries, but got $(length(horiz_column_scenario))"
