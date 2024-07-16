################################################################################
# INPUT_PARAMETERS.jl
# DESCRIPTION: Variables to be modified by the user for each simulation run. 
# 
# Eryn Cangi
# Created April 2024
# Last edited: May 2024
# Currently tested for Julia: 1.8.5
################################################################################

# Set the planet 
# =======================================================================================================
const planet = "Venus"
    # OPTIONS: "Mars", "Venus"

# Input and output files, directory
# =======================================================================================================
const results_dir = code_dir*"../Results_$(planet)/"
const initial_atm_file = "$(planet)-Inputs/" * "INITIAL_GUESS.h5" # File to use to initialize the atmosphere. 
    # OPTIONS: 
    # INITIAL_GUESS_MARS.h5 --> Basic Mars starting file.
    # INITIAL_GUESS_MARS_bxz4YnHk.h5 --> A Mars atmosphere that includes N2O, NO2, and their ions;
    #                                    not particularly motivated by any present-day data.                                         
    # INITIAL_GUESS_VENUS_vGFd5b0a.h5 --> Best Venus initial atmosphere 
const final_atm_file = "final_atmosphere.h5"
const reaction_network_spreadsheet = code_dir*"$(planet)-Inputs/REACTION_NETWORK_$(uppercase(planet)).xlsx"
    # OPTIONS: "REACTION_NETWORK_MIN_IONOSPHERE.xlsx", code_dir*"REACTION_NETWORK_$(uppercase(planet)).xlsx"

# Descriptive attributes of this model run
# =======================================================================================================
const optional_logging_note = "Run after introducing Mike's sign convention and calculation changes in boundaryconditions and fluxcoefs" # Brief summary of simulation goal
const results_version = "v1"  # Helps keep track of attempts if you need to keep changing things

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
    # Vary water, insolation, and temperature at same time: "inclusive-peri" # "inclusive-mean" #"inclusive-ap" 
    
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

# MARS DUST OPTIONS:
const dust_storm_on = false  # This adds a parcel of water at some altitude, sort of like if there was a dust storm. Haven't published on this.
    # OPTIONS: True, false
const H2O_excess = 250  # excess H2O in ppm
    # OPTIONS: Typical amounts are 100-300 ppm.
const HDO_excess = 0.350 # excess HDO in ppm (divide by 1000 to get ppb)
    # OPTIONS: 0.350 is reasonable

# Settings to control species in the model
# -------------------------------------------------------------------
const ions_included = true
const converge_which = "both"
    # OPTIONS: "ions", "neutrals", "both"
const dont_compute_chemistry = [:Ar]
const dont_compute_transport = []
const dont_compute_either_chem_or_transport = []  # Chemical species which should never update their densities, but may participate in chem+transport.
    # OPTIONS: Any species included in the model. 
const assume_photochem_eq = false # whether to turn on photochemical equilibrium for short-lived species

# Turn plots off and on
# =======================================================================================================
const make_P_and_L_plots = true  # Turn off to save several minutes of runtime if you don't need to check for equilibrium.
    # OPTIONS: True, false

# Algorithm tolerances
# =======================================================================================================
const rel_tol = 1e-6
const abs_tol = 1e-12 

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
