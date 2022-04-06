################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Just some standard global constants that need to get used 
# EVERYWHERE. Also the chemical reaction network.
# 
# Eryn Cangi
# Created December 2019
# Last edited: 12 November 2021
# Currently tested for Julia: 1.7.1
################################################################################

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                                                                        !! #
# !!                      !!!!! SUPER IMPORTANT !!!!!                       !! #
# !!     !!! Check the following every time you run the simulation !!!      !! #
# !!                                                                        !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# **************************************************************************** #
#                                                                              #
#                         Main simulation parameters                           #
#                                                                              #
# **************************************************************************** #
# Basic simulation parameters
const simtype = "temp"
const controltemps = [216., 130., 205.]
const ions_included = false
const converge_which = "neutrals"

# Descriptions
const tag = "test_errcheck5" # Optional extra bit for the filename to help indicate what it is
const optional_logging_note = "From neutral baseline: fval relative error is now fval/dndt...checking n values."

# Detailed solver characteristics
const problem_type = "Gear" #"ODE" #  "SS" # 
const gear_timestep_type = "static"
const n_steps = 1000 # Number of logarithmic timesteps to use with gear_timestep_type = 'static'
const dt_incr_factor = 1.1
const dt_decr_factor = 10
const ediff = false
const mdiff = false 
const electron_val = "none"
const error_checking_scheme = "new"#"old"

# Folders and files 
const sim_folder_name = "neutrals_$(problem_type)_$(tag)"
const initial_atm_file = "neutrals_SS_nearconverged.h5"#"converged_neutral_atmosphere.h5"# 
const final_atm_file = "final_atmosphere.h5"
const reaction_network_spreadsheet = code_dir*"REACTION_NETWORK.xlsx"

# Water 
const water_mixing_ratio = 1.38e-4

# Solar conditions
const solarfile = "marssolarphotonflux_solarmean.dat" # you may replace 'mean' with 'max' or 'min'
const SZA = 60 # SZA in degrees 

# Check that float type is properly set
# if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
#     throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in CUSTOMIZATIONS.jl")
# elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
#     println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
# end

# **************************************************************************** #
#                                                                              #
#                       Dust storm excess water options                        #
#                                                                              #
# **************************************************************************** #
const dust_storm_on = false
const H2O_excess = 250 # excess H2O in ppm
const HDO_excess = 0.350 # excess HDO in ppm (divide by 1000 to get ppb)
const excess_peak_alt = 42 # altitude at which to add the extra water 

# **************************************************************************** #
#                                                                              #
#                       Simulation time and tolerances                         #
#                                                                              #
# **************************************************************************** #
const dt_min_and_max = Dict("neutrals"=>[-3, 16], "ions"=>[-4, 6], "both"=>[-4, 14])
const rel_tol = 1e-2 # relative tolerance

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                                                                        !! #
# !!                      !!!!!    END CHECK    !!!!!                       !! #
# !!                                                                        !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# **************************************************************************** #
#                                                                              #
#      Misc. things used when adding new species or changing altitude grid     #
#                                                                              #
# **************************************************************************** #

const make_new_alt_grid = false
const use_nonzero_initial_profiles = true
const do_chem = true 
const do_trans = true 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!! end check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# Sets whether photochemical equilibrium is assumed. Aids in converging ions and neutrals
# together. Generally leave it as is so the code determines it, but you can change it
# if need be
if problem_type == "Gear"
    const assume_photochem_eq = false
else
    const assume_photochem_eq = converge_which == "both" ? true : false
end

# **************************************************************************** #
#                                                                              #
#                    Temperature arrays and water boundaries                   #
#                                                                              #
# **************************************************************************** #

# Temperature 
const T_surf = controltemps[1]
const T_meso = controltemps[2]
const T_exo = controltemps[3]

const Tn_arr = [T_updated(a, controltemps[1], controltemps[2], controltemps[3], "neutral") for a in alt];
const Ti_arr = [T_updated(a, controltemps[1], controltemps[2], controltemps[3], "ion") for a in alt];
const Te_arr = [T_updated(a, controltemps[1], controltemps[2], controltemps[3], "electron") for a in alt];

Temp_keepSVP(z::Float64) = T_updated(z, meanTs, meanTm, meanTe, "neutral") # Needed for boundary conditions.

const Tplasma_arr = Ti_arr .+ Te_arr;
const Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
const Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr)
const fix_SVP = true
const T_top = T_updated(zmax, controltemps[1], controltemps[2], controltemps[3], "neutral")

# **************************************************************************** #
#                                                                              #
#                             Boundary conditions                              #
#                                                                              #
# **************************************************************************** #
const H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # Using this function keeps SVP fixed 
const HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))

const speciesbclist=Dict(
                        :CO2=>["n" 2.1e17; "f" 0.],
                        :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                        :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                        :H2O=>["n" H2Osat[1]; "f" 0.], # bc doesnt matter if H2O fixed
                        :HDO=>["n" HDOsat[1]; "f" 0.],
                        :O=>["f" 0.; "f" 1.2e8],
                        :H2=>["f" 0.; "v" effusion_velocity(T_top, 2.0; zmax)],  # velocities are in cm/s
                        :HD=>["f" 0.; "v" effusion_velocity(T_top, 3.0; zmax)],
                        :H=>["f" 0.; "v" effusion_velocity(T_top, 1.0; zmax)],
                        :D=>["f" 0.; "v" effusion_velocity(T_top, 2.0; zmax)],
                       );

# **************************************************************************** #
#                                                                              #
#                    Species name lists and J rate lists                       #
#                                                                              #
# **************************************************************************** #
# Neutrals --------------------------------------------------------------------
const conv_neutrals = [:Ar, :CO, :CO2, :H, :H2, :H2O, :H2O2, :HO2, :HOCO, :N2, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD];
const new_neutrals = [];

const neutral_species = [];
append!(neutral_species, conv_neutrals)
append!(neutral_species, new_neutrals)

# Ions -------------------------------------------------------------------------
const conv_ions = [:CO2pl];
const new_ions = []; 
const ion_species = [];
append!(ion_species, conv_ions)
append!(ion_species, new_ions)

# Full species list -------------------------------------------------------------
const all_species = [];
append!(all_species, neutral_species)
append!(all_species, ion_species)

# Sort name lists created here -------------------------------------------------
sort!(all_species)
sort!(neutral_species)
sort!(ion_species)

# Photolysis and Photoionization rate symbol lists ----------------------------

# Photolysis and Photoionization rate symbol lists ----------------------------
const conv_Jrates = [# Original neutral photodissociation
                    :JCO2toCOpO,:JCO2toCOpO1D,:JO2toOpO,:JO2toOpO1D,
                    :JO3toO2pO,:JO3toO2pO1D,:JO3toOpOpO,:JH2toHpH,:JOHtoOpH,
                    :JOHtoO1DpH,:JHO2toOHpO,:JH2OtoHpOH,:JH2OtoH2pO1D,:JH2OtoHpHpO,
                    :JH2O2to2OH,:JH2O2toHO2pH,:JH2O2toH2OpO1D,

                    # Original deuterated neutral photodissociation
                    :JHDOtoHpOD, :JHDOtoDpOH, :JHDO2toOHpOD,
                    :JHDOtoHDpO1D, :JHDOtoHpDpO, :JODtoOpD, :JHDtoHpD, :JDO2toODpO,
                    :JHDO2toDO2pH, :JHDO2toHO2pD, :JHDO2toHDOpO1D, :JODtoO1DpD,

                    # New neutral photodissociation (from Roger)
                    # :JCO2toCpOpO, :JCO2toCpO2, :JCOtoCpO, # TODO: Incorporate these to the neutral model.
                  ];
const newJrates = [];
const Jratelist = [];
append!(Jratelist, conv_Jrates)
append!(Jratelist, newJrates)

# This dictionary specifies the species absorbing a photon for each J rate using regular expressions.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])

# To make a dictionary of photolysis products based on the rates, it's easiest to apply regular expressions if 
# we replace all the "p" meaning "plus" with "a", because "p" also appears in "pl" which is part of the spcies name.
# Anyway, her's a dict of photolysis products. 
Jratelist_sub_a_for_p = [Symbol(replace(string(j), r"p(?!l)"=>"a", "JH2O2to2OH"=>"JH2O2toOHaOH")) for j in Jratelist]
const photolysis_products = Dict([x=>[Symbol(m.match) for m in eachmatch(r"(?<=to|a)[A-Z0-9(pl)]+", string(y))] for (x,y) in zip(Jratelist, Jratelist_sub_a_for_p)]);


# **************************************************************************** #
#                                                                              #
#                      Miscellaneous logical groupings                         #
#                                                                              #
# **************************************************************************** #

# Other logical groupings -------------------------------------------------------
const D_bearing_species = [s for s in setdiff(union(neutral_species, ion_species), [:O1D]) if occursin('D', string(s))];
const D_ions = [s for s in ion_species if occursin('D', string(s))];
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];

# Sort name lists created here -------------------------------------------------
sort!(D_bearing_species)
sort!(D_ions)
sort!(N_neutrals)

# **************************************************************************** #
#                                                                              #
#                    Define short- and long-lived species                      #
#                                                                              #
# **************************************************************************** #

# Short lived species, whose chemical lifetime is << diffusion timescale -------
const short_lived_species = [];# technically shortlived but count as longlived: :CH, :HCO, :HO2, :O3, :OH, :O1D, :DO2, :OD...
if assume_photochem_eq
    append!(short_lived_species, [:NO2, :CN, :HNO, :NH, :NH2, :C, :CH])
    append!(short_lived_species, ion_species)
end

# Long lived species ------------------------------------------------------------
const long_lived_species = setdiff(all_species, short_lived_species)

# Sort name lists created here -------------------------------------------------
sort!(short_lived_species)
sort!(long_lived_species)

# **************************************************************************** #
#                                                                              #
#               Species participating in chemistry and transport               #
#                                                                              #
# **************************************************************************** #

# Non-participants -------------------------------------------------------------
const no_chem_species = [:Ar, :CO2pl, :N2];
const no_transport_species = [:Ar, :CO2pl, :N2];

# Participants ------------------------------------------------------------------
const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species ---------------------------------------------------
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)

# Sort name lists created here -------------------------------------------------
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
# its position with the active species vector 
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

