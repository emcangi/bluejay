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

using DataFrames

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
const problem_type = "Gear" #"SS" #"ODE" #   
const timestep_type = "static"# "dynamic"#
const ions_included = false
const converge_which = "neutrals"
const nontherm = false 
const fixed_species = [:Ar, :N2, :CO2pl] # here you may enter any species that you want to be completely fixed (no updates to densities from chemistry or transport)

# Solar conditions
const solarcyc = "mean"
const solarfile = "marssolarphotonflux_solar$(solarcyc).dat" # you may replace 'mean' with 'max' or 'min'
const SZA = 60 # SZA in degrees 

# Descriptions
const tag = "neutrals" # Optional extra bit for the filename to help indicate what it is
const optional_logging_note = "Test if letting water vary leads to a decrease in thermal flux of H"

# Detailed solver characteristics
const n_steps = 500 # Number of logarithmic timesteps to use with gear_timestep_type = 'static'
const dt_incr_factor = 2
const dt_decr_factor = 10
const ediff = false
const mdiff = false 
const electron_val = "none"
const error_checking_scheme = "old"#"new"#

# Temperature and water 
const controltemps = [216., 130., 205.] # mean: 210 # max: 280
const meantemps = [216., 130., 205.] # Used for saturation vapor pressure.
const water_mixing_ratio = 1.38e-4

# Short tags 
const hrshortcode, rshortcode = generate_code(problem_type, ions_included, timestep_type, electron_val, nontherm)

# Folders and files 
const sim_folder_name = "$(hrshortcode)_$(rshortcode)_$(tag)"
const initial_atm_file = "converged_neutral_atmosphere.h5"# "neutrals_SS_nearconverged.h5"#
const final_atm_file = "final_atmosphere.h5"
const reaction_network_spreadsheet = code_dir*"REACTION_NETWORK.xlsx"

# Check that float type is properly set
# if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
#     throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in CUSTOMIZATIONS.jl")
# elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
#     println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
# end

# **************************************************************************** #
#                                                                              #
#                       Simulation time and tolerances                         #
#                                                                              #
# **************************************************************************** #
const dt_min_and_max = Dict("neutrals"=>[-3, 16], "ions"=>[-4, 6], "both"=>[-3, 16])
const rel_tol = 1e-2 # relative tolerance

# **************************************************************************** #
#                                                                              #
#                       Dust storm excess water options                        #
#                                                                              #
# **************************************************************************** #
const dust_storm_on = false
const H2O_excess = 250 # excess H2O in ppm
const HDO_excess = 0.350 # excess HDO in ppm (divide by 1000 to get ppb)
const excess_peak_alt = 42 # altitude at which to add the extra water 

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
const do_chem = true 
const do_trans = true 
const make_new_alt_grid = false
const use_nonzero_initial_profiles = true

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
const T_surf = controltemps[1]
const T_meso = controltemps[2]
const T_exo = controltemps[3]

const Tn_arr = [T(a, controltemps[1], controltemps[2], controltemps[3], "neutral") for a in alt];
const Ti_arr = [T(a, controltemps[1], controltemps[2], controltemps[3], "ion") for a in alt];
const Te_arr = [T(a, controltemps[1], controltemps[2], controltemps[3], "electron") for a in alt];

Temp_keepSVP(z::Float64) = T(z, meantemps..., "neutral") # Needed for boundary conditions.

const Tplasma_arr = Ti_arr .+ Te_arr;
const Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
const Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr)
const fix_SVP = true
const T_top = T(zmax, controltemps[1], controltemps[2], controltemps[3], "neutral")

# **************************************************************************** #
#                                                                              #
#                             Boundary conditions                              #
#                                                                              #
# **************************************************************************** #
const H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # Using this function keeps SVP fixed 
const HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))

const speciesbclist=Dict(:CO2=>Dict("n"=>[2.1e17, NaN], "f"=>[NaN, 0.]),
                        :Ar=>Dict("n"=>[2.0e-2*2.1e17, NaN], "f"=>[NaN, 0.]),
                        :N2=>Dict("n"=>[1.9e-2*2.1e17, NaN], "f"=>[NaN, 0.]),
                        :H2O=>Dict("n"=>[H2Osat[1], NaN], "f"=>[NaN, 0.]), # bc doesnt matter if H2O fixed
                        :HDO=>Dict("n"=>[HDOsat[1], NaN], "f"=>[NaN, 0.]),
                        :O=> Dict("f"=>[0., 1.2e8]),
                        :H2=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(T_top, 2.0; zmax)]),  # velocities are in cm/s
                        :HD=>Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(T_top, 3.0; zmax)]),
                        :H=> Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(T_top, 1.0; zmax)]),
                        :D=> Dict("f"=>[0., NaN], "v"=>[NaN, effusion_velocity(T_top, 2.0; zmax)]), 
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
const neutral_species = [conv_neutrals..., new_neutrals...];

# Ions -------------------------------------------------------------------------
const conv_ions = [:CO2pl];
const new_ions = []; 
const ion_species = [conv_ions..., new_ions...];

# Full species list -------------------------------------------------------------
const all_species = [neutral_species..., ion_species...];

# Sort name lists created here -------------------------------------------------
sort!(all_species)
sort!(neutral_species)
sort!(ion_species)

# Photolysis and Photoionization rate symbol lists ----------------------------

const conv_Jrates, newJrates = format_Jrates(reaction_network_spreadsheet, all_species, "Jratelist"; ions_on=ions_included)
const Jratelist = [conv_Jrates..., newJrates...];

# This dictionary specifies the species absorbing a photon for each J rate using regular expressions.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])

# **************************************************************************** #
#                                                                              #
#                      Miscellaneous logical groupings                         #
#                                                                              #
# **************************************************************************** #
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
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
const no_chem_species = [];
const no_transport_species = [];

# This will append any species that you have picked to be completely fixed. 
for fs in fixed_species
    push!(no_chem_species, fs)
    push!(no_transport_species, fs)
end

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
push!(PARAMETERS_GEN, ("ELECTRON_PROF", electron_val));
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

# This is so ugly because the XLSX package won't write columns of different lengths, so I have to pad all the shorter lists
# with blanks up to the length of the longest list and also transform all the symbols into strings. 
L = max(length(all_species), length(neutral_species), length(ion_species), length(no_chem_species), length(no_transport_species))
PARAMETERS_SPLISTS = DataFrame(AllSpecies=[[string(a) for a in all_species]..., ["" for i in 1:L-length(all_species)]...], 
                               Neutrals=[[string(n) for n in neutral_species]..., ["" for i in 1:L-length(neutral_species)]...], 
                               Ions=[[string(i) for i in ion_species]..., ["" for i in 1:L-length(ion_species)]...],
                               NoChem=[[string(nc) for nc in no_chem_species]..., ["" for i in 1:L-length(no_chem_species)]...],
                               NoTransport=[[string(nt) for nt in no_transport_species]..., ["" for i in 1:L-length(no_transport_species)]...]);

PARAMETERS_SOLVER = DataFrame(Field=[], Value=[]);
PARAMETERS_XSECTS = DataFrame(Species=[], Description=[], Filename=[]);
PARAMETERS_BCS = DataFrame(Species=[], Type=[], Lower=[], Upper=[]);

