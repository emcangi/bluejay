################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Global constants, simulation parameters, reaction networks. 
# 
# Eryn Cangi
# Created December 2019
# Last edited: January 2022
# Currently tested for Julia: 1.6.1
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
const tag = "internally_screaming"#addN(2D)" # Optional extra bit for the filename to help indicate what it is
const problem_type = "SS" # "Gear" #"ODE" #  
const converge_which = "both"
const optional_logging_note = "This is a test run because the stackoverflow error is back and I wanna scream"#This run is to incorporate the species N2D "

# Check that float type is properly set
if problem_type == "Gear" && (ftype_ncur == Float64 || ftype_chem == Float64)
    throw("If problem_type = 'Gear' in PARAMETERS, both ftype_ncur and ftype_chem must = Double64 in CUSTOMIZATIONS.jl")
elseif problem_type != "Gear" && (ftype_ncur == Double64 || ftype_chem == Double64)
    println("problem_type != Gear but using Double64 in CUSTOMIZATIONS.jl")
end

# Folders and files 
const sim_folder_name = "$(simtype)_$(Int64(controltemps[1]))_$(Int64(controltemps[2]))_$(Int64(controltemps[3]))_$(problem_type)_$(tag)"
const initial_atm_file = "converged_gear_20211231.h5"#"converged_full_atmosphere.h5"
# const initial_atm_file = results_dir*"temp_216_130_205_SS_Newest/final_atmosphere.h5"
const final_atm_file = "final_atmosphere.h5"
const reaction_network_spreadsheet = code_dir*"REACTION_NETWORK.xlsx"

# Water 
const water_mixing_ratio = 1.38e-4

# Solar conditions
const solarfile = "marssolarphotonflux_solarmean.dat" # you may replace 'mean' with 'max' or 'min'
const SZA = 60 # SZA in degrees 

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
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-4, 14])
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
const do_chem = true 
const do_trans = true 
const make_new_alt_grid = false
const use_nonzero_initial_profiles = false

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
                        :H2=>["f" 0.; "v" effusion_velocity(T_top, 2.0, zmax)],  # velocities are in cm/s
                        :HD=>["f" 0.; "v" effusion_velocity(T_top, 3.0, zmax)],
                        :H=>["f" 0.; "v" effusion_velocity(T_top, 1.0, zmax)],
                        :D=>["f" 0.; "v" effusion_velocity(T_top, 2.0, zmax)],
                       );

# **************************************************************************** #
#                                                                              #
#                    Species name lists and J rate lists                       #
#                                                                              #
# **************************************************************************** #
# Neutrals --------------------------------------------------------------------
const conv_neutrals = [:Ar, :C, :CH, :CN, :CO, :CO2, :H, :H2, :H2O, :H2O2, 
                       :HCN, :HCO, :HNO, :HO2, :HOCO, 
                       :N, :N2, :NH, :NH2, :N2O, :NO, :NO2, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DCO, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD];
const new_neutrals = []#:Nup2D];

const neutral_species = [];
append!(neutral_species, conv_neutrals)
append!(neutral_species, new_neutrals)

# Ions -------------------------------------------------------------------------
const conv_ions = [:CO2pl, :HCO2pl, :DCO2pl, :Opl, :O2pl, # Nair minimal ionosphere 
                   :Arpl, :ArHpl, :ArDpl, 
                   :Cpl, :CHpl, :CNpl, :COpl, 
                   :Hpl, :Dpl, :H2pl, :HDpl, :H3pl, :H2Dpl, :HD2pl, 
                   :H2Opl,  :HDOpl, :H3Opl, :H2DOpl,
                   :HO2pl, :HCOpl, :DCOpl, :HOCpl, :DOCpl,
                   :HCNpl, :HCNHpl, :HNOpl, :HN2Opl,  
                   :Npl,  :NHpl, :NH2pl, :NH3pl, :N2pl, :N2Hpl, :N2Dpl, :N2Opl, :NOpl, :NO2pl,
                   :OHpl, :ODpl               
                  ];
const new_ions = []; 
const ion_species = [];
append!(ion_species, conv_ions)
append!(ion_species, new_ions)

# Full species list -------------------------------------------------------------
const all_species = [];
append!(all_species, neutral_species)
append!(all_species, ion_species)

# Photolysis and Photoionization rate symbol lists ----------------------------

const conv_Jrates = [# Original neutral photodissociation
                    :JCO2toCOaO,:JCO2toCOaO1D,:JO2toOaO,:JO2toOaO1D,
                    :JO3toO2aO,:JO3toO2aO1D,:JO3toOaOaO,:JH2toHaH,:JOHtoOaH,
                    :JOHtoO1DaH,:JHO2toOHaO,:JH2OtoHaOH,:JH2OtoH2aO1D,:JH2OtoHaHaO,
                    :JH2O2toOHaOH,:JH2O2toHO2aH,:JH2O2toH2OaO1D,

                    # Original deuterated neutral photodissociation
                    :JHDOtoHaOD, :JHDOtoDaOH, :JHDO2toOHaOD,
                    :JHDOtoHDaO1D, :JHDOtoHaDaO, :JODtoOaD, :JHDtoHaD, :JDO2toODaO,
                    :JHDO2toDO2aH, :JHDO2toHO2aD, :JHDO2toHDOaO1D, :JODtoO1DaD,

                    # New neutral photodissociation (from Roger)
                    :JCO2toCaOaO, :JCO2toCaO2, :JCOtoCaO,
                    :JN2OtoN2aO1D, :JNO2toNOaO, :JNOtoNaO,

                    # New photoionization/ion-involved photodissociation (Roger)
                    :JCO2toCO2pl, :JCO2toOplaCO, :JOtoOpl, :JO2toO2pl, # Nair minimal ionosphere
                    :JCO2toCO2plpl, :JCO2toCplplaO2, :JCO2toCOplaOpl,:JCO2toOplaCplaO, :JCO2toCplaO2, :JCO2toCOplaO, 
                    :JCOtoCaOpl, :JCOtoCOpl,  :JCOtoOaCpl, 
                    :JHtoHpl, 
                    :JH2toH2pl, :JH2toHplaH, :JHDtoHDpl, 
                    :JH2OtoH2Opl, 
                    :JH2OtoOplaH2, :JH2OtoHplaOH, :JH2OtoOHplaH, :JHDOtoHDOpl,
                    :JH2O2toH2O2pl, 
                    :JN2toN2pl, :JN2toNplaN, 
                    :JN2OtoN2Opl, :JNO2toNO2pl, :JNOtoNOpl, 
                    :JO3toO3pl,
                  ];
const newJrates = [];
const Jratelist = [];
append!(Jratelist, conv_Jrates)
append!(Jratelist, newJrates)

# These dictionaries specify the species absorbing a photon for each J rate, and the products of the reaction.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])
const photolysis_products = Dict([x=>[Symbol(m.match) for m in eachmatch(r"(?<=to|a)[A-Z0-9(pl)]+", string(x))] for x in Jratelist]);

# **************************************************************************** #
#                                                                              #
#                      Miscellaneous logical groupings                         #
#                                                                              #
# **************************************************************************** #
# List of D bearing species and their H analogue ions.
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
const D_bearing_species = [s for s in setdiff(union(neutral_species, ion_species), [:O1D]) if occursin('D', string(s))];
const D_ions = [s for s in ion_species if occursin('D', string(s))];
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];

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

# **************************************************************************** #
#                                                                              #
#               Species participating in chemistry and transport               #
#                                                                              #
# **************************************************************************** #

# Non-participants -------------------------------------------------------------
const no_chem_species = [];
const no_transport_species = [];

if converge_which == "neutrals"
    append!(no_chem_species, union(conv_ions, N_neutrals)) # This is because the N chemistry is intimiately tied up with the ions.
    append!(no_transport_species, union(conv_ions, N_neutrals, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, setdiff(conv_neutrals, N_neutrals))
    append!(no_transport_species, setdiff(conv_neutrals, N_neutrals))
elseif converge_which == "both"
    append!(no_transport_species, short_lived_species)
end

# Participants ------------------------------------------------------------------
const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species ---------------------------------------------------
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)


# **************************************************************************** #
#                                                                              #
#              Misc. things that depend on things defined above                #
#                                                                              #
# **************************************************************************** #

# Annoyingly, these have to be here because they depend on other things defined above.
# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict(vcat([s=>"--" for s in setdiff(D_bearing_species, [:HD2pl])], [:HD2pl=>":", :Nup2D=>"-."]) )

# Species-specific scale heights - has to be done here instead of in the param file
const Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)], molmass) for sp in all_species])

# To allow water to be active in the upper atmosphere but not the lower atmosphere, we need 
# its position with the active species vector 
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

