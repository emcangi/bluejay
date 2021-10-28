################################################################################
# PARAMETERS-conv3.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Global constants, simulation parameters, reaction networks. 
# USE: Use only when converging an atmosphere after introducing ions, as
# step 3. Adds N-bearing neutrals and all ions to an atmosphere that contains 
# standard neutrals, C, CH, HCO, and the Nair minimal ions.
# 
# Eryn Cangi
# Created December 2019
# Last edited: August 2021
# Currently tested for Julia: 1.6.1
################################################################################

# **************************************************************************** #
# **************************************************************************** #
# **************************************************************************** #
#                      begin modifiable value zone                             #
# **************************************************************************** #
# **************************************************************************** #
# **************************************************************************** #

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                      !!!!! SUPER IMPORTANT !!!!!                       !! #
# !!     !!! Check the following every time you run the simulation !!!      !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
const sim_folder_name = "ions3-QNDF"
const initial_atm_file = "ions2-QNDF-old.h5"
const final_atm_file = "ions3-QNDF-no-quad-NH"
const converge_which = "both"
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-4, 14])
const rel_tol = 1e-4

# This stuff is mutable, but less likely to change.
const make_new_alt_grid = false
const use_nonzero_initial_profiles = true
const do_chem = true 
const do_trans = true 
const solarfile = "marssolarphotonflux_solarmean.dat" # you may replace 'mean' with 'max' or 'min'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!! end check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# this line OVERRIDDEN because in this scenario we need to get the ions and N-neutrals in.
const assume_photochem_eq = converge_which == "both" ? true : false

# Some water stuff to control the general shape of the profile
const upper_lower_bdy = 80e5 # the uppermost layer at which water will be fixed, in cm
const hygropause_alt = 40e5
const MR_mean_water = 1.38e-4

# Species name lists and Jrate symbol lists  ===================================

# Neutrals --------------------------------------------------------------------
const conv_neutrals = [:Ar, :C, :CH, :CO, :CO2, :H, :H2, :H2O, :H2O2, :HCO, :HO2, :HOCO, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :N2, :OD];
const new_neutrals = [:DCO, :N2O, :NO2, :CN, :HCN, :HNO, :N, :NH, :NH2, :NO];

const neutral_species = [];
append!(neutral_species, conv_neutrals)
append!(neutral_species, new_neutrals)

# Ions -------------------------------------------------------------------------
const conv_ions = [:CO2pl, :HCO2pl, :Opl, :O2pl];  # Nair minimal ionosphere  
const new_ions = [:Arpl, :ArHpl, :Cpl, :CHpl, :CNpl, :COpl, 
                  :Hpl, :H2pl, :H2Opl, :H3pl, :H3Opl,
                  :HCNpl, :HCNHpl, :HCOpl, 
                  :HNOpl, :HN2Opl, :HOCpl, :HO2pl, 
                  :Npl,  :NHpl, :NH2pl, :NH3pl, :N2pl, :N2Hpl, :N2Opl, :NOpl, :NO2pl,
                  :OHpl,
                  # Deuterated ions
                   :ArDpl, :Dpl, :DCOpl, :DOCpl, :DCO2pl, :HDpl, :HD2pl, :H2Dpl, :H2DOpl, :HDOpl, :N2Dpl, :ODpl
                 ];

const ion_species = [];
append!(ion_species, conv_ions)
append!(ion_species, new_ions)

# Full species list -------------------------------------------------------------
const all_species = [];
append!(all_species, neutral_species)
append!(all_species, ion_species)

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

                    # New photoionization/ion-involved photodissociation (Roger)
                    :JCO2toCO2pl, :JCO2toOplpCO, :JOtoOpl, :JO2toO2pl, # Nair minimal ionosphere

                    # New neutral photodissociation (from Roger)
                    :JCO2toCpOpO, :JCO2toCpO2, :JCOtoCpO,
                  ];
const newJrates = [# New neutral photodissociation (from Roger)
                    :JN2OtoN2pO1D, :JNO2toNOpO, :JNOtoNpO,

                    # New photoionization/ion-involved photodissociation (Roger)
                    :JCO2toCO2plpl, :JCO2toCplplpO2, :JCO2toCOplpOpl,:JCO2toOplpCplpO, :JCO2toCplpO2, :JCO2toCOplpO, 
                    :JCOtoCpOpl, :JCOtoCOpl,  :JCOtoOpCpl, 
                    :JHtoHpl, 
                    :JH2toH2pl, :JH2toHplpH, :JHDtoHDpl, 
                    :JH2OtoH2Opl, 
                    :JH2OtoOplpH2, :JH2OtoHplpOH, :JH2OtoOHplpH, :JHDOtoHDOpl,
                    :JH2O2toH2O2pl, 
                    :JN2toN2pl, :JN2toNplpN, 
                    :JN2OtoN2Opl, :JNO2toNO2pl, :JNOtoNOpl, 
                    :JO3toO3pl
                ];

const Jratelist = [];
append!(Jratelist, conv_Jrates)
append!(Jratelist, newJrates)

# Other logical groupings -------------------------------------------------------
const D_bearing_species = [s for s in union(neutral_species, ion_species) if occursin('D', string(s))];
const N_species = [s for s in neutral_species if occursin('N', string(s))];

# Short lived species, whose chemical lifetime is << diffusion timescale
const short_lived_species = [];# technically shortlived but count as longlived: :CH, :HCO, :HO2, :O3, :OH, :O1D, :DO2, :OD...
if assume_photochem_eq
    append!(short_lived_species, [:NO2, :CN, :HNO, :NH, :NH2])
    append!(short_lived_species, ion_species)
end

# Long lived species
const long_lived_species = setdiff(all_species, short_lived_species)

# Species participating in chem and transport -----------------------------------

# Non-participants
const no_chem_species = [];
const no_transport_species = [];

if converge_which == "neutrals"
    append!(no_chem_species, union(conv_ions, N_species)) # This is because the N chemistry is intimiately tied up with the ions.
    append!(no_transport_species, union(conv_ions, N_species, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, setdiff(conv_neutrals, N_species))
    append!(no_transport_species, union(short_lived_species, setdiff(conv_neutrals, N_species)))
elseif converge_which == "both"
    # Next two lines are SPECIAL for this file because we are trying to work in ions, N-neutrals. 
    # append!(no_chem_species, setdiff(conv_neutrals, N_species))  # i.e. non-nitrogen bearing neutrals are inactive
    # append!(no_transport_species, setdiff(conv_neutrals, N_species))
    append!(no_transport_species, short_lived_species)
end

# Participants
const chem_species = setdiff(all_species, no_chem_species);
const transport_species = setdiff(all_species, no_transport_species);

# Active and inactive species 
const active_species = union(chem_species, transport_species)
const inactive_species = intersect(no_chem_species, no_transport_species)
const active_longlived = intersect(active_species, long_lived_species)
const active_shortlived = intersect(active_species, short_lived_species)

# To allow water to be active in the upper atmosphere but not the lower atmosphere, we need 
# its position with the active species vector 
const H2Oi = findfirst(x->x==:H2O, active_longlived)
const HDOi = findfirst(x->x==:HDO, active_longlived)

# **************************************************************************** #
# **************************************************************************** #
# **************************************************************************** #
#                        end modifiable value zone                             #
# **************************************************************************** #
# **************************************************************************** #
# **************************************************************************** #

const extra_plots_dir = "/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Auxiliary plots/"
const research_dir = "/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/"
const results_dir = "/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Results/"
const xsecfolder = research_dir * "uvxsect/";

# fundamental constants ========================================================
const kB_MKS = 1.38e-23;        # J/K - needed for saturation vapor pressure empirical equation.
const kB = 1.38e-16;            # erg/K
const bigG = 6.67e-8;           # dyne-cm^2/g^2
const mH = 1.67e-24;            # g 
const marsM = 0.1075*5.972e27;  # g 
const radiusM = 3396e5;         # cm
const q = 4.8032e-10            # statcoulomb (cm^1.5 g^0.5 s^-1)
const DH = 5.5 * 1.6e-4         # SMOW value from Yung 1988

# Altitude grid specifications =================================================
max_alt = 250e5
const alt = convert(Array, (0:2e5:max_alt)) # TODO: Figure out how to get this without being hard-coded
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]; # the altitude grid CELLS but in integers.
const non_bdy_layers = alt[2:end-1]  # all layers, centered on 2 km, 4...248. Excludes the boundary layers which are [-1, 1] and [249, 251].
const num_layers = length(alt) - 2 # there are 124 non-boundary layers.
const plot_grid = non_bdy_layers ./ 1e5;  # for plotting. Points located at atmospheric layer cell centers and in units of km.
const plot_grid_bdys = collect(1:2:((max_alt / 1e5)-1))  # the boundaries; includes the boundary layers at top and bottom.
const upper_lower_bdy_i = Int64(upper_lower_bdy / 2e5) # index of the uppermost layer at which water will be fixed

const zmin = alt[1]
const zmax = alt[end];
const dz = alt[2]-alt[1];
const n_alt_index=Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])

# Mean temperatures and simulation temperature parameters ======================
const meanTs = 216.0
const meanTt = 130.0
const meanTe = 205.0
const meantemps = [meanTs, meanTt, meanTe]

const meanTsint = 216
const meanTtint = 130
const meanTeint = 205

# These are the low and high values for the "standard atmosphere and reasonable 
# climate variations" cases. NOT the full range of temperatures used to make the
# detailed cases stuff, because those include temps up to 350 K.
const lowTs = 160.0
const hiTs = 270.0
const lowTt = 100.0
const hiTt = 160.0
const lowTe = 150.0
const hiTe = 250.0

const highestTe = 350.0  # This is the highest exobase temperature considered.
                          # AFAIK this parameter is only used in making the 3-panel
                          # temperature profile plot.


#  Useful dictionaries ==========================================================
# List of D bearing species and their H analogues. INCOMPLETE, only has ions right now because we don't need the neutral analogue list.
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl,)  

const molmass = Dict(:Ar=>40, 
                     :C=>12, :CH=>13, :CN=>26,
                     :CO=>28, :CO2=>44, 
                     :H=>1, :H2=>2, :H2O=>18, :H2O2=>34, 
                     :HCN=>27, :HCO=>29, :HNO=>31, 
                     :HO2=>33, :HOCO=>45, 
                     :N=>14, :N2=>28,
                     :N2O=>44, :NH=>15, :NH2=>16, :NO=>30, :NO2=>46, 
                     :O=>16, :O1D=>16, :O2=>32, :O3=>48, :OH=>17,

                     # Neutrals .- deuterated
                     :D=>2, :DCO=>30, :DO2=>34, :DOCO=>46, :HD=>3, :HDO=>19, :HDO2=>35, :OD=>18, 

                     # Ions
                     :Arpl=>40, :ArHpl=>41, :ArDpl=>42,
                     :Cpl=>12, :CHpl=>13, :CNpl=>26, 
                     :COpl=>28, :CO2pl=>44,
                     :Hpl=>1, :Dpl=>2, 
                     :H2pl=>2, :HDpl=>3, 
                     :H3pl=>3, :H2Dpl=>4, :HD2pl=>5, 
                     :H2Opl=>18, :HDOpl=>19, 
                     :H3Opl=>19,  :H2DOpl=>20,
                     :HO2pl=>33, 
                     :HCOpl=>29, :HOCpl=>29, :DCOpl=>30, :DOCpl=>30, 
                     :HCO2pl=>45, :DCO2pl=>46, 
                     :HCNpl=>27, :HCNHpl=>28, :HNOpl=>31, :HN2Opl=>45,  
                     :Npl=>14, :N2pl=>28, 
                     :N2Hpl=>29, :N2Dpl=>30, 
                     :N2Opl=>44, 
                     :NHpl=>15, :NH2pl=>16, :NH3pl=>17,
                     :NOpl=>30, :NO2pl=>46,
                     :Opl=>16, :O2pl=>32, :OHpl=>17, :ODpl=>18 
                     )


# Polarizability from NIST. Experimental values from: https://cccbdb.nist.gov/pollistx.asp
# Calculations for species not available ine xperiment from: https://cccbdb.nist.gov/polcalc2x.asp
# Deuterated species not listed in either are estimated by me to be the same as their H-bearing analogue.
# I used the calcualtions that use "Density functional", "aug-cc-PVDZ", and "mPW1PW91" 
# because that was the method that gave the closest answer for HD to the experimental value. 
# I have no idea what any of it means or whether it's reasonable. I'm not a quantum chemist.
# Values are given in cm^3
const polarizability = Dict(# Values available from experiment
                                    :Ar=>1.664e-24, :C=>1.760e-24,  :CO=>1.953e-24,  :CO2=>2.507e-24, 
                                    :H=>0.667e-24,  :H2=>0.787e-24, :H2O=>1.501e-24, :HCN=>2.593e-24, :HD=>0.791e-24, 
                                    :N=>1.1e-24,    :N2=>1.710e-24, :N2O=>2.998e-24, :NO=>1.698e-24, :NO2=>2.910e-24, 
                                    :O=>0.802e-24,  :O2=>1.59e-24,  :O3=>3.079e-24, 

                                    # Values from calculation
                                    :CH=>2.108e-24,   :CN=>3.042e-24,   :D=>0.713e-24, 
                                    :H2O2=>2.143e-24, :HCO=>2.505e-24,  :HDO=>1.358e-24, :HNO=>2.123e-24, 
                                    :HO2=>1.858e-24,  :HOCO=>3.224e-24, :NH=>1.418e-24,  :NH2=>1.752e-24, 
                                    :OH=>1.020e-24,   :OD=>1.020e-24,

                                    # Assumed same as hydrogen analogue
                                    :DCO=>2.505e-24, :DO2=>1.858e-24, :DOCO=>3.224e-24, :HDO2=>2.143e-24, :O1D=>0.802e-24, 
                                    )

const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])

# Common plot specifications =======================================================
const speciescolor = Dict(:H => "#ff0000", :D => "#ff0000", # red
                           :H2 => "#e526d7", :HD =>  "#e526d7", # dark pink/magenta

                           # hydroxides
                           :OH => "#7700d5", :OD => "#7700d5", # purple

                           # water group (roughly, I ain't a chemist)
                           :H2O => "#0083dc", :HDO => "#0083dc", # cornflower blue
                           :H2O2 => "#0000ff", :HDO2 => "#0000ff", # true blue
                           :HO2 => "#046868", :DO2 => "#046868",  # dark teal

                           # O group
                           :O1D => "#808000", # olive
                           :O => "#1a6115",   # forest green
                           :O2 => "#15da09",  # kelly/grass green
                           :O3 => "#269e56",  # light green

                           # CO group
                           :CO2 => "#000000",   # black
                           :CO => "#ff6600",    # orange
                           :HOCO => "#e8ba8c", :DOCO => "#e8ba8c",  #tannish

                           # Other neutrals
                           :Ar=>"#808080", :C=>"#d9c382", :CH=>"#cea3ce", :CN=>"#6d5000", :HCN=>"#479f5e", 
                           :HCO=>"#94c6bf", :DCO=>"#94c6bf",
                           :N=>"#6e748a", :N2=>"#cccccc", :N2O=>"#be8b65", :NCO=>"#633339", :NH=>"#FD8C9B", :NH2=>"#6ceb83", :NO=>"#a27fff", :NO2=>"#fe03cb", 
                           :HNO=>"#76bcfd",

                           # ions
                           :Arpl=>"#808080", :ArHpl=>"#660000", :ArDpl=>"#660000",
                           :Cpl=>"#d9c382", :CHpl=>"#cea3ce", :CNpl=>"#6d5000", 
                           :COpl=>"#ff6600", :CO2pl=>"#000000",  
                           :Dpl=>"#ff0000", 
                           :HDpl=>"#e526d7", :HD2pl=>"#b9675f", :H2Dpl=>"#b9675f",
                           :Hpl=>"#ff0000", :H2pl=>"#e526d7", :H3pl=>"#b9675f",
                           :H2Opl=>"#0083dc", :HDOpl=>"#0083dc",
                           :H3Opl=>"#280041", :H2DOpl=>"#280041", #:HD2Opl=>"#280041",
                           :HCNpl=>"#479f5e", :HCNHpl=>"#50455b", 
                           :HCOpl=>"#3366ff", :DCOpl=>"#3366ff", 
                           :HOCpl=>"#5e90ff", :DOCpl=>"#5e90ff", 
                           :HCO2pl=>"#222222", :DCO2pl=>"#222222",
                           :HNOpl=>"#eb0077", :HN2Opl=>"#a37bb3", :HOCOpl=>"#e8ba8c", :HO2pl=>"#046868", 
                           :Npl=>"#6e748a", :N2pl=>"#cccccc", 
                           :N2Hpl=>"#9a4700", :N2Dpl=>"#9a4700", 
                           :N2Opl=>"#be8b65", 
                           :NHpl=>"#cacdda", :NH2pl=>"#6ceb83", :NH3pl=>"#c8c400", 
                           :NOpl=>"#a27fff", :NO2pl=>"#fe03cb", 
                           :Opl=>"#1a6115", :O2pl=>"#15da09", 
                           :OHpl=>"#7700d5", :ODpl=>"#7700d5"
                           );

# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict([s=>"--" for s in D_bearing_species]);
                
const medgray = "#444444"

# Crosssection filenames ======================================================
const co2file = "CO2.dat"
const co2exfile = "binnedCO2e.csv" # added to shield short λ of sunlight in upper atmo
const h2ofile = "h2oavgtbl.dat"
const hdofile = "HDO.dat"#"HDO_250K.dat"#
const h2o2file = "H2O2.dat"
const hdo2file = "H2O2.dat" #TODO: do HDO2 xsects exist?
const o3file = "O3.dat"
const o3chapfile = "O3Chap.dat"
const o2file = "O2.dat"
const o2_130_190 = "130-190.cf4"
const o2_190_280 = "190-280.cf4"
const o2_280_500 = "280-500.cf4"
const h2file = "binnedH2.csv"
const hdfile = "binnedH2.csv" # TODO: change this to HD file if xsects ever exist
const ohfile = "binnedOH.csv"
const oho1dfile = "binnedOHo1D.csv"
const odfile = "OD.csv"

# Transport network ============================================================
const upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
                    Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
                    for s in transport_species];

const downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
                      Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
                      for s in transport_species];

const local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in transport_species]
                                [Symbol("t"*string(s)*"_down") for s in transport_species]
                                [Symbol("t"*string(s)*"_above_down") for s in transport_species]
                                [Symbol("t"*string(s)*"_below_up") for s in transport_species]]...;];

const transportnet = [[upeqns...;]; [downeqns...;]];

# define names for all the species active in the coupled rates:
const active_longlived_above = [Symbol(string(s)*"_above") for s in active_longlived];
const active_longlived_below = [Symbol(string(s)*"_below") for s in active_longlived];

