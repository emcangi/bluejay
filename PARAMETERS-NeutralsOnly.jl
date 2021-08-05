################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Just some standard global constants that need to get used 
# EVERYWHERE. Also the chemical reaction network.
# 
# Eryn Cangi
# Created December 2019
# Last edited: 14 May 2021
# Currently tested for Julia: 1.4.1
################################################################################

################################################################################
###################### begin modifiable value zone #############################
################################################################################

const max_alt = 250e5
const upper_lower_bdy = 80e5 # the uppermost layer at which water will be fixed, in cm
const hygropause_alt = 40e5
const MR_mean_water = 1.38e-4

# Timesteps and iterations =====================================================
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-4, 14])
const rel_tol = 1e-4

# General species name lists for converged and newly introduced species =======
const conv_neutrals = [:Ar, :CO, :CO2, :H, :H2, :H2O, :H2O2, :HO2, :HOCO, :N2, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD
                      ];
const conv_ions = [:CO2pl];

const new_neutrals = []; # should remain empty but has to be here for code to work
const new_ions = []; # should remain empty but has to be here for code to work

# Photolysis and Photoionization rate symbol lists
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

const nochemspecies = [:Ar, :N2, :CO2pl];
const notransportspecies = [:Ar, :N2, :CO2pl];

################################################################################
############################ end modification zone #############################
################################################################################


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
const DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988                       # a number to append to some plot filenames. so ugly but it's the easiest way

# Altitude grid specifications =================================================
const alt = convert(Array, (0:2e5:max_alt)) # TODO: Figure out how to get this without being hard-coded
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]; # the altitude grid CELLS but in integers.
const non_bdy_layers = alt[2:end-1]  # all layers, centered on 2 km, 4...248. Excludes the boundary layers which are [-1, 1] and [249, 251].
const num_layers = length(alt) - 2 # there are 124 non-boundary layers.
const plot_grid = non_bdy_layers ./ 1e5;  # for plotting. Points located at atmospheric layer cell centers and in units of km.
const plot_grid_bdys = collect(1:2:((max_alt / 1e5)-1))  # the boundaries; includes the boundary layers at top and bottom.
const upper_lower_bdy_i = Int64(upper_lower_bdy / 2e5) # the uppermost layer at which water will be fixed, in cm

const zmin = alt[1]
const zmax = alt[end];
const dz = alt[2]-alt[1];
const n_alt_index=Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])

# Temperatures and water stuff =================================================
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

highestTe = 350.0  # This is the highest exobase temperature considered.
                          # AFAIK this parameter is only used in making the 3-panel
                          # temperature profile plot.

# Set up the master lists of species ==========================================
const Jratelist = [];
append!(Jratelist, conv_Jrates)
append!(Jratelist, newJrates)

const ionlist = [];
append!(ionlist, conv_ions)

const fullspecieslist = [];
append!(fullspecieslist, conv_neutrals)
append!(fullspecieslist, ionlist)

const neutrallist = setdiff(fullspecieslist, ionlist)

const chemspecies = setdiff(fullspecieslist, nochemspecies);
const transportspecies = setdiff(fullspecieslist, notransportspecies);
const activespecies = union(chemspecies, transportspecies)
const inactivespecies = intersect(nochemspecies, notransportspecies)

# NEW - for handling different water behavior in upper/lower atmo
# Get the position of H2O and HDO symbols within active species. This is used so that we can control its behavior differently
# in different parts of the atmosphere.
const H2Oi = findfirst(x->x==:H2O, activespecies)
const HDOi = findfirst(x->x==:HDO, activespecies)

# this gives the indices of inactivespecies within fullspecieslist. Used for 
# constructing tup and tdown in a more efficient way.
const notransport_inds = [findfirst(x->x==ias, fullspecieslist) for ias in inactivespecies]

#  Useful dictionaries ==========================================================
# List of D bearing species and their H analogues. INCOMPLETE, only has ions right now because we don't need the neutral analogue list.
const D_H_analogues = Dict( :ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl)

const speciesmolmasslist = Dict(:Ar=>40, 
                                :C=>12, :CH=>13, :CN=>26,
                                :CO=>28, :CO2=>44, 
                                :H=>1, :H2=>2, :H2O=>18, :H2O2=>34, 
                                :HCN=>27, :HCO=>29, :HNO=>31, 
                                :HO2=>33, :HOCO=>45, 
                                :N=>14, :N2=>28,
                                :N2O=>44, :NH=>15, :NH2=>16, :NO=>30, :NO2=>46, 
                                :O=>16, :O1D=>16, :O2=>32, :O3=>48, :OH=>17,

                                # Neutrals .- deuterated
                                :D=>2, :DO2=>34, :DOCO=>46, :HD=>3, :HDO=>19, :HDO2=>35, :OD=>18, 

                                # Ions
                                :CO2pl=>44
                                )

# Polarizability from NIST. list: https://cccbdb.nist.gov/pollistx.asp
# CH, CN, H2O2, HCO, HNO, HO2, HOCO, NH, NH2, O(1D), OH, D, DO2, DOCO, HDO, HDO2, and OD 
# are not measured by experiment. Calculations from https://cccbdb.nist.gov/polcalc2x.asp
# exist for CH, CN, H2O2, HCO, HNO, HO2, HOCO, NH, NH2, OH, D, HDO, and OD.
# DO2, DOCO, HDO2 and O1D are not found and are estimated by me to be the same as their
# parent species.
# I used the calcualtions that use "Density functional", "aug-cc-PVDZ", and "mPW1PW91" 
# because that was the method that gave the closest answer for HD to the experimental value. 
# I have no idea what any of it means or whether it's reasonable. I'm not a quantum chemist.
# Values are given in cm^3
const species_polarizability = Dict(# Values available from experiment
                                    :Ar=>1.664e-24, :C=>1.760e-24, :CO=>1.953e-24, :CO2=>2.507e-24, 
                                    :H=>0.667e-24, :H2=>0.787e-24, :H2O=>1.501e-24, :HCN=>2.593e-24, :HD=>0.791e-24, 
                                    :N=>1.1e-24, :N2=>1.710e-24, :N2O=>2.998e-24, :NO=>1.698e-24, :NO2=>2.910e-24, 
                                    :O=>0.802e-24,  :O2=>1.59e-24, :O3=>3.079e-24, 

                                    # Values from calculation
                                    :CH=>2.108e-24, :CN=>3.042e-24, :D=>0.713e-24, :DO2=>1.858e-24, :DOCO=>3.224e-24, 
                                    :H2O2=>2.143e-24, :HCO=>2.505, :HDO=>1.358e-24, :HDO2=>2.143e-24, :HNO=>2.123e-24, 
                                    :HO2=>1.858e-24, :HOCO=>3.224e-24, :NH=>1.418e-24, :NH2=>1.752e-24, 
                                    :O1D=>0.802e-24, :OH=>1.020e-24, :OD=>1.020e-24
                                    )

# This dictionary could probably be replaced with some simple regular expression code, but I haven't done it yet.
const absorber = Dict(:JCO2ion =>:CO2,
                :JCO2toCOpO =>:CO2,
                :JCO2toCOpO1D =>:CO2,
                :JO2toOpO =>:O2,
                :JO2toOpO1D =>:O2,
                :JO3toO2pO =>:O3,
                :JO3toO2pO1D =>:O3,
                :JO3toOpOpO =>:O3,
                :JH2toHpH =>:H2,
                :JHDtoHpD => :HD,
                :JOHtoOpH =>:OH,
                :JOHtoO1DpH =>:OH,
                :JODtoOpD =>:OD,
                :JODtoO1DpD => :OD,
                :JHO2toOHpO =>:HO2,
                :JDO2toODpO => :DO2,
                :JH2OtoHpOH =>:H2O,
                :JH2OtoH2pO1D =>:H2O,
                :JH2OtoHpHpO =>:H2O,
                :JH2O2to2OH =>:H2O2,
                :JH2O2toHO2pH =>:H2O2,
                :JH2O2toH2OpO1D =>:H2O2,
                :JHDO2toHDOpO1D => :HDO2,
                :JHDOtoHpOD=>:HDO,
                :JHDO2toOHpOD=>:HDO2,
                :JHDO2toDO2pH => :HDO2,
                :JHDO2toHO2pD => :HDO2,
                :JHDOtoDpOH=>:HDO,
                :JHDOtoHpDpO=>:HDO,
                :JHDOtoHDpO1D=>:HDO,
                # NEW: reactions from Roger's model. 
                :JH2OtoH2Opl=>:H2O,
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
                           :Ar=>"#808080", :C=>"#d9c382", :CH=>"#cea3ce", :CN=>"#6d5000", :HCN=>"#479f5e", :HCO=>"#94c6bf", 
                           :N=>"#6e748a", :N2=>"#cccccc", :N2O=>"#be8b65", :NCO=>"#633339", :NH=>"#FD8C9B", :NH2=>"#6ceb83", :NO=>"#a27fff", :NO2=>"#fe03cb", 
                           :HNO=>"#76bcfd",

                           # ions
                           :Arpl=>"#808080", :ArHpl=>"#660000", :ArDpl=>"#660000",
                           :Cpl=>"#d9c382", :CHpl=>"#cea3ce", :CNpl=>"#6d5000", 
                           :COpl=>"#ff6600", :CO2pl=>"#000000",  
                           :Dpl=>"#ff0000", :DCOpl=>"#3366ff", :HDpl=>"#e526d7", :HD2pl=>"#b9675f", :H2Dpl=>"#b9675f",
                           :Hpl=>"#ff0000", :H2pl=>"#e526d7", :H3pl=>"#b9675f",
                           :H2Opl=>"#0083dc",  :H3Opl=>"#280041",
                           :HCNpl=>"#479f5e", :HCNHpl=>"#50455b", 
                           :HCOpl=>"#3366ff", :HOCpl=>"#5e90ff", :HCO2pl=>"#222222", 
                           :HNOpl=>"#eb0077", :HN2Opl=>"#a37bb3", :HOCOpl=>"#e8ba8c", :HO2pl=>"#046868", 
                           :Npl=>"#6e748a", :N2pl=>"#cccccc", :N2Dpl=>"#9a4700", :N2Hpl=>"#9a4700", :N2Opl=>"#be8b65", 
                           :NHpl=>"#cacdda", :NH2pl=>"#6ceb83", :NH3pl=>"#c8c400", 
                           :NOpl=>"#a27fff", :NO2pl=>"#fe03cb", 
                           :Opl=>"#1a6115", :O2pl=>"#15da09", :OHpl=>"#7700d5", 
                           );

# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict(:D => "--", :HD => "--", :OD => "--", :HDO => "--", :HDO2 => "--", :DO2 => "--", :DOCO => "--", 
                           :ArDpl=>"--", :Dpl=>"--", :DCOpl=>"--", :HDpl=>"--", :HD2pl=>"--", :H2Dpl=>"-.", :N2Dpl=>"--");
                
const medgray = "#444444"

# Crosssection filenames ======================================================
# There's gotta be a better way to do this. Probably a dictionary of Jrates to strings.

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

# Chemistry ====================================================================
# function to replace three body rates with the recommended expression
# threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1))
# threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1))

threebody(k0, kinf) = :($k0 .* M ./ (1 .+ $k0 .* M ./ $kinf).*0.6 .^ ((1 .+ (log10.($k0 .* M ./ $kinf)) .^2).^-1.0))
threebodyca(k0, kinf) = :($k0 ./ (1 .+ $k0 ./ ($kinf ./ M)).*0.6 .^ ((1 .+ (log10.($k0 ./ ($kinf .* M))) .^2).^-1.0))

################################################################################
############################### REACTION NETWORK ###############################
################################################################################

# reactions and multipliers on base rates for deuterium reactions from Yung
# 1988; base rates from this work or Chaffin+ 2017. Note: below, H-ana means 
# the same reaction but with only H-bearing species.
const reactionnet = [   #Photodissociation
             [[:CO2], [:CO, :O], :JCO2toCOpO],
             [[:CO2], [:CO, :O1D], :JCO2toCOpO1D],
             [[:O2], [:O, :O], :JO2toOpO],
             [[:O2], [:O, :O1D], :JO2toOpO1D],
             [[:O3], [:O2, :O], :JO3toO2pO],
             [[:O3], [:O2, :O1D], :JO3toO2pO1D],
             [[:O3], [:O, :O, :O], :JO3toOpOpO],
             [[:H2], [:H, :H], :JH2toHpH],
             [[:HD], [:H, :D], :JHDtoHpD],
             [[:OH], [:O, :H], :JOHtoOpH],
             [[:OH], [:O1D, :H], :JOHtoO1DpH],
             [[:OD], [:O, :D], :JODtoOpD],
             [[:OD], [:O1D, :D], :JODtoO1DpD],
             [[:HO2], [:OH, :O], :JHO2toOHpO], # other branches should be here, but have not been measured
             [[:DO2], [:OD, :O], :JDO2toODpO],
             [[:H2O], [:H, :OH], :JH2OtoHpOH],
             [[:H2O], [:H2, :O1D], :JH2OtoH2pO1D],
             [[:H2O], [:H, :H, :O], :JH2OtoHpHpO],
             [[:HDO], [:H, :OD], :JHDOtoHpOD], 
             [[:HDO], [:D, :OH], :JHDOtoDpOH], 
             [[:HDO], [:HD, :O1D], :JHDOtoHDpO1D], # inspiration from Yung89
             [[:HDO], [:H, :D, :O], :JHDOtoHpDpO], # inspiration from Yung89
             [[:H2O2], [:OH, :OH], :JH2O2to2OH],
             [[:H2O2], [:HO2, :H], :JH2O2toHO2pH],
             [[:H2O2], [:H2O, :O1D], :JH2O2toH2OpO1D],
             [[:HDO2], [:OH, :OD], :JHDO2toOHpOD], # Yung89
             [[:HDO2], [:DO2, :H], :JHDO2toDO2pH],
             [[:HDO2], [:HO2, :D], :JHDO2toHO2pD],
             [[:HDO2], [:HDO, :O1D], :JHDO2toHDOpO1D],
             # NEW: neutral dissociation from Roger Yelle
             # [[:CO2], [:C, :O, :O], :JCO2toCpOpO], # New, non-N neutrals 
             # [[:CO2], [:C, :O2], :JCO2toCpO2],#  New, non-N neutrals 
             # [[:CO], [:C, :O], :JCOtoCpO], # New, non-N neutrals 
             
             # NEW: photoionization from Roger's model
             # [[:CO2], [:CO2pl], :JCO2toCO2pl],  # Nair minimal ionosphere

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13 .* (Tn ./ 300.) .^ 0.6),:(2.1e9 .* (Tn ./ 300.) .^ 6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13 .* (Tn ./ 300.) .^ 0.6),:(2.1e9 .* (Tn ./ 300.) .^ 6.1))], # Yung88: same as H-ana.
             ## to HOCO/DOCO
             [[:CO, :OH], [:HOCO], threebody(:(5.9e-33 .* (Tn ./ 300.) .^ -1.4),:(1.1e-12 .* (Tn ./ 300.) .^ 1.3))], # Sander2011
             [[:CO, :OD], [:DOCO], threebody(:(5.9e-33 .* (Tn ./ 300.) .^ -1.4),:(1.1e-12 .* (Tn ./ 300.) .^ 1.3))],

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 ./ Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5) .* 8.7e-10)],

             # recombination of H. Use EITHER the first line OR the 2nd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32 .* (298 ./ Tn) .^ 2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32 .* (298 ./ Tn) .^ 2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32 .* (298 ./ Tn) .^ 2.27)], # Yung88: rate same as H-ana.

             # Various H reactions 
             [[:H, :HD], [:H2, :D], :(6.31e-11 .* exp.(-4038 ./ Tn))], # TODO:CHECK rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11 .* exp.(-3821 ./ Tn))], # TODO:CHECK NIST (1986, 200-300K): 8.19e-13 .* exp.(-2700/Tn)
             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn) .^ 2)], # NIST database
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn) .^ 2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn) .^ 2)], # not in Yung88. assumed rate
             [[:H, :OD], [:OH, :D], :(3.3e-9 .* (Tn .^ -0.63) ./ (0.72 .* exp.(717 ./ Tn)))], # rate: Yung88. NIST (Howard82): 5.25E-11 .* (Tn/298) .^ -0.63  .- turn off for Case 2
             [[:D, :OH], [:OD, :H], :(3.3e-9 .* Tn .^ -0.63)], # Yung88  .- turn off for Case 2

             [[:H, :O2], [:HO2], threebody(:(2.0 .* 4.4e-32 .* (Tn ./ 300.) .^ -1.3), # Sander2011, 300K+. Yung89: 5.5e-32(Tn ./ 300) .^ -1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11 .* (Tn ./ 300.) .^ 0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0 .* 4.4e-32 .* (Tn ./ 300.) .^ -1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11 .* (Tn ./ 300.) .^ 0.2))],
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Burkholder 2020
             [[:D, :HO2], [:OH, :OD], :(0.71 .* 7.2e-11)], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:H, :DO2], [:OH, :OD], :(7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :HO2], [:H2, :O2], :(0.5 .* 6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:D, :HO2], [:HD, :O2], :(0.71 .* 0.5 .* 6.9e-12)], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:H, :DO2], [:HD, :O2], :(0.5 .* 6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # Burkholder 2020; they do not have O1D though
             [[:D, :HO2], [:HDO, :O1D], :(0.71 .* 1.6e-12)], # Yung88: rate 0.71 .* H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HDO, :O1D], :(1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:H, :DO2], [:HO2, :D], :(1e-10 ./ (0.54 .* exp.(890 ./ Tn)))], # Yung88 (assumed) .- turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 .- turn off for Case 2
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12 .* exp.(-1890 ./ Tn))], # verified NIST 4/3/18. Only valid for Tn>300K. No experiment for lower.
             # [[:H, :HDO2], [:DO2, :H2], :(0)], # Cazaux2010: branching ratio = 0
             # [[:H, :HDO2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:DO2, :H2], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11 .* exp.(-1800 ./ Tn))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO, :OH], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<Tn<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O, :OD], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # see previous line
             [[:D, :H2O2], [:HDO, :OH], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # see previous line
             [[:D, :H2O2], [:H2O, :OD], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # see previous line
             [[:D, :HDO2], [:OD, :HDO], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # added 4/3 with assumed rate from other rxns
             [[:D, :HDO2], [:OH, :D2O], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # sourced from Cazaux et al
             [[:H, :O3], [:OH, :O2], :(1.4e-10 .* exp.(-470 ./ Tn))], # Burkholder 2020
             [[:D, :O3], [:OD, :O2], :(0.71 .* 1.4e-10 .* exp.(-470 ./ Tn))], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89, NIST 4/3/18.

             # H2 and HD loss
             [[:H2, :O], [:OH, :H], :(6.34e-12 .* exp.(-4000 ./ Tn))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12 .* exp.(-4390 ./ Tn))], # NIST. TODO:CHECK
             [[:HD, :O], [:OD, :H], :(1.68e-12 .* exp.(-4400 ./ Tn))], # NIST TODO:CHECK
             [[:H2, :OH], [:H2O, :H], :(2.8e-12 .* exp.(-1800 ./ Tn))], # Burkholder 2020
             [[:H2, :OD], [:HDO, :H], :(2.8e-12 .* exp.(-1800 ./ Tn))], # Yung88: rate same as H-ana (assumed)
             [[:H2, :OD], [:H2O, :D], :(0)], # Yung88 (assumed)
             [[:HD, :OH], [:HDO, :H], :((3 ./ 20.) .* 2.8e-12 .* exp.(-1800 ./ Tn))], # Yung88: rate (3/20) .* H-ana. Sander2011: 5e-12 .* exp.(-2130 ./ Tn)
             [[:HD, :OH], [:H2O, :D], :((3 ./ 20.) .* 2.8e-12 .* exp.(-1800 ./ Tn))], # see prev line
             ### [[:HD, :OD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:HD, :OD], [:D2O, :H], :(???)],  # can't find a rate...?

             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13 .* exp.(460 ./ Tn))], # Burkholder2020. Yung89: 2.3e-13 .* exp.(600/Tn). KIDA 230-420K: 2.2e-13 .* exp.(600/Tn)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13 .* exp.(460 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2 .* 2.1e-33 .* exp.(920 ./ Tn))], # Burkholder2020.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2 .* 2.1e-33 .* exp.(920 ./ Tn))], # added 3/13 with assumed same rate as H analogue
             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11 .* exp.(180 ./ Tn))], # Burkholder 2020 Yung89: 2.2e-11 .* exp.(120/Tn) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11 .* exp.(180 ./ Tn))], # Yung88: rate same as H-ana.
             [[:O, :HO2], [:OH, :O2], :(3.0e-11 .* exp.(200 ./ Tn))], # Burkholder 2020
             [[:O, :DO2], [:OD, :O2], :(3.0e-11 .* exp.(200 ./ Tn))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12 .* exp.(-2000 ./ Tn))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5 .* 1.4e-12 .* exp.(-2000 ./ Tn))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5 .* 1.4e-12 .* exp.(-2000 ./ Tn))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:OH, :OH], [:H2O, :O], :(1.8e-12)], # Burkholder 2020
             [[:OD, :OH], [:HDO, :O], :(1.8e-12)], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3 .* 6.9e-31 .* (Tn ./ 300.) .^ -1.0), :(2.6e-11))], # Burkholder2020. Why 1.3? Reviewer suggestion to adjust for mars?
             [[:OD, :OH], [:HDO2], threebody(:(1.3 .* 6.9e-31 .* (Tn ./ 300.) .^ -1.0), :(2.6e-11))], # Yung88: rate same as H-ana
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12 .* exp.(-940 ./ Tn))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12 .* exp.(-940 ./ Tn))], # Yung88: rate same as H-ana
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # verified NIST 4/3/18. Yung89: 4.6e-11 .* exp.(230/Tn) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # Yung88: same as H-ana.
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12 .* exp.(-160 ./ Tn))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12 .* exp.(-200/Tn). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5 .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: rate 0.5 .* H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5 .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: rate 0.5 .* H-ana.
             [[:O3, :HO2], [:OH, :O2, :O2], :(1.0e-14 .* exp.(-490 ./ Tn))], # Sander2011. Yung89: 1.1e-14 .* exp.(-500/Tn). KIDA 250-340K: 2.03e-16 .* (Tn ./ 300) .^ 4.57 .* exp.(693/Tn). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:O3, :DO2], [:OD, :O2, :O2], :(1.0e-14 .* exp.(-490 ./ Tn))], # Yung88: same as H-ana (assumed)
             
             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8 .* 3.0e-33 .* (300 ./ Tn) .^ 3.25)], # Deighan 2012 
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35 .* exp.(724 ./ Tn))], 
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5 .* 6.0e-34 .* (300 ./ Tn) .^ 2.4)], # Burkholder2020 
             [[:O, :O3], [:O2, :O2], :(8.0e-12 .* exp.(-2060 ./ Tn))],  # Burkholder 2020
             # [[:O, :CO, :M], [:CO2, :M], :(2.2e-33 .* exp.(-1780 ./ Tn))],  # DUPLICATE - commented out in favor of Roger's reaction

             # O1D attack            
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Burkholder 2020
             [[:O1D, :HD], [:H, :OD], :(0.41 .* 1.2e-10)], # Yung88: rate 0.41 .* H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41 .* 1.2e-10)], # Yung88: rate 0.41 .* H-ana (assumed). NIST 1e-10 @298K
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10 .* exp.(60 ./ Tn))], # Burkholder 2020
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10 .* exp.(60 ./ Tn))], # Yung88: rate same as H-ana.
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11 .* exp.(115 ./ Tn))], # Burkholder 2020
             [[:O1D, :O2], [:O, :O2], :(3.3e-11 .* exp.(55 ./ Tn))], # Burkholder 2020 (upd. 31 Dec 2020)
             [[:O1D, :O3], [:O2, :O2], :(2.4e-10)], # Burkholder 2020
             [[:O1D, :O3], [:O, :O, :O2], :(2.4e-10)], # Burkholder 2020

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 / Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)]

              # -----------------

             
             # NEW .- Neutral reactions from Roger Yelle
             # [[:CO2, :H], [:CO, :OH], :(3.38e-10 .* exp.(-13163.0 ./ Tn))],  # KIDA: :(2.51e-10 .* exp.(-13300 ./ Tn))
             # [[:CO2, :O], [:O2, :CO], :(2.46e-11 .* exp.(-26567.0 ./ Tn))], # KIDA agrees
             # [[:H2O, :H], [:OH, :H2], :((1.69e-14 .* Tn .^ 1.2) .* exp.(-9610.0 ./ Tn))], # KIDA :(6.82e-12 .* (Tn ./ 300) .^ 1.6 .* exp.(-9720 ./ Tn))
             # [[:H2O, :O], [:OH, :OH], :((8.20e-14 .* Tn .^ 0.95) .* exp.(-8571.0 ./ Tn))],  # KIDA :(1.85e-11 .* (Tn ./ 300) .^ 0.95 .* exp.(-8570 ./ Tn))
             # [[:HO2, :CO], [:CO2, :OH], :(5.6e-10 .* exp.(-12160.0 ./ Tn))],
             # [[:HO2, :H2], [:H2O2, :H], :(5.0e-11 .* exp.(-13110.0 ./ Tn))],
             # [[:HO2, :H2O], [:H2O2, :OH], :(4.65e-11 .* exp.(-16500.0 ./ Tn))],
             # [[:HOCO, :OH], [:CO2, :H2O], :(1.03e-11)],
             # [[:O, :H], [:OH], :(8.65e-18 .* (Tn .^ -0.38))],
             # [[:O1D, :CO], [:CO, :O], :(4.7e-11 .* exp.(63.0 ./ Tn))],
             # [[:O1D, :CO], [:CO2], :(8.0e-11)],
             # [[:O1D, :H2O2], [:H2O2, :O], :(5.2e-10)],
             # [[:O2, :CO], [:CO2, :O], :(5.99e-12 .* exp.(-24075.0 ./ Tn))],
             # [[:O2, :H], [:OH, :O], :(2.61e-10 .* exp.(-8156.0 ./ Tn))],
             # [[:O2, :H2], [:HO2, :H], :(2.4e-10 .* exp.(-28500.0 ./ Tn))],
             # [[:O2, :H2], [:OH, :OH], :(3.16e-10 .* exp.(-21890.0 ./ Tn))],
             # [[:OH, :H], [:H2, :O], :(8.1e-21 .* (Tn .^ 2.8) .* exp.(-1950.0 ./ Tn))],
             
             # Type 4
             # simpler type, when Troe parameter is 0. Updated 5 Feb 2021 to match Roger's code
             # [[:CO, :H], [:HCO], :(0.0 .+ (2.0e-35 .* (Tn .^ 0.2) .* 1.0 .* (Tn .^ 0.2) .* M) ./ (2.0e-35 .* (Tn .^ 0.2) .* M .+ 1.0 .* (Tn .^ 0.2)))], 

             # More complicated - need the minimum of the two expressions. These updated 5 Feb 2021 to match Roger's code.
             # [[:CO, :O], [:CO2], :(min.($:(1.0 .* exp.(-1509.0 ./ Tn)), $:(0.0 .+ (10 .^ ((log10.(0.4)) ./ (1 .+ ((log10.((1.7e-33 .* exp.(-1509.0 ./ Tn) .* M) ./ (1.0 .* exp.(-1509.0 ./ Tn))) .- 0.4 .- 0.67 .* log10.(0.4)) ./ (0.75 .- 1.27 .* log10.(0.4) .- 0.14 .* (log10.((1.7e-33 .* exp.(-1509.0 ./ Tn) .* M) ./ (1.0 .* exp.(-1509.0 ./ Tn))) .- 0.4 .- 0.67 .* log10.(0.4)))) .^ 2)) .* 1.7e-33 .* exp.(-1509.0 ./ Tn) .* 1.0 .* exp.(-1509.0 ./ Tn) .* M) ./ (1.7e-33 .* exp.(-1509.0 ./ Tn) .* M .+ 1.0 .* exp.(-1509.0 ./ Tn)))))]
             ];
