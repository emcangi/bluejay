################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Just some standard global constants that need to get used 
# EVERYWHERE. Also the chemical reaction network.
# 
# Eryn Cangi
# Created December 2019
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
################################################################################

research_dir = "/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/"
results_dir = "/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Results/"

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Altitude grid discretization =================================================
const alt = convert(Array, (0:2e5:250e5)) # TODO: Figure out how to get this without being hard-coded
const intaltgrid = round.(Int64, alt/1e5)[2:end-1];
const zmin = alt[1]
const zmax = alt[end];
const dz = alt[2]-alt[1];
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])

hygropause_alt = 40e5

# Temperatures and water stuff =================================================
global meanTs = 216.0
global meanTt = 130.0
global meanTe = 205.0
global meantemps = [meanTs, meanTt, meanTe]

global meanTsint = 216
global meanTtint = 130
global meanTeint = 205

# These are the low and high values for the "standard atmosphere and reasonable 
# climate variations" cases. NOT the full range of temperatures used to make the
# detailed cases stuff, because those include temps up to 350 K.
global lowTs = 160.0
global hiTs = 270.0
global lowTt = 100.0
global hiTt = 160.0
global lowTe = 150.0
global hiTe = 250.0

global highestTe = 350.0  # This is the highest exobase temperature considered.
                          # AFAIK this parameter is only used in making the 3-panel
                          # temperature profile plot.

MR_mean_water = 1.38e-4

# Lists ========================================================================
const fullspecieslist = [:Ar, :C, :CH, :CN, :CO, :CO2, 
                         :H, :H2, :H2O, :H2O2, 
                         :HCN, :HCO, :HNO, :HO2, :HOCO, 
                         :N, :N2, :N2O, :NH, :NH2, :NO, :NO2, 
                         :O, :O1D, :O2, :O3, :OH,
                         :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD,   # D-bearing

                         # Ions
                         :Arpl, :ArHpl, :Cpl, :CHpl, :CNpl, 
                         :COpl, 
                         :CO2pl, # Nair minimal ionosphere
                         :Hpl, :H2pl, :H2Opl, :H3pl, :H3Opl,
                         :HCNpl, :HCNHpl, :HCOpl, 
                         :HCO2pl, # Nair minimal ionosphere
                         :HNOpl, :HN2Opl, :HOCpl, :HO2pl, 
                         :Npl,  :NHpl, :NH2pl, :NH3pl, :N2pl, :N2Hpl, :N2Opl, :NOpl, :NO2pl,
                         :Opl, :O2pl, # Nair minimal ionosphere
                         :OHpl, 
                         # Deuterated ions
                         #:ArDpl, :CDpl, :Dpl, 
                         #:DCOpl, :DOCOpl, 
                         #:HDpl, :HD2pl, :H2Dpl, :H2DOpl, 
                         #:N_2Dpl, :NDpl, :ODpl
                         ]; 
# NOTE: These species are produced, but not tracked in the model - not in the fullspecieslist.
# :C2Npl, :C2N2pl, :C2Opl, :C2pl, :CH2pl, :CN2, :H2CNpl, :H2COpl, :H2NOpl, :HCOOH2pl, :HeHpl,
# :HNCOpl, :HNCpl, :N2D, :NCO, :NH3, :NH4pl, :OCNpl, :NHDpl, :NHD2pl, :ND3pl, :HC2Npl,


ionlist = [sp for sp in fullspecieslist if occursin("pl", String(sp))]  # A list of only ions, for calculating :E density.
neutrallist = setdiff(fullspecieslist, ionlist)

specieslist=fullspecieslist;  

# array of species for which photolysis is important. All rates should
# start with J and contain a species in specieslist above, which is used to 
# compute photolysis. 
const Jratelist=[:JCO2ion,:JCO2toCOpO,:JCO2toCOpO1D,:JO2toOpO,:JO2toOpO1D,
                 :JO3toO2pO,:JO3toO2pO1D,:JO3toOpOpO,:JH2toHpH,:JOHtoOpH,
                 :JOHtoO1DpH,:JHO2toOHpO,:JH2OtoHpOH,:JH2OtoH2pO1D,:JH2OtoHpHpO,
                 :JH2O2to2OH,:JH2O2toHO2pH,:JH2O2toH2OpO1D,
                 # deuterated species J rates:
                 :JHDOtoHpOD, :JHDOtoDpOH, :JHDO2toOHpOD,
                 :JHDOtoHDpO1D, :JHDOtoHpDpO, :JODtoOpD, :JHDtoHpD, :JDO2toODpO,
                 :JHDO2toDO2pH, :JHDO2toHO2pD, :JHDO2toHDOpO1D, :JODtoO1DpD,

                 # NEW: photodissociation from Roger
                 :JCO2toCpOpO, :JCO2toCpO2, :JCOtoCpO, :JN2OtoN2pO1D, :JNO2toNOpO,:JNOtoNpO, 
                  
                 # NEW: Photoionization from Roger
                 :JCO2toCO2pl, :JCO2toOplpCO, # Nair minimal ionosphere
                 :JCO2toCO2plpl, :JCO2toCplplpO2, :JCO2toCOplpOpl,:JCO2toOplpCplpO, :JCO2toCplpO2, :JCO2toCOplpO, 
                 :JCOtoCpOpl, :JCOtoCOpl,  :JCOtoOpCpl, 
                 
                 :JHtoHpl, 
                 :JH2toH2pl, :JH2toHplpH, 
                 :JH2OtoH2Opl, :JH2OtoOplpH2, :JH2OtoHplpOH, :JH2OtoOHplpH,
                 :JH2O2toH2O2pl, 
                 
                 :JN2toN2pl, :JN2toNplpN, 
                 :JN2OtoN2Opl, :JNO2toNO2pl, :JNOtoNOpl, 

                 :JOtoOpl, :JO2toO2pl, # Nair minimal ionosphere
                 :JO3toO3pl
                 ];

const nochemspecies = [:H2O, :HDO, :Ar, :N2];
const chemspecies = setdiff(specieslist, nochemspecies);
const notransportspecies = [:H2O, :HDO, :Ar, :N2];
const transportspecies = setdiff(specieslist, notransportspecies);
const activespecies = union(chemspecies, transportspecies)
const inactivespecies = intersect(nochemspecies, notransportspecies)
const speciesmolmasslist = Dict(:Ar=>40, 
                                :C=>12, :CH=>13, :CN=>26,
                                :CO=>28, :CO2=>44, 
                                :H=>1, :H2=>2, :H2O=>18, :H2O2=>34, 
                                :HCN=>27, :HCO=>29, :HNO=>31, 
                                :HO2=>33, :HOCO=>45, 
                                :N=>14, :N2=>28,
                                :N2O=>44, :NH=>15, :NH2=>16, :NO=>30, :NO2=>46, 
                                :O=>16, :O1D=>16, :O2=>32, :O3=>48, :OH=>17,

                                # Neutrals - deuterated
                                :D=>2, :DO2=>34, :DOCO=>46, :HD=>3, :HDO=>19, :HDO2=>35, :OD=>18, 

                                # Ions
                                :Arpl=>40, :ArHpl=>41,    # Ar
                                :Cpl=>12,   # C ions
                                :CHpl=>13,          # -anes
                                :CNpl=>26, # C and N
                                :COpl=>28, :CO2pl=>44, # CO and stuff 
                                :Hpl=>1, :HCNpl=>27, :HCNHpl=>28, :HCOpl=>29, 
                                :HCO2pl=>45, :HNOpl=>31, :HN2Opl=>45, :HOCpl=>29, :HO2pl=>33, 
                                :H2pl=>2, :H2Opl=>18,
                                :H3pl=>3, :H3Opl=>19,   # H ions
                                :Npl=>14, :N2pl=>28, 
                                :N2Opl=>44, :N2Hpl=>29, :NHpl=>15, :NH2pl=>16, :NH3pl=>17,  # ammonia and such
                                :NOpl=>30, :NO2pl=>46,
                                :Opl=>16, :O2pl=>32, :OHpl=>17, # O ions 

                                # Deuterated ions
                                #:ArDpl,
                                #:CDpl, 
                                #:Dpl, 
                                #:DCOpl, :DOCOpl, 
                                #:HDpl, :HD2pl,
                                #:H2Dpl, :H2DOpl, 
                                #:N_2Dpl, :NDpl,
                                # :ODpl

                                # electrons
                                :E=>5.44191266e-4)


# Plotty plot plot stuff =======================================================

rcParams = PyDict(matplotlib."rcParams")

global speciescolor = Dict( # H group
                :H => "#ff0000", :D => "#ff0000", # red
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
                :Arpl=>"#808080", :ArHpl=>"#d70000",
                :Cpl=>"#d9c382", :CHpl=>"#cea3ce", :CNpl=>"#6d5000", 
                :COpl=>"#ff6600", :CO2pl=>"#000000",  
                :Hpl=>"#ff0000", :H2pl=>"#e526d7", :H3pl=>"#dcb3af",
                :H2Opl=>"#0083dc",  :H3Opl=>"#280041",
                :HCNpl=>"#479f5e", :HCNHpl=>"#50455b", 
                :HCOpl=>"#94c6bf", :HOCpl=>"#5e90ff", :HCO2pl=>"#222222", 
                :HNOpl=>"#eb0077", :HN2Opl=>"#a37bb3", :HOCOpl=>"#e8ba8c", :HO2pl=>"#046868", 
                :Npl=>"#6e748a", :N2pl=>"#cccccc", :N2Hpl=>"#9a4700", :N2Opl=>"#be8b65", 
                :NHpl=>"#cacdda", :NH2pl=>"#6ceb83", :NH3pl=>"#c8c400", 
                :NOpl=>"#a27fff", :NO2pl=>"#fe03cb", 
                :Opl=>"#1a6115", :O2pl=>"#15da09", :OHpl=>"#7700d5", 
                );

# D group will have dashed lines; neutrals, solid (default)
global speciesstyle = Dict(:D => "--", :HD => "--", :OD => "--", :HDO => "--", :HDO2 => "--", :DO2 => "--", :DOCO => "--",
                    :CN=>":", :HCN=>":", :HNO=>":", :N=>":", :N2=>":", :N2O=>":", :NH=>":", :NH2=>":", :NO=>":", :NO2=>":",
                    :CNpl=>":", :HCNpl=>":", :HCNHpl=>":", :HNOpl=>":", :HN2Opl=>":", :Npl=>":", :N2pl=>":", :N2Opl=>":", :NHpl=>":", :NH2pl=>":", :NH3pl=>":", 
                    :NOpl=>":", :NO2pl=>":", :N2Hpl=>":",
                    :Ar=>"-.", :C=>"-.", :CO=>"-.", :CO2=>"-.", :O=>"-.", :O1D=>"-.", :O2=>"-.", :O3=>"-.",
                    :Arpl=>"-.", :ArHpl=>"-.", :Cpl=>"-.", :COpl=>"-.", :CO2pl=>"-.", :Opl=>"-.", :O2pl=>"-.", );
                
medgray = "#444444"

# Chemistry ====================================================================
# function to replace three body rates with the recommended expression
threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1))
threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1))

################################################################################
############################### REACTION NETWORK ###############################
################################################################################

# reactions and multipliers on base rates for deuterium reactions from Yung
# 1988; base rates from this work or Chaffin+ 2017. Note: below, H-ana means 
# the same reaction but with only H-bearing species.
reactionnet = [   #Photodissociation
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
             [[:CO], [:C, :O], :JCOtoCpO],
             [[:N2O], [:N2, :O1D], :JN2OtoN2pO1D],
             [[:NO2], [:NO, :O], :JNO2toNOpO],
             [[:NO], [:N, :O], :JNOtoNpO],
             [[:CO2], [:C, :O, :O], :JCO2toCpOpO],
             [[:CO2], [:C, :O2], :JCO2toCpO2],
             # NEW: photoionization from Roger's model
             [[:CO], [:COpl], :JCOtoCOpl],
             [[:CO], [:C, :Opl], :JCOtoCpOpl],
             [[:CO], [:O, :Cpl], :JCOtoOpCpl],
             [[:CO2], [:CO2pl], :JCO2toCO2pl],  # Nair minimal ionosphere
             [[:CO2], [:CO2plpl], :JCO2toCO2plpl],  
             [[:CO2], [:Cplpl, :O2], :JCO2toCplplpO2],
             [[:CO2], [:Cpl, :O2], :JCO2toCplpO2],
             [[:CO2], [:COpl, :Opl], :JCO2toCOplpOpl],
             [[:CO2], [:COpl, :O], :JCO2toCOplpO],
             [[:CO2], [:Opl, :CO], :JCO2toOplpCO], # Nair minimal ionosphere
             [[:CO2], [:Opl, :Cpl, :O], :JCO2toOplpCplpO],
             [[:H2O], [:H2Opl], :JH2OtoH2Opl],
             [[:H2O], [:Opl, :H2], :JH2OtoOplpH2],
             [[:H2O], [:Hpl, :OH], :JH2OtoHplpOH],
             [[:H2O], [:OHpl, :H], :JH2OtoOHplpH],
             [[:N2], [:N2pl], :JN2toN2pl],
             [[:N2], [:Npl, :N], :JN2toNplpN],
             [[:NO2], [:NO2pl], :JNO2toNO2pl],
             [[:NO], [:NOpl], :JNOtoNOpl],
             [[:N2O], [:N2Opl], :JN2OtoN2Opl],
             [[:H], [:Hpl], :JHtoHpl],
             [[:H2], [:H2pl], :JH2toH2pl],
             [[:H2], [:Hpl, :H], :JH2toHplpH],
             [[:H2O2], [:H2O2pl], :JH2O2toH2O2pl], 
             [[:O], [:Opl], :JOtoOpl],   # Nair minimal ionosphere
             [[:O2], [:O2pl], :JO2toO2pl],   # Nair minimal ionosphere
             [[:O3], [:O3pl], :JO3toO3pl],

             # recombination of 
             [[:O, :O, :M], [:O2, :M], :(1.8*3.0e-33*(300 ./ Tn)^3.25)],
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35*exp(724 ./ Tn))],
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5*6.0e-34*(300 ./ Tn)^2.4)],
             [[:O, :O3], [:O2, :O2], :(8.0e-12*exp(-2060 ./ Tn))],  # Burkholder 2020
             [[:O, :CO, :M], [:CO2, :M], :(2.2e-33*exp(-1780 ./ Tn))],

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.3e-11*exp(55 ./ Tn))], # Burkholder 2020 (upd. 31 Dec 2020)
             [[:O1D, :O3], [:O2, :O2], :(2.4e-10)], # Burkholder 2020
             [[:O1D, :O3], [:O, :O, :O2], :(2.4e-10)], # Burkholder 2020
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11*exp(115 ./ Tn))], # Burkholder 2020
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Burkholder 2020
             [[:O1D, :HD], [:H, :OD], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10*exp(60 ./ Tn))], # Burkholder 2020
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp(60 ./ Tn))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12*exp(-4000 ./ Tn))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12*exp(-4390 ./ Tn))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12*exp(-4400 ./ Tn))], # NIST
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11*exp(-4038 ./ Tn))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11*exp(-3821 ./ Tn))], # NIST (1986, 200-300K): 8.19e-13*exp(-2700/Tn)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12*exp(-1800 ./ Tn))], # Burkholder 2020
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ Tn))], # Yung88: rate (3/20)*H-ana. Sander2011: 5e-12*exp(-2130 ./ Tn)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ Tn))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12*exp(-1800 ./ Tn))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:OD, :HD], [:D2O, :H], :(???)],  # can't find a rate...?

             # recombination of H. Use EITHER the first line OR the 2nd and 3rd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32*(298 ./ Tn)^2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32*(298 ./ Tn)^2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32*(298 ./ Tn)^2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9*6.8e-31*(300 ./ Tn)^2)], # Can't find in databases. Mike's rate.
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ Tn)^2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ Tn)^2)], # not in Yung88. assumed rate

             ## H + HO2
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Burkholder 2020
             [[:H, :HO2], [:H2, :O2], :(0.5*6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # Burkholder 2020; they do not have O1D though
             [[:H, :DO2], [:OH, :OD], :(7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HD, :O2], :(0.5*6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HDO, :O1D], :(1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:D, :HO2], [:OH, :OD], :(0.71*7.2e-11)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:D, :HO2], [:HD, :O2], :(0.71*0.5*6.9e-12)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:D, :HO2], [:HDO, :O1D], :(0.71*1.6e-12)], # Yung88: rate 0.71*H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HO2, :D], :(1e-10 ./ (0.54*exp(890 ./ Tn)))], # Yung88 (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12*exp(-1890 ./ Tn))], # verified NIST 4/3/18. Only valid for Tn>300K. No experiment for lower.
             # [[:H, :HDO2], [:DO2, :H2], :(0)], # Cazaux2010: branching ratio = 0
             # [[:H, :HDO2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:DO2, :H2], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11*exp(-1800 ./ Tn))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO, :OH], :(0.5*1.16e-11*exp(-2110 ./ Tn))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<Tn<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O, :OD], :(0.5*1.16e-11*exp(-2110 ./ Tn))], # see previous line
             [[:D, :H2O2], [:HDO, :OH], :(0.5*1.16e-11*exp(-2110 ./ Tn))], # see previous line
             [[:D, :H2O2], [:H2O, :OD], :(0.5*1.16e-11*exp(-2110 ./ Tn))], # see previous line
             [[:D, :HDO2], [:OD, :HDO], :(0.5*1.16e-11*exp(-2110 ./ Tn))], # added 4/3 with assumed rate from other rxns
             [[:D, :HDO2], [:OH, :D2O], :(0.5*1.16e-11*exp(-2110 ./ Tn))], # sourced from Cazaux et al

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0*4.4e-32*(Tn/300.)^-1.3), # Sander2011, 300K+. Yung89: 5.5e-32(Tn/300)^-1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11*(Tn/300.)^0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0*4.4e-32*(Tn/300.)^-1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11*(Tn/300.)^0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10*exp(-470 ./ Tn))], # Burkholder 2020
             [[:D, :O3], [:OD, :O2], :(0.71*1.4e-10*exp(-470 ./ Tn))], # Yung88: rate 0.71*H-ana (assumed). verified Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11*exp(180 ./ Tn))], # Burkholder 2020Yung89: 2.2e-11*exp(120/Tn) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11*exp(180 ./ Tn))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11*exp(200 ./ Tn))], # Burkholder 2020
             [[:O, :DO2], [:OD, :O2], :(3.0e-11*exp(200 ./ Tn))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12*exp(-2000 ./ Tn))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5*1.4e-12*exp(-2000 ./ Tn))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5*1.4e-12*exp(-2000 ./ Tn))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(1.8e-12)], # Burkholder 2020
             [[:OD, :OH], [:HDO, :O], :(1.8e-12)], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3*6.9e-31*(Tn/300.)^-1.0),:(2.6e-11))], # Sander2011. Why 1.3?
             [[:OD, :OH], [:HDO2], threebody(:(1.3*6.9e-31*(Tn/300.)^-1.0),:(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12*exp(-940 ./ Tn))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12*exp(-940 ./ Tn))], # Yung88: rate same as H-ana
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11*exp(250 ./ Tn))], # verified NIST 4/3/18. Yung89: 4.6e-11*exp(230/Tn) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ Tn))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ Tn))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12*exp(-160 ./ Tn))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12*exp(-200/Tn). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12*exp(-160 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5*2.9e-12*exp(-160 ./ Tn))], # Yung88: rate 0.5*H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5*2.9e-12*exp(-160 ./ Tn))], # Yung88: rate 0.5*H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14*exp(-490 ./ Tn))], # Sander2011. Yung89: 1.1e-14*exp(-500/Tn). KIDA 250-340K: 2.03e-16*(Tn/300)^4.57*exp(693/Tn). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(1.0e-14*exp(-490 ./ Tn))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13*exp(460 ./ Tn))], # Sander2011. Yung89: 2.3e-13*exp(600/Tn). KIDA 230-420K: 2.2e-13*exp(600/Tn)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13*exp(460 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2*2.1e-33*exp(920 ./ Tn))], # Sander2011.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2*2.1e-33*exp(920 ./ Tn))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9*(Tn^-0.63)/(0.72*exp(717 ./ Tn)))], # rate: Yung88. NIST (Howard82): 5.25E-11*(Tn/298)^-0.63  - turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9*Tn^-0.63)], # Yung88  - turn off for Case 2

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13*(Tn/300.)^0.6),:(2.1e9*(Tn/300.)^6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13*(Tn/300.)^0.6),:(2.1e9*(Tn/300.)^6.1))], # Yung88: same as H-ana.
             ### possible deuterated analogues below
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33*(Tn/300.)^-1.4),:(1.1e-12*(Tn/300.)^1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33*(Tn/300.)^-1.4),:(1.1e-12*(Tn/300.)^1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 / Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)],

             # NEW - Neutral reactions from Roger Yelle
             # Type 1
             [[:O1D], [:O], :(5.10e-3)],

             # Type 2
             [[:C, :C], [:C2], :(2.16e-11)],
             [[:C, :H], [:CH], :(1.00e-17)],
             [[:C, :N], [:CN], :((6.93e-20*Tn^0.37)*exp(-51.0 ./ Tn))],
             [[:CH, :C], [:C2, :H], :(6.59e-11)],
             [[:CH, :H], [:H2, :C], :(1.31e-10*exp(-80.0 ./ Tn))],
             [[:CH, :H2], [:CH2, :H], :(2.90e-10*exp(-1670.0 ./ Tn))],
             [[:CH, :H2], [:CH3], :((2.92e-16*Tn^-0.71)*exp(-11.6 ./ Tn))],
             [[:CH, :N], [:CN, :H], :(2.77e-10*Tn^-0.09)],
             [[:CH, :O], [:CO, :H], :(6.60e-11)],
             [[:CH, :O], [:HCOP, :e], :(4.20e-13*exp(-850.0 ./ Tn))],
             [[:CH, :O], [:OH, :C], :(2.52e-11*exp(-2381.0 ./ Tn))],
             [[:CN, :H2], [:HCN, :H], :((1.80e-19*Tn^2.6)*exp(-960.0 ./ Tn))],
             [[:CN, :N], [:N2, :C], :(9.80e-10*Tn^-0.4)],
             [[:CN, :NH], [:HCN, :N], :((1.70e-13*Tn^0.5)*exp(-1000.0 ./ Tn))],
             [[:CN, :O], [:CO, :N], :(5.00e-11*exp(-200.0 ./ Tn))],
             [[:CN, :O], [:NO, :C], :(5.37e-11*exp(-13800.0 ./ Tn))],
             [[:CO2, :CH], [:HCO, :CO], :((1.70e-14*Tn^0.5)*exp(-3000.0 ./ Tn))],
             [[:CO2, :H], [:CO, :OH], :(3.38e-10*exp(-13163.0 ./ Tn))],
             [[:CO2, :N], [:NO, :CO], :(3.20e-13*exp(-1710.0 ./ Tn))],
             [[:CO2, :O], [:O2, :CO], :(2.46e-11*exp(-26567.0 ./ Tn))],
             [[:H2, :C], [:CH, :H], :(6.64e-10*exp(-11700.0 ./ Tn))],
             [[:H2O, :H], [:OH, :H2], :((1.69e-14*Tn^1.2)*exp(-9610.0 ./ Tn))],
             [[:H2O, :O], [:OH, :OH], :((8.20e-14*Tn^0.95)*exp(-8571.0 ./ Tn))],
             [[:HCN, :CN], [:NCCN, :H], :((1.45e-21*Tn^1.71)*exp(-770.0 ./ Tn))],
             [[:HCN, :H], [:CN, :H2], :(6.20e-10*exp(-12500.0 ./ Tn))],
             [[:HCN, :O], [:CO, :NH], :((1.09e-15*Tn^1.14)*exp(-3742.0 ./ Tn))],
             [[:HCN, :O], [:NCO, :H], :((5.19e-16*Tn^1.38)*exp(-3693.0 ./ Tn))],
             [[:HCN, :O], [:OH, :CN], :(6.21e-10*exp(-12439.0 ./ Tn))],
             [[:HCN, :OH], [:H2O, :CN], :((3.60e-17*Tn^1.5)*exp(-3887.0 ./ Tn))],
             [[:HCO, :C], [:C2O, :H], :(1.00e-10)],
             [[:HCO, :C], [:CO, :CH], :(1.00e-10)],
             [[:HCO, :CH], [:CH2, :CO], :((5.30e-14*Tn^0.7)*exp(-500.0 ./ Tn))],
             [[:HCO, :CN], [:HCN, :CO], :(1.00e-10)],
             [[:HCO, :H], [:CO, :H2], :(1.50e-10)],
             [[:HCO, :HCO], [:CO, :CO, :H2], :(7.35e-12)],
             [[:HCO, :HCO], [:H2CO, :CO], :(4.26e-11)],
             [[:HCO, :N], [:CO, :NH], :((3.30e-13*Tn^0.5)*exp(-1000.0 ./ Tn))],
             [[:HCO, :N], [:HCN, :O], :(1.70e-10)],
             [[:HCO, :N], [:NCO, :H], :(1.00e-10)],
             [[:HCO, :NO], [:HNO, :CO], :(1.20e-11)],
             [[:HCO, :O], [:CO, :OH], :(5.00e-11)],
             [[:HCO, :O], [:CO2, :H], :(5.00e-11)],
             [[:HCO, :O2], [:CO2, :OH], :(7.60e-13)],
             [[:HCO, :O2], [:HO2, :CO], :(5.20e-12)],
             [[:HCO, :OH], [:H2O, :CO], :(1.80e-10)],
             [[:HNO, :CH], [:CH2, :NO], :(1.73e-11)],
             [[:HNO, :CN], [:HCN, :NO], :(3.00e-11)],
             [[:HNO, :CO], [:CO2, :NH], :(3.32e-12*exp(-6170.0 ./ Tn))],
             [[:HNO, :H], [:NH2, :O], :((5.81e-9*Tn^-0.3)*exp(-14730.0 ./ Tn))],
             [[:HNO, :H], [:NO, :H2], :((7.41e-13*Tn^0.72)*exp(-329.0 ./ Tn))],
             [[:HNO, :HCO], [:H2CO, :NO], :(1.00e-12*exp(-1000.0 ./ Tn))],
             [[:HNO, :N], [:N2O, :H], :((8.26e-14*Tn^0.5)*exp(-1500.0 ./ Tn))],
             [[:HNO, :N], [:NO, :NH], :((1.70e-13*Tn^0.5)*exp(-1000.0 ./ Tn))],
             [[:HNO, :O], [:NO, :OH], :(6.00e-11*Tn^-0.08)],
             [[:HNO, :O], [:NO2, :H], :(1.00e-12)],
             [[:HNO, :O], [:O2, :NH], :((1.70e-13*Tn^0.5)*exp(-3500.0 ./ Tn))],
             [[:HNO, :OH], [:H2O, :NO], :((5.54e-15*Tn^1.23)*exp(44.3 ./ Tn))],
             [[:HO2, :CH], [:CH2, :O2], :((1.70e-14*Tn^0.5)*exp(-7550.0 ./ Tn))],
             [[:HO2, :CH], [:HCO, :OH], :((8.31e-13*Tn^0.5)*exp(-3000.0 ./ Tn))],
             [[:HO2, :CO], [:CO2, :OH], :(5.60e-10*exp(-12160.0 ./ Tn))],
             [[:HO2, :H2], [:H2O2, :H], :(5.00e-11*exp(-13110.0 ./ Tn))],
             [[:HO2, :H2O], [:H2O2, :OH], :(4.65e-11*exp(-16500.0 ./ Tn))],
             [[:HO2, :HCO], [:H2CO, :O2], :(5.00e-11)],
             [[:HO2, :N], [:NO, :OH], :(2.20e-11)],
             [[:HO2, :N], [:O2, :NH], :(1.70e-13)],
             [[:HO2, :NO], [:NO2, :OH], :(3.30e-12*exp(270.0 ./ Tn))],
             [[:HOCO, :OH], [:CO2, :H2O], :(1.03e-11)],
             [[:N, :C], [:CN], :((3.49e-19*Tn^0.14)*exp(-0.18 ./ Tn))],
             [[:N2O, :CO], [:CO2, :N2], :(1.62e-13*exp(-8780.0 ./ Tn))],
             [[:N2O, :H], [:NO, :NH], :((1.11e-1*Tn^-2.16)*exp(-18700.0 ./ Tn))],
             [[:N2O, :H], [:OH, :N2], :((8.08e-22*Tn^3.15)*exp(-3603.0 ./ Tn))],
             [[:N2O, :NO], [:NO2, :N2], :((8.74e-19*Tn^2.23)*exp(-23292.0 ./ Tn))],
             [[:N2O, :O], [:NO, :NO], :(1.15e-10*exp(-13400.0 ./ Tn))],
             [[:N2O, :O], [:O2, :N2], :(1.66e-10*exp(-14100.0 ./ Tn))],
             [[:N2O, :OH], [:HO2, :N2], :(3.70e-13*exp(-2740.0 ./ Tn))],
             [[:NH, :C], [:CH, :N], :((9.99e-13*Tn^0.5)*exp(-4000.0 ./ Tn))],
             [[:NH, :C], [:CN, :H], :(1.20e-10)],
             [[:NH, :H], [:H2, :N], :((9.99e-13*Tn^0.5)*exp(-2400.0 ./ Tn))],
             [[:NH, :N], [:N2, :H], :(4.98e-11)],
             [[:NH, :NH], [:N2, :H, :H], :(1.16e-9)],
             [[:NH, :NH], [:N2, :H2], :(1.70e-11)],
             [[:NH, :NH], [:NH2, :N], :((6.29e-18*Tn^1.8)*exp(70.0 ./ Tn))],
             [[:NH, :O], [:NO, :H], :(1.80e-10*exp(-300.0 ./ Tn))],
             [[:NH, :O], [:OH, :N], :(1.16e-11)],
             [[:NH2, :C], [:HCN, :H], :((5.77e-11*Tn^-0.1)*exp(9.0 ./ Tn))],
             [[:NH2, :C], [:HNC, :H], :((5.77e-11*Tn^-0.1)*exp(9.0 ./ Tn))],
             [[:NH2, :H], [:NH, :H2], :((1.36e-14*Tn^1.02)*exp(-2161.0 ./ Tn))],
             [[:NH2, :H2], [:NH3, :H], :((4.74e-25*Tn^3.89)*exp(-1400.0 ./ Tn))],
             [[:NH2, :NO], [:H2O, :N2], :(3.60e-12*exp(450.0 ./ Tn))],
             [[:NH2, :O], [:HNO, :H], :(1.11e-10*Tn^-0.1)],
             [[:NH2, :O], [:OH, :NH], :(1.24e-11*Tn^-0.1)],
             [[:NH2, :OH], [:H2O, :NH], :((1.08e-15*Tn^1.25)*exp(43.5 ./ Tn))],
             [[:NH2, :OH], [:NH3, :O], :((2.73e-15*Tn^0.76)*exp(-262.0 ./ Tn))],
             [[:NO, :C], [:CN, :O], :(1.49e-10*Tn^-0.16)],
             [[:NO, :C], [:CO, :N], :(2.24e-10*Tn^-0.16)],
             [[:NO, :CH], [:CO, :NH], :(1.52e-11)],
             [[:NO, :CH], [:HCN, :O], :(1.31e-10)],
             [[:NO, :CH], [:HCO, :N], :(1.14e-11)],
             [[:NO, :CH], [:NCO, :H], :(2.47e-11)],
             [[:NO, :CH], [:OH, :CN], :(1.90e-12)],
             [[:NO, :CN], [:CO, :N2], :(1.60e-13)],
             [[:NO, :H], [:NH, :O], :((1.64e-9*Tn^-0.1)*exp(-35220.0 ./ Tn))],
             [[:NO, :N], [:N2, :O], :(2.10e-11*exp(100.0 ./ Tn))],
             [[:NO, :NH], [:N2, :O, :H], :(7.40e-10*exp(-10540.0 ./ Tn))],
             [[:NO, :NH], [:N2O, :H], :((4.56e-9*Tn^-0.78)*exp(-40.0 ./ Tn))],
             [[:NO, :NH], [:OH, :N2], :((1.14e-9*Tn^-0.78)*exp(-40.0 ./ Tn))],
             [[:NO, :O], [:O2, :N], :(1.18e-11*exp(-20413.0 ./ Tn))],
             [[:NO2, :CN], [:CO2, :N2], :((6.12e-11*Tn^-0.752)*exp(-173.0 ./ Tn))],
             [[:NO2, :CN], [:N2O, :CO], :((8.16e-11*Tn^-0.752)*exp(-173.0 ./ Tn))],
             [[:NO2, :CN], [:NCO, :NO], :((8.77e-10*Tn^-0.752)*exp(-173.0 ./ Tn))],
             [[:NO2, :CO], [:CO2, :NO], :(1.48e-10*exp(-17000.0 ./ Tn))],
             [[:NO2, :H], [:NO, :OH], :(4.00e-10*exp(-340.0 ./ Tn))],
             [[:NO2, :N], [:N2, :O, :O], :(2.41e-12)],
             [[:NO2, :N], [:N2O, :O], :(5.80e-12*exp(220.0 ./ Tn))],
             [[:NO2, :N], [:NO, :NO], :(1.00e-12)],
             [[:NO2, :N], [:O2, :N2], :(1.00e-12)],
             [[:NO2, :NH], [:HNO, :NO], :((1.56e-6*Tn^-1.94)*exp(-56.9 ./ Tn))],
             [[:NO2, :NH], [:N2O, :OH], :((1.09e-6*Tn^-1.94)*exp(-56.9 ./ Tn))],
             [[:NO2, :NH2], [:N2O, :H2O], :(2.10e-12*exp(650.0 ./ Tn))],
             [[:NO2, :O], [:O2, :NO], :(5.10e-12*exp(210.0 ./ Tn))],
             [[:O, :C], [:CO], :((1.75e-19*Tn^0.705)*exp(-136.0 ./ Tn))],
             [[:O, :H], [:OH], :(8.65e-18*Tn^-0.38)],
             [[:O1D, :CO], [:CO, :O], :(4.70e-11*exp(63.0 ./ Tn))],
             [[:O1D, :CO], [:CO2], :(8.00e-11)],
             [[:O1D, :H2O2], [:H2O2, :O], :(5.20e-10)],
             [[:O1D, :N2], [:N2, :O], :(2.15e-11*exp(110.0 ./ Tn))],
             [[:O1D, :N2O], [:NO, :NO], :(7.26e-11*exp(20.0 ./ Tn))],
             [[:O1D, :N2O], [:O2, :N2], :(4.64e-11*exp(20.0 ./ Tn))],
             [[:O1D, :NO], [:NO, :O], :(4.00e-11)],
             [[:O1D, :NO2], [:NO2, :O], :(1.13e-10*exp(115.0 ./ Tn))],
             [[:O1D, :NO2], [:O2, :NO], :(2.31e-10)],
             [[:O2, :C], [:CO, :O], :(3.03e-10*Tn^-0.32)],
             [[:O2, :CH], [:CO, :O, :H], :(1.20e-11)],
             [[:O2, :CH], [:CO, :OH], :(8.00e-12)],
             [[:O2, :CH], [:CO2, :H], :(1.20e-11)],
             [[:O2, :CH], [:HCO, :O], :(8.00e-12)],
             [[:O2, :CN], [:CO, :NO], :(3.00e-12*exp(210.0 ./ Tn))],
             [[:O2, :CN], [:NCO, :O], :((5.97e-11*Tn^-0.19)*exp(31.9 ./ Tn))],
             [[:O2, :CO], [:CO2, :O], :(5.99e-12*exp(-24075.0 ./ Tn))],
             [[:O2, :H], [:OH, :O], :(2.61e-10*exp(-8156.0 ./ Tn))],
             [[:O2, :H2], [:HO2, :H], :(2.40e-10*exp(-28500.0 ./ Tn))],
             [[:O2, :H2], [:OH, :OH], :(3.16e-10*exp(-21890.0 ./ Tn))],
             [[:O2, :N], [:NO, :O], :(1.50e-11*exp(-3600.0 ./ Tn))],
             [[:O2, :NH], [:HNO, :O], :(4.00e-11*exp(-6970.0 ./ Tn))],
             [[:O2, :NH], [:NO, :OH], :(1.50e-13*exp(-770.0 ./ Tn))],
             [[:O3, :N], [:O2, :NO], :(1.00e-16)],
             [[:O3, :NO], [:NO2, :O2], :(3.00e-12*exp(-1500.0 ./ Tn))],
             [[:O3, :NO2], [:NO3, :O2], :(1.20e-13*exp(-2450.0 ./ Tn))],
             [[:OH, :C], [:CO, :H], :((7.98e-10*Tn^-0.34)*exp(-0.108 ./ Tn))],
             [[:OH, :CH], [:HCO, :H], :((8.31e-13*Tn^0.5)*exp(-5000.0 ./ Tn))],
             [[:OH, :H], [:H2, :O], :((8.10e-21*Tn^2.8)*exp(-1950.0 ./ Tn))],
             [[:OH, :N], [:NH, :O], :((1.06e-11*Tn^0.1)*exp(-10700.0 ./ Tn))],
             [[:OH, :N], [:NO, :H], :(1.80e-10*Tn^-0.2)],
             [[:OH, :NH], [:H2O, :N], :(3.30e-15*Tn^1.2)],
             [[:OH, :NH], [:HNO, :H], :(3.30e-11)],
             [[:OH, :NH], [:NH2, :O], :((1.66e-12*Tn^0.1)*exp(-5800.0 ./ Tn))],
             [[:OH, :NH], [:NO, :H2], :(4.16e-11)],

             # Type 4
             [[:CH, :H2], [:CH2, :H], :(0.0 + 0.6 * (8.5e-11*(Tn^0.15) - 0.0) * 4.7e-26*(Tn^-1.6) * M / (8.5e-11*(Tn^0.15) - 0.0 + 4.7e-26*(Tn^-1.6) * M))],
             [[:CO, :H], [:HCO], :(0.0 + 0.0 * (1.0*(Tn^0.2) - 0.0) * 2.0e-35*(Tn^0.2) * M / (1.0*(Tn^0.2) - 0.0 + 2.0e-35*(Tn^0.2) * M))],
             [[:CO, :O], [:CO2], :(0.0 + 0.4 * (1.0*exp(-1509.0 ./ Tn) - 0.0) * 1.7e-33*exp(-1509.0 ./ Tn) * M / (1.0*exp(-1509.0 ./ Tn) - 0.0 + 1.7e-33*exp(-1509.0 ./ Tn) * M))],
             [[:N, :H], [:NH], :(0.0 + 0.0 * (1.0 - 0.0) * 5.0e-32 * M / (1.0 - 0.0 + 5.0e-32 * M))],
             [[:NO, :H], [:HNO], :(0.0 + 0.82 * (2.53e-9*(Tn^-0.41) - 0.0) * 9.56e-29*(Tn^-1.17)*exp(-212.0 ./ Tn) * M / (2.53e-9*(Tn^-0.41) - 0.0 + 9.56e-29*(Tn^-1.17)*exp(-212.0 ./ Tn) * M))],
             [[:NO, :O], [:NO2], :(0.0 + 0.8 * (4.9e-10*(Tn^-0.4) - 0.0) * 9.2e-28*(Tn^-1.6) * M / (4.9e-10*(Tn^-0.4) - 0.0 + 9.2e-28*(Tn^-1.6) * M))],
             [[:O, :N], [:NO], :(0.0 + 0.4 * (1.0 - 0.0) * 5.46e-33*exp(155.0 ./ Tn) * M / (1.0 - 0.0 + 5.46e-33*exp(155.0 ./ Tn) * M))],
             [[:O1D, :N2], [:N2O], :(0.0 + 0.0 * (1.0*(Tn^-0.9) - 0.0) * 4.75e-34*(Tn^-0.9) * M / (1.0*(Tn^-0.9) - 0.0 + 4.75e-34*(Tn^-0.9) * M))],

             # Type 5
             # TODO: Check that these are correct. The form of k0 and kinf doesn't seem to agree with the Vuitton
             # appendix Roger sent, but it makes sense according to Sander 2011 which is presumably where that 
             # Vuitton info came from in the first place.
             [[:NO, :OH], [:HONO], threebody(:(1.93e-24*(Tn/300.)^-2.6), :(6.37e-11*(Tn/300.)^-0.1))],
             [[:NO2, :HO2], [:HO2NO2], threebody(:(5.02e-23*(Tn/300.)^-3.4), :(2.21e-11*(Tn/300.)^-0.3))],
             [[:NO2, :O], [:NO3], threebody(:(7.19e-27*(Tn/300.)^-1.8), :(1.19e-9*(Tn/300.)^-0.7))],
             [[:NO2, :OH], [:HONO2], threebody(:(4.86e-23*(Tn/300.)^-3.0), :(2.80e-11))],
             [[:NO2, :OH], [:HOONO], threebody(:(4.17e-22*(Tn/300.)^-3.9), :(7.27e12*(Tn/300.)^-0.5))],

        
             # IONOSPHERE - reactions from Roger Yelle
             # TODO: test; then add D 
             [[:ArHpl, :C], [:CHpl, :Ar], :(1.02e-9)],
             [[:ArHpl, :CO], [:HCOpl, :Ar], :(1.25e-9)],
             [[:ArHpl, :CO2], [:HCO2pl, :Ar], :(1.1e-9)],
             [[:ArHpl, :H2], [:H3pl, :Ar], :(6.3e-10)],
             [[:ArHpl, :N2], [:N2Hpl, :Ar], :(8.0e-10)],
             [[:ArHpl, :O], [:OHpl, :Ar], :(5.9e-10)],
             [[:ArHpl, :O2], [:HO2pl, :Ar], :(5.05e-10)],
             [[:Arpl, :CO], [:COpl, :Ar], :(4.4e-11)],
             [[:Arpl, :CO2], [:CO2pl, :Ar], :(4.8e-10)],
             [[:Arpl, :H2], [:ArHpl, :H], :(8.72e-10)],
             [[:Arpl, :H2], [:H2pl, :Ar], :(1.78e-11)],
             [[:Arpl, :H2O], [:ArHpl, :OH], :(3.24e-10)],
             [[:Arpl, :H2O], [:H2Opl, :Ar], :(1.3e-9)],
             [[:Arpl, :N2], [:N2pl, :Ar], :(1.1e-11)],
             [[:Arpl, :N2O], [:N2Opl, :Ar], :(2.91e-10)],
             [[:Arpl, :N2O], [:N2pl, :Ar, :O], :(3.0e-12)],
             [[:Arpl, :N2O], [:NOpl, :Ar, :N], :(3.0e-12)],
             [[:Arpl, :N2O], [:Opl, :N2, :Ar], :(3.0e-12)],
             [[:Arpl, :NO], [:NOpl, :Ar], :(3.1e-10)],
             [[:Arpl, :NO2], [:NO2pl, :Ar], :(2.76e-11)],
             [[:Arpl, :NO2], [:NOpl, :Ar, :O], :(4.32e-10)],
             [[:Arpl, :O2], [:O2pl, :Ar], :(4.6e-11)],
             [[:CHpl, :C], [:C2pl, :H], :(1.2e-9)],
             [[:CHpl, :CN], [:C2Npl, :H], :(9.53e-9*((300 ./ Ti)^-0.5))],
             [[:CHpl, :CO], [:HCOpl, :C], :(7.0e-12*((300 ./ Ti)^-0.5))],
             [[:CHpl, :CO2], [:HCOpl, :CO], :(1.6e-9)],
             [[:CHpl, :H], [:Cpl, :H2], :(7.5e-10)],
             [[:CHpl, :H2], [:CH2pl, :H], :(1.2e-9)],
             [[:CHpl, :H2O], [:H2COpl, :H], :(1.0e-8*((300 ./ Ti)^-0.5))],
             [[:CHpl, :H2O], [:H3Opl, :C], :(1.45e-9)],
             [[:CHpl, :H2O], [:HCOpl, :H2], :(5.02e-8*((300 ./ Ti)^-0.5))],
             [[:CHpl, :HCN], [:C2Npl, :H2], :(4.2e-10)],
             [[:CHpl, :HCN], [:HC2Npl, :H], :(2.8e-10)],
             [[:CHpl, :HCN], [:HCNHpl, :C], :(2.1e-9)],
             [[:CHpl, :HCO], [:CH2pl, :CO], :(7.97e-9*((300 ./ Ti)^-0.5))],
             [[:CHpl, :HCO], [:HCOpl, :CH], :(7.97e-9*((300 ./ Ti)^-0.5))],
             [[:CHpl, :N], [:CNpl, :H], :(1.9e-10)],
             [[:CHpl, :NH], [:CNpl, :H2], :(1.32e-8*((300 ./ Ti)^-0.5))],
             [[:CHpl, :NO], [:NOpl, :CH], :(7.6e-10)],
             [[:CHpl, :O], [:COpl, :H], :(3.5e-10)],
             [[:CHpl, :O2], [:HCOpl, :O], :(8.73e-10)],
             [[:CHpl, :OH], [:COpl, :H2], :(1.3e-8*((300 ./ Ti)^-0.5))],
             [[:CNpl, :C], [:Cpl, :CN], :(1.1e-10)],
             [[:CNpl, :CH], [:CHpl, :CN], :(1.11e-8*((300 ./ Ti)^-0.5))],
             [[:CNpl, :CO], [:COpl, :CN], :(4.4e-10)],
             [[:CNpl, :CO2], [:C2Opl, :NO], :(3.3e-10)],
             [[:CNpl, :CO2], [:CO2pl, :CN], :(4.4e-10)],
             [[:CNpl, :CO2], [:OCNpl, :CO], :(3.3e-10)],
             [[:CNpl, :H], [:Hpl, :CN], :(6.4e-10)],
             [[:CNpl, :H2], [:HCNpl, :H], :(8.0e-10)],
             [[:CNpl, :H2], [:HNCpl, :H], :(8.0e-10)],
             [[:CNpl, :H2O], [:H2CNpl, :O], :(4.8e-10)],
             [[:CNpl, :H2O], [:H2Opl, :CN], :(3.2e-10)],
             [[:CNpl, :H2O], [:HCNpl, :OH], :(1.6e-9)],
             [[:CNpl, :H2O], [:HCOpl, :NH], :(1.6e-10)],
             [[:CNpl, :H2O], [:HNCOpl, :H], :(6.4e-10)],
             [[:CNpl, :HCN], [:C2N2pl, :H], :(4.59e-10)],
             [[:CNpl, :HCN], [:HCNpl, :CN], :(2.24e-9)],
             [[:CNpl, :HCO], [:HCNpl, :CO], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:CNpl, :HCO], [:HCOpl, :CN], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:CNpl, :N], [:N2pl, :C], :(6.1e-10)],
             [[:CNpl, :N2O], [:N2Opl, :CN], :(4.56e-10)],
             [[:CNpl, :N2O], [:NOpl, :CN2], :(1.52e-10)],
             [[:CNpl, :N2O], [:OCNpl, :N2], :(1.52e-10)],
             [[:CNpl, :NH], [:NHpl, :CN], :(1.13e-8*((300 ./ Ti)^-0.5))],
             [[:CNpl, :NH2], [:NH2pl, :CN], :(1.58e-8*((300 ./ Ti)^-0.5))],
             [[:CNpl, :NO], [:NOpl, :CN], :(5.7e-10)],
             [[:CNpl, :NO], [:OCNpl, :N], :(1.9e-10)],
             [[:CNpl, :O], [:Opl, :CN], :(6.5e-11)],
             [[:CNpl, :O2], [:NOpl, :CO], :(8.6e-11)],
             [[:CNpl, :O2], [:O2pl, :CN], :(2.58e-10)],
             [[:CNpl, :O2], [:OCNpl, :O], :(8.6e-11)],
             [[:CNpl, :OH], [:OHpl, :CN], :(1.11e-8*((300 ./ Ti)^-0.5))],
             [[:CO2pl, :H], [:HCOpl, :O], :(4.47e-10)],
             [[:CO2pl, :H], [:Hpl, :CO2], :(5.53e-11)],
             [[:CO2pl, :H2], [:HCO2pl, :H], :(2.24e-9*((300 ./ Ti)^-0.15))], # Nair minimal ionosphere
             [[:CO2pl, :H2O], [:H2Opl, :CO2], :(1.8e-9)],
             [[:CO2pl, :H2O], [:HCO2pl, :OH], :(6.0e-10)],
             [[:CO2pl, :HCN], [:HCNpl, :CO2], :(8.1e-10)],
             [[:CO2pl, :HCN], [:HCO2pl, :CN], :(9.0e-11)],
             [[:CO2pl, :N], [:COpl, :NO], :(3.4e-10)],
             [[:CO2pl, :NO], [:NOpl, :CO2], :(1.23e-10)],
             [[:CO2pl, :O], [:O2pl, :CO], :(1.6e-10)], # Nair minimal ionosphere
             [[:CO2pl, :O], [:Opl, :CO2], :(1.0e-10)], # Nair minimal ionosphere
             [[:CO2pl, :O2], [:O2pl, :CO2], :(5.5e-11)],
             [[:COpl, :C], [:Cpl, :CO], :(1.1e-10)],
             [[:COpl, :CH], [:CHpl, :CO], :(5.54e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :CH], [:HCOpl, :C], :(5.54e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :CO2], [:CO2pl, :CO], :(1.1e-9)], 
             [[:COpl, :H], [:Hpl, :CO], :(4.0e-10)],
             [[:COpl, :H2], [:HCOpl, :H], :(7.5e-10)],
             [[:COpl, :H2], [:HOCpl, :H], :(7.5e-10)],
             [[:COpl, :H2O], [:H2Opl, :CO], :(1.56e-9)],
             [[:COpl, :H2O], [:HCOpl, :OH], :(8.4e-10)],
             [[:COpl, :HCN], [:HCNpl, :CO], :(3.06e-9)],
             [[:COpl, :HCO], [:HCOpl, :CO], :(1.28e-8*((300 ./ Ti)^-0.5))],
             [[:COpl, :N], [:NOpl, :C], :(8.2e-11)],
             [[:COpl, :NH], [:HCOpl, :N], :(5.54e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :NH], [:NHpl, :CO], :(5.54e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :NH2], [:HCOpl, :NH], :(7.79e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :NH2], [:NH2pl, :CO], :(7.79e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :NO], [:NOpl, :CO], :(4.2e-10)],
             [[:COpl, :O], [:Opl, :CO], :(1.4e-10)], 
             [[:COpl, :O2], [:O2pl, :CO], :(1.5e-10)],
             [[:COpl, :OH], [:HCOpl, :O], :(5.37e-9*((300 ./ Ti)^-0.5))],
             [[:COpl, :OH], [:OHpl, :CO], :(5.37e-9*((300 ./ Ti)^-0.5))],
             [[:Cpl, :C], [:C2pl], :(1.52e-18*((300 ./ Ti)^0.17)*exp.(-101.5 ./ Ti))],
             [[:Cpl, :CH], [:C2pl, :H], :(6.58e-9*((300 ./ Ti)^-0.5))],
             [[:Cpl, :CH], [:CHpl, :C], :(6.58e-9*((300 ./ Ti)^-0.5))],
             [[:Cpl, :CO2], [:CO2pl, :C], :(1.1e-10)],
             [[:Cpl, :CO2], [:COpl, :CO], :(9.9e-10)],
             [[:Cpl, :H], [:CHpl], :(1.7e-17)],
             [[:Cpl, :H2], [:CH2pl], :(3.32e-13*((300 ./ Ti)^-1.3)*exp.(-23.0 ./ Ti))],
             [[:Cpl, :H2], [:CHpl, :H], :(7.4e-10*exp.(-4537.0 ./ Ti))],
             [[:Cpl, :H2O], [:H2Opl, :C], :(2.4e-10)],
             [[:Cpl, :H2O], [:HCOpl, :H], :(1.56e-8*((300 ./ Ti)^-0.5))],
             [[:Cpl, :H2O], [:HOCpl, :H], :(2.16e-9)],
             [[:Cpl, :HCN], [:C2Npl, :H], :(2.95e-9)],
             [[:Cpl, :HCO], [:CHpl, :CO], :(8.31e-9*((300 ./ Ti)^-0.5))],
             [[:Cpl, :HCO], [:HCOpl, :C], :(8.31e-9*((300 ./ Ti)^-0.5))],
             [[:Cpl, :N], [:CNpl], :(7.24e-19*((300 ./ Ti)^0.07)*exp.(-57.5 ./ Ti))],
             [[:Cpl, :N2O], [:NOpl, :CN], :(9.1e-10)],
             [[:Cpl, :NH], [:CNpl, :H], :(1.35e-8*((300 ./ Ti)^-0.5))],
             [[:Cpl, :NH2], [:HCNpl, :H], :(1.91e-8*((300 ./ Ti)^-0.5))],
             [[:Cpl, :NO], [:NOpl, :C], :(7.5e-10)],
             [[:Cpl, :O], [:COpl], :(7.39e-18*((300 ./ Ti)^-0.15)*exp.(-68.0 ./ Ti))],
             [[:Cpl, :O2], [:COpl, :O], :(3.48e-10)],
             [[:Cpl, :O2], [:Opl, :CO], :(5.22e-10)],
             [[:Cpl, :OH], [:COpl, :H], :(1.33e-8*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :C], [:CHpl, :OH], :(1.1e-9)],
             [[:H2Opl, :CH], [:CH2pl, :OH], :(5.89e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :CH], [:CHpl, :H2O], :(5.89e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :CO], [:HCOpl, :OH], :(4.25e-10)],
             [[:H2Opl, :H2], [:H3Opl, :H], :(7.6e-10)],
             [[:H2Opl, :H2O], [:H3Opl, :OH], :(1.85e-9)],
             [[:H2Opl, :HCN], [:HCNHpl, :OH], :(1.05e-9)],
             [[:H2Opl, :HCO], [:H2COpl, :OH], :(4.85e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :HCO], [:H3Opl, :CO], :(4.85e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :HCO], [:HCOpl, :H2O], :(4.85e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :N], [:HNOpl, :H], :(1.12e-10)],
             [[:H2Opl, :N], [:NOpl, :H2], :(2.8e-11)],
             [[:H2Opl, :NH], [:H3Opl, :N], :(1.23e-8*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :NH2], [:NH2pl, :H2O], :(8.49e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :NH2], [:NH3pl, :OH], :(8.49e-9*((300 ./ Ti)^-0.5))],
             [[:H2Opl, :NO], [:NOpl, :H2O], :(4.6e-10)],
             [[:H2Opl, :NO2], [:NO2pl, :H2O], :(1.2e-9)],
             [[:H2Opl, :O], [:O2pl, :H2], :(4.0e-11)],
             [[:H2Opl, :O2], [:O2pl, :H2O], :(3.3e-10)],
             [[:H2Opl, :OH], [:H3Opl, :O], :(1.2e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :Ar], [:ArHpl, :H], :(2.1e-9)],
             [[:H2pl, :C], [:CHpl, :H], :(2.4e-9)],
             [[:H2pl, :CH], [:CH2pl, :H], :(1.23e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :CH], [:CHpl, :H2], :(1.23e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :CN], [:CNpl, :H2], :(2.08e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :CN], [:HCNpl, :H], :(2.08e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :CO], [:COpl, :H2], :(6.44e-10)],
             [[:H2pl, :CO], [:HCOpl, :H], :(2.9e-9)],
             [[:H2pl, :CO2], [:HCO2pl, :H], :(2.35e-9)],
             [[:H2pl, :H], [:Hpl, :H2], :(6.4e-10)],
             [[:H2pl, :H2], [:H3pl, :H], :(2.0e-9)],
             [[:H2pl, :H2O], [:H2Opl, :H2], :(3.87e-9)],
             [[:H2pl, :H2O], [:H3Opl, :H], :(3.43e-9)],
             [[:H2pl, :HCN], [:HCNpl, :H2], :(4.68e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :HCO], [:H3pl, :CO], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :HCO], [:HCOpl, :H2], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :N], [:NHpl, :H], :(1.9e-9)],
             [[:H2pl, :N2], [:N2Hpl, :H], :(2.0e-9)],
             [[:H2pl, :N2O], [:HN2Opl, :H], :(1.32e-9)],
             [[:H2pl, :N2O], [:N2Hpl, :OH], :(7.77e-10)],
             [[:H2pl, :NH], [:NH2pl, :H], :(1.32e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :NH], [:NHpl, :H2], :(1.32e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :NO], [:HNOpl, :H], :(1.1e-9)],
             [[:H2pl, :NO], [:NOpl, :H2], :(1.1e-9)],
             [[:H2pl, :O], [:OHpl, :H], :(1.5e-9)],
             [[:H2pl, :O2], [:HO2pl, :H], :(1.92e-9)],
             [[:H2pl, :O2], [:O2pl, :H2], :(7.83e-10)],
             [[:H2pl, :OH], [:H2Opl, :H], :(1.32e-8*((300 ./ Ti)^-0.5))],
             [[:H2pl, :OH], [:OHpl, :H2], :(1.32e-8*((300 ./ Ti)^-0.5))],
             [[:H3Opl, :C], [:HCOpl, :H2], :(1.0e-11)],
             [[:H3Opl, :CH], [:CH2pl, :H2O], :(1.18e-8*((300 ./ Ti)^-0.5))],
             [[:H3Opl, :HCN], [:HCNHpl, :H2O], :(3.8e-9)],
             [[:H3Opl, :NH2], [:NH3pl, :H2O], :(1.68e-8*((300 ./ Ti)^-0.5))],
             [[:H3pl, :Ar], [:ArHpl, :H2], :(3.65e-10)],
             [[:H3pl, :C], [:CHpl, :H2], :(2.0e-9)],
             [[:H3pl, :CH], [:CH2pl, :H2], :(2.08e-8*((300 ./ Ti)^-0.5))],
             [[:H3pl, :CN], [:HCNpl, :H2], :(3.46e-8*((300 ./ Ti)^-0.5))],
             [[:H3pl, :CO], [:HCOpl, :H2], :(3.06e-9*((300 ./ Ti)^-0.142)*exp.(3.41 ./ Ti))],
             [[:H3pl, :CO], [:HOCpl, :H2], :(5.82e-10*((300 ./ Ti)^0.0661)*exp.(-5.21 ./ Ti))],
             [[:H3pl, :CO2], [:HCO2pl, :H2], :(2.5e-9)],
             [[:H3pl, :H2O], [:H3Opl, :H2], :(5.3e-9)],
             [[:H3pl, :HCN], [:HCNHpl, :H2], :(7.5e-9)],
             [[:H3pl, :HCO], [:H2COpl, :H2], :(2.94e-8*((300 ./ Ti)^-0.5))],
             [[:H3pl, :N], [:NH2pl, :H], :(3.9e-10)],
             [[:H3pl, :N], [:NHpl, :H2], :(2.6e-10)],
             [[:H3pl, :N2], [:N2Hpl, :H2], :(1.63e-9)],
             [[:H3pl, :N2O], [:HN2Opl, :H2], :(2.5e-9)],
             [[:H3pl, :NH], [:NH2pl, :H2], :(2.25e-8*((300 ./ Ti)^-0.5))],
             [[:H3pl, :NO], [:HNOpl, :H2], :(1.94e-9)],
             [[:H3pl, :NO2], [:NO2pl, :H2, :H], :(7.0e-12)],
             [[:H3pl, :NO2], [:NOpl, :OH, :H2], :(6.93e-10)],
             [[:H3pl, :O], [:H2Opl, :H], :(8.33e-10*((300 ./ Ti)^-0.156)*exp.(-1.4 ./ Ti))],
             [[:H3pl, :O], [:OHpl, :H2], :(1.94e-9*((300 ./ Ti)^-0.156)*exp.(-1.4 ./ Ti))],
             [[:H3pl, :O2], [:HO2pl, :H2], :(6.7e-10)],
             [[:H3pl, :OH], [:H2Opl, :H2], :(2.25e-8*((300 ./ Ti)^-0.5))],
             [[:HCNHpl, :CH], [:CH2pl, :HCN], :(5.46e-9*((300 ./ Ti)^-0.5))],
             [[:HCNHpl, :H2O], [:H3Opl, :HCN], :(8.8e-13)],
             [[:HCNpl, :C], [:CHpl, :CN], :(1.1e-9)],
             [[:HCNpl, :CH], [:CH2pl, :CN], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :CO], [:HCOpl, :CN], :(1.38e-10)],
             [[:HCNpl, :CO], [:HNCpl, :CO], :(3.22e-10)],
             [[:HCNpl, :CO2], [:HCO2pl, :CN], :(2.1e-10)],
             [[:HCNpl, :CO2], [:HNCpl, :CO2], :(2.9e-10)],
             [[:HCNpl, :H], [:Hpl, :HCN], :(3.7e-11)],
             [[:HCNpl, :H2], [:HCNHpl, :H], :(8.8e-10)],
             [[:HCNpl, :H2O], [:H2Opl, :HCN], :(3.12e-8*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :H2O], [:H3Opl, :CN], :(3.12e-8*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :HCN], [:HCNHpl, :CN], :(1.45e-9)],
             [[:HCNpl, :HCO], [:H2COpl, :CN], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :HCO], [:HCNHpl, :CO], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :N], [:CHpl, :N2], :(2.2e-10)],
             [[:HCNpl, :N2O], [:N2Opl, :HCN], :(1.08e-9)],
             [[:HCNpl, :NH], [:NH2pl, :CN], :(1.13e-8*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :NH2], [:NH3pl, :CN], :(1.56e-8*((300 ./ Ti)^-0.5))],
             [[:HCNpl, :NO], [:NOpl, :HCN], :(8.1e-10)],
             [[:HCNpl, :O2], [:O2pl, :HCN], :(4.1e-10)],
             [[:HCNpl, :OH], [:H2Opl, :CN], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:HCO2pl, :C], [:CHpl, :CO2], :(1.0e-9)],
             [[:HCO2pl, :CO], [:HCOpl, :CO2], :(7.8e-10)],
             [[:HCO2pl, :H2O], [:H3Opl, :CO2], :(2.65e-9)],
             [[:HCO2pl, :O], [:HCOpl, :O2], :(5.8e-10)],
             [[:HCOpl, :C], [:CHpl, :CO], :(1.1e-9)],
             [[:HCOpl, :CH], [:CH2pl, :CO], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:HCOpl, :H2O], [:H3Opl, :CO], :(2.6e-9)],
             [[:HCOpl, :H2O], [:HCOOH2pl], :(6.64e-10*((300 ./ Ti)^-1.3))],
             [[:HCOpl, :HCN], [:HCNHpl, :CO], :(3.5e-9)],
             [[:HCOpl, :HCO], [:H2COpl, :CO], :(1.26e-8*((300 ./ Ti)^-0.5))],
             [[:HCOpl, :N2O], [:HN2Opl, :CO], :(3.3e-12)],
             [[:HCOpl, :NH], [:NH2pl, :CO], :(1.11e-8*((300 ./ Ti)^-0.5))],
             [[:HCOpl, :NH2], [:NH3pl, :CO], :(1.54e-8*((300 ./ Ti)^-0.5))],
             [[:HCOpl, :OH], [:H2Opl, :CO], :(1.07e-8*((300 ./ Ti)^-0.5))],
             [[:HCOpl, :OH], [:HCO2pl, :H], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:HN2Opl, :CO], [:HCOpl, :N2O], :(5.3e-10)],
             [[:HN2Opl, :H2O], [:H3Opl, :N2O], :(2.83e-9)],
             [[:HNOpl, :C], [:CHpl, :NO], :(1.0e-9)],
             [[:HNOpl, :CH], [:CH2pl, :NO], :(1.07e-8*((300 ./ Ti)^-0.5))],
             [[:HNOpl, :CN], [:HCNpl, :NO], :(1.51e-8*((300 ./ Ti)^-0.5))],
             [[:HNOpl, :CO], [:HCOpl, :NO], :(8.6e-10)],
             [[:HNOpl, :CO2], [:HCO2pl, :NO], :(9.4e-10)],
             [[:HNOpl, :H2O], [:H3Opl, :NO], :(2.3e-9)],
             [[:HNOpl, :HCN], [:HCNHpl, :NO], :(1.71e-8*((300 ./ Ti)^-0.5))],
             [[:HNOpl, :HCO], [:H2COpl, :NO], :(1.25e-8*((300 ./ Ti)^-0.5))],
             [[:HNOpl, :N2], [:N2Hpl, :NO], :(1.0e-11)],
             [[:HNOpl, :NH], [:NH2pl, :NO], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:HNOpl, :NH2], [:NH3pl, :NO], :(1.52e-8*((300 ./ Ti)^-0.5))],
             [[:HNOpl, :NO], [:NOpl, :HNO], :(7.0e-10)],
             [[:HNOpl, :O], [:NO2pl, :H], :(1.0e-12)],
             [[:HNOpl, :OH], [:H2Opl, :NO], :(1.07e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :C], [:CHpl, :O2], :(1.0e-9)],
             [[:HO2pl, :CH], [:CH2pl, :O2], :(1.07e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :CN], [:HCNpl, :O2], :(1.49e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :CO], [:HCOpl, :O2], :(8.4e-10)],
             [[:HO2pl, :CO2], [:HCO2pl, :O2], :(1.1e-9)],
             [[:HO2pl, :H2], [:H3pl, :O2], :(3.3e-10)],
             [[:HO2pl, :H2O], [:H3Opl, :O2], :(1.42e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :HCN], [:HCNHpl, :O2], :(1.68e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :HCO], [:H2COpl, :O2], :(1.23e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :N], [:NO2pl, :H], :(1.0e-12)],
             [[:HO2pl, :N2], [:N2Hpl, :O2], :(8.0e-10)],
             [[:HO2pl, :NH], [:NH2pl, :O2], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :NH2], [:NH3pl, :O2], :(1.51e-8*((300 ./ Ti)^-0.5))],
             [[:HO2pl, :NO], [:HNOpl, :O2], :(7.7e-10)],
             [[:HO2pl, :O], [:OHpl, :O2], :(6.2e-10)],
             [[:HO2pl, :OH], [:H2Opl, :O2], :(1.06e-8*((300 ./ Ti)^-0.5))],
             [[:HOCpl, :CO], [:HCOpl, :CO], :(6.0e-10)],
             [[:HOCpl, :CO2], [:HCO2pl, :CO], :(9.45e-10)],
             [[:HOCpl, :H2], [:H3pl, :CO], :(2.68e-10)],
             [[:HOCpl, :H2], [:HCOpl, :H2], :(3.8e-10)],
             [[:HOCpl, :N2], [:N2Hpl, :CO], :(6.7e-10)],
             [[:HOCpl, :N2O], [:HN2Opl, :CO], :(1.17e-9)],
             [[:HOCpl, :NO], [:HNOpl, :CO], :(7.1e-10)],
             [[:HOCpl, :O2], [:HO2pl, :CO], :(1.9e-10)],
             [[:Hpl, :CH], [:CHpl, :H], :(3.29e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :CO2], [:HCOpl, :O], :(3.8e-9)],
             [[:Hpl, :H], [:H2pl], :(2.34e-22*((300 ./ Ti)^1.49)*exp.(-228.0 ./ Ti))],
             [[:Hpl, :H2O], [:H2Opl, :H], :(8.2e-9)],
             [[:Hpl, :HCN], [:HCNpl, :H], :(1.1e-8)],
             [[:Hpl, :HCO], [:COpl, :H2], :(1.63e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :HCO], [:H2pl, :CO], :(1.63e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :HCO], [:HCOpl, :H], :(1.63e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :HNO], [:NOpl, :H2], :(6.93e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :N2O], [:N2Hpl, :O], :(3.52e-10)],
             [[:Hpl, :N2O], [:N2Opl, :H], :(1.85e-9)],
             [[:Hpl, :NH], [:NHpl, :H], :(3.64e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :NH2], [:NH2pl, :H], :(5.02e-8*((300 ./ Ti)^-0.5))],
             [[:Hpl, :NO], [:NOpl, :H], :(1.9e-9)],
             [[:Hpl, :NO2], [:NOpl, :OH], :(1.9e-9)],
             [[:Hpl, :O], [:Opl, :H], :(3.75e-10)],
             [[:Hpl, :O2], [:O2pl, :H], :(1.17e-9)],
             [[:Hpl, :OH], [:OHpl, :H], :(3.64e-8*((300 ./ Ti)^-0.5))],
             [[:N2Hpl, :C], [:CHpl, :N2], :(1.1e-9)],
             [[:N2Hpl, :CH], [:CH2pl, :N2], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:N2Hpl, :CO], [:HCOpl, :N2], :(8.8e-10)],
             [[:N2Hpl, :CO2], [:HCO2pl, :N2], :(1.07e-9)],
             [[:N2Hpl, :H2], [:H3pl, :N2], :(5.1e-18)],
             [[:N2Hpl, :H2O], [:H3Opl, :N2], :(2.6e-9)],
             [[:N2Hpl, :HCN], [:HCNHpl, :N2], :(3.2e-9)],
             [[:N2Hpl, :HCO], [:H2COpl, :N2], :(1.26e-8*((300 ./ Ti)^-0.5))],
             [[:N2Hpl, :N2O], [:HN2Opl, :N2], :(1.25e-9)],
             [[:N2Hpl, :NH], [:NH2pl, :N2], :(1.11e-8*((300 ./ Ti)^-0.5))],
             [[:N2Hpl, :NH2], [:NH3pl, :N2], :(1.54e-8*((300 ./ Ti)^-0.5))],
             [[:N2Hpl, :NO], [:HNOpl, :N2], :(3.4e-10)],
             [[:N2Hpl, :O], [:OHpl, :N2], :(1.4e-10)],
             [[:N2Hpl, :OH], [:H2Opl, :N2], :(1.07e-8*((300 ./ Ti)^-0.5))],
             [[:N2Opl, :CO], [:CO2pl, :N2], :(1.11e-10)],
             [[:N2Opl, :CO], [:NOpl, :NCO], :(1.89e-10)],
             [[:N2Opl, :H2], [:HN2Opl, :H], :(2.56e-10)],
             [[:N2Opl, :H2], [:N2Hpl, :OH], :(1.04e-10)],
             [[:N2Opl, :H2O], [:H2Opl, :N2O], :(3.27e-8*((300 ./ Ti)^-0.5))],
             [[:N2Opl, :H2O], [:HN2Opl, :OH], :(3.64e-9*((300 ./ Ti)^-0.5))],
             [[:N2Opl, :N2O], [:NOpl, :NO, :N2], :(1.2e-11)],
             [[:N2Opl, :NO], [:NOpl, :N2O], :(2.3e-10)],
             [[:N2Opl, :NO2], [:NO2pl, :N2O], :(2.21e-10)],
             [[:N2Opl, :NO2], [:NOpl, :N2, :O2], :(4.29e-10)],
             [[:N2Opl, :O2], [:NOpl, :NO2], :(4.59e-11)],
             [[:N2Opl, :O2], [:O2pl, :N2O], :(2.24e-10)],
             [[:N2pl, :Ar], [:Arpl, :N2], :(2.0e-13)],
             [[:N2pl, :C], [:Cpl, :N2], :(1.1e-10)],
             [[:N2pl, :CH], [:CHpl, :N2], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:N2pl, :CN], [:CNpl, :N2], :(1.73e-9*((300 ./ Ti)^-0.5))],
             [[:N2pl, :CO], [:COpl, :N2], :(7.3e-11)],
             [[:N2pl, :CO2], [:CO2pl, :N2], :(8.0e-10)],
             [[:N2pl, :H2], [:N2Hpl, :H], :(1.87e-9*exp.(-54.7 ./ Ti))],
             [[:N2pl, :H2O], [:H2Opl, :N2], :(1.9e-9)],
             [[:N2pl, :H2O], [:N2Hpl, :OH], :(5.04e-10)],
             [[:N2pl, :HCN], [:HCNpl, :N2], :(3.9e-10)],
             [[:N2pl, :HCO], [:HCOpl, :N2], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:N2pl, :HCO], [:N2Hpl, :CO], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:N2pl, :N], [:Npl, :N2], :(1.4e-11)],
             [[:N2pl, :N2O], [:N2Opl, :N2], :(6.0e-10)],
             [[:N2pl, :NH], [:NHpl, :N2], :(1.13e-8*((300 ./ Ti)^-0.5))],
             [[:N2pl, :NH2], [:NH2pl, :N2], :(1.54e-8*((300 ./ Ti)^-0.5))],
             [[:N2pl, :NO], [:NOpl, :N2], :(4.4e-10)],
             [[:N2pl, :O], [:NOpl, :N], :(1.3e-10)],
             [[:N2pl, :O], [:Opl, :N2], :(9.8e-12)],
             [[:N2pl, :O2], [:O2pl, :N2], :(5.0e-11)],
             [[:N2pl, :OH], [:OHpl, :N2], :(1.09e-8*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :CH], [:CH2pl, :NH], :(6.06e-9*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :CH], [:CHpl, :NH2], :(6.06e-9*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :CN], [:HCNHpl, :N], :(1.73e-9*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :H2], [:NH3pl, :H], :(1.95e-10)],
             [[:NH2pl, :H2O], [:H3Opl, :NH], :(2.73e-9)],
             [[:NH2pl, :H2O], [:NH3pl, :OH], :(8.7e-11)],
             [[:NH2pl, :H2O], [:NH4pl, :O], :(1.16e-10)],
             [[:NH2pl, :HCN], [:HCNHpl, :NH], :(2.08e-8*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :HCO], [:H2COpl, :NH], :(7.45e-9*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :HCO], [:HCOpl, :NH2], :(7.45e-9*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :N], [:N2Hpl, :H], :(9.1e-11)],
             [[:NH2pl, :NH], [:NH3pl, :N], :(1.26e-8*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :NH2], [:NH3pl, :NH], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:NH2pl, :NO], [:NOpl, :NH2], :(7.0e-10)],
             [[:NH2pl, :O], [:HNOpl, :H], :(7.2e-11)],
             [[:NH2pl, :O2], [:H2NOpl, :O], :(1.19e-10)],
             [[:NH2pl, :O2], [:HNOpl, :OH], :(2.1e-11)],
             [[:NH3pl, :CH], [:NH4pl, :C], :(1.2e-8*((300 ./ Ti)^-0.5))],
             [[:NH3pl, :H2], [:NH4pl, :H], :(4.4e-13)],
             [[:NH3pl, :H2O], [:NH4pl, :OH], :(2.5e-10)],
             [[:NH3pl, :HCO], [:NH4pl, :CO], :(7.27e-9*((300 ./ Ti)^-0.5))],
             [[:NH3pl, :NO], [:NOpl, :NH3], :(7.2e-10)],
             [[:NHpl, :C], [:CHpl, :N], :(1.6e-9)],
             [[:NHpl, :CH], [:CH2pl, :N], :(1.71e-8*((300 ./ Ti)^-0.5))],
             [[:NHpl, :CN], [:HCNpl, :N], :(2.77e-8*((300 ./ Ti)^-0.5))],
             [[:NHpl, :CO], [:HCOpl, :N], :(4.41e-10)],
             [[:NHpl, :CO], [:OCNpl, :H], :(5.39e-10)],
             [[:NHpl, :CO2], [:HCO2pl, :N], :(3.85e-10)],
             [[:NHpl, :CO2], [:HNOpl, :CO], :(3.85e-10)],
             [[:NHpl, :CO2], [:NOpl, :HCO], :(3.3e-10)],
             [[:NHpl, :H2], [:H3pl, :N], :(1.85e-10)],
             [[:NHpl, :H2], [:NH2pl, :H], :(1.05e-9)],
             [[:NHpl, :H2O], [:H2Opl, :NH], :(1.05e-9)],
             [[:NHpl, :H2O], [:H3Opl, :N], :(1.05e-9)],
             [[:NHpl, :H2O], [:HNOpl, :H2], :(3.5e-10)],
             [[:NHpl, :H2O], [:NH2pl, :OH], :(8.75e-10)],
             [[:NHpl, :H2O], [:NH3pl, :O], :(1.75e-10)],
             [[:NHpl, :HCN], [:HCNHpl, :N], :(3.12e-8*((300 ./ Ti)^-0.5))],
             [[:NHpl, :HCO], [:H2COpl, :N], :(2.25e-8*((300 ./ Ti)^-0.5))],
             [[:NHpl, :N], [:N2pl, :H], :(1.3e-9)],
             [[:NHpl, :N2], [:N2Hpl, :N], :(6.5e-10)],
             [[:NHpl, :NH], [:NH2pl, :N], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:NHpl, :NH2], [:NH3pl, :N], :(2.6e-8*((300 ./ Ti)^-0.5))],
             [[:NHpl, :NO], [:N2Hpl, :O], :(1.78e-10)],
             [[:NHpl, :NO], [:NOpl, :NH], :(7.12e-10)],
             [[:NHpl, :O], [:OHpl, :N], :(1.0e-9)],
             [[:NHpl, :O2], [:HO2pl, :N], :(1.64e-10)],
             [[:NHpl, :O2], [:NOpl, :OH], :(2.05e-10)],
             [[:NHpl, :O2], [:O2pl, :NH], :(4.51e-10)],
             [[:NHpl, :OH], [:H2Opl, :N], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:NO2pl, :H], [:NOpl, :OH], :(1.9e-10)],
             [[:NO2pl, :H2], [:NOpl, :H2O], :(1.5e-10)],
             [[:NO2pl, :NO], [:NOpl, :NO2], :(2.75e-10)],
             [[:Npl, :CH], [:CNpl, :H], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:Npl, :CN], [:CNpl, :N], :(1.91e-8*((300 ./ Ti)^-0.5))],
             [[:Npl, :CO], [:COpl, :N], :(4.93e-10)],
             [[:Npl, :CO], [:Cpl, :NO], :(5.6e-12)],
             [[:Npl, :CO], [:NOpl, :C], :(6.16e-11)],
             [[:Npl, :CO2], [:CO2pl, :N], :(9.18e-10)],
             [[:Npl, :CO2], [:COpl, :NO], :(2.02e-10)],
             [[:Npl, :H2], [:NHpl, :H], :(5.0e-10*exp.(-85.0 ./ Ti))],
             [[:Npl, :H2O], [:H2Opl, :N], :(2.7e-9)],
             [[:Npl, :HCN], [:HCNpl, :N], :(3.7e-9)],
             [[:Npl, :HCO], [:HCOpl, :N], :(7.79e-9*((300 ./ Ti)^-0.5))],
             [[:Npl, :HCO], [:NHpl, :CO], :(7.79e-9*((300 ./ Ti)^-0.5))],
             [[:Npl, :N], [:N2pl], :(9.44e-19*((300 ./ Ti)^0.24)*exp.(-26.1 ./ Ti))],
             [[:Npl, :N2O], [:NOpl, :N2], :(5.5e-10)],
             [[:Npl, :NH], [:N2pl, :H], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:Npl, :NH], [:NHpl, :N], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:Npl, :NH2], [:NH2pl, :N], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:Npl, :NO], [:N2pl, :O], :(8.33e-11)],
             [[:Npl, :NO], [:NOpl, :N], :(4.72e-10)],
             [[:Npl, :O2], [:NOpl, :O], :(2.32e-10)],
             [[:Npl, :O2], [:O2pl, :N], :(3.07e-10)],
             [[:Npl, :O2], [:Opl, :NO], :(4.64e-11)],
             [[:Npl, :OH], [:OHpl, :N], :(6.41e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :C], [:COpl, :O], :(5.2e-11)],
             [[:O2pl, :C], [:Cpl, :O2], :(5.2e-11)],
             [[:O2pl, :CH], [:CHpl, :O2], :(5.37e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :CH], [:HCOpl, :O], :(5.37e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :HCO], [:HCOpl, :O2], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :HCO], [:HO2pl, :CO], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :N], [:NOpl, :O], :(1.0e-10)],
             [[:O2pl, :NH], [:HNOpl, :O], :(5.54e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :NH], [:NO2pl, :H], :(5.54e-9*((300 ./ Ti)^-0.5))],
             [[:O2pl, :NH2], [:NH2pl, :O2], :(1.51e-8*((300 ./ Ti)^-0.5))],
             [[:O2pl, :NO], [:NOpl, :O2], :(4.6e-10)],
             [[:O2pl, :NO2], [:NO2pl, :O2], :(6.6e-10)],
             [[:OHpl, :C], [:CHpl, :O], :(1.2e-9)],
             [[:OHpl, :CH], [:CH2pl, :O], :(6.06e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :CH], [:CHpl, :OH], :(6.06e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :CN], [:HCNpl, :O], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:OHpl, :CO], [:HCOpl, :O], :(8.4e-10)],
             [[:OHpl, :CO2], [:HCO2pl, :O], :(1.35e-9)],
             [[:OHpl, :H2], [:H2Opl, :H], :(9.7e-10)],
             [[:OHpl, :H2O], [:H2Opl, :OH], :(1.59e-9)],
             [[:OHpl, :H2O], [:H3Opl, :O], :(1.3e-9)],
             [[:OHpl, :HCN], [:HCNHpl, :O], :(2.08e-8*((300 ./ Ti)^-0.5))],
             [[:OHpl, :HCO], [:H2COpl, :O], :(4.85e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :HCO], [:H2Opl, :CO], :(4.85e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :HCO], [:HCOpl, :OH], :(4.85e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :N], [:NOpl, :H], :(8.9e-10)],
             [[:OHpl, :N2], [:N2Hpl, :O], :(2.4e-10)],
             [[:OHpl, :N2O], [:HN2Opl, :O], :(9.58e-10)],
             [[:OHpl, :N2O], [:N2Opl, :OH], :(2.13e-10)],
             [[:OHpl, :N2O], [:NOpl, :HNO], :(1.46e-10)],
             [[:OHpl, :NH], [:NH2pl, :O], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :NH2], [:NH2pl, :OH], :(8.66e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :NH2], [:NH3pl, :O], :(8.66e-9*((300 ./ Ti)^-0.5))],
             [[:OHpl, :NO], [:HNOpl, :O], :(6.11e-10)],
             [[:OHpl, :NO], [:NOpl, :OH], :(8.15e-10)],
             [[:OHpl, :O], [:O2pl, :H], :(7.1e-10)],
             [[:OHpl, :O2], [:O2pl, :OH], :(3.8e-10)],
             [[:OHpl, :OH], [:H2Opl, :O], :(1.21e-8*((300 ./ Ti)^-0.5))],
             [[:Opl, :CH], [:CHpl, :O], :(6.06e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :CH], [:COpl, :H], :(6.06e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :CN], [:NOpl, :C], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:Opl, :CO2], [:O2pl, :CO], :(1.1e-9)],   # Nair minimal ionosphere
             [[:Opl, :H], [:Hpl, :O], :(6.4e-10)],
             [[:Opl, :H2], [:OHpl, :H], :(1.62e-9)],
             [[:Opl, :H2O], [:H2Opl, :O], :(2.6e-9)],
             [[:Opl, :HCN], [:COpl, :NH], :(1.17e-9)],
             [[:Opl, :HCN], [:HCOpl, :N], :(1.17e-9)],
             [[:Opl, :HCN], [:NOpl, :CH], :(1.17e-9)],
             [[:Opl, :HCO], [:HCOpl, :O], :(7.45e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :HCO], [:OHpl, :CO], :(7.45e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :N2], [:NOpl, :N], :(4.58e-9*((300 ./ Ti)^-1.37)*exp.(-28.592 ./ Ti))],
             [[:Opl, :N2O], [:N2Opl, :O], :(6.3e-10)],
             [[:Opl, :NH], [:NHpl, :O], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :NH], [:NOpl, :H], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :NH2], [:NH2pl, :O], :(1.73e-8*((300 ./ Ti)^-0.5))],
             [[:Opl, :NO], [:NOpl, :O], :(8.0e-13)],
             [[:Opl, :NO2], [:NO2pl, :O], :(1.6e-9)],
             [[:Opl, :NO2], [:NOpl, :O2], :(8.3e-10)],
             [[:Opl, :O2], [:O2pl, :O], :(2.1e-11)],
             [[:Opl, :OH], [:O2pl, :H], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:Opl, :OH], [:OHpl, :O], :(6.24e-9*((300 ./ Ti)^-0.5))],
             [[:ArHpl, :E], [:Ar, :H], :(1.0e-9)],
             [[:Arpl, :E], [:Ar], :(4.0e-12*((300 ./ Te)^0.6))],
             [[:CHpl, :E], [:C, :H], :(1.65e-6*((300 ./ Te)^-0.42))],
             [[:CNpl, :E], [:N, :C], :(3.12e-6*((300 ./ Te)^-0.5))],
             [[:CO2pl, :E], [:CO, :O], :(3.03e-5*((300 ./ Te)^-0.75))],   # Nair minimal ionosphere
             [[:COpl, :E], [:O, :C], :(4.82e-6*((300 ./ Te)^-0.55))],
             [[:COpl, :E], [:O1D, :C], :(2.48e-8*((300 ./ Te)^-0.55))],
             [[:Cpl, :E], [:C], :(6.28e-10*((300 ./ Te)^-0.59))],
             [[:H2Opl, :E], [:O, :H, :H], :(2.08e-5*((300 ./ Te)^-0.74))],
             [[:H2Opl, :E], [:H2, :O], :(2.64e-6*((300 ./ Te)^-0.74))],
             [[:H2Opl, :E], [:OH, :H], :(5.86e-6*((300 ./ Te)^-0.74))],
             [[:H2pl, :E], [:H, :H], :(1.86e-7*((300 ./ Te)^-0.43))],
             [[:H3Opl, :E], [:H2, :O, :H], :(9.68e-8*((300 ./ Te)^-0.5))],
             [[:H3Opl, :E], [:H2O, :H], :(1.86e-6*((300 ./ Te)^-0.5))],
             [[:H3Opl, :E], [:OH, :H, :H], :(4.47e-6*((300 ./ Te)^-0.5))],
             [[:H3Opl, :E], [:OH, :H2], :(1.04e-6*((300 ./ Te)^-0.5))],
             [[:H3pl, :E], [:H, :H, :H], :(8.46e-7*((300 ./ Te)^-0.52))],
             [[:H3pl, :E], [:H2, :H], :(4.54e-7*((300 ./ Te)^-0.52))],
             [[:HCNHpl, :E], [:CN, :H, :H], :(3.79e-6*((300 ./ Te)^-0.65))],
             [[:HCNpl, :E], [:CN, :H], :(3.46e-6*((300 ./ Te)^-0.5))],
             [[:HCO2pl, :E], [:CO, :O, :H], :(2.18e-7)],
             [[:HCO2pl, :E], [:CO, :OH], :(9.18e-8)],
             [[:HCO2pl, :E], [:CO2, :H], :(1.7e-8)],  # Nair minimal ionosphere
             [[:HCOpl, :E], [:CH, :O], :(1.15e-7*((300 ./ Te)^-0.64))],
             [[:HCOpl, :E], [:CO, :H], :(1.06e-5*((300 ./ Te)^-0.64))],
             [[:HCOpl, :E], [:OH, :C], :(8.08e-7*((300 ./ Te)^-0.64))],
             [[:HN2Opl, :E], [:N2, :O, :H], :(3.81e-5*((300 ./ Te)^-0.74))],
             [[:HN2Opl, :E], [:N2, :OH], :(4.38e-5*((300 ./ Te)^-0.74))],
             [[:HNOpl, :E], [:NO, :H], :(5.2e-6*((300 ./ Te)^-0.5))],
             [[:HO2pl, :E], [:O2, :H], :(5.2e-6*((300 ./ Te)^-0.5))],
             [[:HOCpl, :E], [:CH, :O], :(1.7e-9*((300 ./ Te)^1.2))],
             [[:HOCpl, :E], [:CO, :H], :(3.3e-5*((300 ./ Te)^-1.0))],
             [[:HOCpl, :E], [:OH, :C], :(1.19e-8*((300 ./ Te)^1.2))],
             [[:Hpl, :E], [:H], :(6.46e-14*((300 ./ Te)^0.7))],
             [[:N2Hpl, :E], [:N2, :H], :(6.6e-7*((300 ./ Te)^-0.51))],
             [[:N2Hpl, :E], [:NH, :N], :(1.17e-6*((300 ./ Te)^-0.51))],
             [[:N2Opl, :E], [:N, :N, :O], :(1.36e-6*((300 ./ Te)^-0.57))],
             [[:N2Opl, :E], [:N2, :O], :(4.09e-6*((300 ./ Te)^-0.57))],
             [[:N2Opl, :E], [:NO, :N], :(3.07e-6*((300 ./ Te)^-0.57))],
             [[:N2pl, :E], [:N, :N], :(5.09e-7*((300 ./ Te)^-0.39))],
             [[:N2pl, :E], [:N2D, :N2D], :(1.42e-6*((300 ./ Te)^-0.39))],
             [[:NH2pl, :E], [:N, :H, :H], :(1.71e-5*((300 ./ Te)^-0.8)*exp.(-17.1 ./ Te))],
             [[:NH2pl, :E], [:NH, :H], :(8.34e-6*((300 ./ Te)^-0.79)*exp.(-17.1 ./ Te))],
             [[:NH3pl, :E], [:NH, :H, :H], :(2.68e-6*((300 ./ Te)^-0.5))],
             [[:NH3pl, :E], [:NH2, :H], :(2.68e-6*((300 ./ Te)^-0.5))],
             [[:NHpl, :E], [:N, :H], :(7.45e-7*((300 ./ Te)^-0.5))],
             [[:NO2pl, :E], [:NO, :O], :(5.2e-6*((300 ./ Te)^-0.5))],
             [[:NOpl, :E], [:O, :N], :(8.52e-7*((300 ./ Te)^-0.37))],
             [[:NOpl, :E], [:O, :N2D], :(4.53e-9*((300 ./ Te)^0.75))],
             [[:Npl, :E], [:N], :(1.9e-10*((300 ./ Te)^-0.7))],
             [[:O2pl, :E], [:O, :O], :(8.15e-6*((300 ./ Te)^-0.65))],   # Nair minimal ionosphere
             [[:OHpl, :E], [:O, :H], :(6.5e-7*((300 ./ Te)^-0.5))],
             [[:Opl, :E], [:O], :(1.4e-10*((300 ./ Te)^-0.66))]  
             ];
