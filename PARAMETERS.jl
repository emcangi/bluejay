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

research_dir = "/home/emc/GDrive-CU/Research-FF/"
results_dir = research_dir*"Results/"
main_cases_dir = "MainCases/"
det_cases_dir = "DetailedCases/"

# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Altitude grid discretization =================================================
const alt = convert(Array, (0:2e5:250e5))
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

global lowTs = 160.0
global hiTs = 270.0
global lowTt = 100.0
global hiTt = 160.0
global lowTe = 150.0
global hiTe = 250.0

MR_mean_water = 1.38e-4

# Lists ========================================================================
const fullspecieslist = [# neutrals
                         :Ar, :C, :CH, :CN, :CN2, :CO, :CO2, :D, :DO2, :DOCO, :H, :H2, :H2O, :H2O2, :HCN, :HCO, :HD, :HDO, 
                         :HDO2, :HNO, :HO2, :HOCO, :He, :N, :N2, :N2D, :N2O, :NCO, :NH, :NH2, :NH3, :NO, :NO2, :O, :O1D,
                         :O2, :O3, :OD, :OH,
                         # ions
                         :ArHpl, :Arpl, :C2N2pl, :C2Npl, :C2Opl, :C2pl, :CH2pl, 
                         :CHpl, :CNpl, :CO2pl, :COpl, :Cpl, :H2CNpl, 
                         :H2COpl, :H2NOpl, :H2Opl, :H2pl, :H3Opl, :H3pl, :HC2Npl, 
                         :HCNHpl, :HCNpl, :HCO2pl, :HCOOH2pl, :HCOpl, :HN2Opl, 
                         :HNCOpl, :HNCpl, :HNOpl, :HO2pl, :HOCpl, :Hpl, :HeHpl, 
                         :Hepl, :N2Hpl, :N2Opl, :N2pl, :NH2pl, :NH3pl, :NH4pl, 
                         :NHpl, :NO2pl, :NOpl, :Npl, :O2Dpl, :O2pl, :OCNpl, :OHpl, 
                         :Opl,
                         # electrons
                         :E];

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
                 # new March 2018
                 :JHDOtoHDpO1D, :JHDOtoHpDpO, :JODtoOpD, :JHDtoHpD, :JDO2toODpO,
                 :JHDO2toDO2pH, :JHDO2toHO2pD, :JHDO2toHDOpO1D, :JODtoO1DpD];

const nochemspecies = [:H2O, :HDO];
const chemspecies = setdiff(specieslist, nochemspecies);
const notransportspecies = [:H2O, :HDO];
const transportspecies = setdiff(specieslist, notransportspecies);
const speciesmolmasslist = Dict(# neutrals
                                :Ar=>40, :C=>12, :CH=>13, :CN=>26, :CN2=>40, :CO=>28, :CO2=>44, :D=>2, :DO2=>34, 
                                :DOCO=>46, :H=>1, :H2=>2, :H2O=>18, :H2O2=>34, :HCN=>27, :HCO=>29, :HD=>3, :HDO=>19, 
                                :HDO2=>35, :HNO=>31, :HO2=>33, :HOCO=>45, :He=>4, :N=>14, :N2=>28, :N2D, :N2O=>44, 
                                :NCO=>42, :NH=>15, :NH2=>16, :NH3=>17, :NO=>30, :NO2=>46, :O=>16, :O1D=>16, :O2=>32, 
                                :O3=>48, :OD=>18, :OH=>17,
                                # ions
                                :ArHpl=>41, :Arpl=>40, :C2N2pl=>52, :C2Npl=>38, :C2Opl=>40, :C2pl=>24, :CH2pl=>14, 
                                :CHpl=>13, :CNpl=>26, :CO2pl=>44, :COpl=>28, :Cpl=>12, :H2CNpl=>28, 
                                :H2COpl=>30, :H2NOpl=>32, :H2Opl=>18, :H2pl=>2, :H3Opl=>19, :H3pl=>3, :HC2Npl=>39, 
                                :HCNHpl=>28, :HCNpl=>27, :HOCOpl=>45, :HCOOH2pl=>47, :HCOpl=>29, :HN2Opl=>45, 
                                :HNCOpl=>43, :HNCpl=>27, :HNOpl=>31, :HO2pl=>33, :HOCpl=>29, :Hpl=>1, :HeHpl=>5, 
                                :Hepl=>4, :N2Hpl=>29, :N2Opl=>44, :N2pl=>28, :NH2pl=>16, :NH3pl=>17, :NH4pl=>18, 
                                :NHpl=>15, :NO2pl=>46, :NOpl=>30, :Npl=>14, :O2Dpl, :O2pl=>32, :OCNpl=>42, :OHpl=>17, 
                                :Opl=>16,
                                # electrons
                                :E);



# Plotty plot plot stuff =======================================================
speciescolor = Dict( # H group
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
                :CO2 => "#d18564",   # dark peach
                :CO2pl => "#614215", # brown
                :CO => "#ff6600",    # orange
                :HOCO => "#e8ba8c", :DOCO => "#e8ba8c",  #tannish

                # Nitrogen, argon, Helium
                :Ar=>"#808080", :N2=>"#cccccc", :He=>"#c17388", # grays
                :N=>"#6e748a", :N2D=>"#809d03", :N2O=>"#be8b65", :NCO=>"#633339", :NH=>"#cacdda", :NH2=>"#6ceb83", :NH3=>"#224069", :NO=>"#a27fff", :NO2=>"#fe03cb", 
                :HNO=>"#76bcfd",

                # creepy chemistry (mostly organic stuff)
                :C=>"#d9c382", :CH=>"#cea3ce", :CN=>"#6d5000", :CN2=>"#006974", :HCN=>"#479f5e", :HCO=>"#94c6bf", 

                # ions
                :ArHpl=>"#d70000", :Arpl=>"#8c3cff", :C2N2pl=>"#028800", :C2Npl=>"#00acc7", :C2Opl=>"#98ff00", :C2pl=>"#ff7fd1", :CH2pl=>"#6c004f", 
                :CHpl=>"#ffa530", :CNpl=>"#00009d", :CO2pl=>"#867068", :COpl=>"#004942", :Cpl=>"#4f2a00", :H2CNpl=>"#00fdcf", 
                :H2COpl=>"#bcb7ff", :H2NOpl=>"#95b47a", :H2Opl=>"#c004b9", :H2pl=>"#2566a2", :H3Opl=>"#280041", :H3pl=>"#dcb3af", :HC2Npl=>"#fef590", 
                :HCNHpl=>"#50455b", :HCNpl=>"#a47c00", :HOCOpl=>"#ff7166", :HCOOH2pl=>"#3f816e", :HCOpl=>"#82000d", :HN2Opl=>"#a37bb3", 
                :HNCOpl=>"#344e00", :HNCpl=>"#9be4ff", :HNOpl=>"#eb0077", :HO2pl=>"#2d000a", :HOCpl=>"#5e90ff", :Hpl=>"#00c720", :HeHpl=>"#5801aa", 
                :Hepl=>"#001e00", :N2Hpl=>"#9a4700", :N2Opl=>"#969fa6", :N2pl=>"#9b425c", :NH2pl=>"#001f32", :NH3pl=>"#c8c400", :NH4pl=>"#ffd0ff", 
                :NHpl=>"#00be9a", :NO2pl=>"#3715ff", :NOpl=>"#2d2525", :Npl=>"#df58ff", :O2Dpl=>"#bee7c0", :O2pl=>"#7f4598", :OCNpl=>"#524f3c", :OHpl=>"#d86600", 
                :Opl=>"#647438",
                );


speciesstyle = Dict( # H group
                :H => "-", :D => "--",
                :H2 => "-",  :HD => "--",
                # hydroxides
                :OH => "-", :OD => "--",
                # "water group" (roughly, I ain't a chemist)
                :H2O => "-", :HDO => "--",
                :H2O2 => "-", :HDO2 => "--",
                :HO2 => "-", :DO2 => "--",
                # O group
                :O1D => "-", :O => "-", :O2 => "-", :O3 => "-",
                # CO group
                :CO => "-", :CO2 => "-", :CO2pl => "-",
                :HOCO => "-", :DOCO => "--",
                # nonreactants
                :Ar => "-", :N2 => "-",);
                
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
             [[:HO2], [:OH, :O], :JHO2toOHpO], # other branches should be here, but
                                               # have not been measured
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

             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8*3.0e-33*(300 ./ T)^3.25)],
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35*exp(724 ./ T))],
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5*6.0e-34*(300 ./ T)^2.4)],
             [[:O, :O3], [:O2, :O2], :(8.0e-12*exp(-2060 ./ T))],  # Sander 2011
             [[:O, :CO, :M], [:CO2, :M], :(2.2e-33*exp(-1780 ./ T))],

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.2e-11*exp(70 ./ T))], # verified NIST 4/3/18
             [[:O1D, :O3], [:O2, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :O3], [:O, :O, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11*exp(115 ./ T))], # Sander2011. NIST: 7.41e-11*exp(120/T)
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Sander2011. Yung89: 1e-10; NIST 1.1e-10
             [[:O1D, :HD], [:H, :OD], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10*exp(60 ./ T))], # Sander2011. Yung89: 2.2e-10; NIST: 1.62e-10*exp(65/T)
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp(60 ./ T))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12*exp(-4000 ./ T))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12*exp(-4390 ./ T))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12*exp(-4400 ./ T))], # NIST
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11*exp(-4038 ./ T))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11*exp(-3821 ./ T))], # NIST (1986, 200-300K): 8.19e-13*exp(-2700/T)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12*exp(-1800 ./ T))], # Sander2011. Yung89: 5.5e-12*exp(-2000/T). KIDA: 7.7E-12*exp(-2100/T). old rate from Mike: 9.01e-13*exp(-1526/T)
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # Yung88: rate (3/20)*H-ana. Sander2011: 5e-12*exp(-2130 ./ T)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12*exp(-1800 ./ T))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:OD, :HD], [:D2O, :H], :(???)],  # can't find a rate...?

             # recombination of H. Use EITHER the first line OR the 2nd and 3rd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32*(298 ./ T)^2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32*(298 ./ T)^2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32*(298 ./ T)^2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # Can't find in databases. Mike's rate.
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # not in Yung88. assumed rate

             ## H + HO2
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Sander2011. Indep of T for 245<T<300
             [[:H, :HO2], [:H2, :O2], :(0.5*6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # O1D is theoretically mandated
             [[:H, :DO2], [:OH, :OD], :(7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HD, :O2], :(0.5*6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HDO, :O1D], :(1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:D, :HO2], [:OH, :OD], :(0.71*7.2e-11)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:D, :HO2], [:HD, :O2], :(0.71*0.5*6.9e-12)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:D, :HO2], [:HDO, :O1D], :(0.71*1.6e-12)], # Yung88: rate 0.71*H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HO2, :D], :(1e-10/(0.54*exp(890 ./ T)))], # Yung88 (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12*exp(-1890 ./ T))], # verified NIST 4/3/18. Only valid for T>300K. No experiment for lower.
             # [[:H, :HDO2], [:DO2, :H2], :(0)], # Cazaux2010: branching ratio = 0
             # [[:H, :HDO2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:DO2, :H2], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11*exp(-1800 ./ T))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO, :OH], :(0.5*1.16e-11*exp(-2110 ./ T))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<T<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O, :OD], :(0.5*1.16e-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:HDO, :OH], :(0.5*1.16e-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:H2O, :OD], :(0.5*1.16e-11*exp(-2110 ./ T))], # see previous line
             [[:D, :HDO2], [:OD, :HDO], :(0.5*1.16e-11*exp(-2110 ./ T))], # added 4/3 with assumed rate from other rxns
             [[:D, :HDO2], [:OH, :D2O], :(0.5*1.16e-11*exp(-2110/T))], # sourced from Cazaux et al

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0*4.4e-32*(T/300.)^-1.3), # Sander2011, 300K+. Yung89: 5.5e-32(T/300)^-1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11*(T/300.)^0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0*4.4e-32*(T/300.)^-1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11*(T/300.)^0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10*exp(-470 ./ T))], # verified Yung89, NIST 4/3/18
             [[:D, :O3], [:OD, :O2], :(0.71*1.4e-10*exp(-470 ./ T))], # Yung88: rate 0.71*H-ana (assumed). verified Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11*exp(180 ./ T))], # Sander2011. KIDA+NIST 4/3/18 150-500K: 2.4e-11*exp(110 ./ T). Yung89: 2.2e-11*exp(120/T) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11*exp(180 ./ T))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11*exp(200 ./ T))], # Sander2011. KIDA (220-400K): 2.7e-11*exp(224/T)
             [[:O, :DO2], [:OD, :O2], :(3.0e-11*exp(200 ./ T))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12*exp(-2000 ./ T))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5*1.4e-12*exp(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5*1.4e-12*exp(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(4.2e-12*exp(-240 ./ T))], # NIST+KIDA, 200-350K: 6.2e-14*(T/300)^2.62*exp(945 ./ T) changed 4/3/18. Yung89: 4.2e-12*exp(-240/T). old rate w/mystery origin: 1.8e-12.
             [[:OD, :OH], [:HDO, :O], :(4.2e-12*exp(-240 ./ T))], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3*6.9e-31*(T/300.)^-1.0),:(2.6e-11))], # Sander2011. Why 1.3?
             [[:OD, :OH], [:HDO2], threebody(:(1.3*6.9e-31*(T/300.)^-1.0),:(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12*exp(-940 ./ T))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12*exp(-940 ./ T))], # Yung88: rate same as H-ana
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11*exp(250 ./ T))], # verified NIST 4/3/18. Yung89: 4.6e-11*exp(230/T) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ T))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ T))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12*exp(-160 ./ T))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12*exp(-200/T). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12*exp(-160 ./ T))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5*2.9e-12*exp(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5*2.9e-12*exp(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14*exp(-490 ./ T))], # Sander2011. Yung89: 1.1e-14*exp(-500/T). KIDA 250-340K: 2.03e-16*(T/300)^4.57*exp(693/T). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(1.0e-14*exp(-490 ./ T))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13*exp(460 ./ T))], # Sander2011. Yung89: 2.3e-13*exp(600/T). KIDA 230-420K: 2.2e-13*exp(600/T)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13*exp(460 ./ T))], # Yung88: same as H-ana (assumed)
             # *** why do we have he next two reactions? I forgot...
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2*2.1e-33*exp(920 ./ T))], # Sander2011.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2*2.1e-33*exp(920 ./ T))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9*(T^-0.63)/(0.72*exp(717 ./ T)))], # rate: Yung88. NIST (Howard82): 5.25E-11*(T/298)^-0.63  - turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9*T^-0.63)], # Yung88  - turn off for Case 2

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13*(T/300.)^0.6),:(2.1e9*(T/300.)^6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13*(T/300.)^0.6),:(2.1e9*(T/300.)^6.1))], # Yung88: same as H-ana.
             ### possible deuterated analogues below
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33*(T/300.)^-1.4),:(1.1e-12*(T/300.)^1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33*(T/300.)^-1.4),:(1.1e-12*(T/300.)^1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 / Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)],

             # IONOSPHERE
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
             [[:CHpl, :CN], [:C2Npl, :H], :(9.53e-9*((300 ./ T)^-0.5))],
             [[:CHpl, :CO], [:HCOpl, :C], :(7.0e-12*((300 ./ T)^-0.5))],
             [[:CHpl, :CO2], [:HCOpl, :CO], :(1.6e-9)],
             [[:CHpl, :H], [:Cpl, :H2], :(7.5e-10)],
             [[:CHpl, :H2], [:CH2pl, :H], :(1.2e-9)],
             [[:CHpl, :H2O], [:H2COpl, :H], :(1.0e-8*((300 ./ T)^-0.5))],
             [[:CHpl, :H2O], [:H3Opl, :C], :(1.45e-9)],
             [[:CHpl, :H2O], [:HCOpl, :H2], :(5.02e-8*((300 ./ T)^-0.5))],
             [[:CHpl, :HCN], [:C2Npl, :H2], :(4.2e-10)],
             [[:CHpl, :HCN], [:HC2Npl, :H], :(2.8e-10)],
             [[:CHpl, :HCN], [:HCNHpl, :C], :(2.1e-9)],
             [[:CHpl, :HCO], [:CH2pl, :CO], :(7.97e-9*((300 ./ T)^-0.5))],
             [[:CHpl, :HCO], [:HCOpl, :CH], :(7.97e-9*((300 ./ T)^-0.5))],
             [[:CHpl, :N], [:CNpl, :H], :(1.9e-10)],
             [[:CHpl, :NH], [:CNpl, :H2], :(1.32e-8*((300 ./ T)^-0.5))],
             [[:CHpl, :NO], [:NOpl, :CH], :(7.6e-10)],
             [[:CHpl, :O], [:COpl, :H], :(3.5e-10)],
             [[:CHpl, :O2], [:HCOpl, :O], :(8.73e-10)],
             [[:CHpl, :OH], [:COpl, :H2], :(1.3e-8*((300 ./ T)^-0.5))],
             [[:CNpl, :C], [:Cpl, :CN], :(1.1e-10)],
             [[:CNpl, :CH], [:CHpl, :CN], :(1.11e-8*((300 ./ T)^-0.5))],
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
             [[:CNpl, :HCO], [:HCNpl, :CO], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:CNpl, :HCO], [:HCOpl, :CN], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:CNpl, :N], [:N2pl, :C], :(6.1e-10)],
             [[:CNpl, :N2O], [:N2Opl, :CN], :(4.56e-10)],
             [[:CNpl, :N2O], [:NOpl, :CN2], :(1.52e-10)],
             [[:CNpl, :N2O], [:OCNpl, :N2], :(1.52e-10)],
             [[:CNpl, :NH], [:NHpl, :CN], :(1.13e-8*((300 ./ T)^-0.5))],
             [[:CNpl, :NH2], [:NH2pl, :CN], :(1.58e-8*((300 ./ T)^-0.5))],
             [[:CNpl, :NO], [:NOpl, :CN], :(5.7e-10)],
             [[:CNpl, :NO], [:OCNpl, :N], :(1.9e-10)],
             [[:CNpl, :O], [:Opl, :CN], :(6.5e-11)],
             [[:CNpl, :O2], [:NOpl, :CO], :(8.6e-11)],
             [[:CNpl, :O2], [:O2pl, :CN], :(2.58e-10)],
             [[:CNpl, :O2], [:OCNpl, :O], :(8.6e-11)],
             [[:CNpl, :OH], [:OHpl, :CN], :(1.11e-8*((300 ./ T)^-0.5))],
             [[:CO2pl, :H], [:HCOpl, :O], :(4.47e-10)],
             [[:CO2pl, :H], [:Hpl, :CO2], :(5.53e-11)],
             [[:CO2pl, :H2], [:HCO2pl, :H], :(2.24e-9*((300 ./ T)^-0.15))],
             [[:CO2pl, :H2O], [:H2Opl, :CO2], :(1.8e-9)],
             [[:CO2pl, :H2O], [:HCO2pl, :OH], :(6.0e-10)],
             [[:CO2pl, :HCN], [:HCNpl, :CO2], :(8.1e-10)],
             [[:CO2pl, :HCN], [:HCO2pl, :CN], :(9.0e-11)],
             [[:CO2pl, :N], [:COpl, :NO], :(3.4e-10)],
             [[:CO2pl, :NO], [:NOpl, :CO2], :(1.23e-10)],
             [[:CO2pl, :O], [:O2pl, :CO], :(1.6e-10)],
             [[:CO2pl, :O], [:Opl, :CO2], :(1.0e-10)],
             [[:CO2pl, :O2], [:O2pl, :CO2], :(5.5e-11)],
             [[:COpl, :C], [:Cpl, :CO], :(1.1e-10)],
             [[:COpl, :CH], [:CHpl, :CO], :(5.54e-9*((300 ./ T)^-0.5))],
             [[:COpl, :CH], [:HCOpl, :C], :(5.54e-9*((300 ./ T)^-0.5))],
             [[:COpl, :CO2], [:CO2pl, :CO], :(1.1e-9)],
             [[:COpl, :H], [:Hpl, :CO], :(4.0e-10)],
             [[:COpl, :H2], [:HCOpl, :H], :(7.5e-10)],
             [[:COpl, :H2], [:HOCpl, :H], :(7.5e-10)],
             [[:COpl, :H2O], [:H2Opl, :CO], :(1.56e-9)],
             [[:COpl, :H2O], [:HCOpl, :OH], :(8.4e-10)],
             [[:COpl, :HCN], [:HCNpl, :CO], :(3.06e-9)],
             [[:COpl, :HCO], [:HCOpl, :CO], :(1.28e-8*((300 ./ T)^-0.5))],
             [[:COpl, :N], [:NOpl, :C], :(8.2e-11)],
             [[:COpl, :NH], [:HCOpl, :N], :(5.54e-9*((300 ./ T)^-0.5))],
             [[:COpl, :NH], [:NHpl, :CO], :(5.54e-9*((300 ./ T)^-0.5))],
             [[:COpl, :NH2], [:HCOpl, :NH], :(7.79e-9*((300 ./ T)^-0.5))],
             [[:COpl, :NH2], [:NH2pl, :CO], :(7.79e-9*((300 ./ T)^-0.5))],
             [[:COpl, :NO], [:NOpl, :CO], :(4.2e-10)],
             [[:COpl, :O], [:Opl, :CO], :(1.4e-10)],
             [[:COpl, :O2], [:O2pl, :CO], :(1.5e-10)],
             [[:COpl, :OH], [:HCOpl, :O], :(5.37e-9*((300 ./ T)^-0.5))],
             [[:COpl, :OH], [:OHpl, :CO], :(5.37e-9*((300 ./ T)^-0.5))],
             [[:Cpl, :C], [:C2pl], :(1.52e-18*((300 ./ T)^0.17)*expl(-(-101.5) ./ T))],
             [[:Cpl, :CH], [:C2pl, :H], :(6.58e-9*((300 ./ T)^-0.5))],
             [[:Cpl, :CH], [:CHpl, :C], :(6.58e-9*((300 ./ T)^-0.5))],
             [[:Cpl, :CO2], [:CO2pl, :C], :(1.1e-10)],
             [[:Cpl, :CO2], [:COpl, :CO], :(9.9e-10)],
             [[:Cpl, :H], [:CHpl], :(1.7e-17)],
             [[:Cpl, :H2], [:CH2pl], :(3.32e-13*((300 ./ T)^-1.3)*expl(-(-23.0) ./ T))],
             [[:Cpl, :H2], [:CHpl, :H], :(7.4e-10*expl(-(-4537.0) ./ T))],
             [[:Cpl, :H2O], [:H2Opl, :C], :(2.4e-10)],
             [[:Cpl, :H2O], [:HCOpl, :H], :(1.56e-8*((300 ./ T)^-0.5))],
             [[:Cpl, :H2O], [:HOCpl, :H], :(2.16e-9)],
             [[:Cpl, :HCN], [:C2Npl, :H], :(2.95e-9)],
             [[:Cpl, :HCO], [:CHpl, :CO], :(8.31e-9*((300 ./ T)^-0.5))],
             [[:Cpl, :HCO], [:HCOpl, :C], :(8.31e-9*((300 ./ T)^-0.5))],
             [[:Cpl, :N], [:CNpl], :(7.24e-19*((300 ./ T)^0.07)*expl(-(-57.5) ./ T))],
             [[:Cpl, :N2O], [:NOpl, :CN], :(9.1e-10)],
             [[:Cpl, :NH], [:CNpl, :H], :(1.35e-8*((300 ./ T)^-0.5))],
             [[:Cpl, :NH2], [:HCNpl, :H], :(1.91e-8*((300 ./ T)^-0.5))],
             [[:Cpl, :NO], [:NOpl, :C], :(7.5e-10)],
             [[:Cpl, :O], [:COpl], :(7.39e-18*((300 ./ T)^-0.15)*expl(-(-68.0) ./ T))],
             [[:Cpl, :O2], [:COpl, :O], :(3.48e-10)],
             [[:Cpl, :O2], [:Opl, :CO], :(5.22e-10)],
             [[:Cpl, :OH], [:COpl, :H], :(1.33e-8*((300 ./ T)^-0.5))],
             [[:H2Opl, :C], [:CHpl, :OH], :(1.1e-9)],
             [[:H2Opl, :CH], [:CH2pl, :OH], :(5.89e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :CH], [:CHpl, :H2O], :(5.89e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :CO], [:HCOpl, :OH], :(4.25e-10)],
             [[:H2Opl, :H2], [:H3Opl, :H], :(7.6e-10)],
             [[:H2Opl, :H2O], [:H3Opl, :OH], :(1.85e-9)],
             [[:H2Opl, :HCN], [:HCNHpl, :OH], :(1.05e-9)],
             [[:H2Opl, :HCO], [:H2COpl, :OH], :(4.85e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :HCO], [:H3Opl, :CO], :(4.85e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :HCO], [:HCOpl, :H2O], :(4.85e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :N], [:HNOpl, :H], :(1.12e-10)],
             [[:H2Opl, :N], [:NOpl, :H2], :(2.8e-11)],
             [[:H2Opl, :NH], [:H3Opl, :N], :(1.23e-8*((300 ./ T)^-0.5))],
             [[:H2Opl, :NH2], [:NH2pl, :H2O], :(8.49e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :NH2], [:NH3pl, :OH], :(8.49e-9*((300 ./ T)^-0.5))],
             [[:H2Opl, :NO], [:NOpl, :H2O], :(4.6e-10)],
             [[:H2Opl, :NO2], [:NO2pl, :H2O], :(1.2e-9)],
             [[:H2Opl, :O], [:O2pl, :H2], :(4.0e-11)],
             [[:H2Opl, :O2], [:O2pl, :H2O], :(3.3e-10)],
             [[:H2Opl, :OH], [:H3Opl, :O], :(1.2e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :Ar], [:ArHpl, :H], :(2.1e-9)],
             [[:H2pl, :C], [:CHpl, :H], :(2.4e-9)],
             [[:H2pl, :CH], [:CH2pl, :H], :(1.23e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :CH], [:CHpl, :H2], :(1.23e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :CN], [:CNpl, :H2], :(2.08e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :CN], [:HCNpl, :H], :(2.08e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :CO], [:COpl, :H2], :(6.44e-10)],
             [[:H2pl, :CO], [:HCOpl, :H], :(2.9e-9)],
             [[:H2pl, :CO2], [:HCO2pl, :H], :(2.35e-9)],
             [[:H2pl, :H], [:Hpl, :H2], :(6.4e-10)],
             [[:H2pl, :H2], [:H3pl, :H], :(2.0e-9)],
             [[:H2pl, :H2O], [:H2Opl, :H2], :(3.87e-9)],
             [[:H2pl, :H2O], [:H3Opl, :H], :(3.43e-9)],
             [[:H2pl, :HCN], [:HCNpl, :H2], :(4.68e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :HCO], [:H3pl, :CO], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :HCO], [:HCOpl, :H2], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :He], [:HeHpl, :H], :(1.35e-10)],
             [[:H2pl, :N], [:NHpl, :H], :(1.9e-9)],
             [[:H2pl, :N2], [:N2Hpl, :H], :(2.0e-9)],
             [[:H2pl, :N2O], [:HN2Opl, :H], :(1.32e-9)],
             [[:H2pl, :N2O], [:N2Hpl, :OH], :(7.77e-10)],
             [[:H2pl, :NH], [:NH2pl, :H], :(1.32e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :NH], [:NHpl, :H2], :(1.32e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :NO], [:HNOpl, :H], :(1.1e-9)],
             [[:H2pl, :NO], [:NOpl, :H2], :(1.1e-9)],
             [[:H2pl, :O], [:OHpl, :H], :(1.5e-9)],
             [[:H2pl, :O2], [:HO2pl, :H], :(1.92e-9)],
             [[:H2pl, :O2], [:O2pl, :H2], :(7.83e-10)],
             [[:H2pl, :OH], [:H2Opl, :H], :(1.32e-8*((300 ./ T)^-0.5))],
             [[:H2pl, :OH], [:OHpl, :H2], :(1.32e-8*((300 ./ T)^-0.5))],
             [[:H3Opl, :C], [:HCOpl, :H2], :(1.0e-11)],
             [[:H3Opl, :CH], [:CH2pl, :H2O], :(1.18e-8*((300 ./ T)^-0.5))],
             [[:H3Opl, :HCN], [:HCNHpl, :H2O], :(3.8e-9)],
             [[:H3Opl, :NH2], [:NH3pl, :H2O], :(1.68e-8*((300 ./ T)^-0.5))],
             [[:H3pl, :Ar], [:ArHpl, :H2], :(3.65e-10)],
             [[:H3pl, :C], [:CHpl, :H2], :(2.0e-9)],
             [[:H3pl, :CH], [:CH2pl, :H2], :(2.08e-8*((300 ./ T)^-0.5))],
             [[:H3pl, :CN], [:HCNpl, :H2], :(3.46e-8*((300 ./ T)^-0.5))],
             [[:H3pl, :CO], [:HCOpl, :H2], :(3.06e-9*((300 ./ T)^-0.142)*expl(-(3.41) ./ T))],
             [[:H3pl, :CO], [:HOCpl, :H2], :(5.82e-10*((300 ./ T)^0.0661)*expl(-(-5.21) ./ T))],
             [[:H3pl, :CO2], [:HCO2pl, :H2], :(2.5e-9)],
             [[:H3pl, :H2O], [:H3Opl, :H2], :(5.3e-9)],
             [[:H3pl, :HCN], [:HCNHpl, :H2], :(7.5e-9)],
             [[:H3pl, :HCO], [:H2COpl, :H2], :(2.94e-8*((300 ./ T)^-0.5))],
             [[:H3pl, :N], [:NH2pl, :H], :(3.9e-10)],
             [[:H3pl, :N], [:NHpl, :H2], :(2.6e-10)],
             [[:H3pl, :N2], [:N2Hpl, :H2], :(1.63e-9)],
             [[:H3pl, :N2O], [:HN2Opl, :H2], :(2.5e-9)],
             [[:H3pl, :NH], [:NH2pl, :H2], :(2.25e-8*((300 ./ T)^-0.5))],
             [[:H3pl, :NO], [:HNOpl, :H2], :(1.94e-9)],
             [[:H3pl, :NO2], [:NO2pl, :H2, :H], :(7.0e-12)],
             [[:H3pl, :NO2], [:NOpl, :OH, :H2], :(6.93e-10)],
             [[:H3pl, :O], [:H2Opl, :H], :(8.33e-10*((300 ./ T)^-0.156)*expl(-(-1.4) ./ T))],
             [[:H3pl, :O], [:OHpl, :H2], :(1.94e-9*((300 ./ T)^-0.156)*expl(-(-1.4) ./ T))],
             [[:H3pl, :O2], [:HO2pl, :H2], :(6.7e-10)],
             [[:H3pl, :OH], [:H2Opl, :H2], :(2.25e-8*((300 ./ T)^-0.5))],
             [[:HCNHpl, :CH], [:CH2pl, :HCN], :(5.46e-9*((300 ./ T)^-0.5))],
             [[:HCNHpl, :H2O], [:H3Opl, :HCN], :(8.8e-13)],
             [[:HCNpl, :C], [:CHpl, :CN], :(1.1e-9)],
             [[:HCNpl, :CH], [:CH2pl, :CN], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:HCNpl, :CO], [:HCOpl, :CN], :(1.38e-10)],
             [[:HCNpl, :CO], [:HNCpl, :CO], :(3.22e-10)],
             [[:HCNpl, :CO2], [:HCO2pl, :CN], :(2.1e-10)],
             [[:HCNpl, :CO2], [:HNCpl, :CO2], :(2.9e-10)],
             [[:HCNpl, :H], [:Hpl, :HCN], :(3.7e-11)],
             [[:HCNpl, :H2], [:HCNHpl, :H], :(8.8e-10)],
             [[:HCNpl, :H2O], [:H2Opl, :HCN], :(3.12e-8*((300 ./ T)^-0.5))],
             [[:HCNpl, :H2O], [:H3Opl, :CN], :(3.12e-8*((300 ./ T)^-0.5))],
             [[:HCNpl, :HCN], [:HCNHpl, :CN], :(1.45e-9)],
             [[:HCNpl, :HCO], [:H2COpl, :CN], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:HCNpl, :HCO], [:HCNHpl, :CO], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:HCNpl, :N], [:CHpl, :N2], :(2.2e-10)],
             [[:HCNpl, :N2O], [:N2Opl, :HCN], :(1.08e-9)],
             [[:HCNpl, :NH], [:NH2pl, :CN], :(1.13e-8*((300 ./ T)^-0.5))],
             [[:HCNpl, :NH2], [:NH3pl, :CN], :(1.56e-8*((300 ./ T)^-0.5))],
             [[:HCNpl, :NO], [:NOpl, :HCN], :(8.1e-10)],
             [[:HCNpl, :O2], [:O2pl, :HCN], :(4.1e-10)],
             [[:HCNpl, :OH], [:H2Opl, :CN], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:HCO2pl, :C], [:CHpl, :CO2], :(1.0e-9)],
             [[:HCO2pl, :CO], [:HCOpl, :CO2], :(7.8e-10)],
             [[:HCO2pl, :H2O], [:H3Opl, :CO2], :(2.65e-9)],
             [[:HCO2pl, :O], [:HCOpl, :O2], :(5.8e-10)],
             [[:HCOpl, :C], [:CHpl, :CO], :(1.1e-9)],
             [[:HCOpl, :CH], [:CH2pl, :CO], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:HCOpl, :H2O], [:H3Opl, :CO], :(2.6e-9)],
             [[:HCOpl, :H2O], [:HCOOH2pl], :(6.64e-10*((300 ./ T)^-1.3))],
             [[:HCOpl, :HCN], [:HCNHpl, :CO], :(3.5e-9)],
             [[:HCOpl, :HCO], [:H2COpl, :CO], :(1.26e-8*((300 ./ T)^-0.5))],
             [[:HCOpl, :N2O], [:HN2Opl, :CO], :(3.3e-12)],
             [[:HCOpl, :NH], [:NH2pl, :CO], :(1.11e-8*((300 ./ T)^-0.5))],
             [[:HCOpl, :NH2], [:NH3pl, :CO], :(1.54e-8*((300 ./ T)^-0.5))],
             [[:HCOpl, :OH], [:H2Opl, :CO], :(1.07e-8*((300 ./ T)^-0.5))],
             [[:HCOpl, :OH], [:HCO2pl, :H], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :C], [:Cpl, :He], :(8.74e-17*((300 ./ T)^0.75))],
             [[:Hepl, :CH], [:CHpl, :He], :(8.66e-9*((300 ./ T)^-0.5))],
             [[:Hepl, :CH], [:Cpl, :He, :H], :(1.91e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :CN], [:Cpl, :N, :He], :(1.52e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :CN], [:Npl, :C, :He], :(1.52e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :CO], [:Cpl, :O, :He], :(1.6e-9)],
             [[:Hepl, :CO2], [:CO2pl, :He], :(5.0e-11)],
             [[:Hepl, :CO2], [:COpl, :O, :He], :(7.8e-10)],
             [[:Hepl, :CO2], [:Cpl, :O2, :He], :(2.0e-11)],
             [[:Hepl, :CO2], [:O2pl, :C, :He], :(1.1e-11)],
             [[:Hepl, :CO2], [:Opl, :CO, :He], :(1.4e-10)],
             [[:Hepl, :H], [:HeHpl], :(3.43e-15*((300 ./ T)^-0.37))],
             [[:Hepl, :H], [:Hpl, :He], :(2.88e-16*((300 ./ T)^0.25))],
             [[:Hepl, :H2], [:H2pl, :He], :(1.7e-14)],
             [[:Hepl, :H2], [:Hpl, :He, :H], :(8.3e-14)],
             [[:Hepl, :H2O], [:H2Opl, :He], :(5.5e-11*((300 ./ T)^-0.5))],
             [[:Hepl, :H2O], [:Hpl, :OH, :He], :(1.85e-10)],
             [[:Hepl, :H2O], [:OHpl, :He, :H], :(2.6e-10*((300 ./ T)^-0.5))],
             [[:Hepl, :HCN], [:CHpl, :N, :He], :(6.93e-10)],
             [[:Hepl, :HCN], [:CNpl, :He, :H], :(1.55e-9)],
             [[:Hepl, :HCN], [:Cpl, :NH, :He], :(8.25e-10)],
             [[:Hepl, :HCN], [:Npl, :CH, :He], :(2.31e-10)],
             [[:Hepl, :HCO], [:CHpl, :O, :He], :(8.49e-9*((300 ./ T)^-0.5))],
             [[:Hepl, :HCO], [:COpl, :He, :H], :(8.49e-9*((300 ./ T)^-0.5))],
             [[:Hepl, :HCO], [:HeHpl, :CO], :(5.2e-9*((300 ./ T)^-0.5))],
             [[:Hepl, :HNO], [:Hpl, :NO, :He], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :HNO], [:NOpl, :He, :H], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :N2], [:N2pl, :He], :(5.2e-10)],
             [[:Hepl, :N2], [:Npl, :N, :He], :(7.8e-10)],
             [[:Hepl, :N2O], [:N2pl, :O, :He], :(1.24e-9)],
             [[:Hepl, :N2O], [:NOpl, :N, :He], :(4.83e-10)],
             [[:Hepl, :N2O], [:Npl, :NO, :He], :(2.99e-10)],
             [[:Hepl, :N2O], [:Opl, :N2, :He], :(2.76e-10)],
             [[:Hepl, :NH], [:Npl, :He, :H], :(1.91e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :NH2], [:NHpl, :He, :H], :(1.39e-8*((300 ./ T)^-0.5))],
             [[:Hepl, :NO], [:Npl, :O, :He], :(1.35e-9)],
             [[:Hepl, :NO], [:Opl, :N, :He], :(1.02e-10)],
             [[:Hepl, :O2], [:O2pl, :He], :(3.0e-11)],
             [[:Hepl, :O2], [:Opl, :O, :He], :(9.7e-10)],
             [[:Hepl, :OH], [:Opl, :He, :H], :(1.91e-8*((300 ./ T)^-0.5))],
             [[:HN2Opl, :CO], [:HCOpl, :N2O], :(5.3e-10)],
             [[:HN2Opl, :H2O], [:H3Opl, :N2O], :(2.83e-9)],
             [[:HNOpl, :C], [:CHpl, :NO], :(1.0e-9)],
             [[:HNOpl, :CH], [:CH2pl, :NO], :(1.07e-8*((300 ./ T)^-0.5))],
             [[:HNOpl, :CN], [:HCNpl, :NO], :(1.51e-8*((300 ./ T)^-0.5))],
             [[:HNOpl, :CO], [:HCOpl, :NO], :(8.6e-10)],
             [[:HNOpl, :CO2], [:HCO2pl, :NO], :(9.4e-10)],
             [[:HNOpl, :H2O], [:H3Opl, :NO], :(2.3e-9)],
             [[:HNOpl, :HCN], [:HCNHpl, :NO], :(1.71e-8*((300 ./ T)^-0.5))],
             [[:HNOpl, :HCO], [:H2COpl, :NO], :(1.25e-8*((300 ./ T)^-0.5))],
             [[:HNOpl, :N2], [:N2Hpl, :NO], :(1.0e-11)],
             [[:HNOpl, :NH], [:NH2pl, :NO], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:HNOpl, :NH2], [:NH3pl, :NO], :(1.52e-8*((300 ./ T)^-0.5))],
             [[:HNOpl, :NO], [:NOpl, :HNO], :(7.0e-10)],
             [[:HNOpl, :O], [:NO2pl, :H], :(1.0e-12)],
             [[:HNOpl, :OH], [:H2Opl, :NO], :(1.07e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :C], [:CHpl, :O2], :(1.0e-9)],
             [[:HO2pl, :CH], [:CH2pl, :O2], :(1.07e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :CN], [:HCNpl, :O2], :(1.49e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :CO], [:HCOpl, :O2], :(8.4e-10)],
             [[:HO2pl, :CO2], [:HCO2pl, :O2], :(1.1e-9)],
             [[:HO2pl, :H2], [:H3pl, :O2], :(3.3e-10)],
             [[:HO2pl, :H2O], [:H3Opl, :O2], :(1.42e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :HCN], [:HCNHpl, :O2], :(1.68e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :HCO], [:H2COpl, :O2], :(1.23e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :N], [:NO2pl, :H], :(1.0e-12)],
             [[:HO2pl, :N2], [:N2Hpl, :O2], :(8.0e-10)],
             [[:HO2pl, :NH], [:NH2pl, :O2], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :NH2], [:NH3pl, :O2], :(1.51e-8*((300 ./ T)^-0.5))],
             [[:HO2pl, :NO], [:HNOpl, :O2], :(7.7e-10)],
             [[:HO2pl, :O], [:OHpl, :O2], :(6.2e-10)],
             [[:HO2pl, :OH], [:H2Opl, :O2], :(1.06e-8*((300 ./ T)^-0.5))],
             [[:HOCpl, :CO], [:HCOpl, :CO], :(6.0e-10)],
             [[:HOCpl, :CO2], [:HCO2pl, :CO], :(9.45e-10)],
             [[:HOCpl, :H2], [:H3pl, :CO], :(2.68e-10)],
             [[:HOCpl, :H2], [:HCOpl, :H2], :(3.8e-10)],
             [[:HOCpl, :N2], [:N2Hpl, :CO], :(6.7e-10)],
             [[:HOCpl, :N2O], [:HN2Opl, :CO], :(1.17e-9)],
             [[:HOCpl, :NO], [:HNOpl, :CO], :(7.1e-10)],
             [[:HOCpl, :O2], [:HO2pl, :CO], :(1.9e-10)],
             [[:Hpl, :CH], [:CHpl, :H], :(3.29e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :CO2], [:HCOpl, :O], :(3.8e-9)],
             [[:Hpl, :H], [:H2pl], :(2.34e-22*((300 ./ T)^1.49)*expl(-(-228.0) ./ T))],
             [[:Hpl, :H2O], [:H2Opl, :H], :(8.2e-9)],
             [[:Hpl, :HCN], [:HCNpl, :H], :(1.1e-8)],
             [[:Hpl, :HCO], [:COpl, :H2], :(1.63e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :HCO], [:H2pl, :CO], :(1.63e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :HCO], [:HCOpl, :H], :(1.63e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :He], [:HeHpl], :(3.14e-19*((300 ./ T)^-0.24))],
             [[:Hpl, :HNO], [:NOpl, :H2], :(6.93e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :N2O], [:N2Hpl, :O], :(3.52e-10)],
             [[:Hpl, :N2O], [:N2Opl, :H], :(1.85e-9)],
             [[:Hpl, :NH], [:NHpl, :H], :(3.64e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :NH2], [:NH2pl, :H], :(5.02e-8*((300 ./ T)^-0.5))],
             [[:Hpl, :NO], [:NOpl, :H], :(1.9e-9)],
             [[:Hpl, :NO2], [:NOpl, :OH], :(1.9e-9)],
             [[:Hpl, :O], [:Opl, :H], :(3.75e-10)],
             [[:Hpl, :O2], [:O2pl, :H], :(1.17e-9)],
             [[:Hpl, :OH], [:OHpl, :H], :(3.64e-8*((300 ./ T)^-0.5))],
             [[:N2Hpl, :C], [:CHpl, :N2], :(1.1e-9)],
             [[:N2Hpl, :CH], [:CH2pl, :N2], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:N2Hpl, :CO], [:HCOpl, :N2], :(8.8e-10)],
             [[:N2Hpl, :CO2], [:HCO2pl, :N2], :(1.07e-9)],
             [[:N2Hpl, :H2], [:H3pl, :N2], :(5.1e-18)],
             [[:N2Hpl, :H2O], [:H3Opl, :N2], :(2.6e-9)],
             [[:N2Hpl, :HCN], [:HCNHpl, :N2], :(3.2e-9)],
             [[:N2Hpl, :HCO], [:H2COpl, :N2], :(1.26e-8*((300 ./ T)^-0.5))],
             [[:N2Hpl, :N2O], [:HN2Opl, :N2], :(1.25e-9)],
             [[:N2Hpl, :NH], [:NH2pl, :N2], :(1.11e-8*((300 ./ T)^-0.5))],
             [[:N2Hpl, :NH2], [:NH3pl, :N2], :(1.54e-8*((300 ./ T)^-0.5))],
             [[:N2Hpl, :NO], [:HNOpl, :N2], :(3.4e-10)],
             [[:N2Hpl, :O], [:OHpl, :N2], :(1.4e-10)],
             [[:N2Hpl, :OH], [:H2Opl, :N2], :(1.07e-8*((300 ./ T)^-0.5))],
             [[:N2Opl, :CO], [:CO2pl, :N2], :(1.11e-10)],
             [[:N2Opl, :CO], [:NOpl, :NCO], :(1.89e-10)],
             [[:N2Opl, :H2], [:HN2Opl, :H], :(2.56e-10)],
             [[:N2Opl, :H2], [:N2Hpl, :OH], :(1.04e-10)],
             [[:N2Opl, :H2O], [:H2Opl, :N2O], :(3.27e-8*((300 ./ T)^-0.5))],
             [[:N2Opl, :H2O], [:HN2Opl, :OH], :(3.64e-9*((300 ./ T)^-0.5))],
             [[:N2Opl, :N2O], [:NOpl, :NO, :N2], :(1.2e-11)],
             [[:N2Opl, :NO], [:NOpl, :N2O], :(2.3e-10)],
             [[:N2Opl, :NO2], [:NO2pl, :N2O], :(2.21e-10)],
             [[:N2Opl, :NO2], [:NOpl, :N2, :O2], :(4.29e-10)],
             [[:N2Opl, :O2], [:NOpl, :NO2], :(4.59e-11)],
             [[:N2Opl, :O2], [:O2pl, :N2O], :(2.24e-10)],
             [[:N2pl, :Ar], [:Arpl, :N2], :(2.0e-13)],
             [[:N2pl, :C], [:Cpl, :N2], :(1.1e-10)],
             [[:N2pl, :CH], [:CHpl, :N2], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:N2pl, :CN], [:CNpl, :N2], :(1.73e-9*((300 ./ T)^-0.5))],
             [[:N2pl, :CO], [:COpl, :N2], :(7.3e-11)],
             [[:N2pl, :CO2], [:CO2pl, :N2], :(8.0e-10)],
             [[:N2pl, :H2], [:N2Hpl, :H], :(1.87e-9*expl(-(-54.7) ./ T))],
             [[:N2pl, :H2O], [:H2Opl, :N2], :(1.9e-9)],
             [[:N2pl, :H2O], [:N2Hpl, :OH], :(5.04e-10)],
             [[:N2pl, :HCN], [:HCNpl, :N2], :(3.9e-10)],
             [[:N2pl, :HCO], [:HCOpl, :N2], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:N2pl, :HCO], [:N2Hpl, :CO], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:N2pl, :N], [:Npl, :N2], :(1.4e-11)],
             [[:N2pl, :N2O], [:N2Opl, :N2], :(6.0e-10)],
             [[:N2pl, :NH], [:NHpl, :N2], :(1.13e-8*((300 ./ T)^-0.5))],
             [[:N2pl, :NH2], [:NH2pl, :N2], :(1.54e-8*((300 ./ T)^-0.5))],
             [[:N2pl, :NO], [:NOpl, :N2], :(4.4e-10)],
             [[:N2pl, :O], [:NOpl, :N], :(1.3e-10)],
             [[:N2pl, :O], [:Opl, :N2], :(9.8e-12)],
             [[:N2pl, :O2], [:O2pl, :N2], :(5.0e-11)],
             [[:N2pl, :OH], [:OHpl, :N2], :(1.09e-8*((300 ./ T)^-0.5))],
             [[:NH2pl, :CH], [:CH2pl, :NH], :(6.06e-9*((300 ./ T)^-0.5))],
             [[:NH2pl, :CH], [:CHpl, :NH2], :(6.06e-9*((300 ./ T)^-0.5))],
             [[:NH2pl, :CN], [:HCNHpl, :N], :(1.73e-9*((300 ./ T)^-0.5))],
             [[:NH2pl, :H2], [:NH3pl, :H], :(1.95e-10)],
             [[:NH2pl, :H2O], [:H3Opl, :NH], :(2.73e-9)],
             [[:NH2pl, :H2O], [:NH3pl, :OH], :(8.7e-11)],
             [[:NH2pl, :H2O], [:NH4pl, :O], :(1.16e-10)],
             [[:NH2pl, :HCN], [:HCNHpl, :NH], :(2.08e-8*((300 ./ T)^-0.5))],
             [[:NH2pl, :HCO], [:H2COpl, :NH], :(7.45e-9*((300 ./ T)^-0.5))],
             [[:NH2pl, :HCO], [:HCOpl, :NH2], :(7.45e-9*((300 ./ T)^-0.5))],
             [[:NH2pl, :N], [:N2Hpl, :H], :(9.1e-11)],
             [[:NH2pl, :NH], [:NH3pl, :N], :(1.26e-8*((300 ./ T)^-0.5))],
             [[:NH2pl, :NH2], [:NH3pl, :NH], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:NH2pl, :NO], [:NOpl, :NH2], :(7.0e-10)],
             [[:NH2pl, :O], [:HNOpl, :H], :(7.2e-11)],
             [[:NH2pl, :O2], [:H2NOpl, :O], :(1.19e-10)],
             [[:NH2pl, :O2], [:HNOpl, :OH], :(2.1e-11)],
             [[:NH3pl, :CH], [:NH4pl, :C], :(1.2e-8*((300 ./ T)^-0.5))],
             [[:NH3pl, :H2], [:NH4pl, :H], :(4.4e-13)],
             [[:NH3pl, :H2O], [:NH4pl, :OH], :(2.5e-10)],
             [[:NH3pl, :HCO], [:NH4pl, :CO], :(7.27e-9*((300 ./ T)^-0.5))],
             [[:NH3pl, :NO], [:NOpl, :NH3], :(7.2e-10)],
             [[:NHpl, :C], [:CHpl, :N], :(1.6e-9)],
             [[:NHpl, :CH], [:CH2pl, :N], :(1.71e-8*((300 ./ T)^-0.5))],
             [[:NHpl, :CN], [:HCNpl, :N], :(2.77e-8*((300 ./ T)^-0.5))],
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
             [[:NHpl, :HCN], [:HCNHpl, :N], :(3.12e-8*((300 ./ T)^-0.5))],
             [[:NHpl, :HCO], [:H2COpl, :N], :(2.25e-8*((300 ./ T)^-0.5))],
             [[:NHpl, :N], [:N2pl, :H], :(1.3e-9)],
             [[:NHpl, :N2], [:N2Hpl, :N], :(6.5e-10)],
             [[:NHpl, :NH], [:NH2pl, :N], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:NHpl, :NH2], [:NH3pl, :N], :(2.6e-8*((300 ./ T)^-0.5))],
             [[:NHpl, :NO], [:N2Hpl, :O], :(1.78e-10)],
             [[:NHpl, :NO], [:NOpl, :NH], :(7.12e-10)],
             [[:NHpl, :O], [:OHpl, :N], :(1.0e-9)],
             [[:NHpl, :O2], [:HO2pl, :N], :(1.64e-10)],
             [[:NHpl, :O2], [:NOpl, :OH], :(2.05e-10)],
             [[:NHpl, :O2], [:O2pl, :NH], :(4.51e-10)],
             [[:NHpl, :OH], [:H2Opl, :N], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:NO2pl, :H], [:NOpl, :OH], :(1.9e-10)],
             [[:NO2pl, :H2], [:NOpl, :H2O], :(1.5e-10)],
             [[:NO2pl, :NO], [:NOpl, :NO2], :(2.75e-10)],
             [[:Npl, :CH], [:CNpl, :H], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:Npl, :CN], [:CNpl, :N], :(1.91e-8*((300 ./ T)^-0.5))],
             [[:Npl, :CO], [:COpl, :N], :(4.93e-10)],
             [[:Npl, :CO], [:Cpl, :NO], :(5.6e-12)],
             [[:Npl, :CO], [:NOpl, :C], :(6.16e-11)],
             [[:Npl, :CO2], [:CO2pl, :N], :(9.18e-10)],
             [[:Npl, :CO2], [:COpl, :NO], :(2.02e-10)],
             [[:Npl, :H2], [:NHpl, :H], :(5.0e-10*expl(-(-85.0) ./ T))],
             [[:Npl, :H2O], [:H2Opl, :N], :(2.7e-9)],
             [[:Npl, :HCN], [:HCNpl, :N], :(3.7e-9)],
             [[:Npl, :HCO], [:HCOpl, :N], :(7.79e-9*((300 ./ T)^-0.5))],
             [[:Npl, :HCO], [:NHpl, :CO], :(7.79e-9*((300 ./ T)^-0.5))],
             [[:Npl, :N], [:N2pl], :(9.44e-19*((300 ./ T)^0.24)*expl(-(-26.1) ./ T))],
             [[:Npl, :N2O], [:NOpl, :N2], :(5.5e-10)],
             [[:Npl, :NH], [:N2pl, :H], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:Npl, :NH], [:NHpl, :N], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:Npl, :NH2], [:NH2pl, :N], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:Npl, :NO], [:N2pl, :O], :(8.33e-11)],
             [[:Npl, :NO], [:NOpl, :N], :(4.72e-10)],
             [[:Npl, :O2], [:NOpl, :O], :(2.32e-10)],
             [[:Npl, :O2], [:O2pl, :N], :(3.07e-10)],
             [[:Npl, :O2], [:Opl, :NO], :(4.64e-11)],
             [[:Npl, :OH], [:OHpl, :N], :(6.41e-9*((300 ./ T)^-0.5))],
             [[:O2Dpl, :H2], [:H2pl, :O], :(4.4e-8*((300 ./ T)^-0.98)*expl(-(-302.4) ./ T))],
             [[:O2Dpl, :H2], [:Hpl, :OH], :(1.62e-8*((300 ./ T)^-0.95)*expl(-(-335.1) ./ T))],
             [[:O2Dpl, :H2], [:OHpl, :H], :(1.5e-9)],
             [[:O2pl, :C], [:COpl, :O], :(5.2e-11)],
             [[:O2pl, :C], [:Cpl, :O2], :(5.2e-11)],
             [[:O2pl, :CH], [:CHpl, :O2], :(5.37e-9*((300 ./ T)^-0.5))],
             [[:O2pl, :CH], [:HCOpl, :O], :(5.37e-9*((300 ./ T)^-0.5))],
             [[:O2pl, :HCO], [:HCOpl, :O2], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:O2pl, :HCO], [:HO2pl, :CO], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:O2pl, :N], [:NOpl, :O], :(1.0e-10)],
             [[:O2pl, :NH], [:HNOpl, :O], :(5.54e-9*((300 ./ T)^-0.5))],
             [[:O2pl, :NH], [:NO2pl, :H], :(5.54e-9*((300 ./ T)^-0.5))],
             [[:O2pl, :NH2], [:NH2pl, :O2], :(1.51e-8*((300 ./ T)^-0.5))],
             [[:O2pl, :NO], [:NOpl, :O2], :(4.6e-10)],
             [[:O2pl, :NO2], [:NO2pl, :O2], :(6.6e-10)],
             [[:OHpl, :C], [:CHpl, :O], :(1.2e-9)],
             [[:OHpl, :CH], [:CH2pl, :O], :(6.06e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :CH], [:CHpl, :OH], :(6.06e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :CN], [:HCNpl, :O], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:OHpl, :CO], [:HCOpl, :O], :(8.4e-10)],
             [[:OHpl, :CO2], [:HCO2pl, :O], :(1.35e-9)],
             [[:OHpl, :H2], [:H2Opl, :H], :(9.7e-10)],
             [[:OHpl, :H2O], [:H2Opl, :OH], :(1.59e-9)],
             [[:OHpl, :H2O], [:H3Opl, :O], :(1.3e-9)],
             [[:OHpl, :HCN], [:HCNHpl, :O], :(2.08e-8*((300 ./ T)^-0.5))],
             [[:OHpl, :HCO], [:H2COpl, :O], :(4.85e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :HCO], [:H2Opl, :CO], :(4.85e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :HCO], [:HCOpl, :OH], :(4.85e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :N], [:NOpl, :H], :(8.9e-10)],
             [[:OHpl, :N2], [:N2Hpl, :O], :(2.4e-10)],
             [[:OHpl, :N2O], [:HN2Opl, :O], :(9.58e-10)],
             [[:OHpl, :N2O], [:N2Opl, :OH], :(2.13e-10)],
             [[:OHpl, :N2O], [:NOpl, :HNO], :(1.46e-10)],
             [[:OHpl, :NH], [:NH2pl, :O], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :NH2], [:NH2pl, :OH], :(8.66e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :NH2], [:NH3pl, :O], :(8.66e-9*((300 ./ T)^-0.5))],
             [[:OHpl, :NO], [:HNOpl, :O], :(6.11e-10)],
             [[:OHpl, :NO], [:NOpl, :OH], :(8.15e-10)],
             [[:OHpl, :O], [:O2pl, :H], :(7.1e-10)],
             [[:OHpl, :O2], [:O2pl, :OH], :(3.8e-10)],
             [[:OHpl, :OH], [:H2Opl, :O], :(1.21e-8*((300 ./ T)^-0.5))],
             [[:Opl, :CH], [:CHpl, :O], :(6.06e-9*((300 ./ T)^-0.5))],
             [[:Opl, :CH], [:COpl, :H], :(6.06e-9*((300 ./ T)^-0.5))],
             [[:Opl, :CN], [:NOpl, :C], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:Opl, :CO2], [:O2pl, :CO], :(1.1e-9)],
             [[:Opl, :H], [:Hpl, :O], :(6.4e-10)],
             [[:Opl, :H2], [:OHpl, :H], :(1.62e-9)],
             [[:Opl, :H2O], [:H2Opl, :O], :(2.6e-9)],
             [[:Opl, :HCN], [:COpl, :NH], :(1.17e-9)],
             [[:Opl, :HCN], [:HCOpl, :N], :(1.17e-9)],
             [[:Opl, :HCN], [:NOpl, :CH], :(1.17e-9)],
             [[:Opl, :HCO], [:HCOpl, :O], :(7.45e-9*((300 ./ T)^-0.5))],
             [[:Opl, :HCO], [:OHpl, :CO], :(7.45e-9*((300 ./ T)^-0.5))],
             [[:Opl, :N2], [:NOpl, :N], :(4.58e-9*((300 ./ T)^-1.37)*expl(-(-28.592) ./ T))],
             [[:Opl, :N2O], [:N2Opl, :O], :(6.3e-10)],
             [[:Opl, :NH], [:NHpl, :O], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:Opl, :NH], [:NOpl, :H], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:Opl, :NH2], [:NH2pl, :O], :(1.73e-8*((300 ./ T)^-0.5))],
             [[:Opl, :NO], [:NOpl, :O], :(8.0e-13)],
             [[:Opl, :NO2], [:NO2pl, :O], :(1.6e-9)],
             [[:Opl, :NO2], [:NOpl, :O2], :(8.3e-10)],
             [[:Opl, :O2], [:O2pl, :O], :(2.1e-11)],
             [[:Opl, :OH], [:O2pl, :H], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:Opl, :OH], [:OHpl, :O], :(6.24e-9*((300 ./ T)^-0.5))],
             [[:ArHpl, :E], [:Ar, :H], :(1.0e-9)],
             [[:Arpl, :E], [:Ar], :(4.0e-12*((300 ./ T)^0.6))],
             [[:CHpl, :E], [:C, :H], :(1.65e-6*((300 ./ T)^-0.42))],
             [[:CNpl, :E], [:N, :C], :(3.12e-6*((300 ./ T)^-0.5))],
             [[:CO2pl, :E], [:CO, :O], :(3.03e-5*((300 ./ T)^-0.75))],
             [[:COpl, :E], [:O, :C], :(4.82e-6*((300 ./ T)^-0.55))],
             [[:COpl, :E], [:O1D, :C], :(2.48e-8*((300 ./ T)^-0.55))],
             [[:Cpl, :E], [:C], :(6.28e-10*((300 ./ T)^-0.59))],
             [[:H2Opl, :E], [:O, :H, :H], :(2.08e-5*((300 ./ T)^-0.74))],
             [[:H2Opl, :E], [:H2, :O], :(2.64e-6*((300 ./ T)^-0.74))],
             [[:H2Opl, :E], [:OH, :H], :(5.86e-6*((300 ./ T)^-0.74))],
             [[:H2pl, :E], [:H, :H], :(1.86e-7*((300 ./ T)^-0.43))],
             [[:H3Opl, :E], [:H2, :O, :H], :(9.68e-8*((300 ./ T)^-0.5))],
             [[:H3Opl, :E], [:H2O, :H], :(1.86e-6*((300 ./ T)^-0.5))],
             [[:H3Opl, :E], [:OH, :H, :H], :(4.47e-6*((300 ./ T)^-0.5))],
             [[:H3Opl, :E], [:OH, :H2], :(1.04e-6*((300 ./ T)^-0.5))],
             [[:H3pl, :E], [:H, :H, :H], :(8.46e-7*((300 ./ T)^-0.52))],
             [[:H3pl, :E], [:H2, :H], :(4.54e-7*((300 ./ T)^-0.52))],
             [[:HCNHpl, :E], [:CN, :H, :H], :(3.79e-6*((300 ./ T)^-0.65))],
             [[:HCNpl, :E], [:CN, :H], :(3.46e-6*((300 ./ T)^-0.5))],
             [[:HCO2pl, :E], [:CO, :O, :H], :(2.18e-7)],
             [[:HCO2pl, :E], [:CO, :OH], :(9.18e-8)],
             [[:HCO2pl, :E], [:CO2, :H], :(1.7e-8)],
             [[:HCOpl, :E], [:CH, :O], :(1.15e-7*((300 ./ T)^-0.64))],
             [[:HCOpl, :E], [:CO, :H], :(1.06e-5*((300 ./ T)^-0.64))],
             [[:HCOpl, :E], [:OH, :C], :(8.08e-7*((300 ./ T)^-0.64))],
             [[:Hepl, :E], [:He], :(9.28e-11*((300 ./ T)^-0.5))],
             [[:HN2Opl, :E], [:N2, :O, :H], :(3.81e-5*((300 ./ T)^-0.74))],
             [[:HN2Opl, :E], [:N2, :OH], :(4.38e-5*((300 ./ T)^-0.74))],
             [[:HNOpl, :E], [:NO, :H], :(5.2e-6*((300 ./ T)^-0.5))],
             [[:HO2pl, :E], [:O2, :H], :(5.2e-6*((300 ./ T)^-0.5))],
             [[:HOCpl, :E], [:CH, :O], :(1.7e-9*((300 ./ T)^1.2))],
             [[:HOCpl, :E], [:CO, :H], :(3.3e-5*((300 ./ T)^-1.0))],
             [[:HOCpl, :E], [:OH, :C], :(1.19e-8*((300 ./ T)^1.2))],
             [[:Hpl, :E], [:H], :(6.46e-14*((300 ./ T)^0.7))],
             [[:N2Hpl, :E], [:N2, :H], :(6.6e-7*((300 ./ T)^-0.51))],
             [[:N2Hpl, :E], [:NH, :N], :(1.17e-6*((300 ./ T)^-0.51))],
             [[:N2Opl, :E], [:N, :N, :O], :(1.36e-6*((300 ./ T)^-0.57))],
             [[:N2Opl, :E], [:N2, :O], :(4.09e-6*((300 ./ T)^-0.57))],
             [[:N2Opl, :E], [:NO, :N], :(3.07e-6*((300 ./ T)^-0.57))],
             [[:N2pl, :E], [:N, :N], :(5.09e-7*((300 ./ T)^-0.39))],
             [[:N2pl, :E], [:N2D, :N2D], :(1.42e-6*((300 ./ T)^-0.39))],
             [[:NH2pl, :E], [:N, :H, :H], :(1.71e-5*((300 ./ T)^-0.8)*expl(-(-17.1) ./ T))],
             [[:NH2pl, :E], [:NH, :H], :(8.34e-6*((300 ./ T)^-0.79)*expl(-(-17.1) ./ T))],
             [[:NH3pl, :E], [:NH, :H, :H], :(2.68e-6*((300 ./ T)^-0.5))],
             [[:NH3pl, :E], [:NH2, :H], :(2.68e-6*((300 ./ T)^-0.5))],
             [[:NHpl, :E], [:N, :H], :(7.45e-7*((300 ./ T)^-0.5))],
             [[:NO2pl, :E], [:NO, :O], :(5.2e-6*((300 ./ T)^-0.5))],
             [[:NOpl, :E], [:O, :N], :(8.52e-7*((300 ./ T)^-0.37))],
             [[:NOpl, :E], [:O, :N2D], :(4.53e-9*((300 ./ T)^0.75))],
             [[:Npl, :E], [:N], :(1.9e-10*((300 ./ T)^-0.7))],
             [[:O2pl, :E], [:O, :O], :(8.15e-6*((300 ./ T)^-0.65))],
             [[:OHpl, :E], [:O, :H], :(6.5e-7*((300 ./ T)^-0.5))],
             [[:Opl, :E], [:O], :(1.4e-10*((300 ./ T)^-0.66))]
             ]
