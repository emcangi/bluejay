################################################################################
# plot_rxn_rates.jl
# TYPE: Analysis (optional)
# WHICH: Equilibrium experiments
# DESCRIPTION: Plots chemical reaction rates by altitude and reaction, total 
# chemical production or consumption rates, total transport rates, 

# Eryn Cangi
# 5 April 2019
# Last edited: 30 November 2020
# Currently tested for Julia: 1.4.1
################################################################################
using PyPlot
using HDF5
using LaTeXStrings
using PyCall
using PlotUtils
using JLD
using Photochemistry
import Photochemistry.fluxcoefs, Photochemistry.scaleH, Photochemistry.getflux

patches = pyimport("matplotlib.patches")
linez = pyimport("matplotlib.lines")

# Parameters and modified reaction network =====================================
include("/home/emc/GDrive-CU/Research-Modeling/UpperAtmoDH/Code/PARAMETERS.jl")

# NOTE: You cannot use reactionnet as defined in PARAMETERS.jl. It has to be
# this one here because we must change all the operators to be array operators.
# These redefinitions are basically overrides...
# Also, the threebody functions must be as used here for the same reason.
# Net last checked: Sept 2020

threebody(k0, kinf) = :($k0 .* M ./ (1 .+ $k0 .* M ./ $kinf).*0.6 .^ ((1 .+ (log10.($k0 .* M ./ $kinf)) .^2).^-1.0))
threebodyca(k0, kinf) = :($k0 ./ (1 .+ $k0 ./ ($kinf ./ M)).*0.6 .^ ((1 .+ (log10.($k0 ./ ($kinf .* M))) .^2).^-1.0))

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
             # NEW: neutral dissociation from Roger Yelle model
             # [[:CO], [:C, :O], :JCOtoCpO],
             # [[:N2O], [:N2, :O1D], :JN2OtoN2pO1D],
             # [[:NO2], [:NO, :O], :JNO2toNOpO],
             # [[:NO], [:N, :O], :JNOtoNpO],
             # [[:CO2], [:C, :O, :O], :JCO2toCpOpO],
             # [[:CO2], [:C, :O2], :JCO2toCpO2],
             # NEW: photoionization from Roger's model
             # [[:CO2], [:CO2pl], :JCO2toCO2pl], 

             # [[:CO2], [:CO2plpl], :JCO2toCO2plpl],  # turn off if worse problems
             # [[:CO2], [:Cplpl, :O2], :JCO2toCplplpO2],
             # [[:CO2], [:Cpl, :O2], :JCO2toCplpO2],
             # [[:CO2], [:COpl, :Opl], :JCO2toCOplpOpl],
             # [[:CO2], [:COpl, :O], :JCO2toCOplpO],
             # [[:CO2], [:Opl, :CO], :JCO2toOplpCO],
             # [[:CO2], [:Opl, :Cpl, :O], :JCO2toOplpCplpO], # turn off if worse problems

             # [[:H2O], [:H2Opl], :JH2OtoH2Opl],
             # [[:H2O], [:Opl, :H2], :JH2OtoOplpH2],

             # [[:H2O], [:Hpl, :OH], :JH2OtoHplpOH], # turn off if worse problems
             # [[:H2O], [:OHpl, :H], :JH2OtoOHplpH],
             # [[:CO], [:COpl], :JCOtoCOpl],
             # [[:CO], [:C, :Opl], :JCOtoCpOpl],
             # [[:CO], [:O, :Cpl], :JCOtoOpCpl],
             # [[:N2], [:N2pl], :JN2toN2pl],
             # [[:N2], [:Npl, :N], :JN2toNplpN],
             # [[:NO2], [:NO2pl], :JNO2toNO2pl],
             # [[:NO], [:NOpl], :JNOtoNOpl],
             # [[:N2O], [:N2Opl], :JN2OtoN2Opl],
             # [[:H], [:Hpl], :JHtoHpl],
             # [[:H2], [:H2pl], :JH2toH2pl],
             # [[:H2], [:Hpl, :H], :JH2toHplpH],
             # [[:H2O2], [:H2O2pl], :JH2O2toH2O2pl], # turn off if worse problems

             # [[:O], [:Opl], :JOtoOpl],
             # [[:O2], [:O2pl], :JO2toO2pl],
             # [[:O3], [:O3pl], :JO3toO3pl],# turn off if worse problems

             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8 .* 3.0e-33 .* (300 ./ Tn).^3.25)],
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35 .* exp.(724 ./ Tn))],
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5 .* 6.0e-34 .* (300 ./ Tn).^2.4)],
             [[:O, :O3], [:O2, :O2], :(8.0e-12 .* exp.(-2060 ./ Tn))],  # Sander 2011
             [[:O, :CO, :M], [:CO2, :M], :(2.2e-33 .* exp.(-1780 ./ Tn))],

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.2e-11 .* exp.(70 ./ Tn))], # verified NIST 4/3/18
             [[:O1D, :O3], [:O2, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :O3], [:O, :O, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11 .* exp.(115 ./ Tn))], # Sander2011. NIST: 7.41e-11 .* exp.(120/Tn)
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Sander2011. Yung89: 1e-10; NIST 1.1e-10
             [[:O1D, :HD], [:H, :OD], :(0.41 .* 1.2e-10)], # Yung88: rate 0.41 .* H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41 .* 1.2e-10)], # Yung88: rate 0.41 .* H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10 .* exp.(60 ./ Tn))], # Sander2011. Yung89: 2.2e-10; NIST: 1.62e-10 .* exp.(65/Tn)
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10 .* exp.(60 ./ Tn))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12 .* exp.(-4000 ./ Tn))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12 .* exp.(-4390 ./ Tn))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12 .* exp.(-4400 ./ Tn))], # NIST
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11 .* exp.(-4038 ./ Tn))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11 .* exp.(-3821 ./ Tn))], # NIST (1986, 200-300K): 8.19e-13 .* exp.(-2700/Tn)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12 .* exp.(-1800 ./ Tn))], # Sander2011. Yung89: 5.5e-12 .* exp.(-2000/Tn). KIDA: 7.7E-12 .* exp.(-2100/Tn). old rate from Mike: 9.01e-13 .* exp.(-1526/Tn)
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.) .* 2.8e-12 .* exp.(-1800 ./ Tn))], # Yung88: rate (3/20) .* H-ana. Sander2011: 5e-12 .* exp.(-2130 ./ Tn)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.) .* 2.8e-12 .* exp.(-1800 ./ Tn))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12 .* exp.(-1800 ./ Tn))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:OD, :HD], [:D2O, :H], :(???)],  # can't find a rate...?

             # recombination of H. Use EITHER the first line OR the 2nd and 3rd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32 .* (298 ./ Tn).^2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32 .* (298 ./ Tn).^2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32 .* (298 ./ Tn).^2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn).^2)], # Can't find in databases. Mike's rate.
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn).^2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn).^2)], # not in Yung88. assumed rate

             ## H + HO2
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Sander2011. Indep of Tn for 245<Tn<300
             [[:H, :HO2], [:H2, :O2], :(0.5 .* 6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # O1D is theoretically mandated
             [[:H, :DO2], [:OH, :OD], :(7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HD, :O2], :(0.5 .* 6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HDO, :O1D], :(1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:D, :HO2], [:OH, :OD], :(0.71 .* 7.2e-11)], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:D, :HO2], [:HD, :O2], :(0.71 .* 0.5 .* 6.9e-12)], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:D, :HO2], [:HDO, :O1D], :(0.71 .* 1.6e-12)], # Yung88: rate 0.71 .* H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HO2, :D], :(1e-10 ./ (0.54 .* exp.(890 ./ Tn)))], # Yung88 (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12 .* exp.(-1890 ./ Tn))], # verified NIST 4/3/18. Only valid for Tn>300K. No exp.eriment for lower.
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

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0 .* 4.4e-32 .* (Tn/300.).^-1.3), # Sander2011, 300K+. Yung89: 5.5e-32(Tn/300).^-1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11 .* (Tn/300.).^0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0 .* 4.4e-32 .* (Tn/300.).^-1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11 .* (Tn/300.).^0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10 .* exp.(-470 ./ Tn))], # verified Yung89, NIST 4/3/18
             [[:D, :O3], [:OD, :O2], :(0.71 .* 1.4e-10 .* exp.(-470 ./ Tn))], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11 .* exp.(180 ./ Tn))], # Sander2011. KIDA+NIST 4/3/18 150-500K: 2.4e-11 .* exp.(110 ./ Tn). Yung89: 2.2e-11 .* exp.(120/Tn) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11 .* exp.(180 ./ Tn))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11 .* exp.(200 ./ Tn))], # Sander2011. KIDA (220-400K): 2.7e-11 .* exp.(224/Tn)
             [[:O, :DO2], [:OD, :O2], :(3.0e-11 .* exp.(200 ./ Tn))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12 .* exp.(-2000 ./ Tn))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5 .* 1.4e-12 .* exp.(-2000 ./ Tn))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5 .* 1.4e-12 .* exp.(-2000 ./ Tn))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(4.2e-12 .* exp.(-240 ./ Tn))], # NIST+KIDA, 200-350K: 6.2e-14 .* (Tn/300).^2.62 .* exp.(945 ./ Tn) changed 4/3/18. Yung89: 4.2e-12 .* exp.(-240/Tn). old rate w/mystery origin: 1.8e-12.
             [[:OD, :OH], [:HDO, :O], :(4.2e-12 .* exp.(-240 ./ Tn))], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3 .* 6.9e-31 .* (Tn/300.).^-1.0),:(2.6e-11))], # Sander2011. Why 1.3?
             [[:OD, :OH], [:HDO2], threebody(:(1.3 .* 6.9e-31 .* (Tn/300.).^-1.0),:(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12 .* exp.(-940 ./ Tn))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12 .* exp.(-940 ./ Tn))], # Yung88: rate same as H-ana
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # verified NIST 4/3/18. Yung89: 4.6e-11 .* exp.(230/Tn) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12 .* exp.(-160 ./ Tn))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12 .* exp.(-200/Tn). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5 .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: rate 0.5 .* H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5 .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: rate 0.5 .* H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14 .* exp.(-490 ./ Tn))], # Sander2011. Yung89: 1.1e-14 .* exp.(-500/Tn). KIDA 250-340K: 2.03e-16 .* (Tn/300).^4.57 .* exp.(693/Tn). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(1.0e-14 .* exp.(-490 ./ Tn))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13 .* exp.(460 ./ Tn))], # Sander2011. Yung89: 2.3e-13 .* exp.(600/Tn). KIDA 230-420K: 2.2e-13 .* exp.(600/Tn)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13 .* exp.(460 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2 .* 2.1e-33 .* exp.(920 ./ Tn))], # Sander2011.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2 .* 2.1e-33 .* exp.(920 ./ Tn))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9 .* (Tn.^-0.63) ./ (0.72 .* exp.(717 ./ Tn)))], # rate: Yung88. NIST (Howard82): 5.25E-11 .* (Tn/298).^-0.63  - turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9 .* Tn.^-0.63)], # Yung88  - turn off for Case 2

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13 .* (Tn/300.).^0.6),:(2.1e9 .* (Tn/300.).^6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13 .* (Tn/300.).^0.6),:(2.1e9 .* (Tn/300.).^6.1))], # Yung88: same as H-ana.
             ### possible deuterated analogues below
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33 .* (Tn/300.).^-1.4),:(1.1e-12 .* (Tn/300.).^1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33 .* (Tn/300.).^-1.4),:(1.1e-12 .* (Tn/300.).^1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 ./ Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5) .* 8.7e-10)],

             # NEW
             #--deleted--
             ];

# TODO: Write a test here that makes sure all the operators are array operators.


# Functions ====================================================================

# CAUTION: ALL THESE FUNCTIONS ARE PASTED IN FROM converge_new_file.jl AND HEAVILY
# MODIFIED TO WORK IN THIS SCRIPT. 
function getflux(n_current, species, t, exptype)
    #=
    Special overload for this file
    Returns a 1D array of fluxes in and out of a given altitude level for a 
    given species. For looking at vertical distribution of fluxes, but it does 
    not modify the concentrations.

    n_current: Array; species number density by altitude
    dz: Float64; layer thickness in cm
    species: Symbol

    returns: Array of raw flux value (#/cm^2/s) at each altitude layer
    =#

    # BEGIN copy-pasted stuff ----------------------------------------------------------------------------------------
    if exptype=="tropo"
        thetemps = [meanTs, t, meanTe]
    elseif exptype=="exo"
        thetemps = [meanTs, meanTt, t]
    elseif exptype=="surf"
        thetemps = [t, meanTt, meanTe]
    end

    # Boundary conditions are needed to get the flux, so we have to copy this little section over from converge_new_file.
    # We can't make it globally available because it depends on the temperature profile for the effusion velocities.
    # TODO: following 3 lines could be moved to a global declaration
    Temp_keepSVP(z::Float64) = Tpiecewise(z, meanTs, meanTt, meanTe)  
    H2Osat = map(x->Psat(x), map(Temp_keepSVP, alt)) # for holding SVP fixed
    HDOsat = map(x->Psat_HDO(x), map(Temp_keepSVP, alt))  # for holding SVP fixed

    speciesbclist=Dict(
                :CO2=>["n" 2.1e17; "f" 0.],
                :Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
                :N2=>["n" 1.9e-2*2.1e17; "f" 0.],
                :H2O=>["n" H2Osat[1]; "f" 0.],
                :HDO=>["n" HDOsat[1]; "f" 0.],
                :O=>["f" 0.; "f" 1.2e8],
                :H2=>["f" 0.; "v" effusion_velocity(T(zmax, thetemps[1], thetemps[2], thetemps[3], "neutral"), 2.0, zmax)],
                :HD=>["f" 0.; "v" effusion_velocity(T(zmax, thetemps[1], thetemps[2], thetemps[3], "neutral"), 3.0, zmax)],
                :H=>["f" 0.; "v" effusion_velocity(T(zmax, thetemps[1], thetemps[2], thetemps[3], "neutral"), 1.0, zmax)],
                :D=>["f" 0.; "v" effusion_velocity(T(zmax, thetemps[1], thetemps[2], thetemps[3], "neutral"), 2.0, zmax)],
               );
    # END copy-pasted stuff ------------------------------------------------------------------------------------------

    # each element in thesecoefs has the format [downward, upward]
    thesecoefs = [fluxcoefs(a, dz, species, n_current, thetemps) for a in alt[2:end-1]]

    # thesebcs has the format [lower bc; upper bc], where each row contains a 
    # character showing the type of boundary condition, and a number giving its value
    thesebcs = boundaryconditions(species, dz, n_current, thetemps, exptype, speciesbclist)

    thesefluxes = fill(convert(Float64, NaN),length(intaltgrid))

    # in the following line for the lowest layer: 
    # first term is -(influx from layer above - outflux from this layer)
    # second term is (-this layer's lower bc that depends on concentration + bc that doesn't depend on concentration)
    thesefluxes[1] = (-(n_current[species][2]*thesecoefs[2][1]
                        -n_current[species][1]*thesecoefs[1][2]) 
                    +(-n_current[species][1]*thesebcs[1, 1]
                      +thesebcs[1, 2]))/2.0
    for ialt in 2:length(intaltgrid)-1
        thesefluxes[ialt] = (-(n_current[species][ialt+1]*thesecoefs[ialt+1][1]  # coming in from above
                               -n_current[species][ialt]*thesecoefs[ialt][2])    # leaving out to above layer
                             +(-n_current[species][ialt]*thesecoefs[ialt][1]     # leaving to the layer below
                               +n_current[species][ialt-1]*thesecoefs[ialt-1][2]))/2.0  # coming in from below
    end
    thesefluxes[end] = (-(thesebcs[2, 2]
                          - n_current[species][end]*thesebcs[2, 1])
                        + (-n_current[species][end]*thesecoefs[end][1]
                           +n_current[species][end-1]*thesecoefs[end-1][2]))/2.0
    return dz*thesefluxes
end

function fluxcoefs(z, dz, species, n_current, thetemps)
    #=
    Special overload for this file:
    1) generates the coefficients K, D, T, Hs if they are not supplied (most common)
    2) Allows passing in a specific temperature parameter (t) and experiment type

    z: a specific altitude in cm
    dz: thickness of an altitude later (2 km, but in cm)
    species: the species for which to calculate the coefficient. Symbol
    n_current: array of species densities by altitude, the current state of the atmosphere
    t: the particular temperature parameter that is being varied to plot reaction rates

    p: upper layer ("plus")
    0: this layer
    m: lower layer ("minus")
    =#

    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    Tp = T(z+dz, thetemps[1], thetemps[2], thetemps[3], species_type)
    T0 = T(z, thetemps[1], thetemps[2], thetemps[3], species_type)
    Tm = T(z-dz, thetemps[1], thetemps[2], thetemps[3], species_type)

    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)

    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = scaleH(z+dz, species, thetemps)
    Hs0 = scaleH(z, species, thetemps)
    Hsm = scaleH(z-dz, species, thetemps)
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz, [Km , K0, Kp], [Dm , D0, Dp], [Tm , T0, Tp],
                     [Hsm, Hs0, Hsp], [H0m, H00, H0p], species)
end

function scaleH(z, species::Symbol, thetemps)
    #=
    Special overload for this file
    Same as first scaleH, but for a particular atomic/molecular species.
    =#  

    # since we need to pass in temperatures.
    #Th is normally T, but I changed the temperature function to T, so now it's Th. For no reason. T related to scale height I guess
    Th = T(z, thetemps[1], thetemps[2], thetemps[3], charge_type(species))
    mm = speciesmolmasslist[species]
    return boltzmannK*Th/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
end

function boundaryconditions(species, dz, n_current, thetemps, exptype, speciesbclist)
    #= 
    Special overload for this file
    =#

    # TODO: This function shouldn't require changing, but check the calls to fluxcoefs and lower_up and upper_down.

    bcs = speciesbcs(species, speciesbclist)
    if issubset([species],notransportspecies)
        bcs = ["f" 0.; "f" 0.]
    end

    # first element returned corresponds to lower BC, second to upper
    # BC transport rate. Within each element, the two rates correspond
    # to the two equations
    # n_b  -> NULL (first rate, depends on species concentration)
    # NULL -> n_b  (second rate, independent of species concentration 
    bcvec = Float64[0 0;0 0]

    # LOWER
    if bcs[1, 1] == "n"
        bcvec[1,:]=[fluxcoefs(alt[2], dz, species, n_current, thetemps)[1],
                    lower_up(alt[1], dz, species, n_current, thetemps)*bcs[1, 2]]
    elseif bcs[1, 1] == "f"
        bcvec[1,:] = [0.0, bcs[1, 2]/dz]
    elseif bcs[1, 1] == "v"
        bcvec[1,:] = [bcs[1, 2]/dz, 0.0]
    else
        throw("Improper lower boundary condition!")
    end

    # UPPER
    if bcs[2, 1] == "n"
        bcvec[2,:] = [fluxcoefs(alt[end-1],dz, species, n_current, thetemps)[2],
                    upper_down(alt[end],dz, species, n_current, thetemps)*bcs[1, 2]]
    elseif bcs[2, 1] == "f"
            bcvec[2,:] = [0.0,-bcs[2, 2]/dz]
    elseif bcs[2, 1] == "v"
        bcvec[2,:] = [bcs[2, 2]/dz, 0.0]
    else
        throw("Improper upper boundary condition!")
    end

    return bcvec
end

function lower_up(z, dz, species, n_current, thetemps)
    #= 
    Special overload for this file
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one above. 
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus") 

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    n_current: Array; species number density by altitude

    returns: array of fluxcoefs
    =#

    species_type = charge_type(species)

    Tp = T(z+dz, thetemps[1], thetemps[2], thetemps[3], species_type)
    T0 = T(z, thetemps[1], thetemps[2], thetemps[3], species_type)
    Tm = 1

    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = 1
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z,nt0)
    Km = 1

    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = 1
    Hsp = scaleH(z+dz, species, thetemps)
    Hs0 = scaleH(z, species, thetemps)
    Hsm = 1
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
    H0m = 1

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[2]
end

function upper_down(z, dz, species, n_current, thetemps)
    #= 
    Special overload for this file
    define transport coefficients for a given atmospheric layer for
    transport from that layer to the one below. 
    p: layer above ("plus"), 0: layer at altitude z, m: layer below ("minus") 

    z: altitude in cm
    dz: altitude layer thickness ("resolution"), in cm
    species: Symbol; species for which this coefficients are calculated
    n_current: Array; species number density by altitude

    returns: return of fluxcoefs
    =#

    # TODO: Update to match form in converge_atmo_instantiate_ions, but without the "Temp_n" type aliasing
    # since we need to pass in temperatures.
    # set temps of nearby layers; depends on ion/electron/neutral
    species_type = charge_type(species)

    Tp = 1
    T0 = T(z, thetemps[1], thetemps[2], thetemps[3], species_type)
    Tm = T(z-dz, thetemps[1], thetemps[2], thetemps[3], species_type)

    ntp = 1
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = 1
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)

    Dp = 1
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = 1
    Hs0 = scaleH(z, species, thetemps)
    Hsm = scaleH(z-dz, species, thetemps)
    H0p = 1
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz,
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[1]
end

# these functions written specifically for this script
function make_ratexdensity(n_current, rxn_net, t, exptype; species=Nothing, species_role=Nothing)
    #=
    n_current: a given result file for a converged atmosphere
    rxn_net: the reaction network.
    exptype: "surf", "tropo", "exo", just allows for specifiying the temperature profile.
    t: a specified temperature parameter, either T_surf, T_tropo, or T_exo; which is identified by exptype.
    species: only reactions including this species will be plotted. If it has a value, so must species_role.
    species_role: whether to look for the species as a reactant, product, or both.  If it has a value, so must species.
    =#

    # Calculate the sum of all bodies in a layer (M) for third body reactions. 
    # This does it in an array so we can easily plot.
    M_by_alt = sum([n_current[sp] for sp in fullspecieslist]) 
    E_by_alt = sum([n_current[sp] for sp in ionlist])

    # Empty dictionary in which to store reaction rates by altitude.
    rxn_dat =  Dict{String,Array{Float64, 1}}()

    # TODO: Put in some logic to catch problems where species or species_role is specified but not the other.

    for rxn in rxn_net
        reactants = rxn[1]
        products = rxn[2]  # vector of the product symbols

        # Construct temperature profile (by altitude)
        temps_neutrals = Array{Float64}(undef, length(alt)-2)
        temps_ions = Array{Float64}(undef, length(alt)-2)
        temps_electrons = Array{Float64}(undef, length(alt)-2)

        three_temps = Dict("surf"=>[t, meanTt, meanTe],
                          "tropo"=>[meanTs, meanTt, meanTe],
                          "exo"=>[meanTs, meanTt, t])

        i = 1
        for i in range(1, length=length(alt)-2)
            temps_neutrals[i] = T(alt[i], three_temps[exptype][1], three_temps[exptype][2], three_temps[exptype][3], "neutral")
            temps_ions[i] = T(alt[i], three_temps[exptype][1], three_temps[exptype][2], three_temps[exptype][3], "ion")
            temps_electrons[i] = T(alt[i], three_temps[exptype][1], three_temps[exptype][2], three_temps[exptype][3], "electron")
            i += 1
        end


        # This skips reactions that don't invovle the species in question.
        if species != Nothing
            role = Dict("reactant"=>reactants, "product"=>products,
                        "both"=>[reactants, products])
            if ~in(species, role[species_role])
                continue
            end
        end

        # get the reactants and products in string form for use in plot labels
        reacts = join(rxn[1], " + ")
        prods = join(rxn[2], " + ")
        rxn_str = string(reacts) * " --> " * string(prods)

        # calculate the reaction strength, rate coefficient * species density.
        # in standard notation, rate = k[A][B]...
        if typeof(rxn[3]) == Symbol # for photodissociation
            alt_arr = n_current[rxn[1][1]] # gets reactant density by altitude
            rate_arr = n_current[rxn[3]]
            ratexdens = rate_arr .* alt_arr
            # put ratexdens in dictionary - photodissociation has only 1 reactant
            rxn_dat[rxn_str] = ratexdens
        else
            # for reactions with more than one reactant
            ratexdens = ones(length(alt)-2)
            j = 1
            for r in rxn[1]
                if r != :M && r != :E
                    # species densities by altitude
                    ratexdens .*= n_current[r]  # multiply by each reactant density
                elseif r == :M
                    ratexdens .*= M_by_alt
                elseif r == :E
                    ratexdens .*= E_by_alt
                else
                    throw("Got an unknown symbol in a reaction rate: $(r)")
                end
            end
            rate = rxn[3]

            # WARNING: This works. I don't really know how. At some point I did. WITCHCRAFT
            @eval ratefunc(Tn, Ti, Te, M, E) = $rate
            rate_arr = Base.invokelatest(ratefunc, temps_neutrals, temps_ions, temps_electrons, M_by_alt, E_by_alt)
            # println("The reaction is: $(rxn_str)")
            # println("The reaction rate is: $(rxn[3])")
            # println("Size of rate_arr: $(size(rate_arr))")
            # println("Size of ratexdens: $(size(ratexdens))")
            # println()
            ratexdens .*= rate_arr
            rxn_dat[rxn_str] = ratexdens
        end
    end
    return rxn_dat
end

function troubleshoot_rates(sp, prob_file, t, exptype; plot_indiv_rxns=false, thresh=1e-6, xlims=Nothing)
    #=
    sp: species to care about
    prob_file: the file that won't converge
    t: the non-mean value of the temperature for T_exptype, i.e. might be 150 for T_surf instead of the mean 216.
    exptype: whether the "surf", "tropo" or "exo" temperature is the parameter being examined
    plot_indiv_rxns: whether to plot lines for the fastest chemical reactions. if false, only the total will be plotted.
                     I don't recommend turning this parameter on
    thresh: a threshhold of reaction rate. plot will only plot the reactions that have values above this threshhold
    xlims: x-axis limits for the main axis that plots chemical reaction rate, and the 
           secondary axis that plots flux. Order: [axis1min, axis1max, axis2min, axis2max]

    This function is for troubleshooting imbalances in the chemical reaction network. It should be used with 
    plot_indiv_rxns=true (yes really. )
    =#
    # ------------------------------------------------------------------------------
    # Basic setup: set up the pathname for the readfile
    readfile = basepath * prob_file

    # --------------------------------------------------------------------------------
    # Load the simulation's atmosphere array
    ncur = get_ncurrent(readfile)

    # --------------------------------------------------------------------------------
    # calculate reaction rates x density of the species at each level of the atmosphere.
    rxd_prod = make_ratexdensity(ncur, reactionnet, t, exptype, species=sp, species_role="product")
    rxd_consume = make_ratexdensity(ncur, reactionnet, t, exptype, species=sp, species_role="reactant")

    # ---------------------------------------------------------------------------------
    # Calculate the fluxes for the species
    fluxarray = getflux(ncur, sp, t, exptype)
    
    # ---------------------------------------------------------------------------------
    # Plot reaction rates and transport rates by altitude
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

    fig, ax = subplots(figsize=(8,6))
    
    subplots_adjust(wspace=0, bottom=0.25)
    plot_bg(ax)

    # Calculate the total reactions per second for this species of interest
    total_prod_rate = Array{Float64}(undef, length(alt)-2)
    total_prod_rate .= 0
    total_loss_rate = Array{Float64}(undef, length(alt)-2)
    total_loss_rate .= 0

    minx = 1e10  # will store the minimum value found in reaction rates
    maxx = 0

    # Collect chem production equations and total 
    for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
        lbl = "$(kv[1])"
        if plot_indiv_rxns == true
            if any(x->x>thresh, kv[2])
                ax.semilogx(kv[2], alt[1:length(alt)-2]./1e5, linestyle="-", linewidth=1, label=lbl)
            end
            if xlims == Nothing
                if minimum(kv[2]) <= minx
                    minx = minimum(kv[2])
                end
                if maximum(kv[2]) >= maxx
                    maxx = maximum(kv[2])
                end
            end
        end
        total_prod_rate += kv[2]
    end

    # Collect chem consumption equations and total 
    for kv in rxd_consume  # loop through the dict of format reaction => [rates by altitude]
        lbl = "$(kv[1])"
        if plot_indiv_rxns == true
            if any(x->x>thresh, kv[2])
                ax.semilogx(kv[2], alt[1:length(alt)-2]./1e5, linestyle="--", linewidth=1, label=lbl)
            end
            if xlims == Nothing
                if minimum(kv[2]) <= minx
                    minx = minimum(kv[2])
                end
                if maximum(kv[2]) >= maxx
                    maxx = maximum(kv[2])
                end
            end
        end
        total_loss_rate += kv[2]
    end

    # check to see if the total rates are going to move the xmin and xmax
    if minimum(total_prod_rate) <= minx
        minx = minimum(total_prod_rate)
    end
    if maximum(total_prod_rate) >= maxx
        maxx = maximum(total_prod_rate)
    end

    if minimum(total_loss_rate) <= minx
        minx = minimum(total_loss_rate)
    end
    if maximum(total_loss_rate) >= maxx
        maxx = maximum(total_loss_rate)
    end

    ax.semilogx(total_prod_rate, alt[1:length(alt)-2]./1e5, color="blue", linewidth=3, label="Total production")
    ax.semilogx(total_loss_rate, alt[1:length(alt)-2]./1e5, color="red", linewidth=3, label="Total loss")

    # make the legend
    blue_line = linez.Line2D([0], [0], color="blue", linewidth=1, linestyle="-")
    red_line = linez.Line2D([0], [0], color="red", linewidth=1, linestyle="-")
    handles = [blue_line, red_line]
    labels = [ "Total chem. prod.", "Total chem. loss"]
    ax.legend()

    # labels and such 
    if xlims == Nothing
        if minx < 1e-20
            minx = 1e-20
        end
        maxx = 10^(ceil(log10(maxx)))
        xlims = [minx, maxx]
    end
    ax.set_title("Reaction rates $(string(sp))", fontsize=20)
    ax.set_xlim(xlims[1], xlims[2])
    ax.set_ylabel("Altitude (km)")
    ax.set_xlabel("Chemical reaction rate ("*L"cm^{-3}s^{-1})")
    savefig(results_dir*"/oldnetwork/chem_rates_$(sp).png", bbox_inches="tight", dpi=300)
    close(fig)
end

# do the stuff ===================================================================
basepath = results_dir

const oldspecieslist = [# neutrals
                         :Ar, :CO, :CO2, 
                         :H, :H2, :H2O, :H2O2, :HO2, :HOCO, 
                         :N2, :O, :O1D, :O2, :O3, :OH,

                         # Neutrals - deuterated
                         :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD, 

                         # Ions
                         :CO2pl, ]; 


for s in oldspecieslist
    troubleshoot_rates(s, "ncurrent_0.01_4.h5", 216.0, "surf"; plot_indiv_rxns=false, thresh=0)
end
