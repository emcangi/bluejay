
################################################################################
############################### REACTION NETWORK ###############################
################################################################################

# function to replace three body rates with the recommended expression

threebody(k0, kinf) = :($k0 .* M ./ (1 .+ $k0 .* M ./ $kinf).*0.6 .^ ((1 .+ (log10.($k0 .* M ./ $kinf)) .^2).^-1.0))
threebodyca(k0, kinf) = :($k0 ./ (1 .+ $k0 ./ ($kinf ./ M)).*0.6 .^ ((1 .+ (log10.($k0 ./ ($kinf .* M))) .^2).^-1.0))

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
             [[:CO2], [:C, :O, :O], :JCO2toCpOpO], 
             [[:CO2], [:C, :O2], :JCO2toCpO2],
             [[:CO], [:C, :O], :JCOtoCpO], 
             [[:N2O], [:N2, :O1D], :JN2OtoN2pO1D],
             [[:NO2], [:NO, :O], :JNO2toNOpO],
             [[:NO], [:N, :O], :JNOtoNpO],
             
             # NEW: photoionization from Roger's model
             [[:CO2], [:CO2pl], :JCO2toCO2pl],  # Nair minimal ionosphere
             [[:CO2], [:Opl, :CO], :JCO2toOplpCO], # Nair minimal ionosphere
             [[:CO2], [:CO2plpl], :JCO2toCO2plpl],  
             [[:CO2], [:Cplpl, :O2], :JCO2toCplplpO2],
             [[:CO2], [:Cpl, :O2], :JCO2toCplpO2],
             [[:CO2], [:COpl, :Opl], :JCO2toCOplpOpl],
             [[:CO2], [:COpl, :O], :JCO2toCOplpO],             
             [[:CO2], [:Opl, :Cpl, :O], :JCO2toOplpCplpO],
             [[:CO], [:COpl], :JCOtoCOpl],
             [[:CO], [:C, :Opl], :JCOtoCpOpl],
             [[:CO], [:O, :Cpl], :JCOtoOpCpl],
             [[:H2O], [:H2Opl], :JH2OtoH2Opl],
             [[:HDO], [:HDOpl], :JHDOtoHDOpl], # NEW
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
             [[:HD], [:HDpl], :JHDtoHDpl], # NEW
             [[:H2], [:Hpl, :H], :JH2toHplpH],
             [[:H2O2], [:H2O2pl], :JH2O2toH2O2pl], 
             [[:O], [:Opl], :JOtoOpl],   # Nair minimal ionosphere
             [[:O2], [:O2pl], :JO2toO2pl],   # Nair minimal ionosphere
             [[:O3], [:O3pl], :JO3toO3pl],

             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8 .* 3.0e-33 .* (300 ./ Tn) .^ 3.25)], # Deighan 2012 # Checked no dups 
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35 .* exp.(724 ./ Tn))], # Checked no dups 
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5 .* 6.0e-34 .* (300 ./ Tn) .^ 2.4)], # Burkholder2020 # Checked no dups 
             [[:O, :O3], [:O2, :O2], :(8.0e-12 .* exp.(-2060 ./ Tn))],  # Burkholder 2020
             # [[:O, :CO, :M], [:CO2, :M], :(2.2e-33 .* exp.(-1780 ./ Tn))],  # DUPLICATE - commented out in favor of Roger's reaction

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.3e-11 .* exp.(55 ./ Tn))], # Burkholder 2020 (upd. 31 Dec 2020)
             [[:O1D, :O3], [:O2, :O2], :(2.4e-10)], # Burkholder 2020
             [[:O1D, :O3], [:O, :O, :O2], :(2.4e-10)], # Burkholder 2020
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11 .* exp.(115 ./ Tn))], # Burkholder 2020
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Burkholder 2020
             [[:O1D, :HD], [:H, :OD], :(0.41 .* 1.2e-10)], # Yung88: rate 0.41 .* H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41 .* 1.2e-10)], # Yung88: rate 0.41 .* H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10 .* exp.(60 ./ Tn))], # Burkholder 2020
             [[:O1D, :HDO], [:OD, :OH], :(((19 ./ 18) .^ -0.5) .* 1.63e-10 .* exp.(60 ./ Tn))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12 .* exp.(-4000 ./ Tn))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12 .* exp.(-4390 ./ Tn))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12 .* exp.(-4400 ./ Tn))], # NIST 
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11 .* exp.(-4038 ./ Tn))], # TODO:CHECK rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11 .* exp.(-3821 ./ Tn))], # TODO:CHECK NIST (1986, 200-300K): 8.19e-13 .* exp.(-2700/Tn)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12 .* exp.(-1800 ./ Tn))], # Burkholder 2020
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.) .* 2.8e-12 .* exp.(-1800 ./ Tn))], # Yung88: rate (3/20) .* H-ana. Sander2011: 5e-12 .* exp.(-2130 ./ Tn)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.) .* 2.8e-12 .* exp.(-1800 ./ Tn))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12 .* exp.(-1800 ./ Tn))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)

             # recombination of H. Use EITHER the first line OR the 2nd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32 .* (298 ./ Tn) .^ 2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32 .* (298 ./ Tn) .^ 2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(((2 ./ 1) .^ -0.5) .* 1.6e-32 .* (298 ./ Tn) .^ 2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9 .* 6.8e-31 .* (300 ./ Tn) .^ 2)], # NIST database
             [[:H, :OD, :CO2], [:HDO, :CO2], :(((18 ./ 17) .^ -0.5) .* 1.9 .* 6.8e-31 .* (300 ./ Tn) .^ 2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(((2 ./ 1) .^ -0.5) .* 1.9 .* 6.8e-31 .* (300 ./ Tn) .^ 2)], # not in Yung88. assumed rate

             ## H + HO2
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Burkholder 2020
             [[:H, :HO2], [:H2, :O2], :(0.5 .* 6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # Burkholder 2020; they do not have O1D though
             [[:H, :DO2], [:OH, :OD], :(((34 ./ 33) .^ -0.5) .* 7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HD, :O2], :(((34 ./ 33) .^ -0.5) .* 0.5 .* 6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HDO, :O1D], :(((34 ./ 33) .^ -0.5) .* 1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:D, :HO2], [:OH, :OD], :(0.71 .* 7.2e-11)], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:D, :HO2], [:HD, :O2], :(0.71 .* 0.5 .* 6.9e-12)], # Yung88: rate 0.71 .* H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:D, :HO2], [:HDO, :O1D], :(0.71 .* 1.6e-12)], # Yung88: rate 0.71 .* H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HO2, :D], :(1e-10 ./ (0.54 .* exp.(890 ./ Tn)))], # Yung88 (assumed) .- turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 .- turn off for Case 2

             ## H + H2O2
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
             [[:D, :HDO2], [:OD, :HDO], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # TODO: Mass scaling for doubly deuterated? I previously added this and assumed rate.
             [[:D, :HDO2], [:OH, :D2O], :(0.5 .* 1.16e-11 .* exp.(-2110 ./ Tn))], # see previous line

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0 .* 4.4e-32 .* (Tn ./ 300.) .^ -1.3), # Sander2011, 300K+. Yung89: 5.5e-32(Tn ./ 300) .^ -1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11 .* (Tn ./ 300.) .^ 0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0 .* 4.4e-32 .* (Tn ./ 300.) .^ -1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11 .* (Tn ./ 300.) .^ 0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10 .* exp.(-470 ./ Tn))], # Burkholder 2020
             [[:D, :O3], [:OD, :O2], :(0.71 .* 1.4e-10 .* exp.(-470 ./ Tn))], # Yung88, Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11 .* exp.(180 ./ Tn))], # Burkholder 2020 Yung89: 2.2e-11 .* exp.(120/Tn) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(((18 ./ 17) .^ -0.5) .* 1.8e-11 .* exp.(180 ./ Tn))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11 .* exp.(200 ./ Tn))], # Burkholder 2020
             [[:O, :DO2], [:OD, :O2], :(((34 ./ 33) .^ -0.5) .* 3.0e-11 .* exp.(200 ./ Tn))], # Yung88: rate same as H-ana. 
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12 .* exp.(-2000 ./ Tn))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5 .* ((35 ./ 34) .^ -0.5) .* 1.4e-12 .* exp.(-2000 ./ Tn))], # Yung88, Yung89
             [[:O, :HDO2], [:OH, :DO2], :(0.5 .* ((35 ./ 34) .^ -0.5) .* 1.4e-12 .* exp.(-2000 ./ Tn))], # Modifier added by me
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(1.8e-12)], # Burkholder 2020
             [[:OD, :OH], [:HDO, :O], :(((18 ./ 17) .^ -0.5) .* 1.8e-12)], 
             [[:OH, :OH], [:H2O2], threebody(:(1.3 .* 6.9e-31 .* (Tn ./ 300.) .^ -1.0), :(2.6e-11))], # Burkholder2020. Why 1.3? Reviewer suggestion to adjust for mars?
             [[:OD, :OH], [:HDO2], threebody(:(1.3 .* 6.9e-31 .* (Tn ./ 300.) .^ -1.0), :(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12 .* exp.(-940 ./ Tn))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(((18 ./ 17) .^ -0.5) .* 1.7e-12 .* exp.(-940 ./ Tn))],
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11 .* exp.(250 ./ Tn))], # verified NIST 4/3/18. Yung89: 4.6e-11 .* exp.(230/Tn) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(((34 ./ 33) .^ -0.5) .* 4.8e-11 .* exp.(250 ./ Tn))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(((18 ./ 17) .^ -0.5) .* 4.8e-11 .* exp.(250 ./ Tn))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12 .* exp.(-160 ./ Tn))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12 .* exp.(-200/Tn). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(((18 ./ 17) .^ -0.5) .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5 .* ((35 ./ 34) .^ -0.5) .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: rate 0.5 .* H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5 .* ((35 ./ 34) .^ -0.5) .* 2.9e-12 .* exp.(-160 ./ Tn))], # Yung88: rate 0.5 .* H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14 .* exp.(-490 ./ Tn))], # Sander2011. Yung89: 1.1e-14 .* exp.(-500/Tn). KIDA 250-340K: 2.03e-16 .* (Tn ./ 300) .^ 4.57 .* exp.(693/Tn). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(((34 ./ 33) .^ -0.5) .* 1.0e-14 .* exp.(-490 ./ Tn))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13 .* exp.(460 ./ Tn))], # Burkholder2020. Yung89: 2.3e-13 .* exp.(600/Tn). KIDA 230-420K: 2.2e-13 .* exp.(600/Tn)
             [[:DO2, :HO2], [:HDO2, :O2], :(((34 ./ 33) .^ -0.5) .* 3.0e-13 .* exp.(460 ./ Tn))], # Yung88: same as H-ana (assumed)
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2 .* 2.1e-33 .* exp.(920 ./ Tn))], # Burkholder2020.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(((34 ./ 33) .^ -0.5) .* 2 .* 2.1e-33 .* exp.(920 ./ Tn))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9 .* (Tn .^ -0.63) ./ (0.72 .* exp.(717 ./ Tn)))], # rate: Yung88. NIST (Howard82): 5.25E-11 .* (Tn/298) .^ -0.63  .- turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9 .* Tn .^ -0.63)], # Yung88  .- turn off for Case 2

             # CO2 recombination due to odd H
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13 .* (Tn ./ 300.) .^ 0.6), :(2.1e9 .* (Tn ./ 300.) .^ 6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13 .* (Tn ./ 300.) .^ 0.6), :(2.1e9 .* (Tn ./ 300.) .^ 6.1))], # Yung88: same as H-ana.
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33 .* (Tn ./ 300.) .^ -1.4), :(1.1e-12 .* (Tn ./ 300.) .^ 1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33 .* (Tn ./ 300.) .^ -1.4), :(1.1e-12 .* (Tn ./ 300.) .^ 1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(((46 ./ 45) .^ -0.5) .* 2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 ./ Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5) .* 8.7e-10)],

             # NEW .- Neutral reactions from Roger Yelle
             # Type 1
             [[:O1D], [:O], :(5.10e-3)],

             # Type 2
             # [[:C, :C], [:C2], :(2.16e-11)],  # UMIST: 4.36e-18 .* ((Tn ./ 300) .^ 0.35) .* exp.(-161.30 ./ Tn)
             [[:C, :H], [:CH], :(1.00e-17)],  # UMIST agrees
             [[:C, :N], [:CN], :((6.93e-20 .* Tn .^ 0.37) .* exp.(-51.0 ./ Tn))], # KIDA: :(7.76e-19 .* ((Tn ./ 300) .^ 0.14) .* exp.(-0.18 ./ Tn) ) for T=40-400K
             [[:CH, :C], [:C2, :H], :(6.59e-11)], # KIDA: 2.4e-10
             [[:CH, :H], [:H2, :C], :(1.31e-10 .* exp.(-80.0 ./ Tn))], # KIDA: :(1.24e-10 .* (Tn ./ 300) .^ 0.26)
             [[:CH, :H2], [:CH2, :H], :(2.90e-10 .* exp.(-1670.0 ./ Tn))], # KIDA agrees
             [[:CH, :H2], [:CH3], :((2.92e-16 .* Tn .^ -0.71) .* exp.(-11.6 ./ Tn))],  # KIDA: :(3.25e-17 .* (Tn ./ 300) .^ -0.6)
             [[:CH, :N], [:CN, :H], :(2.77e-10 .* Tn .^ -0.09)], # KIDA: :(1.4e-10 .8 (Tn ./ 300) .^ 0.41)
             [[:CH, :O], [:CO, :H], :(6.60e-11)], # KIDA: 7e-11
             [[:CH, :O], [:HCOpl], :(4.20e-13 .* exp.(-850.0 ./ Tn))], # KIDA: :(2e-11 .* (Tn ./ 300) .^ 0.44)
             [[:CH, :O], [:OH, :C], :(2.52e-11 .* exp.(-2381.0 ./ Tn))],  # KIDA agrees
             [[:CN, :H2], [:HCN, :H], :((1.80e-19 .* Tn .^ 2.6) .* exp.(-960.0 ./ Tn))],  # KIDA: :(4.96e-13 .* ((Tn ./ 300) .^ 2.6) .* exp.(-960  ./ Tn) )
             [[:CN, :N], [:N2, :C], :(9.80e-10 .* Tn .^ -0.4)],  # KIDA: :(8.8e-11 .* (Tn ./ 300) .^ 0.42)
             [[:CN, :NH], [:HCN, :N], :((1.70e-13 .* Tn .^ 0.5) .* exp.(-1000.0 ./ Tn))], # KIDA: :(2.94e-12 .* (Tn ./ 300) .^ 0.5 .* exp.(-1000 ./ Tn))
             [[:CN, :O], [:CO, :N], :(5.00e-11 .* exp.(-200.0 ./ Tn))], # KIDA: just 5e-11 
             [[:CN, :O], [:NO, :C], :(5.37e-11 .* exp.(-13800.0 ./ Tn))],  # KIDA: :(3.81e-11 .* (Tn ./ 300) .^ 0.5 .* exp.(-14500 ./ Tn))
             [[:CO2, :CH], [:HCO, :CO], :((1.70e-14 .* Tn .^ 0.5) .* exp.(-3000.0 ./ Tn))],  # KIDA: :(2.94e-13 .* (Tn ./ 300) .^ 0.5 .* exp.(3000 ./ Tn))
             [[:CO2, :H], [:CO, :OH], :(3.38e-10 .* exp.(-13163.0 ./ Tn))],  # KIDA: :(2.51e-10 .* exp.(-13300 ./ Tn))
             [[:CO2, :N], [:NO, :CO], :(3.20e-13 .* exp.(-1710.0 ./ Tn))],  # KIDA agrees 
             [[:CO2, :O], [:O2, :CO], :(2.46e-11 .* exp.(-26567.0 ./ Tn))], # KIDA agrees
             [[:H2, :C], [:CH, :H], :(6.64e-10 .* exp.(-11700.0 ./ Tn))],  # KIDA agrees
             [[:H2O, :H], [:OH, :H2], :((1.69e-14 .* Tn .^ 1.2) .* exp.(-9610.0 ./ Tn))], # KIDA :(6.82e-12 .* (Tn ./ 300) .^ 1.6 .* exp.(-9720 ./ Tn))
             [[:H2O, :O], [:OH, :OH], :((8.20e-14 .* Tn .^ 0.95) .* exp.(-8571.0 ./ Tn))],  # KIDA :(1.85e-11 .* (Tn ./ 300) .^ 0.95 .* exp.(-8570 ./ Tn))
             [[:HCN, :CN], [:NCCN, :H], :(1.45e-21 .* (Tn .^ 1.71) .* exp.(-770.0 ./ Tn))],
             [[:HCN, :H], [:CN, :H2], :(6.2e-10 .* exp.(-12500.0 ./ Tn))],
             [[:HCN, :O], [:CO, :NH], :(1.09e-15 .* (Tn .^ 1.14) .* exp.(-3742.0 ./ Tn))],
             [[:HCN, :O], [:NCO, :H], :(5.19e-16 .* (Tn .^ 1.38) .* exp.(-3693.0 ./ Tn))],
             [[:HCN, :O], [:OH, :CN], :(6.21e-10 .* exp.(-12439.0 ./ Tn))],
             [[:HCN, :OH], [:H2O, :CN], :(3.6e-17 .* (Tn .^ 1.5) .* exp.(-3887.0 ./ Tn))],
             [[:HCO, :C], [:C2O, :H], :(1.0e-10)],
             [[:HCO, :C], [:CO, :CH], :(1.0e-10)],
             [[:HCO, :CH], [:CH2, :CO], :(5.3e-14 .* (Tn .^ 0.7) .* exp.(-500.0 ./ Tn))],
             [[:HCO, :CN], [:HCN, :CO], :(1.0e-10)],

             [[:HCO, :H], [:CO, :H2], :(1.5e-10)],
             # NEW
             [[:DCO, :H], [:CO, :HD], :(((30 ./ 29) .^ -0.5) .* 1.5e-10)],  
             [[:HCO, :D], [:CO, :HD], :(((2 ./ 1) .^ -0.5) .* 1.5e-10)], 
             # [[:DCO, :D], [:CO, :HD], :(((2 ./ 1) .^ -0.5) .* 1.5e-10)],  # TODO: Doubly deuterated mass scaling factor?

             [[:HCO, :HCO], [:CO, :CO, :H2], :(7.35e-12)],
             [[:HCO, :HCO], [:H2CO, :CO], :(4.26e-11)],
             [[:HCO, :N], [:CO, :NH], :(3.3e-13 .* (Tn .^ 0.5) .* exp.(-1000.0 ./ Tn))],
             [[:HCO, :N], [:HCN, :O], :(1.7e-10)],
             [[:HCO, :N], [:NCO, :H], :(1.0e-10)],
             [[:HCO, :NO], [:HNO, :CO], :(1.2e-11)],

             [[:HCO, :O], [:CO, :OH], :(5.0e-11)],
             # NEW
             [[:DCO, :O], [:CO, :OD], :(((30 ./ 29) .^ -0.5) .* 5.0e-11)],

             [[:HCO, :O], [:CO2, :H], :(5.0e-11)],
             # NEW
             [[:DCO, :O], [:CO2, :D], :(((30 ./ 29) .^ -0.5) .* 5.0e-11)],

             [[:HCO, :O2], [:CO2, :OH], :(7.6e-13)],
             # NEW 
             [[:DCO, :O2], [:CO2, :OD], :(((30 ./ 29) .^ -0.5) .* 7.6e-13)],

             [[:HCO, :O2], [:HO2, :CO], :(5.2e-12)],
             # NEW
             [[:DCO, :O2], [:DO2, :CO], :(((30 ./ 29) .^ -0.5) .* 5.2e-12)],

             [[:HCO, :OH], [:H2O, :CO], :(1.8e-10)],
             # NEW
             [[:DCO, :OH], [:HDO, :CO], :(0.5 .* ((30 ./ 29) .^ -0.5) .* 1.8e-10)],
             [[:HCO, :OD], [:HDO, :CO], :(0.5 .* ((30 ./ 29) .^ -0.5) .* 1.8e-10)],

             [[:HNO, :CH], [:CH2, :NO], :(1.73e-11)],
             [[:HNO, :CN], [:HCN, :NO], :(3.0e-11)],
             [[:HNO, :CO], [:CO2, :NH], :(3.32e-12 .* exp.(-6170.0 ./ Tn))],
             [[:HNO, :H], [:NH2, :O], :(5.81e-9 .* (Tn .^ -0.3) .* exp.(-14730.0 ./ Tn))],
             [[:HNO, :H], [:NO, :H2], :(7.41e-13 .* (Tn .^ 0.72) .* exp.(-329.0 ./ Tn))],
             [[:HNO, :HCO], [:H2CO, :NO], :(1.0e-12 .* exp.(-1000.0 ./ Tn))],
             [[:HNO, :N], [:N2O, :H], :(8.26e-14 .* (Tn .^ 0.5) .* exp.(-1500.0 ./ Tn))],
             [[:HNO, :N], [:NO, :NH], :(1.7e-13 .* (Tn .^ 0.5) .* exp.(-1000.0 ./ Tn))],
             [[:HNO, :O], [:NO, :OH], :(6.0e-11 .* (Tn .^ -0.08))],
             [[:HNO, :O], [:NO2, :H], :(1.0e-12)],
             [[:HNO, :O], [:O2, :NH], :(1.7e-13 .* (Tn .^ 0.5) .* exp.(-3500.0 ./ Tn))],
             [[:HNO, :OH], [:H2O, :NO], :(5.54e-15 .* (Tn .^ 1.23) .* exp.(44.3 ./ Tn))],
             [[:HO2, :CH], [:CH2, :O2], :(1.7e-14 .* (Tn .^ 0.5) .* exp.(-7550.0 ./ Tn))],
             [[:HO2, :CH], [:HCO, :OH], :(8.31e-13 .* (Tn .^ 0.5) .* exp.(-3000.0 ./ Tn))],
             [[:HO2, :CO], [:CO2, :OH], :(5.6e-10 .* exp.(-12160.0 ./ Tn))],
             [[:HO2, :H2], [:H2O2, :H], :(5.0e-11 .* exp.(-13110.0 ./ Tn))],
             [[:HO2, :H2O], [:H2O2, :OH], :(4.65e-11 .* exp.(-16500.0 ./ Tn))],

             [[:HO2, :HCO], [:H2CO, :O2], :(5.0e-11)],
             # NEW
             [[:DO2, :HCO], [:HDCO, :O2], :(((34 ./ 33) .^ -0.5) .* 5.0e-11)],
             [[:HO2, :DCO], [:HDCO, :O2], :(((30 ./ 29) .^ -0.5) .* 5.0e-11)],

             [[:HO2, :N], [:NO, :OH], :(2.2e-11)],
             # NEW
             [[:DO2, :N], [:NO, :OD], :(((34 ./ 33) .^ -0.5) .* 2.2e-11)],

             [[:HO2, :N], [:O2, :NH], :(1.7e-13)],

             [[:HO2, :NO], [:NO2, :OH], :(3.3e-12 .* exp.(270.0 ./ Tn))],
             # NEW 
             [[:DO2, :NO], [:NO2, :OD], :(((34 ./ 33) .^ -0.5) .* 3.3e-12 .* exp.(270.0 ./ Tn))],

             [[:HOCO, :OH], [:CO2, :H2O], :(1.03e-11)],
             # NEW
             [[:DOCO, :OH], [:CO2, :HDO], :(((46 ./ 45) .^ -0.5) .* 1.03e-11)],
             [[:HOCO, :OD], [:CO2, :HDO], :(((18 ./ 17) .^ -0.5) .* 1.03e-11)],
             # [[:DOCO, :OD], [:CO2, :HDO], :(1.03e-11)], # TODO: Doubly deuterated mass scaling factor?

             # [[:N, :C], [:CN], :(3.49e-19 .* (Tn .^ 0.14) .* exp.(-0.18 ./ Tn))], # ROGER'S DUPLICATE 
             [[:N2O, :CO], [:CO2, :N2], :(1.62e-13 .* exp.(-8780.0 ./ Tn))],
             [[:N2O, :H], [:NO, :NH], :(0.111 .* (Tn .^ -2.16) .* exp.(-18700.0 ./ Tn))],
             [[:N2O, :H], [:OH, :N2], :(8.08e-22 .* (Tn .^ 3.15) .* exp.(-3603.0 ./ Tn))],
             [[:N2O, :NO], [:NO2, :N2], :(8.74e-19 .* (Tn .^ 2.23) .* exp.(-23292.0 ./ Tn))],
             [[:N2O, :O], [:NO, :NO], :(1.15e-10 .* exp.(-13400.0 ./ Tn))],
             [[:N2O, :O], [:O2, :N2], :(1.66e-10 .* exp.(-14100.0 ./ Tn))],
             [[:N2O, :OH], [:HO2, :N2], :(3.7e-13 .* exp.(-2740.0 ./ Tn))],
             [[:NH, :C], [:CH, :N], :(9.99e-13 .* (Tn .^ 0.5) .* exp.(-4000.0 ./ Tn))],
             [[:NH, :C], [:CN, :H], :(1.2e-10)],
             [[:NH, :H], [:H2, :N], :(9.99e-13 .* (Tn .^ 0.5) .* exp.(-2400.0 ./ Tn))],
             [[:NH, :N], [:N2, :H], :(4.98e-11)],
             # [[:NH, :NH], [:N2, :H, :H], :(1.16e-9)], # Quadratic NH 
             # [[:NH, :NH], [:N2, :H2], :(1.7e-11)],
             # [[:NH, :NH], [:NH2, :N], :(6.29e-18 .* (Tn .^ 1.8) .* exp.(70.0 ./ Tn))],
             [[:NH, :O], [:NO, :H], :(1.8e-10 .* exp.(-300.0 ./ Tn))],
             [[:NH, :O], [:OH, :N], :(1.16e-11)],
             [[:NH2, :C], [:HCN, :H], :(5.77e-11 .* (Tn .^ -0.1) .* exp.(9.0 ./ Tn))],
             [[:NH2, :C], [:HNC, :H], :(5.77e-11 .* (Tn .^ -0.1) .* exp.(9.0 ./ Tn))],
             [[:NH2, :H], [:NH, :H2], :(1.36e-14 .* (Tn .^ 1.02) .* exp.(-2161.0 ./ Tn))],
             [[:NH2, :H2], [:NH3, :H], :(4.74e-25 .* (Tn .^ 3.89) .* exp.(-1400.0 ./ Tn))],
             [[:NH2, :NO], [:H2O, :N2], :(3.6e-12 .* exp.(450.0 ./ Tn))],
             [[:NH2, :O], [:HNO, :H], :(1.11e-10 .* (Tn .^ -0.1))],
             [[:NH2, :O], [:OH, :NH], :(1.24e-11 .* (Tn .^ -0.1))],
             [[:NH2, :OH], [:H2O, :NH], :(1.08e-15 .* (Tn .^ 1.25) .* exp.(43.5 ./ Tn))],
             [[:NH2, :OH], [:NH3, :O], :(2.73e-15 .* (Tn .^ 0.76) .* exp.(-262.0 ./ Tn))],
             [[:NO, :C], [:CN, :O], :(1.49e-10 .* (Tn .^ -0.16))],
             [[:NO, :C], [:CO, :N], :(2.24e-10 .* (Tn .^ -0.16))],
             [[:NO, :CH], [:CO, :NH], :(1.52e-11)],
             [[:NO, :CH], [:HCN, :O], :(1.31e-10)],
             [[:NO, :CH], [:HCO, :N], :(1.14e-11)],
             [[:NO, :CH], [:NCO, :H], :(2.47e-11)],
             [[:NO, :CH], [:OH, :CN], :(1.9e-12)],
             [[:NO, :CN], [:CO, :N2], :(1.6e-13)],
             [[:NO, :H], [:NH, :O], :(1.64e-9 .* (Tn .^ -0.1) .* exp.(-35220.0 ./ Tn))],
             [[:NO, :N], [:N2, :O], :(2.1e-11 .* exp.(100.0 ./ Tn))],
             [[:NO, :NH], [:N2, :O, :H], :(7.4e-10 .* exp.(-10540.0 ./ Tn))],
             [[:NO, :NH], [:N2O, :H], :(4.56e-9 .* (Tn .^ -0.78) .* exp.(-40.0 ./ Tn))],
             [[:NO, :NH], [:OH, :N2], :(1.14e-9 .* (Tn .^ -0.78) .* exp.(-40.0 ./ Tn))],
             [[:NO, :O], [:O2, :N], :(1.18e-11 .* exp.(-20413.0 ./ Tn))],
             [[:NO2, :CN], [:CO2, :N2], :(6.12e-11 .* (Tn .^ -0.752) .* exp.(-173.0 ./ Tn))],
             [[:NO2, :CN], [:N2O, :CO], :(8.16e-11 .* (Tn .^ -0.752) .* exp.(-173.0 ./ Tn))],
             [[:NO2, :CN], [:NCO, :NO], :(8.77e-10 .* (Tn .^ -0.752) .* exp.(-173.0 ./ Tn))],
             [[:NO2, :CO], [:CO2, :NO], :(1.48e-10 .* exp.(-17000.0 ./ Tn))],

             [[:NO2, :H], [:NO, :OH], :(4.0e-10 .* exp.(-340.0 ./ Tn))],
             # NEW 
             [[:NO2, :D], [:NO, :OD], :(((2 ./ 1) .^ -0.5) .* 4.0e-10 .* exp.(-340.0 ./ Tn))],

             [[:NO2, :N], [:N2, :O, :O], :(2.41e-12)],
             [[:NO2, :N], [:N2O, :O], :(5.8e-12 .* exp.(220.0 ./ Tn))],
             [[:NO2, :N], [:NO, :NO], :(1.0e-12)],
             [[:NO2, :N], [:O2, :N2], :(1.0e-12)],
             [[:NO2, :NH], [:HNO, :NO], :(1.56e-6 .* (Tn .^ -1.94) .* exp.(-56.9 ./ Tn))],
             [[:NO2, :NH], [:N2O, :OH], :(1.09e-6 .* (Tn .^ -1.94) .* exp.(-56.9 ./ Tn))],
             [[:NO2, :NH2], [:N2O, :H2O], :(2.1e-12 .* exp.(650.0 ./ Tn))],
             [[:NO2, :O], [:O2, :NO], :(5.1e-12 .* exp.(210.0 ./ Tn))],
             [[:O, :C], [:CO], :(1.75e-19 .* (Tn .^ 0.705) .* exp.(-136.0 ./ Tn))],

             [[:O, :H], [:OH], :(8.65e-18 .* (Tn .^ -0.38))],
             # NEW 
             [[:O, :D], [:OD], :(((2 ./ 1) .^ -0.5) .* 8.65e-18 .* (Tn .^ -0.38))],

             [[:O1D, :CO], [:CO, :O], :(4.7e-11 .* exp.(63.0 ./ Tn))],
             [[:O1D, :CO], [:CO2], :(8.0e-11)],
             [[:O1D, :H2O2], [:H2O2, :O], :(5.2e-10)],
             [[:O1D, :N2], [:N2, :O], :(2.15e-11 .* exp.(110.0 ./ Tn))],
             [[:O1D, :N2O], [:NO, :NO], :(7.26e-11 .* exp.(20.0 ./ Tn))],
             [[:O1D, :N2O], [:O2, :N2], :(4.64e-11 .* exp.(20.0 ./ Tn))],
             [[:O1D, :NO], [:NO, :O], :(4.0e-11)],
             [[:O1D, :NO2], [:NO2, :O], :(1.13e-10 .* exp.(115.0 ./ Tn))],
             [[:O1D, :NO2], [:O2, :NO], :(2.31e-10)],
             [[:O2, :C], [:CO, :O], :(3.03e-10 .* (Tn .^ -0.32))],
             [[:O2, :CH], [:CO, :O, :H], :(1.2e-11)],
             [[:O2, :CH], [:CO, :OH], :(8.0e-12)],
             [[:O2, :CH], [:CO2, :H], :(1.2e-11)],
             [[:O2, :CH], [:HCO, :O], :(8.0e-12)],
             [[:O2, :CN], [:CO, :NO], :(3.0e-12 .* exp.(210.0 ./ Tn))],
             [[:O2, :CN], [:NCO, :O], :(5.97e-11 .* (Tn .^ -0.19) .* exp.(31.9 ./ Tn))],
             [[:O2, :CO], [:CO2, :O], :(5.99e-12 .* exp.(-24075.0 ./ Tn))],
             [[:O2, :H], [:OH, :O], :(2.61e-10 .* exp.(-8156.0 ./ Tn))],
             [[:O2, :H2], [:HO2, :H], :(2.4e-10 .* exp.(-28500.0 ./ Tn))],
             [[:O2, :H2], [:OH, :OH], :(3.16e-10 .* exp.(-21890.0 ./ Tn))],
             [[:O2, :N], [:NO, :O], :(1.5e-11 .* exp.(-3600.0 ./ Tn))],
             [[:O2, :NH], [:HNO, :O], :(4.0e-11 .* exp.(-6970.0 ./ Tn))],
             [[:O2, :NH], [:NO, :OH], :(1.5e-13 .* exp.(-770.0 ./ Tn))],
             [[:O3, :N], [:O2, :NO], :(1.0e-16)],
             [[:O3, :NO], [:NO2, :O2], :(3.0e-12 .* exp.(-1500.0 ./ Tn))],
             [[:O3, :NO2], [:NO3, :O2], :(1.2e-13 .* exp.(-2450.0 ./ Tn))],
             [[:OH, :C], [:CO, :H], :(7.98e-10 .* (Tn .^ -0.34) .* exp.(-0.108 ./ Tn))],
             [[:OH, :CH], [:HCO, :H], :(8.31e-13 .* (Tn .^ 0.5) .* exp.(-5000.0 ./ Tn))],
             [[:OH, :H], [:H2, :O], :(8.1e-21 .* (Tn .^ 2.8) .* exp.(-1950.0 ./ Tn))],
             [[:OH, :N], [:NH, :O], :(1.06e-11 .* (Tn .^ 0.1) .* exp.(-10700.0 ./ Tn))],
             [[:OH, :N], [:NO, :H], :(1.8e-10 .* (Tn .^ -0.2))],
             [[:OH, :NH], [:H2O, :N], :(3.3e-15 .* (Tn .^ 1.2))],
             [[:OH, :NH], [:HNO, :H], :(3.3e-11)],
             [[:OH, :NH], [:NH2, :O], :(1.66e-12 .* (Tn .^ 0.1) .* exp.(-5800.0 ./ Tn))],
             [[:OH, :NH], [:NO, :H2], :(4.16e-11)],

             # Type 4
             # simpler type, when Troe parameter is 0. Updated 5 Feb 2021 to match Roger's code
             [[:N, :H], [:NH], :(0.0 .+ (5.0e-32 .* 1.0 .* M) ./ (5.0e-32 .* M .+ 1.0))],  
             [[:CO, :H], [:HCO], :(0.0 .+ (2.0e-35 .* (Tn .^ 0.2) .* 1.0 .* (Tn .^ 0.2) .* M) ./ (2.0e-35 .* (Tn .^ 0.2) .* M .+ 1.0 .* (Tn .^ 0.2)))], 
             # NEW TODO: CO + H is complicated, need to treat this one individually to get D rxn.
             [[:CO, :D], [:DCO], :(0.0 .+ (2.0e-35 .* (Tn .^ 0.2) .* 1.0 .* (Tn .^ 0.2) .* M) ./ (2.0e-35 .* (Tn .^ 0.2) .* M .+ 1.0 .* (Tn .^ 0.2)))], 

             [[:O1D, :N2], [:N2O], :(0.0 .+ (4.75e-34 .* (Tn .^ -0.9) .* 1.0 .* (Tn .^ -0.9) .* M) ./ (4.75e-34 .* (Tn .^ -0.9) .* M .+ 1.0 .* (Tn .^ -0.9)))], 
                        
             # More complicated - need the minimum of the two expressions. These updated 5 Feb 2021 to match Roger's code.
             # ROGER'S DUPLICATE #[[:CH, :H2], [:CH2, :H], :(min.($:(8.5e-11 .* (Tn .^ 0.15)), $:(0.0 .+ (10 .^ ((log10.(0.6)) ./ (1 .+ ((log10.((4.7e-26 .* (Tn .^ -1.6) .* M) ./ (8.5e-11 .* (Tn .^ 0.15))) .- 0.4 .- 0.67 .* log10.(0.6)) ./ (0.75 .- 1.27 .* log10.(0.6) .- 0.14 .* (log10.((4.7e-26 .* (Tn .^ -1.6) .* M) ./ (8.5e-11 .* (Tn .^ 0.15))) .- 0.4 .- 0.67 .* log10.(0.6)))) .^ 2)) .* 4.7e-26 .* (Tn .^ -1.6) .* 8.5e-11 .* (Tn .^ 0.15) .* M) ./ (4.7e-26 .* (Tn .^ -1.6) .* M .+ 8.5e-11 .* (Tn .^ 0.15)))))],
             [[:CO, :O], [:CO2], :(min.($:(1.0 .* exp.(-1509.0 ./ Tn)), $:(0.0 .+ (10 .^ ((log10.(0.4)) ./ (1 .+ ((log10.((1.7e-33 .* exp.(-1509.0 ./ Tn) .* M) ./ (1.0 .* exp.(-1509.0 ./ Tn))) .- 0.4 .- 0.67 .* log10.(0.4)) ./ (0.75 .- 1.27 .* log10.(0.4) .- 0.14 .* (log10.((1.7e-33 .* exp.(-1509.0 ./ Tn) .* M) ./ (1.0 .* exp.(-1509.0 ./ Tn))) .- 0.4 .- 0.67 .* log10.(0.4)))) .^ 2)) .* 1.7e-33 .* exp.(-1509.0 ./ Tn) .* 1.0 .* exp.(-1509.0 ./ Tn) .* M) ./ (1.7e-33 .* exp.(-1509.0 ./ Tn) .* M .+ 1.0 .* exp.(-1509.0 ./ Tn)))))],
             
             [[:NO, :H], [:HNO], :(min.($:(2.53e-9 .* (Tn .^ -0.41)), $:(0.0 .+ (10 .^ ((log10.(0.82)) ./ (1 .+ ((log10.((9.56e-29 .* (Tn .^ -1.17) .* exp.(-212.0 ./ Tn) .* M) ./ (2.53e-9 .* (Tn .^ -0.41))) .- 0.4 .- 0.67 .* log10.(0.82)) ./ (0.75 .- 1.27 .* log10.(0.82) .- 0.14 .* (log10.((9.56e-29 .* (Tn .^ -1.17) .* exp.(-212.0 ./ Tn) .* M) ./ (2.53e-9 .* (Tn .^ -0.41))) .- 0.4 .- 0.67 .* log10.(0.82)))) .^ 2)) .* 9.56e-29 .* (Tn .^ -1.17) .* exp.(-212.0 ./ Tn) .* 2.53e-9 .* (Tn .^ -0.41) .* M) ./ (9.56e-29 .* (Tn .^ -1.17) .* exp.(-212.0 ./ Tn) .* M .+ 2.53e-9 .* (Tn .^ -0.41)))))],
             # NEW: complicated rate. requires special calculation. or assumption that it's the same.

             [[:NO, :O], [:NO2], :(min.($:(4.9e-10 .* (Tn .^ -0.4)), $:(0.0 .+ (10 .^ ((log10.(0.8)) ./ (1 .+ ((log10.((9.2e-28 .* (Tn .^ -1.6) .* M) ./ (4.9e-10 .* (Tn .^ -0.4))) .- 0.4 .- 0.67 .* log10.(0.8)) ./ (0.75 .- 1.27 .* log10.(0.8) .- 0.14 .* (log10.((9.2e-28 .* (Tn .^ -1.6) .* M) ./ (4.9e-10 .* (Tn .^ -0.4))) .- 0.4 .- 0.67 .* log10.(0.8)))) .^ 2)) .* 9.2e-28 .* (Tn .^ -1.6) .* 4.9e-10 .* (Tn .^ -0.4) .* M) ./ (9.2e-28 .* (Tn .^ -1.6) .* M .+ 4.9e-10 .* (Tn .^ -0.4)))))],
             [[:O, :N], [:NO], :(min.($:(1.0 .* Tn .^ 0), $:(0.0 .+ (10 .^ ((log10.(0.4)) ./ (1 .+ ((log10.((5.46e-33 .* exp.(155.0 ./ Tn) .* M) ./ (1.0)) .- 0.4 .- 0.67 .* log10.(0.4)) ./ (0.75 .- 1.27 .* log10.(0.4) .- 0.14 .* (log10.((5.46e-33 .* exp.(155.0 ./ Tn) .* M) ./ (1.0)) .- 0.4 .- 0.67 .* log10.(0.4)))) .^ 2)) .* 5.46e-33 .* exp.(155.0 ./ Tn) .* 1.0 .* M) ./ (5.46e-33 .* exp.(155.0 ./ Tn) .* M .+ 1.0))))],
             # Note: one of the values to compare for O+N->NO is just 1.0, but the other value to compare can be an array. To make 1.0 into an array or a float as needed,
             # we just multiply by T .^ 0, so when T is an array, so will be the value, and same for when its a float.

             # Type 5
             [[:NO, :OH], [:HONO], threebody(:(1.93e-24 .* (Tn .^ -2.6)), :(6.37e-11 .* (Tn .^ -0.1)))],
             [[:NO2, :HO2], [:HO2NO2], threebody(:(5.02e-23 .* (Tn .^ -3.4)), :(2.21e-11 .* (Tn .^ -0.3)))],
             [[:NO2, :O], [:NO3], threebody(:(7.19e-27 .* (Tn .^ -1.8)), :(1.19e-9 .* (Tn .^ -0.7)))],
             [[:NO2, :OH], [:HONO2], threebody(:(4.86e-23 .* (Tn .^ -3.0)), :(2.8e-11))],
             [[:NO2, :OH], [:HOONO], threebody(:(4.17e-22 .* (Tn .^ -3.9)), :(7.27e12 .* (Tn .^ -0.5)))],
                     
             # IONOSPHERE .- reactions from Roger Yelle. Updated 5 Feb 2021 to match Roger's code
             [[:ArHpl, :C], [:CHpl, :Ar], :(2 .* 1.02e-9)],

             [[:ArHpl, :CO], [:HCOpl, :Ar], :(2 .* 1.25e-9)],
             # NEW
             # [[:ArDpl, :CO], [:DCOpl, :Ar], :(((42 ./ 41) .^ -0.5) .* 2 .* 1.25e-9)], # DUPLICATE

             [[:ArHpl, :CO2], [:HCO2pl, :Ar], :(2 .* 1.1e-9)],
             [[:ArHpl, :H2], [:H3pl, :Ar], :(2 .* 6.3e-10)],
             [[:ArHpl, :N2], [:N2Hpl, :Ar], :(2 .* 8.0e-10)],
             [[:ArHpl, :O], [:OHpl, :Ar], :(2 .* 5.9e-10)],
             [[:ArHpl, :O2], [:HO2pl, :Ar], :(2 .* 5.05e-10)],
             [[:Arpl, :CO], [:COpl, :Ar], :(2 .* 4.4e-11)],
             [[:Arpl, :CO2], [:CO2pl, :Ar], :(2 .* 4.8e-10)],
             [[:Arpl, :H2], [:ArHpl, :H], :(2 .* 8.72e-10)],
             [[:Arpl, :H2], [:H2pl, :Ar], :(2 .* 1.78e-11)],
             [[:Arpl, :H2O], [:ArHpl, :OH], :(2 .* 3.24e-10)],
             [[:Arpl, :H2O], [:H2Opl, :Ar], :(2 .* 1.3e-9)],
             [[:Arpl, :N2], [:N2pl, :Ar], :(2 .* 1.1e-11)],
             [[:Arpl, :N2O], [:N2Opl, :Ar], :(2 .* 2.91e-10)],
             [[:Arpl, :N2O], [:N2pl, :Ar, :O], :(2 .* 3.0e-12)],
             [[:Arpl, :N2O], [:NOpl, :Ar, :N], :(2 .* 3.0e-12)],
             [[:Arpl, :N2O], [:Opl, :N2, :Ar], :(2 .* 3.0e-12)],
             [[:Arpl, :NO], [:NOpl, :Ar], :(2 .* 3.1e-10)],
             [[:Arpl, :NO2], [:NO2pl, :Ar], :(2 .* 2.76e-11)],
             [[:Arpl, :NO2], [:NOpl, :Ar, :O], :(2 .* 4.32e-10)],
             [[:Arpl, :O2], [:O2pl, :Ar], :(2 .* 4.6e-11)],
             [[:CHpl, :C], [:C2pl, :H], :(2 .* 1.2e-9)],
             [[:CHpl, :CN], [:C2Npl, :H], :(2 .* 9.53e-9 .* (Ti .^ -0.5))],
             [[:CHpl, :CO], [:HCOpl, :C], :(2 .* 7.0e-12 .* (Ti .^ -0.5))],
             [[:CHpl, :CO2], [:HCOpl, :CO], :(2 .* 1.6e-9)],
             [[:CHpl, :H], [:Cpl, :H2], :(2 .* 7.5e-10)],
             [[:CHpl, :H2], [:CH2pl, :H], :(2 .* 1.2e-9)],
             [[:CHpl, :H2O], [:H2COpl, :H], :(2 .* 1.0e-8 .* (Ti .^ -0.5))],
             [[:CHpl, :H2O], [:H3Opl, :C], :(2 .* 1.45e-9)],
             [[:CHpl, :H2O], [:HCOpl, :H2], :(2 .* 5.02e-8 .* (Ti .^ -0.5))],
             [[:CHpl, :HCN], [:C2Npl, :H2], :(2 .* 4.2e-10)],
             [[:CHpl, :HCN], [:HC2Npl, :H], :(2 .* 2.8e-10)],
             [[:CHpl, :HCN], [:HCNHpl, :C], :(2 .* 2.1e-9)],
             [[:CHpl, :HCO], [:CH2pl, :CO], :(2 .* 7.97e-9 .* (Ti .^ -0.5))],
             [[:CHpl, :HCO], [:HCOpl, :CH], :(2 .* 7.97e-9 .* (Ti .^ -0.5))],
             [[:CHpl, :N], [:CNpl, :H], :(2 .* 1.9e-10)],
             [[:CHpl, :NH], [:CNpl, :H2], :(2 .* 1.32e-8 .* (Ti .^ -0.5))],
             [[:CHpl, :NO], [:NOpl, :CH], :(2 .* 7.6e-10)],
             [[:CHpl, :O], [:COpl, :H], :(2 .* 3.5e-10)],
             [[:CHpl, :O2], [:HCOpl, :O], :(2 .* 8.73e-10)],
             [[:CHpl, :OH], [:COpl, :H2], :(2 .* 1.3e-8 .* (Ti .^ -0.5))],
             [[:CNpl, :C], [:Cpl, :CN], :(2 .* 1.1e-10)],
             [[:CNpl, :CH], [:CHpl, :CN], :(2 .* 1.11e-8 .* (Ti .^ -0.5))],
             [[:CNpl, :CO], [:COpl, :CN], :(2 .* 4.4e-10)],
             [[:CNpl, :CO2], [:C2Opl, :NO], :(2 .* 3.3e-10)],
             [[:CNpl, :CO2], [:CO2pl, :CN], :(2 .* 4.4e-10)],
             [[:CNpl, :CO2], [:OCNpl, :CO], :(2 .* 3.3e-10)],
             [[:CNpl, :H], [:Hpl, :CN], :(2 .* 6.4e-10)],
             [[:CNpl, :H2], [:HCNpl, :H], :(2 .* 8.0e-10)],
             [[:CNpl, :H2], [:HNCpl, :H], :(2 .* 8.0e-10)],
             [[:CNpl, :H2O], [:H2CNpl, :O], :(2 .* 4.8e-10)],
             [[:CNpl, :H2O], [:H2Opl, :CN], :(2 .* 3.2e-10)],
             [[:CNpl, :H2O], [:HCNpl, :OH], :(2 .* 1.6e-9)],
             [[:CNpl, :H2O], [:HCOpl, :NH], :(2 .* 1.6e-10)],
             [[:CNpl, :H2O], [:HNCOpl, :H], :(2 .* 6.4e-10)],
             [[:CNpl, :HCN], [:C2N2pl, :H], :(2 .* 4.59e-10)],
             [[:CNpl, :HCN], [:HCNpl, :CN], :(2 .* 2.24e-9)],
             [[:CNpl, :HCO], [:HCNpl, :CO], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:CNpl, :HCO], [:HCOpl, :CN], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:CNpl, :N], [:N2pl, :C], :(2 .* 6.1e-10)],
             [[:CNpl, :N2O], [:N2Opl, :CN], :(2 .* 4.56e-10)],
             [[:CNpl, :N2O], [:NOpl, :CN2], :(2 .* 1.52e-10)],
             [[:CNpl, :N2O], [:OCNpl, :N2], :(2 .* 1.52e-10)],
             [[:CNpl, :NH], [:NHpl, :CN], :(2 .* 1.13e-8 .* (Ti .^ -0.5))],
             [[:CNpl, :NH2], [:NH2pl, :CN], :(2 .* 1.58e-8 .* (Ti .^ -0.5))],
             [[:CNpl, :NO], [:NOpl, :CN], :(2 .* 5.7e-10)],
             [[:CNpl, :NO], [:OCNpl, :N], :(2 .* 1.9e-10)],
             [[:CNpl, :O], [:Opl, :CN], :(2 .* 6.5e-11)],
             [[:CNpl, :O2], [:NOpl, :CO], :(2 .* 8.6e-11)],
             [[:CNpl, :O2], [:O2pl, :CN], :(2 .* 2.58e-10)],
             [[:CNpl, :O2], [:OCNpl, :O], :(2 .* 8.6e-11)],
             [[:CNpl, :OH], [:OHpl, :CN], :(2 .* 1.11e-8 .* (Ti .^ -0.5))],

             [[:CO2pl, :H], [:HCOpl, :O], :(2 .* 4.47e-10)],
             [[:CO2pl, :H], [:Hpl, :CO2], :(2 .* 5.53e-11)],
             # NEW 
             # [[:CO2pl, :D], [:DCOpl, :O], :(((2 ./ 1) .^ -0.5) .* 2 .* 4.47e-10)], # DUPLICATE
             # [[:CO2pl, :D], [:Dpl, :CO2], :(((2 ./ 1) .^ -0.5) .* 2 .* 5.53e-11)], # DUPLICATE

             [[:CO2pl, :H2], [:HCO2pl, :H], :(4.7e-10)], # Nair minimal ionosphere. # Borodi2009: :(9.5e-10 .* ((Ti ./ 300) .^ -0.15)). Roger: :(2 .* 2.24e-9 .* (Ti .^ -0.15)). Nair: in use. BAD RATE? 2.24e-9 .* ((300 ./ Ti) .^ -0.15)
             # NEW
             [[:CO2pl, :HD], [:DCO2pl, :H], :(0.5 .* ((3 ./ 2) .^ -0.5) .* 4.7e-10)], # factor of 0.5 to account for two branches.
             [[:CO2pl, :HD], [:HCO2pl, :D], :(0.5 .* ((3 ./ 2) .^ -0.5) .* 4.7e-10)], # factor of 0.5 to account for two branches.

             [[:CO2pl, :H2O], [:H2Opl, :CO2], :(2 .* 1.8e-9)],
             [[:CO2pl, :H2O], [:HCO2pl, :OH], :(2 .* 6.0e-10)],
             # NEW 
             [[:CO2pl, :HDO], [:HDOpl, :CO2], :(((19 ./ 18) .^ -0.5) .* 2 .* 1.8e-9)],
             [[:CO2pl, :HDO], [:DCO2pl, :OH], :(0.5 .* ((19 ./ 18) .^ -0.5) .* 2 .* 6.0e-10)],
             [[:CO2pl, :HDO], [:HCO2pl, :OD], :(0.5 .* ((19 ./ 18) .^ -0.5) .* 2 .* 6.0e-10)],

             [[:CO2pl, :HCN], [:HCNpl, :CO2], :(2 .* 8.1e-10)],
             [[:CO2pl, :HCN], [:HCO2pl, :CN], :(2 .* 9.0e-11)],
             [[:CO2pl, :N], [:COpl, :NO], :(2 .* 3.4e-10)],
             [[:CO2pl, :NO], [:NOpl, :CO2], :(2 .* 1.23e-10)],
             [[:CO2pl, :O], [:O2pl, :CO], :(1.6e-10)], # Nair minimal ionosphere # Fehsenfeld1970: -in use-. Tenewitz 2018, sect 3b. 2e-11 * 0.98: :(1.96e-11). Roger: :(2 .* 1.6e-10). 
             [[:CO2pl, :O], [:Opl, :CO2], :(9.6e-11)], # Nair minimal ionosphere. # Fehsenfeld1970: -in use- Tenewitz 2018, sect 3b. 2e-11 * 0.02: :(4e-13). Nair: :(1.0e-10) Roger: :(2 .* 1.0e-10).
             [[:CO2pl, :O2], [:O2pl, :CO2], :(2 .* 5.5e-11)],
             [[:COpl, :C], [:Cpl, :CO], :(2 .* 1.1e-10)],
             [[:COpl, :CH], [:CHpl, :CO], :(2 .* 5.54e-9 .* (Ti .^ -0.5))],
             [[:COpl, :CH], [:HCOpl, :C], :(2 .* 5.54e-9 .* (Ti .^ -0.5))],
             [[:COpl, :CO2], [:CO2pl, :CO], :(2 .* 1.1e-9)],

             [[:COpl, :H], [:Hpl, :CO], :(2 .* 4.0e-10)],
             # NEW
             # [[:COpl, :D], [:Dpl, :CO], :(((2 ./ 1) .^ -0.5) .* 2 .* 4.0e-10)], # DUPLICATE

             [[:COpl, :H2], [:HCOpl, :H], :(2 .* 7.5e-10)],
             [[:COpl, :H2], [:HOCpl, :H], :(2 .* 7.5e-10)],
             # NEW
             [[:COpl, :HD], [:DCOpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 7.5e-10)],
             [[:COpl, :HD], [:HCOpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 7.5e-10)],
             [[:COpl, :HD], [:DOCpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 7.5e-10)],
             [[:COpl, :HD], [:HOCpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 7.5e-10)],

             [[:COpl, :H2O], [:H2Opl, :CO], :(2 .* 1.56e-9)],
             # NEW
             [[:COpl, :HDO], [:HDOpl, :CO], :(((19 ./ 18) .^ -0.5) .* 2 .* 1.56e-9)],

             [[:COpl, :H2O], [:HCOpl, :OH], :(2 .* 8.4e-10)],
             # NEW
             [[:COpl, :HDO], [:DCOpl, :OH], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 8.4e-10)],
             [[:COpl, :HDO], [:HCOpl, :OD], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 8.4e-10)],

             [[:COpl, :HCN], [:HCNpl, :CO], :(2 .* 3.06e-9)],
             [[:COpl, :HCO], [:HCOpl, :CO], :(2 .* 1.28e-8 .* (Ti .^ -0.5))],
             [[:COpl, :N], [:NOpl, :C], :(2 .* 8.2e-11)],
             [[:COpl, :NH], [:HCOpl, :N], :(2 .* 5.54e-9 .* (Ti .^ -0.5))],
             [[:COpl, :NH], [:NHpl, :CO], :(2 .* 5.54e-9 .* (Ti .^ -0.5))],
             [[:COpl, :NH2], [:HCOpl, :NH], :(2 .* 7.79e-9 .* (Ti .^ -0.5))],
             [[:COpl, :NH2], [:NH2pl, :CO], :(2 .* 7.79e-9 .* (Ti .^ -0.5))],
             [[:COpl, :NO], [:NOpl, :CO], :(2 .* 4.2e-10)],
             [[:COpl, :O], [:Opl, :CO], :(2 .* 1.4e-10)],
             [[:COpl, :O2], [:O2pl, :CO], :(2 .* 1.5e-10)],
             [[:COpl, :OH], [:HCOpl, :O], :(2 .* 5.37e-9 .* (Ti .^ -0.5))],
             [[:COpl, :OH], [:OHpl, :CO], :(2 .* 5.37e-9 .* (Ti .^ -0.5))],
             [[:Cpl, :C], [:C2pl], :(2 .* 1.52e-18 .* (Ti .^ 0.17) .* exp.(-101.5 ./ Ti))],
             [[:Cpl, :CH], [:C2pl, :H], :(2 .* 6.58e-9 .* (Ti .^ -0.5))],
             [[:Cpl, :CH], [:CHpl, :C], :(2 .* 6.58e-9 .* (Ti .^ -0.5))],
             [[:Cpl, :CO2], [:CO2pl, :C], :(2 .* 1.1e-10)],
             [[:Cpl, :CO2], [:COpl, :CO], :(2 .* 9.9e-10)],
             [[:Cpl, :H], [:CHpl], :(2 .* 1.7e-17)],
             [[:Cpl, :H2], [:CH2pl], :(2 .* 3.32e-13 .* (Ti .^ -1.3) .* exp.(-23.0 ./ Ti))],
             [[:Cpl, :H2], [:CHpl, :H], :(2 .* 7.4e-10 .* exp.(-4537.0 ./ Ti))],

             [[:Cpl, :H2O], [:H2Opl, :C], :(2 .* 2.4e-10)],
             # NEW
             [[:Cpl, :HDO], [:HDOpl, :C], :(((19 ./ 18) .^ -0.5) .* 2 .* 2.4e-10)],

             [[:Cpl, :H2O], [:HCOpl, :H], :(2 .* 1.56e-8 .* (Ti .^ -0.5))],
             # NEW
             [[:Cpl, :HDO], [:DCOpl, :H], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 1.56e-8 .* (Ti .^ -0.5))],
             [[:Cpl, :HDO], [:HCOpl, :D], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 1.56e-8 .* (Ti .^ -0.5))],

             [[:Cpl, :H2O], [:HOCpl, :H], :(2 .* 2.16e-9)],
             # NEW
             [[:Cpl, :HDO], [:DOCpl, :H], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 2.16e-9)],
             [[:Cpl, :HDO], [:HOCpl, :D], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 2.16e-9)],

             [[:Cpl, :HCN], [:C2Npl, :H], :(2 .* 2.95e-9)],
             [[:Cpl, :HCO], [:CHpl, :CO], :(2 .* 8.31e-9 .* (Ti .^ -0.5))],
             [[:Cpl, :HCO], [:HCOpl, :C], :(2 .* 8.31e-9 .* (Ti .^ -0.5))],
             [[:Cpl, :N], [:CNpl], :(2 .* 7.24e-19 .* (Ti .^ 0.07) .* exp.(-57.5 ./ Ti))],
             [[:Cpl, :N2O], [:NOpl, :CN], :(2 .* 9.1e-10)],
             [[:Cpl, :NH], [:CNpl, :H], :(2 .* 1.35e-8 .* (Ti .^ -0.5))],
             [[:Cpl, :NH2], [:HCNpl, :H], :(2 .* 1.91e-8 .* (Ti .^ -0.5))],
             [[:Cpl, :NO], [:NOpl, :C], :(2 .* 7.5e-10)],
             [[:Cpl, :O], [:COpl], :(2 .* 7.39e-18 .* (Ti .^ -0.15) .* exp.(-68.0 ./ Ti))],
             [[:Cpl, :O2], [:COpl, :O], :(2 .* 3.48e-10)],
             [[:Cpl, :O2], [:Opl, :CO], :(2 .* 5.22e-10)],
             [[:Cpl, :OH], [:COpl, :H], :(2 .* 1.33e-8 .* (Ti .^ -0.5))],
             [[:H2Opl, :C], [:CHpl, :OH], :(2 .* 1.1e-9)],
             [[:H2Opl, :CH], [:CH2pl, :OH], :(2 .* 5.89e-9 .* (Ti .^ -0.5))],
             [[:H2Opl, :CH], [:CHpl, :H2O], :(2 .* 5.89e-9 .* (Ti .^ -0.5))],

             [[:H2Opl, :CO], [:HCOpl, :OH], :(2 .* 4.25e-10)],
             # NEW 
             [[:HDOpl, :CO], [:DCOpl, :OH], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 4.25e-10)],
             [[:HDOpl, :CO], [:HCOpl, :OD], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 4.25e-10)],

             [[:H2Opl, :H2], [:H3Opl, :H], :(2 .* 7.6e-10)],
             # NEW
             [[:HDOpl, :H2], [:H2DOpl, :H], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 7.6e-10)],
             [[:HDOpl, :H2], [:H3Opl, :D], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 7.6e-10)],
             [[:H2Opl, :HD], [:H2DOpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 7.6e-10)],
             [[:H2Opl, :HD], [:H3Opl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 7.6e-10)],
             # [[:HDOpl, :HD], [:HD2Opl, :H], :(2 .* 7.6e-10)], # TODO: doubly deuterated?
             # [[:HDOpl, :HD], [:H2DOpl, :D], :(2 .* 7.6e-10)],

             [[:H2Opl, :H2O], [:H3Opl, :OH], :(2 .* 1.85e-9)],
             [[:H2Opl, :HCN], [:HCNHpl, :OH], :(2 .* 1.05e-9)],
             [[:H2Opl, :HCO], [:H2COpl, :OH], :(2 .* 4.85e-9 .* (Ti .^ -0.5))],
             [[:H2Opl, :HCO], [:H3Opl, :CO], :(2 .* 4.85e-9 .* (Ti .^ -0.5))],
             [[:H2Opl, :HCO], [:HCOpl, :H2O], :(2 .* 4.85e-9 .* (Ti .^ -0.5))],
             
             [[:H2Opl, :N], [:HNOpl, :H], :(2 .* 1.12e-10)],
             # NEW
             [[:HDOpl, :N], [:DNOpl, :H], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 1.12e-10)],
             [[:HDOpl, :N], [:HNOpl, :D], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 1.12e-10)],

             [[:H2Opl, :N], [:NOpl, :H2], :(2 .* 2.8e-11)],
             # NEW 
             [[:HDOpl, :N], [:NOpl, :HD], :(((19 ./ 18) .^ -0.5) .* 2 .* 2.8e-11)],

             [[:H2Opl, :NH], [:H3Opl, :N], :(2 .* 1.23e-8 .* (Ti .^ -0.5))],
             [[:H2Opl, :NH2], [:NH2pl, :H2O], :(2 .* 8.49e-9 .* (Ti .^ -0.5))],
             [[:H2Opl, :NH2], [:NH3pl, :OH], :(2 .* 8.49e-9 .* (Ti .^ -0.5))],

             [[:H2Opl, :NO], [:NOpl, :H2O], :(2 .* 4.6e-10)],
             # NEW 
             [[:HDOpl, :NO], [:NOpl, :HDO], :(((19 ./ 18) .^ -0.5) .* 2 .* 4.6e-10)],

             [[:H2Opl, :NO2], [:NO2pl, :H2O], :(2 .* 1.2e-9)],

             [[:H2Opl, :O], [:O2pl, :H2], :(2 .* 4.0e-11)],
             # NEW
             [[:HDOpl, :O], [:O2pl, :HD], :(((19 ./ 18) .^ -0.5) .* 2 .* 4.0e-11)],

             [[:H2Opl, :O2], [:O2pl, :H2O], :(2 .* 3.3e-10)],
             # NEW
             [[:HDOpl, :O2], [:O2pl, :HDO], :(((19 ./ 18) .^ -0.5) .* 2 .* 3.3e-10)],

             [[:H2Opl, :OH], [:H3Opl, :O], :(2 .* 1.2e-8 .* (Ti .^ -0.5))],
             
             [[:H2pl, :Ar], [:ArHpl, :H], :(2 .* 2.1e-9)],
             # NEW
             # [[:HDpl, :Ar], [:ArDpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.1e-9)], # DUPLICATE
             [[:HDpl, :Ar], [:ArHpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.1e-9)],

             [[:H2pl, :C], [:CHpl, :H], :(2 .* 2.4e-9)],
             [[:H2pl, :CH], [:CH2pl, :H], :(2 .* 1.23e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :CH], [:CHpl, :H2], :(2 .* 1.23e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :CN], [:CNpl, :H2], :(2 .* 2.08e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :CN], [:HCNpl, :H], :(2 .* 2.08e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :CO], [:COpl, :H2], :(2 .* 6.44e-10)],

             [[:H2pl, :CO], [:HCOpl, :H], :(2 .* 2.9e-9)],
             # NEW
             [[:HDpl, :CO], [:DCOpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.9e-9)],
             [[:HDpl, :CO], [:HCOpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.9e-9)],

             [[:H2pl, :CO2], [:HCO2pl, :H], :(2 .* 2.35e-9)],
             # NEW 
             [[:HDpl, :CO2], [:DCO2pl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.35e-9)],
             [[:HDpl, :CO2], [:HCO2pl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.35e-9)],

             [[:H2pl, :H], [:Hpl, :H2], :(2 .* 6.4e-10)],
             [[:H2pl, :H2], [:H3pl, :H], :(2 .* 2.0e-9)],
             [[:H2pl, :H2O], [:H2Opl, :H2], :(2 .* 3.87e-9)],
             [[:H2pl, :H2O], [:H3Opl, :H], :(2 .* 3.43e-9)],
             [[:H2pl, :HCN], [:HCNpl, :H2], :(2 .* 4.68e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :HCO], [:H3pl, :CO], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :HCO], [:HCOpl, :H2], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :N], [:NHpl, :H], :(2 .* 1.9e-9)],

             [[:H2pl, :N2], [:N2Hpl, :H], :(2 .* 2.0e-9)],
             # NEW
             [[:HDpl, :N2], [:N2Dpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.0e-9)],
             [[:HDpl, :N2], [:N2Hpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 2.0e-9)],

             [[:H2pl, :N2O], [:HN2Opl, :H], :(2 .* 1.32e-9)],
             [[:H2pl, :N2O], [:N2Hpl, :OH], :(2 .* 7.77e-10)],
             [[:H2pl, :NH], [:NH2pl, :H], :(2 .* 1.32e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :NH], [:NHpl, :H2], :(2 .* 1.32e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :NO], [:HNOpl, :H], :(2 .* 1.1e-9)],
             [[:H2pl, :NO], [:NOpl, :H2], :(2 .* 1.1e-9)],

             [[:H2pl, :O], [:OHpl, :H], :(2 .* 1.5e-9)],
             # NEW 
             [[:HDpl, :O], [:ODpl, :H], :(0.5 .* ((3 ./ 2) .^ -0.5) .* 2 .* 1.5e-9)],
             [[:HDpl, :O], [:OHpl, :D], :(0.5 .* ((3 ./ 2) .^ -0.5) .* 2 .* 1.5e-9)],

             [[:H2pl, :O2], [:HO2pl, :H], :(2 .* 1.92e-9)],
             # NEW
             [[:HDpl, :O2], [:DO2pl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 1.92e-9)],
             [[:HDpl, :O2], [:HO2pl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 1.92e-9)],

             [[:H2pl, :O2], [:O2pl, :H2], :(2 .* 7.83e-10)],
             [[:H2pl, :OH], [:H2Opl, :H], :(2 .* 1.32e-8 .* (Ti .^ -0.5))],
             [[:H2pl, :OH], [:OHpl, :H2], :(2 .* 1.32e-8 .* (Ti .^ -0.5))],
             [[:H3Opl, :C], [:HCOpl, :H2], :(2 .* 1.0e-11)],
             [[:H3Opl, :CH], [:CH2pl, :H2O], :(2 .* 1.18e-8 .* (Ti .^ -0.5))],
             [[:H3Opl, :HCN], [:HCNHpl, :H2O], :(2 .* 3.8e-9)],
             [[:H3Opl, :NH2], [:NH3pl, :H2O], :(2 .* 1.68e-8 .* (Ti .^ -0.5))],
             [[:H3pl, :Ar], [:ArHpl, :H2], :(2 .* 3.65e-10)],
             [[:H3pl, :C], [:CHpl, :H2], :(2 .* 2.0e-9)],
             [[:H3pl, :CH], [:CH2pl, :H2], :(2 .* 2.08e-8 .* (Ti .^ -0.5))],
             [[:H3pl, :CN], [:HCNpl, :H2], :(2 .* 3.46e-8 .* (Ti .^ -0.5))],
             [[:H3pl, :CO], [:HCOpl, :H2], :(2 .* 3.06e-9 .* (Ti .^ -0.142) .* exp.(3.41 ./ Ti))],
             [[:H3pl, :CO], [:HOCpl, :H2], :(2 .* 5.82e-10 .* (Ti .^ 0.0661) .* exp.(-5.21 ./ Ti))],
             [[:H3pl, :CO2], [:HCO2pl, :H2], :(2 .* 2.5e-9)],
             [[:H3pl, :H2O], [:H3Opl, :H2], :(2 .* 5.3e-9)],
             [[:H3pl, :HCN], [:HCNHpl, :H2], :(2 .* 7.5e-9)],
             [[:H3pl, :HCO], [:H2COpl, :H2], :(2 .* 2.94e-8 .* (Ti .^ -0.5))],
             [[:H3pl, :N], [:NH2pl, :H], :(2 .* 3.9e-10)],
             [[:H3pl, :N], [:NHpl, :H2], :(2 .* 2.6e-10)],
             [[:H3pl, :N2], [:N2Hpl, :H2], :(2 .* 1.63e-9)],
             [[:H3pl, :N2O], [:HN2Opl, :H2], :(2 .* 2.5e-9)],
             [[:H3pl, :NH], [:NH2pl, :H2], :(2 .* 2.25e-8 .* (Ti .^ -0.5))],
             [[:H3pl, :NO], [:HNOpl, :H2], :(2 .* 1.94e-9)],
             [[:H3pl, :NO2], [:NO2pl, :H2, :H], :(2 .* 7.0e-12)],
             [[:H3pl, :NO2], [:NOpl, :OH, :H2], :(2 .* 6.93e-10)],
             [[:H3pl, :O], [:H2Opl, :H], :(2 .* 8.33e-10 .* (Ti .^ -0.156) .* exp.(-1.4 ./ Ti))],
             [[:H3pl, :O], [:OHpl, :H2], :(2 .* 1.94e-9 .* (Ti .^ -0.156) .* exp.(-1.4 ./ Ti))],
             [[:H3pl, :O2], [:HO2pl, :H2], :(2 .* 6.7e-10)],
             [[:H3pl, :OH], [:H2Opl, :H2], :(2 .* 2.25e-8 .* (Ti .^ -0.5))],
             [[:HCNHpl, :CH], [:CH2pl, :HCN], :(2 .* 5.46e-9 .* (Ti .^ -0.5))],
             [[:HCNHpl, :H2O], [:H3Opl, :HCN], :(2 .* 8.8e-13)],
             [[:HCNpl, :C], [:CHpl, :CN], :(2 .* 1.1e-9)],
             [[:HCNpl, :CH], [:CH2pl, :CN], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             [[:HCNpl, :CO], [:HCOpl, :CN], :(2 .* 1.38e-10)],
             [[:HCNpl, :CO], [:HNCpl, :CO], :(2 .* 3.22e-10)],
             [[:HCNpl, :CO2], [:HCO2pl, :CN], :(2 .* 2.1e-10)],
             [[:HCNpl, :CO2], [:HNCpl, :CO2], :(2 .* 2.9e-10)],
             [[:HCNpl, :H], [:Hpl, :HCN], :(2 .* 3.7e-11)],
             [[:HCNpl, :H2], [:HCNHpl, :H], :(2 .* 8.8e-10)],
             [[:HCNpl, :H2O], [:H2Opl, :HCN], :(2 .* 3.12e-8 .* (Ti .^ -0.5))],
             [[:HCNpl, :H2O], [:H3Opl, :CN], :(2 .* 3.12e-8 .* (Ti .^ -0.5))],
             [[:HCNpl, :HCN], [:HCNHpl, :CN], :(2 .* 1.45e-9)],
             [[:HCNpl, :HCO], [:H2COpl, :CN], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:HCNpl, :HCO], [:HCNHpl, :CO], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:HCNpl, :N], [:CHpl, :N2], :(2 .* 2.2e-10)],
             [[:HCNpl, :N2O], [:N2Opl, :HCN], :(2 .* 1.08e-9)],
             [[:HCNpl, :NH], [:NH2pl, :CN], :(2 .* 1.13e-8 .* (Ti .^ -0.5))],
             [[:HCNpl, :NH2], [:NH3pl, :CN], :(2 .* 1.56e-8 .* (Ti .^ -0.5))],
             [[:HCNpl, :NO], [:NOpl, :HCN], :(2 .* 8.1e-10)],
             [[:HCNpl, :O2], [:O2pl, :HCN], :(2 .* 4.1e-10)],
             [[:HCNpl, :OH], [:H2Opl, :CN], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             [[:HCO2pl, :C], [:CHpl, :CO2], :(2 .* 1.0e-9)],
             
             [[:HCO2pl, :CO], [:HCOpl, :CO2], :(2 .* 7.8e-10)],
             # NEW
             [[:DCO2pl, :CO], [:DCOpl, :CO2], :(((46 ./ 45) .^ -0.5) .* 2 .* 7.8e-10)],

             [[:HCO2pl, :H2O], [:H3Opl, :CO2], :(2 .* 2.65e-9)],
             # NEW
             [[:DCO2pl, :H2O], [:H2DOpl, :CO2], :(((46 ./ 45) .^ -0.5) .* 2 .* 2.65e-9)],
             [[:HCO2pl, :HDO], [:H2DOpl, :CO2], :(((19 ./ 18) .^ -0.5) .* 2 .* 2.65e-9)],

             [[:HCO2pl, :O], [:HCOpl, :O2], :(2 .* 5.8e-10)],
             # NEW
             [[:DCO2pl, :O], [:DCOpl, :O2], :(((46 ./ 45) .^ -0.5) .* 2 .* 5.8e-10)],

             [[:HCOpl, :C], [:CHpl, :CO], :(2 .* 1.1e-9)],
             # NEW
             [[:DCOpl, :C], [:CDpl, :CO], :(((30 ./ 29) .^ -0.5) .* 2 .* 1.1e-9)],

             [[:HCOpl, :CH], [:CH2pl, :CO], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],

             [[:HCOpl, :H2O], [:H3Opl, :CO], :(2 .* 2.6e-9)],
             # NEW
             [[:DCOpl, :H2O], [:H2DOpl, :CO], :(((30 ./ 29) .^ -0.5) .* 2 .* 2.6e-9)],
             [[:HCOpl, :HDO], [:H2DOpl, :CO], :(((30 ./ 29) .^ -0.5) .* 2 .* 2.6e-9)],
             # [[:DCOpl, :HDO], [:HD2Opl, :CO], :(2 .* 2.6e-9)],  # TODO: Doubly deuterated mass scaling factor?

             [[:HCOpl, :H2O], [:HCOOH2pl], :(2 .* 6.64e-10 .* (Ti .^ -1.3))],
             [[:HCOpl, :HCN], [:HCNHpl, :CO], :(2 .* 3.5e-9)],
             [[:HCOpl, :HCO], [:H2COpl, :CO], :(2 .* 1.26e-8 .* (Ti .^ -0.5))],
             [[:HCOpl, :N2O], [:HN2Opl, :CO], :(2 .* 3.3e-12)],
             [[:HCOpl, :NH], [:NH2pl, :CO], :(2 .* 1.11e-8 .* (Ti .^ -0.5))],
             [[:HCOpl, :NH2], [:NH3pl, :CO], :(2 .* 1.54e-8 .* (Ti .^ -0.5))],
             [[:HCOpl, :OH], [:H2Opl, :CO], :(2 .* 1.07e-8 .* (Ti .^ -0.5))],
             [[:HCOpl, :OH], [:HCO2pl, :H], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:HN2Opl, :CO], [:HCOpl, :N2O], :(2 .* 5.3e-10)],
             [[:HN2Opl, :H2O], [:H3Opl, :N2O], :(2 .* 2.83e-9)],
             [[:HNOpl, :C], [:CHpl, :NO], :(2 .* 1.0e-9)],
             [[:HNOpl, :CH], [:CH2pl, :NO], :(2 .* 1.07e-8 .* (Ti .^ -0.5))],
             [[:HNOpl, :CN], [:HCNpl, :NO], :(2 .* 1.51e-8 .* (Ti .^ -0.5))],
             [[:HNOpl, :CO], [:HCOpl, :NO], :(2 .* 8.6e-10)],
             [[:HNOpl, :CO2], [:HCO2pl, :NO], :(2 .* 9.4e-10)],
             [[:HNOpl, :H2O], [:H3Opl, :NO], :(2 .* 2.3e-9)],
             [[:HNOpl, :HCN], [:HCNHpl, :NO], :(2 .* 1.71e-8 .* (Ti .^ -0.5))],
             [[:HNOpl, :HCO], [:H2COpl, :NO], :(2 .* 1.25e-8 .* (Ti .^ -0.5))],
             [[:HNOpl, :N2], [:N2Hpl, :NO], :(2 .* 1.0e-11)],
             [[:HNOpl, :NH], [:NH2pl, :NO], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             [[:HNOpl, :NH2], [:NH3pl, :NO], :(2 .* 1.52e-8 .* (Ti .^ -0.5))],
             [[:HNOpl, :NO], [:NOpl, :HNO], :(2 .* 7.0e-10)],
             [[:HNOpl, :O], [:NO2pl, :H], :(2 .* 1.0e-12)],
             [[:HNOpl, :OH], [:H2Opl, :NO], :(2 .* 1.07e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :C], [:CHpl, :O2], :(2 .* 1.0e-9)],
             [[:HO2pl, :CH], [:CH2pl, :O2], :(2 .* 1.07e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :CN], [:HCNpl, :O2], :(2 .* 1.49e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :CO], [:HCOpl, :O2], :(2 .* 8.4e-10)],
             [[:HO2pl, :CO2], [:HCO2pl, :O2], :(2 .* 1.1e-9)],
             [[:HO2pl, :H2], [:H3pl, :O2], :(2 .* 3.3e-10)],
             [[:HO2pl, :H2O], [:H3Opl, :O2], :(2 .* 1.42e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :HCN], [:HCNHpl, :O2], :(2 .* 1.68e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :HCO], [:H2COpl, :O2], :(2 .* 1.23e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :N], [:NO2pl, :H], :(2 .* 1.0e-12)],
             [[:HO2pl, :N2], [:N2Hpl, :O2], :(2 .* 8.0e-10)],
             [[:HO2pl, :NH], [:NH2pl, :O2], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :NH2], [:NH3pl, :O2], :(2 .* 1.51e-8 .* (Ti .^ -0.5))],
             [[:HO2pl, :NO], [:HNOpl, :O2], :(2 .* 7.7e-10)],
             [[:HO2pl, :O], [:OHpl, :O2], :(2 .* 6.2e-10)],
             [[:HO2pl, :OH], [:H2Opl, :O2], :(2 .* 1.06e-8 .* (Ti .^ -0.5))],

             [[:HOCpl, :CO], [:HCOpl, :CO], :(2 .* 6.0e-10)],
             # NEW
             [[:DOCpl, :CO], [:DCOpl, :CO], :(((30 ./ 29) .^ -0.5) .* 2 .* 6.0e-10)],

             [[:HOCpl, :CO2], [:HCO2pl, :CO], :(2 .* 9.45e-10)],
             [[:HOCpl, :H2], [:H3pl, :CO], :(2 .* 2.68e-10)],
             [[:HOCpl, :H2], [:HCOpl, :H2], :(2 .* 3.8e-10)],
             [[:HOCpl, :N2], [:N2Hpl, :CO], :(2 .* 6.7e-10)],
             [[:HOCpl, :N2O], [:HN2Opl, :CO], :(2 .* 1.17e-9)],
             [[:HOCpl, :NO], [:HNOpl, :CO], :(2 .* 7.1e-10)],
             [[:HOCpl, :O2], [:HO2pl, :CO], :(2 .* 1.9e-10)],
             [[:Hpl, :CH], [:CHpl, :H], :(2 .* 3.29e-8 .* (Ti .^ -0.5))],
             
             [[:Hpl, :CO2], [:HCOpl, :O], :(2 .* 3.8e-9)],
             # NEW
             # [[:Dpl, :CO2], [:DCOpl, :O], :(((2 ./ 1) .^ -0.5) .* 2 .* 3.8e-9)], # DUPLICATE

             [[:Hpl, :H], [:H2pl], :(2 .* 2.34e-22 .* (Ti .^ 1.49) .* exp.(-228.0 ./ Ti))],

             [[:Hpl, :H2O], [:H2Opl, :H], :(2 .* 8.2e-9)],
             # NEW
             [[:Dpl, :H2O], [:HDOpl, :H], :(((2 ./ 1) .^ -0.5) .* 0.5 .* 2 .* 8.2e-9)],
             # [[:Dpl, :H2O], [:H2Opl, :D], :(((2 ./ 1) .^ -0.5) .* 0.5 .* 2 .* 8.2e-9)], # DUPLICATE
             [[:Hpl, :HDO], [:HDOpl, :H], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 8.2e-9)],
             [[:Hpl, :HDO], [:H2Opl, :D], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 8.2e-9)],
             # [[:Dpl, :HDO], [:HDOpl, :D], :(((?? ./ ??) .^ -0.5) .* 2 .* 8.2e-9)], # TODO: doubly deuterated?

             [[:Hpl, :HCN], [:HCNpl, :H], :(2 .* 1.1e-8)],
             [[:Hpl, :HCO], [:COpl, :H2], :(2 .* 1.63e-8 .* (Ti .^ -0.5))],
             [[:Hpl, :HCO], [:H2pl, :CO], :(2 .* 1.63e-8 .* (Ti .^ -0.5))],
             [[:Hpl, :HCO], [:HCOpl, :H], :(2 .* 1.63e-8 .* (Ti .^ -0.5))],
             [[:Hpl, :HNO], [:NOpl, :H2], :(2 .* 6.93e-8 .* (Ti .^ -0.5))],
             [[:Hpl, :N2O], [:N2Hpl, :O], :(2 .* 3.52e-10)],
             [[:Hpl, :N2O], [:N2Opl, :H], :(2 .* 1.85e-9)],
             [[:Hpl, :NH], [:NHpl, :H], :(2 .* 3.64e-8 .* (Ti .^ -0.5))],
             [[:Hpl, :NH2], [:NH2pl, :H], :(2 .* 5.02e-8 .* (Ti .^ -0.5))],
             [[:Hpl, :NO], [:NOpl, :H], :(2 .* 1.9e-9)],
             [[:Hpl, :NO2], [:NOpl, :OH], :(2 .* 1.9e-9)],

             

             [[:Hpl, :O2], [:O2pl, :H], :(2 .* 1.17e-9)],
             # NEW 
             # [[:Dpl, :O2], [:O2pl, :D], :(((2 ./ 1) .^ -0.5) .* 2 .* 1.17e-9)], # DUPLICATE

             [[:Hpl, :OH], [:OHpl, :H], :(2 .* 3.64e-8 .* (Ti .^ -0.5))],
             [[:N2Hpl, :C], [:CHpl, :N2], :(2 .* 1.1e-9)],
             [[:N2Hpl, :CH], [:CH2pl, :N2], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             
             [[:N2Hpl, :CO], [:HCOpl, :N2], :(2 .* 8.8e-10)],
             # NEW
             [[:N2Dpl, :CO], [:DCOpl, :N2], :(((30 ./ 29) .^ -0.5) .* 2 .* 8.8e-10)], 

             [[:N2Hpl, :CO2], [:HCO2pl, :N2], :(2 .* 1.07e-9)],
             [[:N2Hpl, :H2], [:H3pl, :N2], :(2 .* 5.1e-18)],
             [[:N2Hpl, :H2O], [:H3Opl, :N2], :(2 .* 2.6e-9)],
             [[:N2Hpl, :HCN], [:HCNHpl, :N2], :(2 .* 3.2e-9)],
             [[:N2Hpl, :HCO], [:H2COpl, :N2], :(2 .* 1.26e-8 .* (Ti .^ -0.5))],
             [[:N2Hpl, :N2O], [:HN2Opl, :N2], :(2 .* 1.25e-9)],
             [[:N2Hpl, :NH], [:NH2pl, :N2], :(2 .* 1.11e-8 .* (Ti .^ -0.5))],
             [[:N2Hpl, :NH2], [:NH3pl, :N2], :(2 .* 1.54e-8 .* (Ti .^ -0.5))],
             [[:N2Hpl, :NO], [:HNOpl, :N2], :(2 .* 3.4e-10)],

             [[:N2Hpl, :O], [:OHpl, :N2], :(2 .* 1.4e-10)],
             # NEW
             [[:N2Dpl, :O], [:ODpl, :N2], :(((30 ./ 29) .^ -0.5) .* 2 .* 1.4e-10)],

             [[:N2Hpl, :OH], [:H2Opl, :N2], :(2 .* 1.07e-8 .* (Ti .^ -0.5))],
             [[:N2Opl, :CO], [:CO2pl, :N2], :(2 .* 1.11e-10)],
             [[:N2Opl, :CO], [:NOpl, :NCO], :(2 .* 1.89e-10)],
             [[:N2Opl, :H2], [:HN2Opl, :H], :(2 .* 2.56e-10)],
             [[:N2Opl, :H2], [:N2Hpl, :OH], :(2 .* 1.04e-10)],
             [[:N2Opl, :H2O], [:H2Opl, :N2O], :(2 .* 3.27e-8 .* (Ti .^ -0.5))],
             [[:N2Opl, :H2O], [:HN2Opl, :OH], :(2 .* 3.64e-9 .* (Ti .^ -0.5))],
             [[:N2Opl, :N2O], [:NOpl, :NO, :N2], :(2 .* 1.2e-11)],
             [[:N2Opl, :NO], [:NOpl, :N2O], :(2 .* 2.3e-10)],
             [[:N2Opl, :NO2], [:NO2pl, :N2O], :(2 .* 2.21e-10)],
             [[:N2Opl, :NO2], [:NOpl, :N2, :O2], :(2 .* 4.29e-10)],
             [[:N2Opl, :O2], [:NOpl, :NO2], :(2 .* 4.59e-11)],
             [[:N2Opl, :O2], [:O2pl, :N2O], :(2 .* 2.24e-10)],
             [[:N2pl, :Ar], [:Arpl, :N2], :(2 .* 2.0e-13)],
             [[:N2pl, :C], [:Cpl, :N2], :(2 .* 1.1e-10)],
             [[:N2pl, :CH], [:CHpl, :N2], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             [[:N2pl, :CN], [:CNpl, :N2], :(2 .* 1.73e-9 .* (Ti .^ -0.5))],
             [[:N2pl, :CO], [:COpl, :N2], :(2 .* 7.3e-11)],
             [[:N2pl, :CO2], [:CO2pl, :N2], :(2 .* 8.0e-10)],

             [[:N2pl, :H2], [:N2Hpl, :H], :(2 .* 1.87e-9 .* exp.(-54.7 ./ Ti))],
             # NEW
             # [[:N2pl, :HD], [:N2Hpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 1.87e-9 .* exp.(-54.7 ./ Ti))], # DUPLICATE
             # [[:N2pl, :HD], [:N2Dpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 1.87e-9 .* exp.(-54.7 ./ Ti))], # DUPLICATE

             [[:N2pl, :H2O], [:H2Opl, :N2], :(2 .* 1.9e-9)],
             # NEW 
             [[:N2pl, :HDO], [:HDOpl, :N2], :(((19 ./ 18) .^ -0.5) .* 2 .* 1.9e-9)],

             [[:N2pl, :H2O], [:N2Hpl, :OH], :(2 .* 5.04e-10)],
             # NEW
             [[:N2pl, :HDO], [:N2Dpl, :OH], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 5.04e-10)],
             [[:N2pl, :HDO], [:N2Hpl, :OD], :(((19 ./ 18) .^ -0.5) .* 0.5 .* 2 .* 5.04e-10)],

             [[:N2pl, :HCN], [:HCNpl, :N2], :(2 .* 3.9e-10)],
             [[:N2pl, :HCO], [:HCOpl, :N2], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:N2pl, :HCO], [:N2Hpl, :CO], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:N2pl, :N], [:Npl, :N2], :(2 .* 1.4e-11)],
             [[:N2pl, :N2O], [:N2Opl, :N2], :(2 .* 6.0e-10)],
             [[:N2pl, :NH], [:NHpl, :N2], :(2 .* 1.13e-8 .* (Ti .^ -0.5))],
             [[:N2pl, :NH2], [:NH2pl, :N2], :(2 .* 1.54e-8 .* (Ti .^ -0.5))],
             [[:N2pl, :NO], [:NOpl, :N2], :(2 .* 4.4e-10)],
             [[:N2pl, :O], [:NOpl, :N], :(2 .* 1.3e-10)],
             [[:N2pl, :O], [:Opl, :N2], :(2 .* 9.8e-12)],
             [[:N2pl, :O2], [:O2pl, :N2], :(2 .* 5.0e-11)],
             [[:N2pl, :OH], [:OHpl, :N2], :(2 .* 1.09e-8 .* (Ti .^ -0.5))],
             [[:NH2pl, :CH], [:CH2pl, :NH], :(2 .* 6.06e-9 .* (Ti .^ -0.5))],
             [[:NH2pl, :CH], [:CHpl, :NH2], :(2 .* 6.06e-9 .* (Ti .^ -0.5))],
             [[:NH2pl, :CN], [:HCNHpl, :N], :(2 .* 1.73e-9 .* (Ti .^ -0.5))],
             [[:NH2pl, :H2], [:NH3pl, :H], :(2 .* 1.95e-10)],
             [[:NH2pl, :H2O], [:H3Opl, :NH], :(2 .* 2.73e-9)],
             [[:NH2pl, :H2O], [:NH3pl, :OH], :(2 .* 8.7e-11)],
             [[:NH2pl, :H2O], [:NH4pl, :O], :(2 .* 1.16e-10)],
             [[:NH2pl, :HCN], [:HCNHpl, :NH], :(2 .* 2.08e-8 .* (Ti .^ -0.5))],
             [[:NH2pl, :HCO], [:H2COpl, :NH], :(2 .* 7.45e-9 .* (Ti .^ -0.5))],
             [[:NH2pl, :HCO], [:HCOpl, :NH2], :(2 .* 7.45e-9 .* (Ti .^ -0.5))],
             [[:NH2pl, :N], [:N2Hpl, :H], :(2 .* 9.1e-11)],
             [[:NH2pl, :NH], [:NH3pl, :N], :(2 .* 1.26e-8 .* (Ti .^ -0.5))],
             [[:NH2pl, :NH2], [:NH3pl, :NH], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:NH2pl, :NO], [:NOpl, :NH2], :(2 .* 7.0e-10)],
             [[:NH2pl, :O], [:HNOpl, :H], :(2 .* 7.2e-11)],
             [[:NH2pl, :O2], [:H2NOpl, :O], :(2 .* 1.19e-10)],
             [[:NH2pl, :O2], [:HNOpl, :OH], :(2 .* 2.1e-11)],
             [[:NH3pl, :CH], [:NH4pl, :C], :(2 .* 1.2e-8 .* (Ti .^ -0.5))],
             [[:NH3pl, :H2], [:NH4pl, :H], :(2 .* 4.4e-13)],
             [[:NH3pl, :H2O], [:NH4pl, :OH], :(2 .* 2.5e-10)],
             [[:NH3pl, :HCO], [:NH4pl, :CO], :(2 .* 7.27e-9 .* (Ti .^ -0.5))],
             [[:NH3pl, :NO], [:NOpl, :NH3], :(2 .* 7.2e-10)],
             [[:NHpl, :C], [:CHpl, :N], :(2 .* 1.6e-9)],
             [[:NHpl, :CH], [:CH2pl, :N], :(2 .* 1.71e-8 .* (Ti .^ -0.5))],
             [[:NHpl, :CN], [:HCNpl, :N], :(2 .* 2.77e-8 .* (Ti .^ -0.5))],
             [[:NHpl, :CO], [:HCOpl, :N], :(2 .* 4.41e-10)],
             [[:NHpl, :CO], [:OCNpl, :H], :(2 .* 5.39e-10)],
             [[:NHpl, :CO2], [:HCO2pl, :N], :(2 .* 3.85e-10)],
             [[:NHpl, :CO2], [:HNOpl, :CO], :(2 .* 3.85e-10)],
             [[:NHpl, :CO2], [:NOpl, :HCO], :(2 .* 3.3e-10)],
             [[:NHpl, :H2], [:H3pl, :N], :(2 .* 1.85e-10)],
             [[:NHpl, :H2], [:NH2pl, :H], :(2 .* 1.05e-9)],
             [[:NHpl, :H2O], [:H2Opl, :NH], :(2 .* 1.05e-9)],
             [[:NHpl, :H2O], [:H3Opl, :N], :(2 .* 1.05e-9)],
             [[:NHpl, :H2O], [:HNOpl, :H2], :(2 .* 3.5e-10)],
             [[:NHpl, :H2O], [:NH2pl, :OH], :(2 .* 8.75e-10)],
             [[:NHpl, :H2O], [:NH3pl, :O], :(2 .* 1.75e-10)],
             [[:NHpl, :HCN], [:HCNHpl, :N], :(2 .* 3.12e-8 .* (Ti .^ -0.5))],
             [[:NHpl, :HCO], [:H2COpl, :N], :(2 .* 2.25e-8 .* (Ti .^ -0.5))],
             [[:NHpl, :N], [:N2pl, :H], :(2 .* 1.3e-9)],
             [[:NHpl, :N2], [:N2Hpl, :N], :(2 .* 6.5e-10)],
             [[:NHpl, :NH], [:NH2pl, :N], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:NHpl, :NH2], [:NH3pl, :N], :(2 .* 2.6e-8 .* (Ti .^ -0.5))],
             [[:NHpl, :NO], [:N2Hpl, :O], :(2 .* 1.78e-10)],
             [[:NHpl, :NO], [:NOpl, :NH], :(2 .* 7.12e-10)],
             [[:NHpl, :O], [:OHpl, :N], :(2 .* 1.0e-9)],
             [[:NHpl, :O2], [:HO2pl, :N], :(2 .* 1.64e-10)],
             [[:NHpl, :O2], [:NOpl, :OH], :(2 .* 2.05e-10)],
             [[:NHpl, :O2], [:O2pl, :NH], :(2 .* 4.51e-10)],
             [[:NHpl, :OH], [:H2Opl, :N], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:NO2pl, :H], [:NOpl, :OH], :(2 .* 1.9e-10)],
             [[:NO2pl, :H2], [:NOpl, :H2O], :(2 .* 1.5e-10)],
             [[:NO2pl, :NO], [:NOpl, :NO2], :(2 .* 2.75e-10)],
             [[:Npl, :CH], [:CNpl, :H], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:Npl, :CN], [:CNpl, :N], :(2 .* 1.91e-8 .* (Ti .^ -0.5))],
             [[:Npl, :CO], [:COpl, :N], :(2 .* 4.93e-10)],
             [[:Npl, :CO], [:Cpl, :NO], :(2 .* 5.6e-12)],
             [[:Npl, :CO], [:NOpl, :C], :(2 .* 6.16e-11)],
             [[:Npl, :CO2], [:CO2pl, :N], :(2 .* 9.18e-10)],
             [[:Npl, :CO2], [:COpl, :NO], :(2 .* 2.02e-10)],

             [[:Npl, :H2], [:NHpl, :H], :(2 .* 5.0e-10 .* exp.(-85.0 ./ Ti))],
             # NEW
             # [[:Npl, :HD], [:NDpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 5.0e-10 .* exp.(-85.0 ./ Ti))], # DUPLICATE
             # [[:Npl, :HD], [:NHpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 5.0e-10 .* exp.(-85.0 ./ Ti))], # DUPLICATE

             [[:Npl, :H2O], [:H2Opl, :N], :(2 .* 2.7e-9)],
             [[:Npl, :HCN], [:HCNpl, :N], :(2 .* 3.7e-9)],
             [[:Npl, :HCO], [:HCOpl, :N], :(2 .* 7.79e-9 .* (Ti .^ -0.5))],
             [[:Npl, :HCO], [:NHpl, :CO], :(2 .* 7.79e-9 .* (Ti .^ -0.5))],
             [[:Npl, :N], [:N2pl], :(2 .* 9.44e-19 .* (Ti .^ 0.24) .* exp.(-26.1 ./ Ti))],
             [[:Npl, :N2O], [:NOpl, :N2], :(2 .* 5.5e-10)],
             [[:Npl, :NH], [:N2pl, :H], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:Npl, :NH], [:NHpl, :N], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:Npl, :NH2], [:NH2pl, :N], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:Npl, :NO], [:N2pl, :O], :(2 .* 8.33e-11)],
             [[:Npl, :NO], [:NOpl, :N], :(2 .* 4.72e-10)],
             [[:Npl, :O2], [:NOpl, :O], :(2 .* 2.32e-10)],
             [[:Npl, :O2], [:O2pl, :N], :(2 .* 3.07e-10)],
             [[:Npl, :O2], [:Opl, :NO], :(2 .* 4.64e-11)],
             [[:Npl, :OH], [:OHpl, :N], :(2 .* 6.41e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :C], [:COpl, :O], :(2 .* 5.2e-11)],
             [[:O2pl, :C], [:Cpl, :O2], :(2 .* 5.2e-11)],
             [[:O2pl, :CH], [:CHpl, :O2], :(2 .* 5.37e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :CH], [:HCOpl, :O], :(2 .* 5.37e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :HCO], [:HCOpl, :O2], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :HCO], [:HO2pl, :CO], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :N], [:NOpl, :O], :(2 .* 1.0e-10)],
             [[:O2pl, :NH], [:HNOpl, :O], :(2 .* 5.54e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :NH], [:NO2pl, :H], :(2 .* 5.54e-9 .* (Ti .^ -0.5))],
             [[:O2pl, :NH2], [:NH2pl, :O2], :(2 .* 1.51e-8 .* (Ti .^ -0.5))],
             [[:O2pl, :NO], [:NOpl, :O2], :(2 .* 4.6e-10)],
             [[:O2pl, :NO2], [:NO2pl, :O2], :(2 .* 6.6e-10)],
             [[:OHpl, :C], [:CHpl, :O], :(2 .* 1.2e-9)],
             [[:OHpl, :CH], [:CH2pl, :O], :(2 .* 6.06e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :CH], [:CHpl, :OH], :(2 .* 6.06e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :CN], [:HCNpl, :O], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],

             [[:OHpl, :CO], [:HCOpl, :O], :(2 .* 8.4e-10)],
             # NEW
             [[:ODpl, :CO], [:DCOpl, :O], :(((18 ./ 17) .^ -0.5) .* 2 .* 8.4e-10)],

             [[:OHpl, :CO2], [:HCO2pl, :O], :(2 .* 1.35e-9)],
             # NEW
             [[:ODpl, :CO2], [:DCO2pl, :O], :(((18 ./ 17) .^ -0.5) .* 2 .* 1.35e-9)],

             [[:OHpl, :H2], [:H2Opl, :H], :(2 .* 9.7e-10)],
             # NEW 
             [[:ODpl, :H2], [:HDOpl, :H], :(0.5 .* ((18 ./ 17) .^ -0.5) .* 2 .* 9.7e-10)],
             [[:ODpl, :H2], [:H2Opl, :D], :(0.5 .* ((18 ./ 17) .^ -0.5) .* 2 .* 9.7e-10)],
             [[:OHpl, :HD], [:HDOpl, :H], :(((3 ./ 2) .^ -0.5) .* 2 .* 9.7e-10)],
             [[:OHpl, :HD], [:H2Opl, :D], :(((3 ./ 2) .^ -0.5) .* 2 .* 9.7e-10)],
             # [[:ODpl, :HD], [:HDOpl, :D], :(2 .* 9.7e-10)],  # TODO: Doubly deuterated mass scaling factor?


             [[:OHpl, :H2O], [:H2Opl, :OH], :(2 .* 1.59e-9)],
             [[:OHpl, :H2O], [:H3Opl, :O], :(2 .* 1.3e-9)],
             [[:OHpl, :HCN], [:HCNHpl, :O], :(2 .* 2.08e-8 .* (Ti .^ -0.5))],
             [[:OHpl, :HCO], [:H2COpl, :O], :(2 .* 4.85e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :HCO], [:H2Opl, :CO], :(2 .* 4.85e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :HCO], [:HCOpl, :OH], :(2 .* 4.85e-9 .* (Ti .^ -0.5))],

             [[:OHpl, :N], [:NOpl, :H], :(2 .* 8.9e-10)],
             # NEW 
             [[:ODpl, :N], [:NOpl, :D], :(((18 ./ 17) .^ -0.5) .* 2 .* 8.9e-10)],

             [[:OHpl, :N2], [:N2Hpl, :O], :(2 .* 2.4e-10)],
             # NEW
             [[:ODpl, :N2], [:N2Dpl, :O], :(((18 ./ 17) .^ -0.5) .* 2 .* 2.4e-10)],

             [[:OHpl, :N2O], [:HN2Opl, :O], :(2 .* 9.58e-10)],
             [[:OHpl, :N2O], [:N2Opl, :OH], :(2 .* 2.13e-10)],
             [[:OHpl, :N2O], [:NOpl, :HNO], :(2 .* 1.46e-10)],
             [[:OHpl, :NH], [:NH2pl, :O], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :NH2], [:NH2pl, :OH], :(2 .* 8.66e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :NH2], [:NH3pl, :O], :(2 .* 8.66e-9 .* (Ti .^ -0.5))],
             [[:OHpl, :NO], [:HNOpl, :O], :(2 .* 6.11e-10)],
             [[:OHpl, :NO], [:NOpl, :OH], :(2 .* 8.15e-10)],

             [[:OHpl, :O], [:O2pl, :H], :(2 .* 7.1e-10)],
             # NEW
             [[:ODpl, :O], [:O2pl, :D], :(((18 ./ 17) .^ -0.5) .* 2 .* 7.1e-10)],

             [[:OHpl, :O2], [:O2pl, :OH], :(2 .* 3.8e-10)],
             # NEW 
             [[:ODpl, :O2], [:O2pl, :OD], :(((18 ./ 17) .^ -0.5) .* 2 .* 3.8e-10)],

             [[:OHpl, :OH], [:H2Opl, :O], :(2 .* 1.21e-8 .* (Ti .^ -0.5))],
             [[:Opl, :CH], [:CHpl, :O], :(2 .* 6.06e-9 .* (Ti .^ -0.5))],
             [[:Opl, :CH], [:COpl, :H], :(2 .* 6.06e-9 .* (Ti .^ -0.5))],
             [[:Opl, :CN], [:NOpl, :C], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:Opl, :CO2], [:O2pl, :CO], :(9.6e-10)], # Nair minimal ionosphere # Roger: :(2 .* 1.1e-9) # Nair: -in use-  
             
             [[:Opl, :H], [:Hpl, :O], :(2 .* 6.4e-10)],
             # NEW
             [[:Opl, :D], [:Dpl, :O], :(((2 ./ 1) .^ -0.5) .* 2 .* 6.4e-10)],


             [[:Opl, :H2], [:OHpl, :H], :(2 .* 1.62e-9)],
             # NEW
             # [[:Opl, :HD], [:OHpl, :D], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 1.62e-9)], # DUPLICATE
             # [[:Opl, :HD], [:ODpl, :H], :(((3 ./ 2) .^ -0.5) .* 0.5 .* 2 .* 1.62e-9)], # DUPLICATE

             [[:Opl, :H2O], [:H2Opl, :O], :(2 .* 2.6e-9)],
             # NEW
             [[:Opl, :HDO], [:HDOpl, :O], :(((19 ./ 18) .^ -0.5) .* 2 .* 2.6e-9)],

             [[:Opl, :HCN], [:COpl, :NH], :(2 .* 1.17e-9)],
             [[:Opl, :HCN], [:HCOpl, :N], :(2 .* 1.17e-9)],
             [[:Opl, :HCN], [:NOpl, :CH], :(2 .* 1.17e-9)],
             [[:Opl, :HCO], [:HCOpl, :O], :(2 .* 7.45e-9 .* (Ti .^ -0.5))],
             [[:Opl, :HCO], [:OHpl, :CO], :(2 .* 7.45e-9 .* (Ti .^ -0.5))],
             [[:Opl, :N2], [:NOpl, :N], :(2 .* 4.58e-9 .* (Ti .^ -1.37) .* exp.(-28.592 ./ Ti))],
             [[:Opl, :N2O], [:N2Opl, :O], :(2 .* 6.3e-10)],
             [[:Opl, :NH], [:NHpl, :O], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:Opl, :NH], [:NOpl, :H], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:Opl, :NH2], [:NH2pl, :O], :(2 .* 1.73e-8 .* (Ti .^ -0.5))],
             [[:Opl, :NO], [:NOpl, :O], :(2 .* 8.0e-13)],
             [[:Opl, :NO2], [:NO2pl, :O], :(2 .* 1.6e-9)],
             [[:Opl, :NO2], [:NOpl, :O2], :(2 .* 8.3e-10)],
             [[:Opl, :O2], [:O2pl, :O], :(2 .* 2.1e-11)],
             [[:Opl, :OH], [:O2pl, :H], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],
             [[:Opl, :OH], [:OHpl, :O], :(2 .* 6.24e-9 .* (Ti .^ -0.5))],

             # Recombination 
             [[:ArHpl, :E], [:Ar, :H], :(2 .* 1.0e-9)],
             [[:Arpl, :E], [:Ar], :(2 .* 4.0e-12 .* (Te .^ 0.6))],
             [[:CHpl, :E], [:C, :H], :(2 .* 1.65e-6 .* (Te .^ -0.42))],
             [[:CNpl, :E], [:N, :C], :(2 .* 3.12e-6 .* (Te .^ -0.5))],
             [[:CO2pl, :E], [:CO, :O], :(3.8e-7)], # Nair minimal ionosphere # Vuitton: :(4.2e-7 .* (Te ./ 300) .^ -0.75) Roger: :(2 .* 3.03e-5 .* (Te .^ -0.75)) Nair: in use. BAD RATE?: 3.03e-5 .* ((300 ./ Te) .^ -0.75)
             [[:COpl, :E], [:O, :C], :(2 .* 4.82e-6 .* (Te .^ -0.55))],
             [[:COpl, :E], [:O1D, :C], :(2 .* 2.48e-8 .* (Te .^ -0.55))],
             [[:Cpl, :E], [:C], :(2 .* 6.28e-10 .* (Te .^ -0.59))],

             [[:H2Opl, :E], [:O, :H, :H], :(2 .* 2.08e-5 .* (Te .^ -0.74))],
             # NEW 
             [[:HDOpl, :E], [:O, :D, :H], :(((19 ./ 18) .^ -0.5) .* 2 .* 2.08e-5 .* (Te .^ -0.74))],

             [[:H2Opl, :E], [:H2, :O], :(2 .* 2.64e-6 .* (Te .^ -0.74))],
             # NEW 
             [[:HDOpl, :E], [:HD, :O], :(((19 ./ 18) .^ -0.5) .* 2 .* 2.64e-6 .* (Te .^ -0.74))],

             [[:H2Opl, :E], [:OH, :H], :(2 .* 5.86e-6 .* (Te .^ -0.74))],
             # NEW 
             [[:HDOpl, :E], [:OH, :D], :(0.5 .* ((19 ./ 18) .^ -0.5) .* 2 .* 5.86e-6 .* (Te .^ -0.74))],
             [[:HDOpl, :E], [:OD, :H], :(0.5 .* ((19 ./ 18) .^ -0.5) .* 2 .* 5.86e-6 .* (Te .^ -0.74))],

             
             [[:H3Opl, :E], [:H2, :O, :H], :(2 .* 9.68e-8 .* (Te .^ -0.5))],
             [[:H2DOpl, :E], [:HD, :O, :H], :(0.5 .* ((20 ./ 19) .^ -0.5) .* 2 .* 9.68e-8 .* (Te .^ -0.5))], # TODO: H2DO+ and HD2O+ network???
             [[:H2DOpl, :E], [:H2, :O, :D], :(0.5 .* ((20 ./ 19) .^ -0.5) .* 2 .* 9.68e-8 .* (Te .^ -0.5))],
             # [[:HD2Opl, :E], [:HD, :O, :D], :(((21 ./ 19) .^ -0.5) .* 2 .* 9.68e-8 .* (Te .^ -0.5))],
             [[:H3Opl, :E], [:H2O, :H], :(2 .* 1.86e-6 .* (Te .^ -0.5))],
             [[:H2DOpl, :E], [:HDO, :H], :(0.5 .* ((20 ./ 19) .^ -0.5) .* 2 .* 1.86e-6 .* (Te .^ -0.5))],
             [[:H2DOpl, :E], [:H2O, :D], :(0.5 .* ((20 ./ 19) .^ -0.5) .* 2 .* 1.86e-6 .* (Te .^ -0.5))],
             # [[:HD2Opl, :E], [:HDO, :D], :(((21 ./ 19) .^ -0.5) .* 2 .* 1.86e-6 .* (Te .^ -0.5))],

             [[:H3Opl, :E], [:OH, :H, :H], :(2 .* 4.47e-6 .* (Te .^ -0.5))],
             # NEW
             [[:H2DOpl, :E], [:OD, :H, :H], :(0.5 .* ((20 ./ 19) .^ -0.5) .* 2 .* 4.47e-6 .* (Te .^ -0.5))],
             [[:H2DOpl, :E], [:OH, :D, :H], :(0.5 .* ((20 ./ 19) .^ -0.5) .* 2 .* 4.47e-6 .* (Te .^ -0.5))],
             # [[:HD2Opl, :E], [:OD, :D, :H], :(0.5 .* ((21 ./ 19) .^ -0.5) .* 2 .* 4.47e-6 .* (Te .^ -0.5))],
             # [[:HD2Opl, :E], [:OH, :D, :D], :(0.5 .* ((21 ./ 19) .^ -0.5) .* 2 .* 4.47e-6 .* (Te .^ -0.5))],
             [[:H3Opl, :E], [:OH, :H2], :(2 .* 1.04e-6 .* (Te .^ -0.5))],
             # NEW
             [[:H2DOpl, :E], [:OD, :H2], :(((20 ./ 19) .^ -0.5) .* 0.5 .* 2 .* 1.04e-6 .* (Te .^ -0.5))],
             [[:H2DOpl, :E], [:OH, :HD], :(((20 ./ 19) .^ -0.5) .* 0.5 .* 2 .* 1.04e-6 .* (Te .^ -0.5))],
             # [[:HD2Opl, :E], [:OH, :H2], :(((21 ./ 19) .^ -0.5) .* 0.5 .* 2 .* 1.04e-6 .* (Te .^ -0.5))],
             # [[:HD2Opl, :E], [:OD, :HD], :(((21 ./ 19) .^ -0.5) .* 0.5 .* 2 .* 1.04e-6 .* (Te .^ -0.5))],

             [[:H3pl, :E], [:H, :H, :H], :(2 .* 8.46e-7 .* (Te .^ -0.52))],
             [[:H3pl, :E], [:H2, :H], :(2 .* 4.54e-7 .* (Te .^ -0.52))],
             [[:HCNHpl, :E], [:CN, :H, :H], :(2 .* 3.79e-6 .* (Te .^ -0.65))],
             [[:HCNpl, :E], [:CN, :H], :(2 .* 3.46e-6 .* (Te .^ -0.5))],

             [[:HCO2pl, :E], [:CO, :O, :H], :(8.1e-7 .* (Te ./ 300) .^ -0.64)], # Roger: :(2 .* 2.18e-7)
             # NEW
             [[:DCO2pl, :E], [:CO, :O, :D], :(((46 ./ 45) .^ -0.5) .* 8.1e-7 .* (Te ./ 300) .^ -0.64)],
             
             [[:HCO2pl, :E], [:CO, :OH], :(3.2e-7 .* (Te ./ 300) .^ -0.64)], # Roger: :(2 .* 9.18e-8)
             # NEW
             [[:DCO2pl, :E], [:CO, :OD], :(((46 ./ 45) .^ -0.5) .* 3.2e-7 .* (Te ./ 300) .^ -0.64)], # Roger: :(2 .* 9.18e-8)

             [[:HCO2pl, :E], [:CO2, :H], :(3.0e-7)], # Nair minimal ionosphere # New: :(6e-8 .* (Te ./ 300) .^ -0.64) (where did I get this? probably Vuitton?) Roger: :(2 .* 1.7e-8) Nair: in use 
             # NEW
             [[:DCO2pl, :E], [:CO2, :D], :(((46 ./ 45) .^ -0.5) .* 3.0e-7)],

             [[:HCOpl, :E], [:CH, :O], :(2 .* 1.15e-7 .* (Te .^ -0.64))],
             [[:HCOpl, :E], [:CO, :H], :(2 .* 1.06e-5 .* (Te .^ -0.64))],
             [[:HCOpl, :E], [:OH, :C], :(2 .* 8.08e-7 .* (Te .^ -0.64))],
             # NEW
             [[:DCOpl, :E], [:CD, :O], :(((30 ./ 29) .^ -0.5) .* 2 .* 1.15e-7 .* (Te .^ -0.64))], # No branching ratio because each rate is different.
             [[:DCOpl, :E], [:CO, :D], :(((30 ./ 29) .^ -0.5) .* 2 .* 1.06e-5 .* (Te .^ -0.64))],
             [[:DCOpl, :E], [:OD, :C], :(((30 ./ 29) .^ -0.5) .* 2 .* 8.08e-7 .* (Te .^ -0.64))],

             [[:HN2Opl, :E], [:N2, :O, :H], :(2 .* 3.81e-5 .* (Te .^ -0.74))],
             [[:HN2Opl, :E], [:N2, :OH], :(2 .* 4.38e-5 .* (Te .^ -0.74))],
             [[:HNOpl, :E], [:NO, :H], :(2 .* 5.2e-6 .* (Te .^ -0.5))],
             [[:HO2pl, :E], [:O2, :H], :(2 .* 5.2e-6 .* (Te .^ -0.5))],
             [[:HOCpl, :E], [:CH, :O], :(2 .* 1.7e-9 .* (Te .^ 1.2))],
             [[:HOCpl, :E], [:CO, :H], :(2 .* 3.3e-5 .* (Te .^ -1.0))],

             [[:HOCpl, :E], [:OH, :C], :(2 .* 1.19e-8 .* (Te .^ 1.2))],
             # NEW
             [[:DOCpl, :E], [:OD, :C], :(((30 ./ 29) .^ -0.5) .* 2 .* 1.19e-8 .* (Te .^ 1.2))],

             [[:Hpl, :E], [:H], :(2 .* 6.46e-14 .* (Te .^ 0.7))],

             [[:N2Hpl, :E], [:N2, :H], :(2 .* 6.6e-7 .* (Te .^ -0.51))],
             # NEW 
             [[:N2Dpl, :E], [:N2, :D], :(((30 ./ 29) .^ -0.5) .* 2 .* 6.6e-7 .* (Te .^ -0.51))],

             [[:N2Hpl, :E], [:NH, :N], :(2 .* 1.17e-6 .* (Te .^ -0.51))],
             [[:N2Opl, :E], [:N, :N, :O], :(2 .* 1.36e-6 .* (Te .^ -0.57))],
             [[:N2Opl, :E], [:N2, :O], :(2 .* 4.09e-6 .* (Te .^ -0.57))],
             [[:N2Opl, :E], [:NO, :N], :(2 .* 3.07e-6 .* (Te .^ -0.57))],
             [[:N2pl, :E], [:N, :N], :(2 .* 5.09e-7 .* (Te .^ -0.39))],
             [[:N2pl, :E], [:N, :N], :(2 .* 1.42e-6 .* (Te .^ -0.39))], # These are secretly N(²D) but I don't track that species. 
             [[:NH2pl, :E], [:N, :H, :H], :(2 .* 1.71e-5 .* (Te .^ -0.8) .* exp.(-17.1 ./ Te))],
             [[:NH2pl, :E], [:NH, :H], :(2 .* 8.34e-6 .* (Te .^ -0.79) .* exp.(-17.1 ./ Te))],
             [[:NH3pl, :E], [:NH, :H, :H], :(2 .* 2.68e-6 .* (Te .^ -0.5))],
             [[:NH3pl, :E], [:NH2, :H], :(2 .* 2.68e-6 .* (Te .^ -0.5))],
             [[:NHpl, :E], [:N, :H], :(2 .* 7.45e-7 .* (Te .^ -0.5))],
             [[:NO2pl, :E], [:NO, :O], :(2 .* 5.2e-6 .* (Te .^ -0.5))],
             [[:NOpl, :E], [:O, :N], :(6.93e-6 .* (Te .^ -0.5))], # NEW rate from Fox+ 2015 ("Chemistry of protonated..."). Really, 95% BR to N(2D) and 5% to else. Roger's rate rom a '90 and '77 pub: 2 .* 8.52e-7 .* (Te .^ -0.37)
             [[:Npl, :E], [:N], :(2 .* 1.9e-10 .* (Te .^ -0.7))],
             [[:O2pl, :E], [:O, :O], :(6.6e-5 .* Te .^ -1.0)], # Nair minimal ionosphere.  # Vuitton: :(1.95e-7 .* (Te ./ 300) .^ (-0.7)) Roger: :(2 .* 8.15e-6 .* (Te .^ -0.65)) Nair: in use. BAD RATE?: 8.15e-6 .* ((300 ./ Te) .^ -0.65)
             
             [[:OHpl, :E], [:O, :H], :(2 .* 6.5e-7 .* (Te .^ -0.5))],
             # NEW 
             [[:ODpl, :E], [:O, :D], :(((18 ./ 17) .^ -0.5) .* 2 .* 6.5e-7 .* (Te .^ -0.5))],

             [[:Opl, :E], [:O], :(2 .* 1.4e-10 .* (Te .^ -0.66))],

             # D IONS!
             [[:Arpl, :HD], [:HDpl, :Ar], :(0.06*8.0e-10)],
             [[:Arpl, :HD], [:ArHpl, :D], :(0.46*8.0e-10)],
             [[:Arpl, :HD], [:ArDpl, :H], :(0.48*8.0e-10)],
             [[:ArHpl, :HD], [:H2Dpl, :Ar], :(8.6e-10)],
             [[:ArDpl, :H2], [:H2Dpl, :Ar], :(8.8e-10)],
             [[:ArDpl, :H2], [:ArHpl, :HD], :(4.5e-10)],
             [[:ArDpl, :N2], [:N2Dpl, :Ar], :(6e-10)],
             [[:ArDpl, :CO], [:DCOpl, :Ar], :(1.25e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:ArDpl, :CO2], [:OCODpl, :Ar], :(1.1e-9)],
             [[:Cpl, :HD], [:CHpl, :D], :(0.17*1.2e-16)],
             [[:Cpl, :HD], [:CDpl, :H], :(0.83*1.2e-16)],
             [[:COpl, :D], [:Dpl, :CO], :(9e-11)], # HAS DUPLICATE (COMPARE RATES)
             [[:CO2pl, :D], [:DCOpl, :O], :(0.76*8.4e-11)], # HAS DUPLICATE (COMPARE RATES)
             [[:CO2pl, :D], [:Dpl, :CO2], :(0.24*8.4e-11)], # HAS DUPLICATE (COMPARE RATES)
             [[:Dpl, :H2], [:Hpl, :HD], :(2.2e-9)],
             
             [[:Dpl, :H2O], [:H2Opl, :D], :(5.2e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:Dpl, :O2], [:O2pl, :D], :(1.6e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:Dpl, :CO2], [:DCOpl, :O], :(2.6e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:Dpl, :CO2], [:CO2pl, :D], :(2.6e-9)],
             [[:Dpl, :N2O], [:N2Opl, :D], :(0.85 * 1.7e-9)],
             [[:Dpl, :N2O], [:N2Dpl, :O], :(0.15 * 1.7e-9)],
             [[:Dpl, :NO], [:NOpl, :D], :(1.8e-9)],
             [[:Dpl, :NO2], [:NOpl, :OD], :(0.95*1.6e-9)],
             [[:Dpl, :NO2], [:NO2pl, :D], :(0.05*1.6e-9)],
             [[:DCOpl, :H], [:HCOpl, :D], :(1.5e-11)],
             [[:Hpl, :HD], [:Dpl, :H2], :(1.1e-10)],
             [[:HCNpl, :D], [:Dpl, :HCN], :(3.7e-11)],
             [[:HCOpl, :D], [:DCOpl, :H], :(4.25e-11)],
             [[:H2Dpl, :H2], [:H3pl, :HD], :(5.3e-10)],
             [[:H2Dpl, :HD], [:H3pl, :D2], :(0.1*5.0e-10)],
             [[:H2Dpl, :HD], [:HD2pl, :H2], :(0.9*5.0e-10)],
             [[:HDpl, :HD], [:H2Dpl, :D], :(0.55*1.53e-9)],
             [[:HDpl, :HD], [:HD2pl, :H], :(0.45*1.53e-9)],
             [[:HDpl, :Ar], [:ArDpl, :H], :(1.53e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:HD2pl, :H2], [:H3pl, :D2], :(0.25*7.6e-10)],
             [[:HD2pl, :H2], [:H2Dpl, :HD], :(0.75*7.6e-10)],
             [[:HD2pl, :HD], [:H2Dpl, :D2], :(0.25*4.5e-10)],
             [[:HD2pl, :HD], [:D3pl, :H2], :(0.75*4.5e-10)],
             [[:Npl, :HD], [:NHpl, :D], :(0.25*3.1e-10)],# HAS DUPLICATE (COMPARE RATES)
             [[:Npl, :HD], [:NDpl, :H], :(0.75*3.1e-10)], # HAS DUPLICATE (COMPARE RATES)
             [[:N2pl, :HD], [:N2Hpl, :D], :(0.49*1.34e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:N2pl, :HD], [:N2Dpl, :H], :(0.51*1.34e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:N2Hpl, :D], [:N2Dpl, :H], :(8e-11)],
             [[:N2Dpl, :H], [:N2Hpl, :D], :(2.5e-11)],
             [[:Opl, :HD], [:OHpl, :D], :(0.54*1.25e-9)], # HAS DUPLICATE (COMPARE RATES)
             [[:Opl, :HD], [:ODpl, :H], :(0.46*1.25e-9)], # HAS DUPLICATE (COMPARE RATES)

             # KEY NON-THERMAL MECHANISM REACTIONS

             ## Molecular hydrogen/hydrogen deuteride ion recombination
             [[:H2pl, :E], [:H, :H], :(2 .* 1.86e-7 .* (Te .^ -0.43))],
             [[:HDpl, :E], [:H, :D], :(((3 ./ 2) .^ -0.5) .* 2 .* 1.86e-7 .* (Te .^ -0.43))], 

             ## Charge exchange,
             [[:Hpl, :H], [:H, :Hpl], :(6.5e-11 .* (Ti .^ 0.5))], # NEW - Yung 1989
             [[:Hpl, :O], [:H, :Opl], :(2 .* 3.75e-10)], # Anicich 1993/Roger (factor of 2). Yung 1989 rate: :(0.8 .* 6.5e-11 .* (Ti .^ 0.5)
             [[:Dpl, :H], [:D, :Hpl], :(0.87 .* 6.5e-11 .* (Ti .^ 0.5))], # NEW - Yung 1989
             [[:Dpl, :O], [:D, :Opl], :(2.8e-10)], # Anicich 2003. Guesstimate by mass-scaling: :(((2 ./ 1) .^ -0.5) .* 2 .* 3.75e-10)

             ];
