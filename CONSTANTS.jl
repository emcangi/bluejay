################################################################################
# CONSTANTS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: General constants. 
# 
# Eryn Cangi
# Created November 2021
# Currently tested for Julia: 1.6.1
################################################################################

const kB_MKS = 1.38e-23;        # J/K - needed for saturation vapor pressure empirical equation.
const kB = 1.38e-16;            # erg/K
const bigG = 6.67e-8;           # dyne-cm^2/g^2
const mH = 1.67e-24;            # g 
const marsM = 0.1075*5.972e27;  # g 
const radiusM = 3396e5;         # cm
const q = 4.8032e-10            # statcoulomb (cm^1.5 g^0.5 s^-1)
const DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Polarizability from NIST. Experimental values from: https://cccbdb.nist.gov/pollistx.asp
# Calculations for species not available in experiment from: https://cccbdb.nist.gov/polcalc2x.asp
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