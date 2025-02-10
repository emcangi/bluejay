################################################################################
# CONSTANTS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: General constants. 
# 
# Eryn Cangi
# Created November 2021
# Currently tested for Julia: 1.8.5
################################################################################

#                              Universal constants
# ===============================================================================
const kB_MKS = 1.38e-23;        # J/K - needed for saturation vapor pressure empirical equation.
const kB = 1.38e-16;            # erg/K
const bigG = 6.67e-8;           # dyne-cm^2/g^2
const mH = 1.67e-24;            # g 
const q = -4.8032e-10            # statcoulomb (cm^1.5 g^0.5 s^-1)
const density_water_cgs = 0.997408  # used mostly for conversions
const s_per_yr = 3.154e7
const SMOW = 1.6e-4             # Standard mean ocean water

# Polarizability from NIST. Experimental values from: https://cccbdb.nist.gov/pollistx.asp
# Calculations for species not available in experiment from: https://cccbdb.nist.gov/polcalc2x.asp
# Deuterated species not listed in either are estimated by me to be the same as their H-bearing analogue.
# I used the calcualtions that use "Density functional", "aug-cc-PVDZ", and "mPW1PW91" 
# because that was the method that gave the closest answer for HD to the experimental value. 
# I have no idea what any of it means or whether it's reasonable. I'm not a quantum chemist.
# Values are given in cm^3
const polarizability = Dict(# Values available from experiment
                            :Ar=>1.664e-24, :C=>1.760e-24,  :CO=>1.953e-24,  :CO2=>2.507e-24, 
                            :H=>0.667e-24,  :H2=>0.8023e-24, # H2 source: Yoon 2010 <-- Kusakabe+ 2004, Phys. Rev. A
                            :H2O=>1.501e-24, :HCN=>2.593e-24, :HD=>0.7976e-24, # HD source: Yoon 2010 <-- Kusakabe+ 2004, Phys. Rev. A
                            :N=>1.1e-24,    :N2=>1.710e-24, :N2O=>2.998e-24, :NO=>1.698e-24, :NO2=>2.910e-24, 
                            :O=>0.802e-24,  :O2=>1.59e-24,  :O3=>3.079e-24, 

                            # Values from calculation
                            :CH=>2.108e-24,   :CN=>3.042e-24,   :D=>0.713e-24, 
                            :H2O2=>2.143e-24, :HCO=>2.505e-24,  :HDO=>1.358e-24, :HNO=>2.123e-24, 
                            :HO2=>1.858e-24,  :HOCO=>3.224e-24, :NH=>1.418e-24,  :NH2=>1.752e-24, 
                            :OH=>1.020e-24,   :OD=>1.020e-24, :HO2NO2=>4.734e-24,

    # I used density functional mPW1PW91, and basis set aug-cc-pVDZ, but I dont think this was the best choice since aug-cc-pVQZ is generally considered to be better (DZ is double zeta, and QZ is quadupal zeta), and I don't know much about what are actually considered the best desnity functionals/hamaltonians and basis sets.

                            #Chlorine species experiment
                            :HCl=> 2.515e-24, :Cl=> 2.180e-24, :Cl2=> 4.610e-24, :COCl2=>6.790e-24,
    
                            # chlorine species calculated
                            :DCl=> 2.409e-24, :ClO=> 2.907e-24, :ClCO=>4.474e-24, :ClCO3=>4.474e-24, :ClNO=>4.790e-24, :ClO2=>4.948e-24,

                            #Sulfur species experiment
                            :S=> 2.900e-24,  :SO2=> 3.882e-24, :SO3=>4.297e-24, :OCS=>5.090e-24,

                            #Sulfur species calculation 
                            :SO=> 3.351e-24, :H2SO4=>5.533e-24, :S2=>5.882e-24, :S3=>8.378e-24, :S2O=>4.297e-24, :HSO3=>5.121e-24, 

                            # Chlorine and Sulfur species calculated
                            :SCl=>4.297e-24, :SCl2=>7.448e-24, :ClS2=>8.788e-24, :SO2Cl2=>8.719e-24, :OSCl=>5.896e-24, :S2Cl2=>11.307e-24,
    
                            # Assumed same as hydrogen analogue
                            :DCO=>2.505e-24, :DO2=>1.858e-24, :DOCO=>3.224e-24, :HDO2=>2.143e-24, :O1D=>0.802e-24, 
                            :HDSO4=>5.533e-24, :DSO3=>5.121e-24, :DO2NO2=>4.734e-24,
                            # Assumed same as non-excited version
                            :Nup2D=>1.710e-24,

                            # placeholder for new Cl and S speices post REU that didn't have calculations
                            :ClCO3=>4.297e-24, :ClSO2=>4.297e-24, :SNO=>4.297e-24, :S2O2=>4.297e-24, 
                            )

const molmass = Dict(:H=>1, :Hpl=>1, 
                     :H2=>2, :H2pl=>2, :D=>2, :Dpl=>2, 
                     :HD=>3, :HDpl=>3, :H3pl=>3, 
                     :H2Dpl=>4, 
                     :HD2pl=>5, 
                     :C=>12, :Cpl=>12,  
                     :CH=>13, :CHpl=>13, 
                     :N=>14, :Npl=>14,
                     :NH=>15, :NHpl=>15, 
                     :NH2=>16, :NH2pl=>16, :O=>16, :O1D=>16, :Opl=>16, 
                     :NH3pl=>17, :OH=>17, :OHpl=>17, 
                     :H2O=>18, :OD=>18, :ODpl=>18, :H2Opl=>18, 
                     :HDO=>19, :HDOpl=>19, :H3Opl=>19,
                     :H2DOpl=>20,
                     :CN=>26, :CNpl=>26, 
                     :HCN=>27, :HCNpl=>27, 
                     :CO=>28, :COpl=>28, :HCNHpl=>28, :N2=>28, :Nup2D=>28, :N2pl=>28, 
                     :HCO=>29, :HCOpl=>29, :HOCpl=>29, :N2Hpl=>29, 
                     :DCO=>30, :DCOpl=>30, :DOCpl=>30, :NO=>30,  :NOpl=>30, :N2Dpl=>30, 
                     :HNO=>31, :HNOpl=>31, 
                     :O2=>32, :O2pl=>32, 
                     :HO2=>33, :HO2pl=>33, 
                     :DO2=>34, :H2O2=>34, 
                     :HDO2=>35,
                     :Ar=>40, :Arpl=>40, 
                     :ArHpl=>41, 
                     :ArDpl=>42,
                     :CO2=>44, :CO2pl=>44, :N2O=>44, :N2Opl=>44, 
                     :HOCO=>45, :HCO2pl=>45, :HN2Opl=>45,  
                     :DOCO=>46, :DCO2pl=>46, :NO2=>46, :NO2pl=>46,
                     :O3=>48, :HCl=>36, :Cl=>35, :ClO=>51, :ClCO=>63, :Cl2=>71, :DCl=>37,
                     :S=>32, :SO=>48, :SO2=>64, :SO3=>80, :H2SO4=>98, :HDSO4=>99,
    
#an ugly way to put in the Post REU Cl and S species except for Cl2 which was alreaddy writen above
    :ClNO=>67, :COCl2=>99, :ClCO3=>95, :ClO2=>67, :SCl=>68, :SCl2=>103, :SO2Cl2=>135, :OSCl=>84, :ClSO2=>100, :SNO=>62, :S2=>64, :S3=>96, :S2O=>80,
    :S2O2=>96, :OCS=>60, :HSO3=>81, :DSO3=>82, :HO2NO2=>79, :DO2NO2=>80, :S2Cl2=>135, :ClS2=>100
                     )

#= Some links for molecule diamters 
    https://pubs.acs.org/doi/full/10.1021/jp412588f?casa_token=fK7ezVqJuxUAAAAA%3AzrZ_oV-cDO5_XfoiZVB9mvF3arIfkANbBuVdcarx62ZJkyP-mpBPs8QwlIQeS_kBzj-6JCoyHirOK1m_ 
    Quantum Mechanical Basis for Kinetic Diameters of Small Gaseous Molecules.Nada Mehio†Sheng Dai†‡De-en Jiang*‡. ACS publications
=#
const diameters = Dict(# I think I want these in cm

                # Simple molecules from Breck, D. W. Zeolite Molecular Sieves: Structure, Chemistry and Use. John Wiley & Sons, Inc.: New York, 1974.
                :H2=> 2.89e-8, :O2=>3.46e-8, :N2=>3.64e-8, :CO=>3.76e-8, :CO2 => 3.30e-8,

    )

const collision_xsect = Dict(:H=>4e-15, # Zhang 2009
                             :D=>4.5e-15, 
                             :H2=>4.5e-15, # assume same as D since they have the same mass...
                             :HD=>5e-15 # assumption that adding a proton or neutron adds 0.5e-15 to the cross section...
                            ) # Units of cm^2; Bohr radius 8.79e-17


#                                 Float types
# ===============================================================================
# This section was introduced by Mike to attempt to solve issues with convergence. 
# However, the model is now running without use of doubles, so this may be removable,
# but has been left just in case it's ever needed again. Use of Doubles increases
# model run time but potentially delivers higher stability in complex simulations. 
# This is probably NOT the best place for it, but is the minimally annoying thing to 
# do after a major reorganization April 2024 by Eryn. This way it can be shared between
# the photochemistry module and the converge_new_file.jl code. 

ftype_ncur = Float64  # used to store n_current values
    # OPTIONS: Float64, Double64   
ftype_chem = Float64 #Double64 #  used to compute chemical reaction rates and chemical jacobian
    # OPTIONS: Float64, Double64   