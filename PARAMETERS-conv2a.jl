################################################################################
# PARAMETERS-conv2a.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Global constants, simulation parameters, reaction networks. 
# USE: This takes a file where the minimal ionosphere has converged, then adds in
# all the new neutrals and the rest of the ions and converges them all at once. 
#
# Eryn Cangi
# Created December 2019
# Last edited: August 2021
# Currently tested for Julia: 1.6.1
################################################################################

# **************************************************************************** #

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!                      !!!!! SUPER IMPORTANT !!!!!                       !! #
# !!     !!! Check the following every time you run the simulation !!!      !! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# NOTE: Special settings for this one!
const sim_folder_name = "ions2a-QNDF"
const initial_atm_file = "ions1-QNDF.h5"
const final_atm_file = "ions2a-QNDF"
const converge_which = "both"
const dt_min_and_max = Dict("neutrals"=>[-3, 14], "ions"=>[-4, 6], "both"=>[-4, 7])
const rel_tol = 1e-4

# This stuff is mutable, but less likely to change.
const make_new_alt_grid = false
const use_nonzero_initial_profiles = true
const do_chem = true 
const do_trans = true 
const solarfile = "marssolarphotonflux_solarmean.dat" # you may replace 'mean' with 'max' or 'min'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!! end check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# Sets whether photochemical equilibrium is assumed. Aids in converging ions and neutrals
# together. Generally leave it as is so the code determines it, but you can change it
# if need be. 

# Here set to false because we want all species active. 
const assume_photochem_eq = converge_which == false#"both" ? true : false

# Some water stuff to control the general shape of the profile
const upper_lower_bdy = 80e5 # the uppermost layer at which water will be fixed, in cm
const hygropause_alt = 40e5
const MR_mean_water = 1.38e-4

# Species name lists and Jrate symbol lists  ===================================

# Neutrals --------------------------------------------------------------------
const conv_neutrals = [:Ar, :CO, :CO2, :H, :H2, :H2O, :H2O2, :HO2, :HOCO, :N2, 
                       :O, :O1D, :O2, :O3, :OH,
                       :D, :DO2, :DOCO, :HD, :HDO, :HDO2, :OD,];
const new_neutrals = [:C, :CH, :HCO, :N2O, :NO2, :CN, :HCN, :HNO, :N, :NH, :NH2, :NO];
const neutral_species = [];
append!(neutral_species, conv_neutrals)
append!(neutral_species, new_neutrals)


# Ions -------------------------------------------------------------------------
const conv_ions = [:CO2pl, :HCO2pl, :Opl, :O2pl]; # Nair minimal ionosphere 
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
                  ];
const newJrates = [# New neutral photodissociation (from Roger)
                    :JCO2toCpOpO, :JCO2toCpO2, :JCOtoCpO,
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
                    :JO3toO3pl];

const Jratelist = [];
append!(Jratelist, conv_Jrates)
append!(Jratelist, newJrates)

# This dictionary specifies the species absorbing a photon for each J rate using regular expressions.
const absorber = Dict([x=>Symbol(match(r"(?<=J).+(?=to)", string(x)).match) for x in Jratelist])

# Other logical groupings -------------------------------------------------------
# List of D bearing species and their H analogue ions.
const D_H_analogues = Dict(:ArDpl=>:ArHpl, :Dpl=>:Hpl, :DCOpl=>:HCOpl, :HDpl=>:H2pl, :HD2pl=>:H3pl, :H2Dpl=>:H3pl, :N2Dpl=>:N2Hpl,
                           :DCO2pl=>:HCO2pl, :DOCpl=>:HOCpl, :H2DOpl=>:H3Opl, :HDOpl=>:H2Opl, :ODpl=>:OHpl)  
const D_bearing_species = [s for s in union(neutral_species, ion_species) if occursin('D', string(s))];
const D_ions = [s for s in ion_species if occursin('D', string(s))];
const N_neutrals = [s for s in neutral_species if occursin('N', string(s))];

# Short lived species, whose chemical lifetime is << diffusion timescale
const short_lived_species = [];
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
    append!(no_chem_species, union(conv_ions, N_neutrals)) # This is because the N chemistry is intimiately tied up with the ions.
    append!(no_transport_species, union(conv_ions, N_neutrals, short_lived_species))
elseif converge_which == "ions"
    append!(no_chem_species, setdiff(conv_neutrals, N_neutrals))
    append!(no_transport_species, union(short_lived_species, setdiff(conv_neutrals, N_neutrals)))
elseif converge_which == "both"
    # Next two lines are SPECIAL for this file because we are trying to work in ions, N-neutrals. 
    append!(no_chem_species, setdiff([:C, :CH, :HCO], setdiff(conv_neutrals, N_neutrals)))  # i.e. non-nitrogen bearing neutrals are inactive
    append!(no_transport_species, setdiff([:C, :CH, :HCO], setdiff(conv_neutrals, N_neutrals)))
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

# Annoyingly, this has to be here because D_bearing_species is defined in this file.
# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict([s=>"--" for s in D_bearing_species]);
