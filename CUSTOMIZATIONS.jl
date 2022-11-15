################################################################################
# CUSTOMIZATIONS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Inputs for functions in the photochemistry module that could
# ostensibly change or need to be called from outside the module. Includes,
# for instance, the names of files containing photochemical cross sections,
# colors to use for species density plots, and so forth.
#
# Eryn Cangi
# Updated: January 2022
# Currently tested for Julia: 1.6.1
################################################################################

using DoubleFloats

const code_dir = "$(@__DIR__)/"
# const extra_plots_dir = code_dir*"../Auxiliary plots/"
const results_dir = code_dir*"../Results/"#-Paper3/"
const xsecfolder = code_dir*"uvxsect/";

# Float types for calculations =================================================
# needed by both Photochemistry.jl and converge_new_file so it has to go here
ftype_ncur = Float64 #Double64 #  used to store n_current values
ftype_chem = Float64 #Double64 #  used to compute chemical reaction rates and chemical jacobian


# Altitude grid specifications =================================================
const max_alt = 250e5
const dz = 2e5
const alt = convert(Array, (0:dz:max_alt))
const n_all_layers = length(alt)
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]; # the altitude grid CELLS but in integers.
const non_bdy_layers = alt[2:end-1]  # all layers, centered on 2 km, 4...248. Excludes the boundary layers which are [-1, 1] and [249, 251].
const num_layers = length(non_bdy_layers) # there are 124 non-boundary layers.
const plot_grid = non_bdy_layers ./ 1e5;  # for plotting. Points located at atmospheric layer cell centers and in units of km.

const zmin = alt[1]
const zmax = alt[end];
const n_alt_index=Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])

# water altitude stuff 
const hygropause_alt = 40e5

# Common plot specifications =======================================================
const speciescolor = Dict( # PRIMARY NEUTRALS + IONS
                    :CO2 =>"#000000", :CO2pl=>"#000000",
                    :CO =>"#ff6600", :COpl=>"#ff6600",
                    :N2=>"#cccccc", :N2pl=>"#cccccc", :Nup2D=>"#cccccc",
                    :Ar=>"#808080", :Arpl=>"#808080", :ArHpl=>"#069668", :ArDpl=>"#069668",
                    :HCO=> "#33bbf9", :HCOpl=>"#33bbf9", :DCO=> "#33bbf9", :DCOpl=>"#33bbf9",
                    :HOCpl=>"#c71a85", :DOCpl=>"#c71a85",
                    :HOCO =>"#667522", :DOCO =>"#667522", :HCO2pl=>"#667522", :DCO2pl=>"#667522",
                    :O1D =>"#7c3b6b",
                    :O => "#7922b4", :Opl=>"#7922b4",
                    :O2 => "#38e278", :O2pl=>"#38e278",
                    :O3 =>"#1c4bb4",

                    # ODD HYDROGEN + IONS
                    :H => "#ed3e7e", :Hpl=>"#ed3e7e", :D => "#ed3e7e", :Dpl=>"#ed3e7e", 
                    :H2 => "#964550", :H2pl=>"#964550", :HD => "#964550", :HDpl=>"#964550", 
                    :H3pl=> "#02531d", :HD2pl=>"#02531d", :H2Dpl=>"#02531d", 
                    :H2O => "#0258ff", :H2Opl=>"#0258ff", :HDO => "#0258ff", :HDOpl=>"#0258ff", 
                    :H3Opl=> "#4f2381", :H2DOpl=>"#4f2381", 
                    :H2O2 =>"#d48f4d", :HDO2 =>"#d48f4d", 
                    :HO2 => "#609111", :HO2pl=>"#609111", :DO2 => "#609111", 
                    
                    :OH => "#a968d2", :OHpl=>"#a968d2", :OD => "#a968d2", :ODpl=>"#a968d2", 


                    # NITROGEN NEUTRALS + IONS
                    :C=>"#d9c382",:Cpl=>"#d9c382",
                    :CH=>"#cea3ce",:CHpl=>"#cea3ce",
                    :CN=>"#6d5000",:CNpl=>"#6d5000",
                    :HCN=>"#519169",:HCNpl=>"#519169",     
                    :HCNHpl=>"#1d3971",
                    :HNO=>"#ff1c5d",:HNOpl=>"#ff1c5d",
                    :HN2Opl=>"#fb7562",
                    :N=>"#6e748a",:Npl=>"#6e748a",
                    :N2O=>"#ca9260",:N2Opl=>"#ca9260",
                    :NH=>"#059dc5",:NHpl=>"#059dc5", 
                    :NH2=>"#0b522e",:NH2pl=>"#0b522e", 
                    :NH3pl=>"#4ba40b", 
                    :NO=>"#e639b1",:NOpl=>"#e639b1",
                    :NO2=>"#a492e5", :NO2pl=>"#a492e5",  
                    :N2Hpl=>"#611115",:N2Dpl=>"#611115", 
                    );

const medgray = "#444444"

# Crosssection filenames ======================================================

const photochem_data_files = Dict(:CO2=>Dict("main"=>"CO2.dat"), 
                                   :H2O=>Dict("main"=>"h2oavgtbl.dat"), 
                                   :HDO=>Dict("main"=>"HDO.dat"), 
                                   :H2O2=>Dict("main"=>"H2O2.dat"), 
                                   :HDO2=>Dict("main"=>"H2O2.dat"), 
                                   :O3=>Dict("main"=>"O3.dat", "chapman"=>"O3Chap.dat"), 
                                   :O2=>Dict("main"=>"O2.dat", "schr_short"=>"130-190.cf4", "schr_mid"=>"190-280.cf4", "schr_long"=>"280-500.cf4"), 
                                   :H2=>Dict("main"=>"binnedH2.csv"), 
                                   :HD=>Dict("main"=>"binnedH2.csv"), 
                                   :OH=>Dict("main"=>"binnedOH.csv", "O1D+H"=>"binnedOHo1D.csv"), 
                                   :OD=>Dict("main"=>"OD.csv"))