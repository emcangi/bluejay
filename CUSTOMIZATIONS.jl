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
const extra_plots_dir = code_dir*"../Auxiliary plots/"
const results_dir = code_dir*"../Results/"
const xsecfolder = code_dir*"uvxsect/";

# Float types for calculations =================================================
# needed by both Photochemistry.jl and converge_new_file so it has to go here
ftype_ncur = Double64 #Float64 #  used to store n_current values
ftype_chem = Double64 #Float64 #  used to compute chemical reaction rates and chemical jacobian


# Altitude grid specifications =================================================
const max_alt = 250e5
const dz = 2e5
const alt = convert(Array, (0:dz:max_alt))
const intaltgrid = round.(Int64, alt/1e5)[2:end-1]; # the altitude grid CELLS but in integers.
const non_bdy_layers = alt[2:end-1]  # all layers, centered on 2 km, 4...248. Excludes the boundary layers which are [-1, 1] and [249, 251].
const num_layers = length(non_bdy_layers) # there are 124 non-boundary layers.
const plot_grid = non_bdy_layers ./ 1e5;  # for plotting. Points located at atmospheric layer cell centers and in units of km.
const plot_grid_bdys = collect(1:2:((max_alt / 1e5)-1))  # the boundaries; includes the boundary layers at top and bottom.

const zmin = alt[1]
const zmax = alt[end];
const n_alt_index=Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])

# water altitude stuff 
const upper_lower_bdy = 80e5 # the uppermost layer at which water will be fixed, in cm
const upper_lower_bdy_i = Int64(upper_lower_bdy / dz) # the uppermost layer at which water will be fixed, in cm
const hygropause_alt = 40e5

# Mean temperatures and simulation temperature parameters ======================
const meanTs = 216.0
const meanTm = 130.0
const meanTe = 205.0
const meantemps = [meanTs, meanTm, meanTe]

const meanTsint = 216
const meanTmint = 130
const meanTeint = 205

# These are the low and high values for the "standard atmosphere and reasonable 
# climate variations" cases. NOT the full range of temperatures used to make the
# detailed cases stuff, because those include temps up to 350 K.
const lowTs = 160.0
const hiTs = 270.0
const lowTm = 100.0
const hiTm = 160.0
const lowTe = 150.0
const hiTe = 250.0

const highestTe = 350.0  # This is the highest exobase temperature considered.
                          # AFAIK this parameter is only used in making the 3-panel
                          # temperature profile plot.

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
