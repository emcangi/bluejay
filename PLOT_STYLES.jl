################################################################################
# PLOT_STYLES.jl
# DESCRIPTION: Sets some style items for plotting. This is a separate file from 
# MODEL_SETUP.jl because this file also has to be frequently loaded separately 
# in analysis notebooks.
# 
# Eryn Cangi
# Created April 2024
# Last edited: May 2024
# Currently tested for Julia: 1.8.5
################################################################################

# ***************************************************************************************************** #
#                                                                                                       #
#                                         Plot styling                                                  #
#                                                                                                       #
# ***************************************************************************************************** #

# Sans-serif font to use 
sansserif_choice = "Arial"

# Monospace font to use
monospace_choice = "FreeMono"

const speciescolor = Dict( # PRIMARY NEUTRALS + IONS
                    :CO2 =>"#333", :CO2pl=>"#333",
                    :CO =>"#ff6600", :COpl=>"#ff6600",
                    :N2=>"#aaaaaa", :N2pl=>"#aaaaaa", :Nup2D=>"#aaa",
                    :Ar=>"#808080", :Arpl=>"#808080", :ArHpl=>"#956979", :ArDpl=>"#956979",
                    :HCO=> "#33bbf9", :HCOpl=>"#33bbf9", :DCO=> "#33bbf9", :DCOpl=>"#33bbf9",
                    :HOCpl=>"#c71a85", :DOCpl=>"#c71a85",
                    :HOCO =>"#667522", :DOCO =>"#667522", :HCO2pl=>"#667522", :DCO2pl=>"#667522",
                    :O1D =>"#7c3b6b",
                    :O => "#7922b4", :Opl=>"#7922b4",
                    :O2 => "#1BAD53", :O2pl=>"#1BAD53",
                    :O3 =>"#1c4bb4",

                    # ODD HYDROGEN + IONS
                    :H => "#ed3e7e", :Hpl=>"#ed3e7e", :D => "#ed3e7e", :Dpl=>"#ed3e7e", 
                    :H2 => "#964550", :H2pl=>"#964550", :HD => "#964550", :HDpl=>"#964550", 
                    :H3pl=> "#A81047", :H2Dpl=>"#A81047", # :HD2pl=>"#A81047",
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
                    #:HCNHpl=>"#1d3971",
                    :HNO=>"#ff1c5d",:HNOpl=>"#ff1c5d",
                    #:HN2Opl=>"#fb7562",
                    :N=>"#6e748a",:Npl=>"#6e748a",
                    :N2O=>"#ca9260",:N2Opl=>"#ca9260",
                    :NH=>"#059dc5",:NHpl=>"#059dc5", 
                    #:NH2=>"#0b522e",:NH2pl=>"#0b522e", 
                    #:NH3pl=>"#4ba40b", 
                    :NO=>"#e639b1",:NOpl=>"#e639b1",
                    :NO2=>"#a492e5", :NO2pl=>"#a492e5",  
                    :N2Hpl=>"#611115",:N2Dpl=>"#611115",
                    :HO2NO2=>"#FF3232", :DO2NO2=>"#FF3232",

                    # many post REU Cl and S species still have place holder numbers
                    :Cl=>"#2EE4EC", :Cl2=>"#89FAD6",
                    :HCl=>"#1F540F", :DCl=>"#DC8181",  
                    :ClO=>"#EEEA0C", :ClCO=>"#EDAC0D", :ClO2=>"#FF3232", :ClCO3=>"#FF3232",
                    :ClNO=>"#FF3232", :COCl2=>"#FF3232",
                    :S=>"#FF3232", :S2=>"#FF3232", :S3=>"#FF3232",
                    :SO=>"#E8C3FF", :SO2=>"#74B18C", :SO3=>"#C9FF55", 
                    :S2O=>"#FF3232",:S2O2=>"#FF3232",
                    :H2SO4=>"#542D5B", :HDSO4=>"#ACFF93", :HSO3=>"#FF3232", :DSO3=>"#FF3232",
                    :SCl=>"#FF3232", :SCl2=>"#FF3232", :S2Cl2=>"#FF3232", :ClS2=>"#FF3232",
                    :SO2Cl2=>"#FF3232", :OSCl=>"#FF3232", :ClSO2=>"#FF3232", 
                    :SNO=>"#FF3232",   :OCS=>"#FF3232", 
                    );

# NOTE: Some code is repeated here below, also occurring in get_deuterated, to figure out which species are deuterated. 
# This is because this file needs to be called by Photochemistry.jl, and I couldn't bear to not have a plot style variable
# somewhere other than PLOT_STYLES.jl. --Eryn
# TODO: This needs to be modified if you add new deuterated species! But, it's only plot styles, so I figured it was ok to 
# be a little lax about this.
known_species = keys(speciescolor)
Dspc = [s for s in setdiff(known_species, [:Nup2D, :O1D]) if occursin('D', string(s))]

# D group will have dashed lines; neutrals, solid (default)
const speciesstyle = Dict(vcat([s=>"--" for s in Dspc], [s=>"-" for s in setdiff(keys(speciescolor), Dspc)], [:HD2pl=>":", :Nup2D=>"-."]) )


const medgray = "#444444"

