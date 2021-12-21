This module contains all the functions needed to run the 1D photochemical simulation in converge_new_file.jl. It has become very large. Each function contains its own description, and this overview file contains a brief description of the categories.

# Basic utility functions

These functions are basic functions usedmore than once to operate on other objects. They have little to do the specific science addressed by this project and are often used for reading command line input, formatting inputs, looking for nans, etc. 

create_folder
deletefirst
find_nonfinites
format_chemistry_string
format_sec_or_min
fluxsymbol
getpos
input
nans_present
next_in_loop
searchdir
search_subfolders
subtract_difflength
       
# Plotting functions

Functions that make plots and do plotting-adjacent tasks.

get_colors
get_grad_colors
plot_atm
plot_bg
plot_extinction
plot_Jrates
plot_rxns
plot_temp_prof
plot_water_profile

                    
# Reaction rate functions

Functions that help calculate reaction rates or column reaction rates.
       
get_column_rates
make_ratexdensity
rxn_chem
rxn_photo
       
# Atmosphere array manipulation

These functions manipulate the array of densities by altitude that describes the atmosphere.

flatten_atm
get_ncurrent
n_tot
unflatten_atm
write_ncurrent

# Boundary condition functions

These functions set up the boundary conditions for the simulation. 

boundaryconditions
effusion_velocity

# transport functions

These functions define molecular diffusion, flux, eddy diffusion, and scale heights--anything to do with physical transport. 

Dcoef
Dcoef!
fluxcoefs
flux_param_arrays
flux_pos_and_neg
get_flux
Keddy
scaleH                                 

# Chemistry functions

Functions that define and set up the basic chemistry in the simulation. 

chemical_jacobian
getrate
loss_equations
loss_rate
make_chemjac_key
make_net_change_expr
meanmass
production_equations
production_rate

# Photochemical equilibrium functions

Functions used to set up expressions for calculating densities of species that are assumed to be in phoeochemical equilibrium.

choose_solutions
construct_quadratic
group_terms
loss_coef
linear_in_species_density

# Photochemistry functions

Functions that process data inputs for photochemical cross sections and build a dictionary of cross sections.

binupO2
co2xsect
h2o2xsect_l
h2o2xsect
hdo2xsect
ho2xsect_l
o2xsect
O3O1Dquantumyield
padtosolar
populate_xsect_dict
quantumyield 

# Temperature functions

Atmospheric temperatures by altitude.

T_all
Tpiecewise    
                                                                                                
# Water profile functions  

Functions that give the saturation vapor pressure for H2O and HDO vapor. Currently, they return the same values, but are written as separate to allow for updates in the future.

Psat Psat_HDO                                                                                                                                                                               

