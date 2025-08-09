# Bluejay expansion from 1-D to 2-D: OLD UPDATES

Note that this is a simplified version with 4 species and 7 altitude bins.

## Nomenclature
- I added a new variable, n_horiz, which is the number of vertical columns or slices. It is currently set to 3, but can be set to 1, which should produce the same output as the original 1-D version, or to any other number.

- In many places in the code, I set things to iterate through each vertical column. In these cases, I always named the index of the vertical column ihoriz, e.g., for ihoriz in 1:n_horiz…

## What I changed:
- Where the atmospheric state is saved, values for each vertical column can be saved. For example:
    - Dictionary atm_dict now has structure {species1:[[densities with altitude for column 1][densities with altitude for col. 2][densities with alt. for col. 3]…], species2:[[…}

- For functions chemJmat and ratefn, a loop over n_horiz was added, so that sparse matrices (dimensions n_horiz*length(active_longlived)*num_layers by ihoriz*length(active_longlived)*num_layers) are produced. When this was implemented, the outputs for each vertical column were identical, as expected. Note that this was only tested with four species in the model.

- Horizontal transport:
    - A new function update_horiz_transport_coefficients, calculates horizontal transport coefficients at the edges and between columns.
    - A new function boundaryconditions_horiz in Core.jl allows edge boundary conditions (i.e., into and out of the first and last vertical columns) to be set and saved in bc_dict_horiz, which now has structure [[0.0 0.0; 0.0 0.0] for ialt in 1:GV.num_layers] for each species. Currently, these are set to zero flux. If density or deposition velocity boundary conditions are required instead, this will need some further work.
    - Dictionary fluxcoefs_horiz_all has structure Dict{Symbol, Array{Float64}}(s=>zeros(9,2) for s in GV.all_species), and has not yet been expanded to include different values for each column (because at the moment all values are zero for each column). This needs to be expanded. Eventually fluxcoefs_horiz_all and fluxcoefs_all should have the same dimensions as one another.

- Vertical transport for each column:
    - The dictionary that stores the lower and upper boundary conditions, bc_dict, now has structure [[X X; X X] for ihoriz in 1:n_horiz] for each species (defined in Core.jl).
    - The vertical transport variables fluxcoefs_all and fluxcoef_dict have been expanded to include values for each vertical column. They now have structures Dict{Symbol, Vector{Array{ftype_ncur}}} ([s=>[fill(0., GV.n_all_layers, 2) for ihoriz in 1:n_horiz] for s in species_list]). Note that the values have not yet been changed to be different for each column, and vertical transport is currently the same for each column.
    - The diffusion coefficient dictionary, Dcoef_dict, has been expanded to include multiple columns. It now has structure: Dict{Symbol, Vector{Vector{ftype_ncur}}}([s=>[[Dcoefs with altitude] for ihoriz in 1:n_horiz] for s in species_list])
    - The eddy diffusion structure, K, has been expanded to include different eddy diffusion profiles for each vertical column. It now has structure: [[K with altitude] for ihoriz in 1:n_horiz].
    - The dictionary with scale height values (H0_dict) has been expanded to take values for all vertical columns. It now has the structure: Dict{String, Vector{Vector{ftype_ncur}}}("neutral"=>[[values with altitude] for ihoriz in 1:n_horiz], "ion"=>[values with altitude] for ihoriz in 1:n_horiz]).
    - Different lower and upper boundary conditions for each vertical column can be set in MODEL_SETUP.jl. Dictionary speciesbclist now has a vector of n_horiz vectors for each species and BC type.

## My next planned steps were:
- I was about to expand MODEL_SETUP.jl to incorporate multiple temperature profiles for each vertical column.
- Currently, horizontal transport is set to zero. The model is currently not set up for non-zero horizontal fluxes. Need to add the terms for the horizontal transport to the ∂ni/∂t equation (see Mike’s summary sheet with equations). I had just added horiz_wind_v in MODEL_SETUP.jl to include some horizontal velocities (set to 0 currently).
- All of the places that might need to be looked at/checked/adapted in the future are indicated with the comment "#MULTICOL WARNING"

## My medium-term planned steps were (not necessarily in a linear order):
- Test that the model is doing the right thing.
    - Make sure that if 1) there are zero flux edge boundary conditions and horizontal transport is set to zero, 2) the lower and upper boundary conditions are the same for each vertical column, and 3) that the other inputs (e.g., temperature, vertical mixing) are the same for each vertical column, the output densities in each vertical column are identical.
    - Make sure that if n_horiz is 1, the model produces exactly the same output as the original 1-D version.
    - Set up a test model with two altitude bins by 2 vertical slices and without chemistry. Compare to analytical solution.
- Connect the front and back edge such that the model is cyclic rather than linear.
- I reduced the species list to four species (O, O2, O3 and Opl); put all the species back in.
- I reduced the altitude grid to 7 bins (by reducing zmax in MODEL_SETUP.jl); put zmax back to 250e5 (line commented out).

## Small things to note:
- Some further edits will need to be made to function diffusion_timescale in AnalyzeChemAndTransport if this is ever called (it isn’t set to be called in the current set-up).
- In some places where function n_tot is called, I have hardcoded the argument ihoriz as 1 for the time being, which will need to be changed for flexibility going forward. There is a '#MULTICOL WARNING' comment in each instance of this.

# Bluejay expansion from 1-D to 2-D: NEW UPDATES

- Continuing working on a simplified version with 4 species (O, O2, O3, and Opl) and 7 altitude bins (90 to 106 km, dz = 2 km) and 3 horizontal columns
- Two-dimensional temperature arrays are defined with horizontal index first and altitude second. MODEL_SETUP.jl initializes Tn_temp, Ti_temp, and Te_temp with shape (n_horiz, num_layers+2) and fills them using Tn_temp[ihoriz, :] etc., then stores them as global constants Tn_arr, Ti_arr, and Te_arr
- Dcoef! operates column-by-column by extracting T_arr_2D[ihoriz, :]; the indexing matches the (n_horiz, num_layers+2) layout used throughout the codebase
- Solar absorption arrays from optical_depth are built as a vector over horizontal columns containing altitude arrays (solarabs[ihoriz][ialt]), explicitly documented as “structured as [n_horiz, num_layers]” and accessed with that ordering when computing Jrates
- Vertical boundary conditions use a dictionary keyed by species and column: bc_dict[sp][ihoriz][row, col]
- Horizontal winds and boundary conditions
    - The variable `horiz_wind_v` in `MODEL_SETUP.jl` provides a horizontal wind profile for each column. By default a constant wind speed (see `horiz_wind_speed` in `INPUT_PARAMETERS.jl`) in cm/s is applied at all altitudes so that horizontal transport is active (set to zero for no horizontal transport). `update_horiz_transport_coefficients` uses these velocities to create forward and backward transport rates. 
    - Edge fluxes are set through the `speciesbclist_horiz` dictionary; specify altitude profiles for each species to impose non-zero flux at the back and front edges. Passing `cyclic=true` treats the domain as periodic so that flux leaving one edge enters from the opposite side and horizontal coefficients wrap between the first and last columns.
    - Horizontal advection employs an upwind scheme that averages the local and neighbouring wind speeds so that flux leaving one column exactly enters the next.
- Verification of Jacobian and Rate Functions: To verify that the multicolumn model keeps a consistent ordering of species, altitude, and horizontal columns when building transport matrices.
- Horizontal boundary conditions are stored separately per altitude with bc_dict_horiz[sp][ialt][row, col] and populated consistently inside the loop over ialt
- Make sure the horizontal transport coefficients are placed correctly in the transport matrices (chemical Jacobian and rate functions)

## Horizontal advection timescale

- Use appropriate wind speeds and column widths, and advection time scales.

## Merged multicolumn with master

- I verified that 2-D multicolumn gives the same final atmosphere as the original 1-D master for both Mars and Venus. I turned off the horizontal transport and all columns are identical as the 1-D; this confirms the 2-D multicolumn is working with all the species and altitudes as in the 1-D master.

- When `n_horiz=1`, the 2-D multicolumn model falls back into a 1-D model.

- Merged multicolumn branch with master branch.

## How to create diff files in my local computer

- From multicolumn branch
<!-- (master in red and multicolumn in green) -->
git diff master multicolumn -- Photochemistry/src/Core.jl > Core.diff
git diff master multicolumn -- Photochemistry/src/AnalyzeChemAndTransport.jl > AnalyzeChemAndTransport.diff
git diff master multicolumn -- Photochemistry/src/BasicUtilities.jl > BasicUtilities.diff
git diff master multicolumn -- Photochemistry/src/Crosssections.jl > Crosssections.diff
git diff master multicolumn -- Photochemistry/src/FileIO.jl > FileIO.diff
git diff master multicolumn -- Photochemistry/src/JuliaODEsolver.jl > JuliaODEsolver.diff
git diff master multicolumn -- Photochemistry/src/Photochemistry.jl > Photochemistry.diff
git diff master multicolumn -- Photochemistry/src/Plotting.jl > Plotting.diff
git diff master multicolumn -- Photochemistry/src/ReactionNetwork.jl > ReactionNetwork.diff
git diff master multicolumn -- Photochemistry/src/UnitConversions.jl > UnitConversions.diff
git diff master multicolumn -- CONSTANTS.jl > CONSTANTS.diff
git diff master multicolumn -- MODEL_SETUP.jl > MODEL_SETUP.diff
git diff master multicolumn -- converge_new_file.jl > converge_new_file.diff
git diff master multicolumn -- PLOT_STYLES.jl > PLOT_STYLES.diff
git diff master multicolumn -- INPUT_PARAMETERS.jl > INPUT_PARAMETERS.diff
git diff master multicolumn -- SeasonsAnalysis/src/MakePlots.jl > MakePlots.diff
git diff master multicolumn -- SeasonsAnalysis/src/Support.jl > Support.diff
git diff master multicolumn -- Solar\ spectra/scale_solar_spectrum.py > scale_solar_spectrum.diff
<!-- git diff master multicolumn -- SeasonsAnalysis/src/SeasonsAnalysis.jl > SeasonsAnalysis.diff -->
<!-- git diff master multicolumn -- scale_solar_spectrum.py > scale_solar_spectrum.diff -->

- From master branch

<!-- OPTION 1 (master in red and multicolumn in green) -->
git diff master..multicolumn -- Photochemistry/src/Core.jl > Core.diff
git diff master..multicolumn -- Photochemistry/src/AnalyzeChemAndTransport.jl > AnalyzeChemAndTransport.diff
git diff master..multicolumn -- Photochemistry/src/BasicUtilities.jl > BasicUtilities.diff
git diff master..multicolumn -- Photochemistry/src/Crosssections.jl > Crosssections.diff
git diff master..multicolumn -- Photochemistry/src/FileIO.jl > FileIO.diff
git diff master..multicolumn -- Photochemistry/src/JuliaODEsolver.jl > JuliaODEsolver.diff
git diff master..multicolumn -- Photochemistry/src/Photochemistry.jl > Photochemistry.diff
git diff master..multicolumn -- Photochemistry/src/Plotting.jl > Plotting.diff
git diff master..multicolumn -- Photochemistry/src/ReactionNetwork.jl > ReactionNetwork.diff
git diff master..multicolumn -- Photochemistry/src/UnitConversions.jl > UnitConversions.diff
git diff master..multicolumn -- CONSTANTS.jl > CONSTANTS.diff
git diff master..multicolumn -- MODEL_SETUP.jl > MODEL_SETUP.diff
git diff master..multicolumn -- converge_new_file.jl > converge_new_file.diff
git diff master..multicolumn -- PLOT_STYLES.jl > PLOT_STYLES.diff
git diff master..multicolumn -- INPUT_PARAMETERS.jl > INPUT_PARAMETERS.diff
git diff master..multicolumn -- SeasonsAnalysis/src/MakePlots.jl > MakePlots.diff
git diff master..multicolumn -- SeasonsAnalysis/src/Support.jl > Support.diff
git diff master..multicolumn -- "Solar spectra/scale_solar_spectrum.py" > scale_solar_spectrum.diff

<!-- OPTION 2 (master in green and multicolumn in red, useful after resolving merge conflicts) -->
git diff multicolumn..master -- Photochemistry/src/Core.jl > Core.diff
git diff multicolumn..master -- Photochemistry/src/AnalyzeChemAndTransport.jl > AnalyzeChemAndTransport.diff
git diff multicolumn..master -- Photochemistry/src/BasicUtilities.jl > BasicUtilities.diff
git diff multicolumn..master -- Photochemistry/src/Crosssections.jl > Crosssections.diff
git diff multicolumn..master -- Photochemistry/src/FileIO.jl > FileIO.diff
git diff multicolumn..master -- Photochemistry/src/JuliaODEsolver.jl > JuliaODEsolver.diff
git diff multicolumn..master -- Photochemistry/src/Photochemistry.jl > Photochemistry.diff
git diff multicolumn..master -- Photochemistry/src/Plotting.jl > Plotting.diff
git diff multicolumn..master -- Photochemistry/src/ReactionNetwork.jl > ReactionNetwork.diff
git diff multicolumn..master -- Photochemistry/src/UnitConversions.jl > UnitConversions.diff
git diff multicolumn..master -- CONSTANTS.jl > CONSTANTS.diff
git diff multicolumn..master -- MODEL_SETUP.jl > MODEL_SETUP.diff
git diff multicolumn..master -- converge_new_file.jl > converge_new_file.diff
git diff multicolumn..master -- PLOT_STYLES.jl > PLOT_STYLES.diff
git diff multicolumn..master -- INPUT_PARAMETERS.jl > INPUT_PARAMETERS.diff
git diff multicolumn..master -- SeasonsAnalysis/src/MakePlots.jl > MakePlots.diff
git diff multicolumn..master -- SeasonsAnalysis/src/Support.jl > Support.diff
git diff multicolumn..master -- "Solar spectra/scale_solar_spectrum.py" > scale_solar_spectrum.diff
