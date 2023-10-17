# bluejay: atmospheric photochemistry in Julia

## Installation 

bluejay is currently up-to-date with Julia 1.8.5. It will likely work with newer versions, but it has not yet been tested. 

At this time, there are no compiled binaries. The model is provided as a collection of modules (under the Photochemistry module) and related scripts. To install, simply ensure that Julia is installed in the location of your choice and added to your environment $PATH variable, and fork the repo. The root directory of the model must contain:

- The Photochemistry module folder
- The `uvxsect` folder
- All scripts needed to run the model, described later in Usage

### Requirements

Julia packages the model depends on:

- DataFrames
- Dates
- DelimitedFiles
- DifferentialEquations
- GeneralizedGenerated
- HDF5
- IncompleteLU
- IterativeSolvers
- JLD
- LaTeXStrings
- LinearAlgebra
- Printf
- PyPlot
- PyCall
- PlotUtils
- Revise
- SparseArrays
- XLSX
    
### Photochemistry module contents

The Photochemistry module contains several submodules:

- **AnalyzeChemAndTransport.jl**: Functions which help in analyzing the chemistry and transport of a given saved atmosphere file. Here you can compute various volume rates, timescales, the total escape, the fractionation factor, and limiting flux.
- **Atmosphere.jl**: Functions that manipulate the atmosphere matrix in some way, or perform some calculation. Includes column densities, electron densities, find_exobase, calculation of total density, mean mass, optical depth, and scale height.
- **BasicUtilities.jl**: non-photochemistry or atmosphere specific functions. Mostly very short helper functions and string manipulation.
- **Core.jl**: The photochemical model core. Functions to manipulate the atmosphere matrix, set up the chemical jacobian and chemistry networks, functions which calculate escape, functions to utilize photochemical equilibrium, basic transport functions.
- **Crosssections.jl**: Functions which load cross section data and populate the dictionary of cross sections. Note, the J rates are hard-coded in as lists of reactants and products at this time due to the many different ways they are handled. This may be improvable/automatable in the future.
- **FileIO.jl**: Anything that involves getting info out of or writing info into a file, including saving the model results and parameter logs.
- **JuliaODEsolver.jl**: Optional submodule with functions to utilize the Julia ODE solvers. Not presently used in our work as the Gear solver (included in the `convege_new_file` script) works better.
- **Plotting.jl**: Functions which make plots for showing model inputs, model evolution, and results.
- **ReactionNetwork.jl**: Everything needed to ingest an Excel spreadsheet of reaction rate data and turn it into a symbolic chemical network for Julia to read. Also includes functions which calculate enthapies of reaction and modify a reaction rate spreadsheet. These functions are NOT called by the model directly and MUST be run by the user when establishing a new simulation, say, for a new planet, because escape energy changes with different gravities.
- **Temperature.jl**: Functions which set up the temperature profile of the atmosphere. Much of this relies on fits to specific data and should be modified for new models of new planets or Mars in other climate scenarios.
- **Water.jl**: Functions that set up the water profile and manipulate it in various ways. Similar to Temperature, these functions should be carefully manipulated by the user based on the model goals.

    
## Usage

**Initial steps**:
- Ensure that the root folder contains the `Photochemistry` folder, the `uvxsect` folder, and the running scripts:
  - `converge_new_file.jl`
  - `CUSTOMIZATIONS.jl`: Some basic model parameters that don't change much. Importantly, the paths to folders where output will be stored are here.
  - `PARAMETERS.jl`: Parameters that the user can change before running the model. 
  - `CONSTANTS.jl`: Global constants. Only needs to be changed if changing the planet.
- Ensure the root folder also contains critical inputs:
  - `REACTION_NETWORK.jl`, reaction rate data for the chemical network
  - `INITIAL_GUESS.jl`, an initial guess for species densities by altitude
- Set your individual file path and font preferences in CUSTOIMZATIONS.jl:
  - Change the value of `results_dir` to a path you would like to use for the results output
  - Set fontnames you would like to use for sans-serif and monospace fonts. Serif fonts are not modified because no code within the module uses a serif font.

**If modeling a new planet**:
1. Change constants in `CONSTANTS.jl` for the planet in question. If you change any variable names, make sure to Ctrl-F replace them all submodules of Photochemistry. 
2. Use `scale_solar_spectrum.py` (Python, not Julia) to scale the solar min/mean/max spectra to the orbital AU appropriate to your chosen planet. If using for exoplanets, please provide your own solar spectrum in final units of photons/s/cm^2/nm.
3. Run `modify_rxn_spreadsheet` (from the `ReactionNetwork.jl` submodule) at the command line so that enthalpies of reaction will be re-calculated and saved to a new REACTION_NETWORK spreadsheet. WARNING: If you don't do this, non-thermal escape fluxes will NOT be reliable!

**Running the model**:
1. Modify `PARAMETERS.jl` to your chosen conditions for the simulation. Note: You may need to modify the `Temperature.jl` and `Water.jl` modules in nontrivial ways. 
2. At the command line, navigate to the model root folder and run `julia converge_new_file.jl`. 

## Credits

The core of bluejay was originally written by <a href="https://github.com/planetarymike">Michael Chaffin</a> under <a href="https://github.com/planetarymike/chaffin_natgeo_mars_photochemistry">this repo</a>. It was expanded to include deuterium chemistry by Eryn Cangi in <a href="https://github.com/emcangi/dh_fractionation">this iteration of the model</a>. 

## Citing and usage

You may fork and use or expand this model for your own research purposes. All usages of bluejay or results derived from it should be acknowledged with at the very least, a citation; if your use of bluejay constitutes a non-negligible or substantial part of a paper, please extend an offer of coauthorship to Eryn Cangi and Mike Chaffin. 

For citations, please use the "Cite this repository" button on Github, but make sure to replace the DOI with the most current. This badge will display the most current release's DOI and link to the Zenodo repository: [![DOI](https://zenodo.org/badge/285645653.svg)](https://zenodo.org/badge/latestdoi/285645653)
  
## Getting help and contributing

This model was developed as part of my (Eryn's) PhD, so some inefficiencies likely remain. If you require any help in using the model, encounter bugs, or would like to contribute, please open a Github issue.

## License

This software is distributed under the <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GNU General Public License (GPL) 3.0</a>. Please see file COPYING.txt for the full license description.

## Name

Mike thought that a bird-themed naming scheme would be fun. He suggested "Sparrow" for this model, combining "[chemical] species" + "arrow" to signify the motion of species vertically through the atmosphere. I pointed out that "Jay" would be a much punnier name (the reaction rates of photochemical processes are referred to as "J rates", also the model is written in Julia, which starts with J), and Mike pointed out that Jay is also a general name so it may be confusing. The blue jay is well recognized as a bird, and has a short and easy to remember name. I hope to make a pretty logo soon.
