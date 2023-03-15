# bluejay: atmospheric photochemistry in Julia

## Installation 

Jay is currently up-to-date with Julia 1.8.5. It will likely work with newer versions, but it has not yet been tested. 

At this time, there are no compiled binaries. The model is provided as a collection of modules (under the Photochemistry module) and related scripts. To install, simply ensure that Julia is installed in the location of your choice, and fork the repo. The root directory of the model must contain:

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
- **Atmosphere.jl**: Functions that manipulate the atmosphere matrix in some way, or perform some calculation. Includes column densities, electron densities, find_exobase, calculation of total density, optical depth, and scale height.
- **BasicUtilities.jl**: non-photochemistry or atmosphere specific functions. Mostly very short helper functions and string manipulation.
- **Core.jl**: The photochemical model core. Functions to manipulate the atmosphere matrix, set up the chemical jacobian and chemistry networks, functions which calculate escape, functions to utilize photochemical equilibrium, basic transport functions.
- **Crosssections.jl**: Functions which load cross section data and populate the dictionary of cross sections. Note, the used J rates are hard-coded in at this time. This could stand to be improved in the future.
- **FileIO.jl**: Anything that involves getting info out of or writing info into a file, including saving the model results and parameter logs.
- **JuliaODEsolver.jl**: Optional submodule with functions to utilize the Julia ODE solvers. Not presently used in our work as the Gear solver (included in the `convege_new_file` script) works better.
- **Plotting.jl**: Functions which make plots for showing results!
- **ReactionNetwork.jl**: Everything needed to ingest an Excel spreadsheet of reaction rate data and turn it into a symbolic chemical network for Julia to read. Also includes functions which calculate enthapies of reaction and modify a reaction rate spreadsheet. These functions are NOT called by the model directly and MUST be run by the user when establishing a new simulation, say, for a new planet, because escape energy changes with different gravities.
- **Temperature.jl**: Functions which set up the temperature profile of the atmosphere. Much of this relies on fits to specific data and should be modified for new models of new planets or Mars in other scenarios.
- **Water.jl**: Functions that set up the water profile and manipulate it in various ways. Similar to Temperature, these functions should be carefully manipulated by the user based on the model goals.

    
## Usage

**Initial steps**:
- Ensure that the root folder contains the `Photochemistry` folder, the `uvxsect` folder, and the running scripts:
  - `converge_new_file.jl`
  - `CUSTOMIZATIONS.jl`: Some basic model parameters that don't change much, as well as names of photochemical cross section data files.
  - `PARAMETERS.jl`: Parameters that the user can change before running the model. 
  - `CONSTANTS.jl`: Global constants. Only needs to be changed if changing the planet.
- Ensure the root folder also contains critical inputs:
  - `REACTION_NETWORK.jl`, reaction rate data for the chemical network
  - `INITIAL_GUESS.jl`, an initial guess for species densities by altitude
- Change the value of `results_dir` to the names of folders you would like to use for the results output.

**If modeling a new planet**:
1. Change constants in `CONSTANTS.jl`
2. Run `modify_rxn_spreadsheet` (from the `ReactionNetwork.jl` submodule) at the command line

**Running the model**:
1. Modify `PARAMETERS.jl` to your chosen conditions for the simulation. Note: You may need to modify the `Temperature.jl` and `Water.jl` modules in nontrivial ways.
2. At the command line, navigate to the model root folder and run `julia converge_new_file.jl`. 

## Credits

The core of Jay was originally written by <a href="https://github.com/planetarymike">Michael Chaffin</a> under <a href="https://github.com/planetarymike/chaffin_natgeo_mars_photochemistry">this repo</a>. It was expanded to include deuterium chemistry by Eryn Cangi in <a href="https://github.com/emcangi/dh_fractionation">this iteration of the model</a>. 

## Citing and usage

You may fork and use or expand this model for your own research purposes. Any usage of the model should be acknowledged with a citation. You may use this Bibtex entry:

```
@misc{Cangi2022,
  author       = {Eryn M Cangi and
  Mike Chaffin},
  title        = {emcangi/dh\_ions: Jay},
  month        = dec,
  year         = 2022, 
  publisher    = {Zenodo},
  version      = {v1.2}, #<----Fill in with the most current version, please
  doi          = {10.5281/zenodo.7662300}, #<----Fill in with the most current DOI, shown in the badge below.
  url          = {https://doi.org/10.5281/zenodo.7662300} #<----Fill in with the most current DOI, shown in the badge below.
}
```

Click here to access the Zenodo DOI. This badge will display the most current release's DOI: [![DOI](https://zenodo.org/badge/285645653.svg)](https://zenodo.org/badge/latestdoi/285645653)
  
If you do not use LaTeX, please use the above information and adapt it to your preferred style.

## Getting help and contributing

This model was developed as part of Eryn's PhD, and I am not a professional software developer. If you require any help in using the model, encounter bugs, or would like to contribute, please open a Github issue.

## License

This software is distributed under the <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GNU General Public License (GPL) 3.0</a>. Please see file COPYING.txt for the full license description.

## Name

Mike thought that a bird-themed naming scheme would be fun. He suggested "Sparrow" for this model, combining "[chemical] species" + "arrow" to signify the motion of species vertically through the atmosphere. Eryn pointed out that "Jay" would be a much punnier name (the reaction rates of photochemical processes are referred to as "J rates", also the model is written in Julia, which starts with J), and Mike pointed out that Jay is also a general name. THe blue jay is well recognized as a bird, and has a short and easy to remember name. Eryn hopes to make a pretty logo soon.
