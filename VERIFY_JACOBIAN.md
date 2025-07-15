# Verification of Jacobian and Rate Functions

This note documents that the multicolumn model keeps a consistent ordering of species, altitude, and horizontal columns when building transport matrices.

## Indexing order
- `flatten_atm` builds vectors by species → altitude → horizontal.
- `unflatten_atm` reshapes back with `(length(species_list), GV.num_layers, n_horiz)`.

## Horizontal offsets
- In `chemJmat`, the base index for a given column/altitude uses `(ihoriz-1)*(length(GV.active_longlived)*GV.num_layers)`.
- When connecting neighboring columns, the Jacobian indices shift by `GV.num_layers*length(GV.active_longlived)`.

## Reshaping in rate functions
- `ratefn` reshapes inputs using `(species, altitude, horizontal)` ordering before looping over columns and altitudes.

## Temperature arrays
- `MODEL_SETUP.jl` creates constant 2‑D temperature profiles `Tn_arr`, `Ti_arr`,
  and `Te_arr` with shape `(n_horiz, num_layers+2)`.
- Each array is accessed as `Tn_arr[ihoriz, ialt]` so the horizontal index comes
  first followed by altitude.
- `Tplasma_arr` stores `Ti_arr .+ Te_arr` for computing ion diffusion.

## Solar fluxes
- `optical_depth` builds `solarabs` with shape `[n_horiz][num_layers]`, where
  each element is a vector of optical depths across all wavelengths.
- `update_Jrates!` loops over columns and altitudes, replacing
  `solarabs[ihoriz][ialt]` with the local actinic flux
  `GV.solarflux[:, 2] .* exp.(-solarabs[ihoriz][ialt])` before integrating with
  `GV.crosssection[j][ihoriz][ialt+1]` to compute Jrates.
<!-- - `optical_depth` returns `solarabs`, a vector over columns containing altitude
  arrays, i.e. `solarabs[ihoriz][ialt]` is an array of wavelength points.
- `update_Jrates!` multiplies `GV.solarflux[:, 2]` by `exp.(-solarabs[ihoriz][ialt])`
  for each column and altitude before integrating with column‑specific cross
  sections. -->

These observations verify that the original Julia code consistently follows the species → altitude → horizontal convention across the major files.