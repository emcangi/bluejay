# Implementation of Horizontal Transport in bluejay (2-D Multicolumn Photochemical Model for Venus and Mars)

## Overview

Horizontal transport is crucial in accurately modeling atmospheric chemistry and dynamics on planets such as Venus and Mars. The updated code implements horizontal transport by combining advection (transport driven by horizontal winds) and diffusion (transport driven by concentration gradients). This document details the implementation approach, equations, and routines involved in achieving horizontal transport within our 2-D multicolumn photochemical model.

## Key Components and Descriptions

### 1. `speciesbclist_horiz`

Defined in `MODEL_SETUP.jl`, this dictionary specifies boundary flux conditions for each chemical species at the domain edges. By default, both front and back edges have zero-flux conditions:

```julia
for sp in all_species
    auto_speciesbclist_horiz[sp] = Dict("f"=>[[0.0 for ialt in 1:num_layers] for c in 1:2])
end
```

This structure allows for easy modification to impose non-zero flux conditions:

- Positive values represent influx into the domain.
- Negative values represent outflux from the domain.

### 2. `boundaryconditions_horiz`

This routine builds the altitude-dependent boundary condition dictionary (`bc_dict_horiz`) from `speciesbclist_horiz`. For periodic domains (`cyclic=true`), all edge fluxes are set to zero to simulate continuity at domain edges:

```julia
bc_dict_horiz[sp][ialt] .= 0.0
```

For non-periodic domains, the imposed boundary fluxes at each altitude (`ialt`) are:

```julia
f_backedge  = [0, -back_flux / GV.dx]
f_frontedge = [0, front_flux / GV.dx]
```

This ensures that specified edge fluxes correctly inject or remove material.

### 3. `fluxcoefs_horiz`

The core computational routine that calculates horizontal transport coefficients (`fluxcoef_dict`) based on diffusion and wind-driven advection:

**Horizontal Diffusion:** Diffusion coefficients between adjacent columns are averaged:

$$
K_{avg} = \frac{K_{i} + K_{i \pm 1}}{2}, \quad D_{avg} = \frac{D_{i} + D_{i \pm 1}}{2}
$$

Horizontal diffusion flux coefficient:

$$
\text{diffusion flux} = \frac{K_{avg} + D_{avg}}{\Delta x^2}
$$

**Horizontal Advection (Upwind Scheme):** Upwind advection ensures numerical stability by considering the wind direction:

- Forward advection coefficient:

$$
\text{adv front} = \frac{\max(v_{local},0) + \max(-v_{front},0)}{\Delta x}
$$

- Backward advection coefficient:

$$
\text{adv back} = \frac{\max(-v_{local},0) + \max(v_{back},0)}{\Delta x}
$$

### Scaling of Horizontal Transport Coefficients

Horizontal transport coefficients are constructed using standard finite‑difference formulas. Diffusion terms arise from discretizing the second spatial derivative in Fick’s law: \$\partial n/\partial t \propto (K+D),\partial^{2}n/\partial x^{2}\$. A centered second derivative contributes a factor of \$\Delta x^{2}\$ in the denominator, so the discrete diffusion coefficient becomes \$(K+D)/\Delta x^{2}\$.

Advection represents the first‑derivative term, \$\partial n/\partial t = -v,\partial n/\partial x\$. Its finite‑difference form introduces only one power of \$\Delta x\$ in the denominator, leading to an advection coefficient of \$v/\Delta x\$.

Consequently **diffusion scales with \$1/\Delta x^{2}\$ whereas advection scales with \$1/\Delta x\$.**

### 4. `update_horiz_transport_coefficients`

This high-level function:

- Computes eddy and molecular diffusion coefficients via `update_diffusion_and_scaleH`.
- Calls `fluxcoefs_horiz` to calculate horizontal transport coefficients (`fluxcoefs_horiz_all`).
- Assembles forward (`tforwards`) and backward (`tbackwards`) transport arrays.
- Applies boundary conditions via `boundaryconditions_horiz`, assembling edge boundary coefficients (`tbackedge`, `tfrontedge`).

These assembled arrays (`tforwards`, `tbackwards`, `tbackedge`, `tfrontedge`) directly integrate with the solver for modeling horizontal transport between columns.

## Integration and Workflow

The horizontal transport is integrated as follows:

- The horizontal wind velocity profiles (`horiz_wind_v`) specified in `MODEL_SETUP.jl` determine advection.
- Edge fluxes defined in `speciesbclist_horiz` inject/remove material at boundaries.
- `update_horiz_transport_coefficients` combines these inputs to produce horizontal transport coefficients required by the solver.

## Numerical Verification

The implementation was validated using a test script (`horizontal_transport_test.jl`) with an analytic solution. This test ensures accurate transport coefficients and validates the correctness of indexing and numerical schemes employed.

**Example:**

- Two altitude bins and two horizontal columns.
- Horizontal cell width $(\Delta x = 1\,\text{cm})$, wind speed $(v = 10\,\text{cm/s})$:

$$
\text{Advection rate} = \frac{v}{\Delta x} = 10\,\text{s}^{-1}
$$

Transport matrix comparison verifies agreement with the analytic solution.

## Strengths of the Current Approach

The model implements horizontal transport by introducing separate horizontal columns coupled through diffusion and wind‑driven advection.

- Upwind Scheme for Advection: This is a robust numerical method for ensuring stability, especially when handling sharp concentration gradients and high wind speeds.
- Averaging Diffusion Coefficients: Taking the average of diffusion coefficients between adjacent columns is an effective way to represent continuous horizontal mixing accurately.
- Modular Implementation: Clearly defined routines (`speciesbclist_horiz`, `fluxcoefs_horiz`, etc.) make the approach adaptable and easy to debug or modify.
- Flexible Boundary Conditions: Allowing for periodic and non-periodic boundary conditions, including customizable fluxes, adds significant versatility.

## Potential Improvements or Alternatives

Implementing more advanced advection schemes or coupling with external wind data could enhance accuracy for cases with strong horizontal variability.

- Higher-order advection: Upwind schemes are robust but only first-order accurate. Using higher-order methods (e.g., flux limiters, TVD schemes) could reduce numerical diffusion, especially when strong gradients occur across columns.
- Variable or spatially dependent horizontal winds: The current implementation allows a wind profile per column. More sophisticated models often include spatial variability in the advection term (e.g., using winds from a GCM or prescribing shear with altitude). Extending the wind arrays or reading them from external data would improve realism.
- Mass-conserving diffusion: The code already averages diffusion coefficients, but ensuring strict global mass conservation (especially with cyclic boundaries) might require verifying that flux leaving one column is exactly balanced by the flux entering the next. Checks or constraints in `fluxcoefs_horiz` could reinforce this.
- Two-way coupling with vertical transport: Horizontal transport is computed separately from vertical flux coefficients. For atmospheres where horizontal transport interacts with vertical mixing (e.g., along isentropes), a more integrated solver might be beneficial.

## Summary

The implemented horizontal transport scheme in the 2-D multicolumn photochemical model effectively combines advection (using an upwind numerical scheme) and diffusion (via averaged diffusion coefficients). The modular design, flexible boundary conditions, and comprehensive verification provide a robust and adaptable framework suitable for modeling planetary atmospheres such as those of Venus and Mars.

## References

Horizontal transport is implemented using an upwind numerical scheme for advection, consistent with methods described by LeVeque (2002) and Hundsdorfer & Verwer (2003), with averaged diffusion coefficients following Brasseur & Solomon (2005) and Jacobson (2005). Similar approaches have been successfully applied in planetary atmospheric models (e.g., Lefèvre et al., 2004).

1. **General Photochemical and Transport Modeling References**
   - Yung, Y. L., & DeMore, W. B. (1999). *Photochemistry of Planetary Atmospheres*. Oxford University Press. (Provides foundational theory and examples of photochemical models for planetary atmospheres, including basic transport equations.)
   - Jacobson, M. Z. (2005). *Fundamentals of Atmospheric Modeling* (2nd ed.). Cambridge University Press. (Detailed coverage of numerical transport schemes, including horizontal advection and diffusion, with explicit numerical formulations.)

2. **Upwind Numerical Advection Schemes**
   - Hundsdorfer, W., & Verwer, J. G. (2003). *Numerical Solution of Time-Dependent Advection-Diffusion-Reaction Equations*. Springer Series in Computational Mathematics, Vol. 33. (Clear and thorough description of upwind advection schemes and stability considerations.)
   - LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press. (Comprehensive reference for numerical flux schemes, including detailed discussion of upwind and related numerical methods.)

3. **Diffusion Coefficients and Numerical Averaging**
   - Brasseur, G. P., & Solomon, S. (2005). *Aeronomy of the Middle Atmosphere* (3rd ed.). Springer Netherlands. (Discusses diffusion in planetary atmospheres, averaging techniques for diffusivity between adjacent grid cells, and practical considerations in numerical models.)
   - Seinfeld, J. H., & Pandis, S. N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change* (2nd ed.). John Wiley & Sons. (Detailed treatments of numerical diffusion approaches, stability, and accuracy concerns in horizontal transport.)

4. **Planetary Atmosphere Modeling Examples**
   - Krasnopolsky, V. A. (2019). "Spectroscopy and Photochemistry of Planetary Atmospheres and Ionospheres: Mars, Venus, Titan, Triton and Pluto." *Cambridge Planetary Science, Series Number 23* (Reviews photochemical modeling frameworks used specifically for Venus, Mars, and other planetary atmospheres, including transport considerations.)
   - Lefèvre, F., Lebonnois, S., Montmessin, F., & Forget, F. (2004). "Three-dimensional modeling of ozone on Mars." *Journal of Geophysical Research: Planets, 109*(E7), E07004. <https://doi.org/10.1029/2004JE002268>
   <!-- (A practical example of implementing horizontal transport with advection-diffusion schemes on Mars.) -->
   <!-- - Montmessin, F. et al., (2011). "A layer of ozone detected in the nightside upper atmosphere of Venus." *Icarus, 216*, 82-85. <https://doi.org/10.1016/j.icarus.2011.08.010>
   - Montmessin, F., & Lefèvre, F. (2013). "Transport-driven formation of a polar ozone layer on Mars." *Nature Geoscience, 6*, 930–933. <https://doi.org/10.1038/ngeo1957> (Detailed application of horizontal transport modeling including wind-driven advection and diffusion on Mars.) -->

# Summary of Tsai et al. (2024): Global Chemical Transport on Hot Jupiters

## Overview
Tsai *et al.* (2024) introduce the 2‑dimensional (horizontal and vertical) photochemical transport model, **VULCAN 2D**, designed to analyze global chemical transport dynamics in tidally locked hot‑Jupiter atmospheres. Their model addresses shortcomings in previous pseudo‑2D models by accommodating non‑uniform horizontal (zonal) winds derived from General Circulation Models (GCMs).

## Model Formulation
### Continuity Equation
VULCAN 2D solves the continuity equation in 2-D Cartesian coordinates (horizontal $`x`$ and vertical $`z`$):

$$
\frac{\partial n(x,z,t)}{\partial t}=P-L-\frac{\partial \phi_z}{\partial z}-\frac{\partial \phi_x}{\partial x}.
$$

where

* $n$ – species number density  
* $P$ – chemical production rate  
* $L$ – chemical loss rate  
* $\phi_z$, $\phi_x$ – vertical and horizontal fluxes, respectively

### Vertical Transport Flux ($\phi_z$)
Vertical flux comprises

* vertical advection  
* eddy diffusion  
* molecular diffusion  

### Horizontal Advection Flux (upwind scheme)
The first‑order upwind differencing used for horizontal advection is

$$
\phi_{x-\frac12} =
\begin{cases}
  v_{k-\frac12}\,n_{k-1}, & v_{k-\frac12}>0 \\
  v_{k-\frac12}\,n_k, & v_{k-\frac12}<0
\end{cases}
$$

$$
\phi_{x+\frac12} =
\begin{cases}
  -\,v_{k+\frac12}\,n_k, & v_{k+\frac12}>0 \\
  -\,v_{k+\frac12}\,n_{k+1}, & v_{k+\frac12}<0
\end{cases}
$$

where $v_{k\pm\frac12}$ is the horizontal wind speed at the cell interface.

### Horizontal Diffusion
Horizontal diffusion is implemented but primarily relevant for meridional (north–south) transport on giant planets. In the zonal (east–west) direction of hot‑Jupiter atmospheres, strong advection usually dominates.

## Key Results and Insights
* **Validation** – model outputs are consistent with analytical solutions and previous pseudo‑2D/3D‑GCM studies.  
* **Mixing regimes** – horizontal transport dominates below $`\sim0.1\,\text{mbar}`$, whereas vertical mixing dominates above this level.  
* **Chemical gradients** – photochemically active species ($`\mathrm{CH}_4`$, $`\mathrm{NH}_3`$, $`\mathrm{HCN}`$) show strong longitudinal gradients; relatively inert species ($`\mathrm{H_2O}`$, $`\mathrm{CO}`$) remain nearly uniform.
* **Limb asymmetry** – horizontal transport drives significant morning–evening limb composition differences, especially in carbon‑rich atmospheres (e.g. HD 209458 b with supersolar C/O).  

## Limitations of Pseudo‑2D Models
Assuming uniform zonal winds can introduce large errors when day–night flows or wind variability are strong. Tsai *et al.* quantify these errors and demonstrate the need for true 2‑D treatment.

## Lessons for *bluejay* Model
The current horizontal‑transport implementation is already robust, but the following refinements could improve clarity and maintainability:

1. **Explicit interface‑flux notation**  
   $$\phi_{x-\frac12}=v_{x-\frac12}\,n_{\text{upwind}}$$  

2. **Isobaric vertical coordinate**  
   $$z=-H\ln\!\left(\frac{p}{P_s}\right)$$  
   ensures fluxes are computed on constant‑pressure levels, matching common practice.

3. **Reduced explicit horizontal diffusion** – rely on advection unless meridional mixing requires otherwise.

4. **CFL‑limited timestep**  
   $$\Delta t<\min\!\left(\frac{\Delta x}{v_x}\right)$$  
   to preserve numerical stability.

## Conclusion
The present framework is sound. Minor tweaks—explicit flux notation and clear isobaric coordinates—will bring it fully in line with the methodology of Tsai *et al.* (2024).

## Reference
Tsai, S.‑M. *et al.* (2024), “Global Chemical Transport on Hot Jupiters: Insights from the 2‑D VULCAN Photochemical Model,” *ApJ*, **963**, 41. <https://doi.org/10.3847/1538-4357/ad1600>
