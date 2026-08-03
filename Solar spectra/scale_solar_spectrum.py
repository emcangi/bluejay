# coding=unicode-escape
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as interp
import datetime
from scipy.optimize import curve_fit

# for converting to photons
h = 6.626e-34
c = 3e8 # m/s

# Conversion functions:
def convert_W_m2_to_photons_per_s_cm2(Wm2, wavelength):
    """
    Converts from W/m^2/nm to photons/s/cm^2.
    multiply the irradiance by λ/hc to get photons; by 1/10000 to convert to 1/cm^2.

    Parameters
    ----------
    Wm2 : array
          Irradiance data in units of W/m2/nm. 
    wavelength : array
                 array of wavelengths at which each flux is recorded, in nm.

    Returns
    ----------
    An array of fluxes in photons/s/cm^2
    """
    return Wm2 * ((wavelength * 10**(-9))/(h*c)) * (1/10000)


def convert_photons_per_s_cm2_to_W_m2(photcm2s, wavelength):
    """
    Presumably this can be used to convert the photon flux back into a W/m^2/nm, 
    but I don't remember why I wrote it. Probably to verify the conversion.
    """
    return photcm2s * (h*c/(wavelength*1e-9)) * (10000/1)


def find_nearest(array, value):
    """Just fines the nearest entry in array closest to value."""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
   
def interpolate_solar_spectrum(spec, planet, show_plots=True, extrap_tail=True,
                               interp_start=0.5, interp_end=2389.51, dl=1, desctag=""):
    """
    Parameters
    ---------- 
    spec : Pandas dataframe
           Contains solar spectrum data from 0.1-2390 nm or so. The binning may be arbitrary.
           Columns:
               λ: wavelength in nm
               irradiance: irradiance in W/m^2/nm
    planet : string
             planet name
    show_plots : bool
                 whether to print out associated plots to the window manager
    extrap_tail : bool
                  whether to add an extrapolation from 2390-2400 (some datasets cut off around 2390)
    interp_start, interp_end : floats
                  Wavelengths in between which, the code will interpolate. I.e., this is the limit of
                  the original source spectrum.
    dl : int
         interpolation spacing in nm
    desctag : string
              A short, descriptive explanation of the spectrum's applicability.
    Returns
    ----------
    whole_spectrum : Pandas DataFrame
                     Wavelengths from 0.5-2399.5 nm, with bin centers at integer wavelengths in nm
                     Rebinned irradiances, units unchanged
    """

    col = [spec.columns[0], spec.columns[1]]
    
    f_irr = interp.interp1d(spec[col[0]], spec[col[1]])

    newx = np.arange(interp_start, interp_end, dl)
    new_irr = f_irr(newx)

    # check to make sure it looks right
    if show_plots:
        plt.figure(figsize=(10,5))
        plt.title("Check that spectrum looks okay after interpolation")
        plt.plot(spec[col[0]], spec[col[1]], color="blue", linewidth=5, label="data", alpha=0.5)
        plt.plot(newx, new_irr, color="black", label="interpolated")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Irradiance (W m^-2 nm^-1)")
        plt.legend()
        plt.show()

    # make a new dataframe containing the interpolated data
    interp_data = pd.DataFrame(np.column_stack((newx, new_irr)), columns=[col[0], col[1]])

    # put the interpolated data together with the data that started out fine for one new data frame 
    i_interp = find_nearest(spec[col[0]], interp_start)
    newsolardata = pd.concat([spec[:i_interp+1], interp_data[1:]], ignore_index="true")
    
    if extrap_tail:
        # Model wants wavelengths up to 2399.5 nm, but using the base interpolation limit of 2389.51,
        # some wavelengths are missing. Here we fit the end of the tail so we can get it out to the same wavelength
        def func(x, a, b, c):
            return a * np.exp(b * x) + c

        datatofit = newsolardata[1000:] # WARNING: this 1000 is hard coded and arbitrary
        print("WARNING: Fitting the red-ward irradiance tail using data points from position 1000 to the end. This is arbitrary and works for the Sun.")
        print("Change the code near this line number if it looks bad.")
        popt, pcov = curve_fit(func, datatofit[col[0]], datatofit[col[1]], p0=(1.5, -0.0001, 1))
        if show_plots:
            plt.figure(figsize=(10,5))
            plt.title("Fitting the long irradiance tail")
            plt.plot(newsolardata[col[0]], newsolardata[col[1]], label="data + interpolated values")
            plt.plot(datatofit[col[0]], func(datatofit[col[0]], popt[0], popt[1], popt[2]), color="purple", label="fit to irradiance")
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Irradiance (W m^-2 nm^-1)")
            plt.legend()
            plt.show()
        
        # make a dataframe containing extrapolated data
        extrap_wavelengths = np.arange(2390.5, 2400.5, 1)
        extrap_irradiance = func(extrap_wavelengths, popt[0], popt[1], popt[2])
        tail = pd.DataFrame(np.column_stack((extrap_wavelengths, extrap_irradiance)), columns=[col[0], col[1]])
        if show_plots:
            plt.figure(figsize=(10,5))
            plt.title("Shows final extrapolated tail")
            plt.plot(newsolardata[col[0]], newsolardata[col[1]], color="green", label="data + interpolated values")
            plt.plot(tail[col[0]], tail[col[1]], color="purple", label="extrapolated tail")
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Irradiance (W m^-2 nm^-1)")
            plt.legend()
            plt.show()

        whole_spectrum = pd.concat([newsolardata, tail], ignore_index=True)
    else:
        whole_spectrum = newsolardata
    
    return whole_spectrum

    
def scale_to_new_distance(spectrum_df, star_to_target_distance=1, star_to_instrument_distance=1, non_scaled_wavelengths=None):
    """
    Given some spectrum in photon/s/nm, scale the spectrum to a new AU to allow simulation of 
    orbital scenarios for which we don't have data.

    FOR EXAMPLE: If your spectrum was originally recorded by an Earth instrument and you want 
    a Mars spectrum for average conditions, you would enter star_to_target_distance=1.524 and
    star_to_instrument_distance=1. 

    If both star_to_instrument_distance and star_to_target_distance are 1, no scaling occurs.

    Parameters
    ----------
    spectrum_df : Pandas DataFrame with columns:
                  wavelength (nm)
                  photon flux (phot/s/cm^2/nm)
    star_to_target_distance : float
                              Distance (in AU) between the stellar source and the notional planet, i.e., 
                              the situation you want to simulate that you don't have data for.
                              Default = 1. 
    star_to_instrument_distance : float
                                  Distance (in AU) between the stellar source and the original instrument 
                                  that measured it. Default is 1, assuming an Earth orbit.
    non_scaled_wavelengths : list of floats
                           Wavelengths (in nm) between which we will NOT apply the scaling factor. 
                           This allows you to provide data at your target from a real instrument, thus
                           without having to scale it from a data set from an instrument somewhere else. This was
                           originally included because I was using spectra that were a composite of data 
                           from MAVEN/EUVM (so, taken at Mars' orbital distance)  and SORCE/SOLSTICE 
                           (taken at Earth's orbital distance) to make solar spectra for Mars. So I didn't need 
                           to scale the Mars data, just the Earth-based data.
                           IF YOU AREN'T USING PLANETARY MISSION DATA - LEAVE THIS AS NONE!

    Returns
    ----------
    spectrum_df : Pandas DataFrame
                  The same dataframe, but now with the photons/s/nm scaled to the AU you asked for. 


    """

    sq_distance_ratio = (star_to_instrument_distance**2)/(star_to_target_distance**2) 

    if non_scaled_wavelengths is not None:
        scale_rows = ~spectrum_df["wavelength (nm)"].between(non_scaled_wavelengths[0], non_scaled_wavelengths[1])
        spectrum_df.loc[scale_rows, "photon flux (phot/s/cm^2/nm)"] = spectrum_df.loc[scale_rows, "photon flux (phot/s/cm^2/nm)"] * sq_distance_ratio
    else:
        spectrum_df.loc[:, "photon flux (phot/s/cm^2/nm)"] = spectrum_df.loc[:, "photon flux (phot/s/cm^2/nm)"] * sq_distance_ratio

    return spectrum_df

# MAIN =======================================================================================================================
# Get user input -------------------------------------------------------------------------------------------------------------
planet_name = input("Enter planet name (the planet for which you want to create a spectrum): ") # This just ends up in the final filenames.
solarfile = input("Enter the file with solar spectrum data in W/m^2/nm: ") # Full path name.

# Distances
theAU = float(input("Enter AU at which you'd like the output: ")) # This is the AU you want to scale to, for which we don't have data
data_AU = input("Enter the AU of the distance between source and observer for the base spectrum (just press enter to use default of 1 AU): ")
if data_AU == "":
    data_AU = 1
else:
    data_AU = float(data_AU)

descriptive_tag = input("Enter a descriptive tag for this spectrum (please use underscores): ") # goes into filename
in_situ_data = input(f"Is any of the data taken by a mission that orbits {planet_name}? (Please enter True or False): ")
if (in_situ_data=="t") or (in_situ_data=="True") or (in_situ_data=="T"):
    min_data_wavelength = float(input("Please enter the minimum wavelength (nm) of your in-situ data: "))
    max_data_wavelength = float(input("Please enter the maximum wavelength (nm) of your in-situ data: "))
    non_scaled_wavelengths = [min_data_wavelength, max_data_wavelength]
    print(f"In situ data used between {min_data_wavelength} and {max_data_wavelength} nm. Those data will not be scaled by AU. "
           + "In all other regions, data will be scaled by (1/{theAU})^2.")
else:
    non_scaled_wavelengths = None
    print(f"Scaling all data to a distance of {theAU} AU.")

show_plots = input("Do you want to visually inspect the plots and save the final result? (Please enter True or False): ")
if (show_plots=="t") or (show_plots=="True") or (show_plots=="T"):
    show_plots = True
else:
    show_plots = False

# Count header rows in the source file, hopefully they start with # or this will silently fail
hdrcount = 0
with open(solarfile) as f: 
    for line in f:
        if line.startswith("#"):
            hdrcount += 1

# Load data, interpolate ----------------------------------------------------------------------------------------------------
print(f"loading {solarfile}")
solarspec = np.loadtxt(solarfile, skiprows=hdrcount)
solarspec_df = pd.DataFrame(solarspec, columns=["wavelength (nm)", "irradiance (W/m^2/nm)"])
solarspec_df_interp = interpolate_solar_spectrum(solarspec_df, planet_name, show_plots=show_plots,
                                                 desctag=descriptive_tag)

# Convert the units 
solarspec_df_interp["photon flux (phot/s/cm^2/nm)"] = np.asarray(list(map(convert_W_m2_to_photons_per_s_cm2, 
                                                                         solarspec_df_interp["irradiance (W/m^2/nm)"], # Flux
                                                                         solarspec_df_interp["wavelength (nm)"] # Wavelengths
                                                                         )
                                                                    )
                                                                )

# Scale to a new AU
solarspec_df_interp = scale_to_new_distance(solarspec_df_interp, star_to_target_distance=theAU, star_to_instrument_distance=data_AU, 
                                            non_scaled_wavelengths=non_scaled_wavelengths)


# Show the final spectrum
if show_plots:
    plt.figure(figsize=(10,5))
    plt.title("Finalized spectrum, in photon flux, with interpolation and extrapolation")
    plt.plot(solarspec_df_interp["wavelength (nm)"], solarspec_df_interp["photon flux (phot/s/cm^2/nm)"])
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Photons (ph cm^-2 s^-1)")
    plt.savefig(f"final_spectrum_{planet_name}_{descriptive_tag}.png", bbox_inches="tight")
    plt.show()
    print()


# Write out the result to a file --------------------------------------------------------------------------------------------

# Build a custom header and writeout the file, with custom header comments before the column names
header = [f"# AU: {theAU}"+'\\n',
          f"# Source data: {solarfile}"+'\\n',
          f"# File creation date: {datetime.datetime.now().strftime('%B %d, %Y')}"+'\\n\\n']
          
outputfile = f"{planet_name}solarphotonflux_{descriptive_tag}.dat"

with open(outputfile, "w") as f:  
    for line in header:
        f.write(line)
    f.write("# ")
    solarspec_df_interp.to_csv(f, sep='\t', float_format="%.2f", columns=["wavelength (nm)", "photon flux (phot/s/cm^2/nm)"], index=False)

