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

def convert_W_m2_to_photons_per_s_cm2(Wm2, wavelength):
    return Wm2 * ((wavelength * 10**(-9))/(h*c)) * (1/10000)

def convert_photons_per_s_cm2_to_W_m2(photcm2s, wavelength):
    return photcm2s * (h*c/(wavelength*1e-9)) * (10000/1)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
   
def interpolate_solar_spectrum(spec, AU, planet, show_plots=True, extrap_tail=True, scale_below=0, scale_above=0, interp_start=0.5, interp_end=2389.51, dl=1, desctag=""):
    """
    Input: 
        spec: Pandas dataframe
              Contains solar spectrum data from 0.1-2390 nm or so. The binning may be arbitrary.
              Columns:
                  λ: wavelength in nm
                  irradiance: irradiance in W/m^2/nm
        AU: float
            Orbital semimajor axis at which to scale the spectrum.
        planet: string
                planet name
        show_plots: bool
                    whether to print out associated plots to the window manager
        extrap_tail: bool
                     whether to add an extrapolation from 2390-2400 (some datasets cut off around 2390)
        scale_below, scale_below: floats
                                  At wavelengths inside [scale_below, scale_above], the data are assumed
                                  to be taken in situ at the planet and to not need to be scaled. 
                                  OUTSIDE this range, data will be scaled by 1/(AU^2). 
                                  e.g., EUVM data from MAVEN at Mars are taken from 0.5--189.5 nm, so when 
                                  EUVM data is included in the input file, one should set scale_below=0.5
                                  and scale_above=189.5.
        interp_start: float
                      a wavelength in nm above which wavelengths will be interpolated to half wavelengths.
        dl: int
            interpolation spacing in nm
        desctag: string
                 A short, descriptive explanation of the spectrum's applicability.
    Output:
        A Pandas dataframe whole_spectrum from 0.5-2399.5 nm, properly binned every half nm, a CSV, and a plot.
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
        # This data is still missing 10 or so datapoints near the end (Mike's solar photon flux goes out to 2399.5).
        # fit the end so we can get it out to the same wavelength
        def func(x, a, b, c):
            return a * np.exp(b * x) + c

        # do the same fit for irradiance for completeness, even though I don't really need the irradiance data
        datatofit = newsolardata[1000:]
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
    
    #make the FINAL dataframe with ALL data for solar max;
    # First convert from W/m^2/nm to the units below, and yes, I have checkd it multiple times lol including on 11/3/22. 
    # multiply the irradiance by λ/hc to get photons; by 1/10000 to convert to cm^2; and 1/AU^2 to convert to planet orbit.
    # whole_spectrum["photon flux (phot/s/cm^2/nm)"] = whole_spectrum[col[1]] * ((whole_spectrum[col[0]] * 10**(-9))/(h*c)) * 1/(AU**2) * (1/10000)
    
    whole_spectrum["photon flux (phot/s/cm^2/nm)"] = np.asarray(list(map(convert_W_m2_to_photons_per_s_cm2, whole_spectrum[col[1]], whole_spectrum[col[0]])))

    # Now scale to the right AU.
    # You can select whether to only scale above a certain wavelength if you 
    # are using MAVEN EUVM data for the short wavelengths!
    if scale_above==0:
        i = 0
    else:
        i = find_nearest(whole_spectrum[col[0]], scale_above)
    
    whole_spectrum.loc[i:, "photon flux (phot/s/cm^2/nm)"] = whole_spectrum.loc[i:, "photon flux (phot/s/cm^2/nm)"] * (1**2)/(AU**2) 

    if show_plots:
        plt.figure(figsize=(10,5))
        plt.title("Finalized spectrum, in photon flux, with interpolation and extrapolation")
        plt.plot(whole_spectrum[col[0]], whole_spectrum["photon flux (phot/s/cm^2/nm)"])
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Photons (ph cm^-2 s^-1)")
        plt.savefig(f"final_spectrum_{planet}_{desctag}.png", bbox_inches="tight")
        plt.show()
        print()
    
    return whole_spectrum

# MAIN =======================================================================================================================
# Get user input -------------------------------------------------------------------------------------------------------------
planet_name = input("Enter planet name: ")
solarfile = input("Enter the file with solar spectrum data in W/m^2/nm: ")
theAU = float(input("Enter AU at which you'd like the output: "))
descriptive_tag = input("Enter a descriptive tag for this spectrum (please use underscores): ")
in_situ_data = input(f"Is any of the data taken by a mission that orbits {planet_name}? (Please enter True or False): ")

if in_situ_data == "True": 
    scale_below = float(input("Please enter the minimum wavelength (nm) of your in-situ data: "))
    scale_above = float(input("Please enter the maximum wavelength (nm) of your in-situ data: "))
    print(f"In situ data used between {scale_below} and {scale_above} nm. Those data will not be scaled by AU. In all other regions, data will be scaled by (1/{theAU})^2.")
else:
    scale_below = 0
    scale_above = 0
    print(f"Scaling all data to a distance of {theAU} AU, using SORCE and TIMED (Earth-based) dataonly.")

# Count header rows in the source file, hopefully they start with # or this will silently fail
hdrcount = 0
with open(solarfile) as f: 
    for line in f:
        if line.startswith("#"):
            hdrcount += 1

# Load data, interpolate, and scale -----------------------------------------------------------------------------------------
print(f"loading {solarfile}")
solarspec = np.loadtxt(solarfile, skiprows=hdrcount)
solarspec_df = pd.DataFrame(solarspec, columns=["wavelength (nm)", "irradiance (W/m^2/nm)"])
solarspec_df_tidy = interpolate_solar_spectrum(solarspec_df, theAU, planet_name, show_plots=True, scale_below=scale_below, scale_above=scale_above, desctag=descriptive_tag)

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
    solarspec_df_tidy.to_csv(f, sep='\t', float_format="%.2f", columns=["wavelength (nm)", "photon flux (phot/s/cm^2/nm)"], index=False)

