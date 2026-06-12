# coding=unicode-escape
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as interp
from copy import deepcopy
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
                                  OUTSIDE this range, data will be scaled by 1/(AU^2), assuming an Earth-
                                  orbit origin for the data, e.g. SORCE.
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
    
    
    return whole_spectrum

# MAIN =======================================================================================================================
# Get user input -------------------------------------------------------------------------------------------------------------
show_plots=False
planet_name = "Mars" # input("Enter planet name: ")
solarfiles = ["solarmin_20081029_SORCE_and_TIMED.dat", 
             "solarmean_2004-02-07_SORCE_and_TIMED.dat", 
             "solarmax_20150604_and_20020322_SORCE_and_TIMED.dat"] # input("Enter the file with solar spectrum data (taken in Earth orbit) in W/m^2/nm: ")
theAU = [1.624, 1.656, 1.537, 1.411, 1.387, 1.483] # float(input("Enter AU at which you'd like the output: "))

use_data = True # input(f"Should we apply data from another instrument at {planet_name}? (Please enter True or False): ")

datafiles = ["euvm_60deg_bins/mvn_euv_l3_Ls30_2023-2-27.csv",
            "euvm_60deg_bins/mvn_euv_l3_Ls90_2023-7-13.csv",
            "euvm_60deg_bins/mvn_euv_l3_Ls150_2023-11-18.csv",
            "euvm_60deg_bins/mvn_euv_l3_Ls210_2024-3-4.csv",
            "euvm_60deg_bins/mvn_euv_l3_Ls270_2022-7-22.csv",
            "euvm_60deg_bins/mvn_euv_l3_Ls330_2022-10-31.csv",] # input(f"Please enter the filename with data from {planet_name}: ")

scale_below = 0
scale_above = 0
datspec = None

if use_data == True: 
    tempdat = pd.read_csv(datafiles[0], delimiter=",", header=0, usecols=["wavelength (nm)", "irradiance (W/m^2)"])
    scale_below = tempdat["wavelength (nm)"][0]
    scale_above = tempdat["wavelength (nm)"][len(tempdat)-1]
    print(f"In situ data used between {scale_below} and {scale_above} nm. Those data will not be scaled by AU. In all other regions, data will be scaled by the specified AU.")

descriptive_tag = ["mean-Ls0-60", 
                   "mean-Ls60-120", 
                   "mean-Ls120-180", 
                   "mean-Ls180-240", 
                   "mean-Ls240-300", 
                   "mean-Ls310-360"] # input("Enter a descriptive tag for this spectrum (please use underscores): ")

for solarfile in solarfiles:
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

    for mdf, AU, desc in zip(datafiles, theAU, descriptive_tag):

        # Open mars data file
        datspec = pd.read_csv(mdf, delimiter=",", header=0, usecols=["wavelength (nm)", "irradiance (W/m^2)"])

        # make a copy
        solarspec_thisiter = deepcopy(solarspec_df)

        # Apply data from in situ
        paste_over_these = solarspec_thisiter.index[solarspec_thisiter["wavelength (nm)"].isin(datspec["wavelength (nm)"])].tolist()
        solarspec_thisiter.loc[paste_over_these, "irradiance (W/m^2/nm)"] = datspec["irradiance (W/m^2)"]

        # Interpolate
        interp_spec = interpolate_solar_spectrum(solarspec_thisiter, AU, planet_name, show_plots=show_plots, 
                                                 desctag=desc)

        # Convert from W/m^2/nm to the units below, and yes, I have checked it multiple times lol including on 11/3/22. 
        # multiply the irradiance by λ/hc to get photons; by 1/10000 to convert to cm^2; and 1/AU^2 to convert to planet orbit.
        # interp_spec["photon flux (phot/s/cm^2/nm)"] = interp_spec[col[1]] * ((interp_spec[col[0]] * 10**(-9))/(h*c)) * 1/(AU**2) * (1/10000)
        interp_spec["photon flux (phot/s/cm^2/nm)"] = np.asarray(list(map(convert_W_m2_to_photons_per_s_cm2, interp_spec[interp_spec.columns[1]], interp_spec[interp_spec.columns[0]])))

        # Now scale to the right AU, avoiding regions where we grabbed data from the planet iteself.
        # Above some value
        if scale_above==0:
            i = 0
        else:
            i = find_nearest(interp_spec[interp_spec.columns[0]], scale_above)
        interp_spec.loc[i:, "photon flux (phot/s/cm^2/nm)"] = interp_spec.loc[i:, "photon flux (phot/s/cm^2/nm)"] * (1**2)/(AU**2)

        # Below some value
        if scale_below==0:
            i = 0
        else:
            i = find_nearest(interp_spec[interp_spec.columns[0]], scale_above)
        interp_spec.loc[:i, "photon flux (phot/s/cm^2/nm)"] = interp_spec.loc[:i, "photon flux (phot/s/cm^2/nm)"] * (1**2)/(AU**2)


        if show_plots:
            plt.figure(figsize=(10,5))
            plt.title("Finalized spectrum, in photon flux, with interpolation and extrapolation")
            plt.plot(interp_spec[interp_spec.columns[0]], interp_spec["photon flux (phot/s/cm^2/nm)"])
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Photons (ph cm^-2 s^-1)")
            plt.savefig(f"final_spectrum_{planet_name}_{descriptive_tag}.png", bbox_inches="tight")
            plt.show()
            print()


        # Write out the result to a file --------------------------------------------------------------------------------------------

        # Build a custom header and writeout the file, with custom header comments before the column names
        sourceinfo = f"# Solar source: {solarfile}"
        if use_data==True:
            sourceinfo += f", data at {planet_name} from {scale_below} to {scale_above}: {mdf}"
        sourceinfo += '\\n'

        header = [f"# AU: {AU}"+'\\n',
                  sourceinfo,
                  f"# File creation date: {datetime.datetime.now().strftime('%B %d, %Y')}"+'\\n\\n']
                  
        outputfile = f"{planet_name}solarphotonflux_{desc}.dat"

        with open(outputfile, "w") as f:  
            for line in header:
                f.write(line)
            f.write("# ")
            interp_spec.to_csv(f, sep='\t', float_format="%.2f", columns=["wavelength (nm)", "photon flux (phot/s/cm^2/nm)"], index=False)

