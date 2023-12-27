# coding=unicode-escape
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as interp
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
   
def interpolate_solar_spectrum(spec, AU, show_plots=True, extrap_tail=False, scale_above=0, interp_start=0.5, interp_end=2399.51, dl=1, desctag=""):
    """
    Input: 
        spec: a Pandas dataframe of SOLAR MEAN spectrum data from 0.1-2390 nm or so, with bins centered
              on every half nm wavelength until 310 nm, and arbitrarily binned thereafter.
              Columns:
                  λ: wavelength in nm
                  irradiance: irradiance in W/m^2/nm
        AU: Orbital radius for which to scale the spectrum. 1.524 is suggested as mean Mars orbital radius.
        show_plots: whether to print out associated plots
        extrap_tail: whether to add an extrapolation from 2390-2400
        scale_above: a wavelength in nanometers above which the result will be scaled by 1/(AU^2). 
                     if using EUV data, you can enter the max wavelength here to not scale that data.
        interp_start: a wavelength in nm above which wavelengths will be interpolated to half wavelengths.
        dl: interpolation spacing in nm.
        desctag: to append to the plot
    Output:
        A Pandas dataframe whole_spectrum from 0.5-2399.5 nm, properly binned every half nm, and a CSV.
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

        whole_spectrum = newsolardata.append(tail, ignore_index=True)
    else:
        whole_spectrum = newsolardata
    
    #make the FINAL dataframe with ALL data for solar max;
    # First convert from W/m^2/nm to the units below, and yes, I have checkd it multiple times lol including on 11/3/22. 
    # multiply the irradiance by λ/hc to get photons; by 1/10000 to convert to cm^2; and 1/AU^2 to convert to Mars orbit.
    # whole_spectrum["photon flux (phot/s/cm^2/nm)"] = whole_spectrum[col[1]] * ((whole_spectrum[col[0]] * 10**(-9))/(h*c)) * 1/(AU**2) * (1/10000)
    
    whole_spectrum["photon flux (phot/s/cm^2/nm)"] = np.asarray(list(map(convert_W_m2_to_photons_per_s_cm2, whole_spectrum[col[1]], whole_spectrum[col[0]])))

    # Now scale to the right AU.
    # You can select whether to only scale above a certain wavelength if you 
    # are using MAVEN EUVM data for the short wavelengths!
    if scale_above==0:
        i = 0
    else:
        i = find_nearest(whole_spectrum[col[0]], scale_above)

    # print(f"Scaling to {AU} AU above {whole_spectrum[i]} nm")
    
    whole_spectrum["photon flux (phot/s/cm^2/nm)"][i:] = whole_spectrum["photon flux (phot/s/cm^2/nm)"][i:] * (1**2)/(AU**2) 

    if show_plots:
        plt.figure(figsize=(10,5))
        plt.title("Finalized spectrum, in photon flux, with interpolation and extrapolation")
        plt.plot(whole_spectrum[col[0]], whole_spectrum["photon flux (phot/s/cm^2/nm)"])
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Photons (ph cm^-2 s^-1)")
        plt.savefig(f"final_spectrum_{desctag}.png", bbox_inches="tight")
        plt.show()
        print()
    
    return whole_spectrum


# Get user input
solarfile = input("Enter the file with solar spectrum data in W/m^2/nm or press enter to use default (composite_SSI_20221122.dat, EUVM up to 189.5 and TSIS above): ")
# solarfile = input() if input()!="" else "composite_SSI_20221122.dat"
print("Enter number of header rows: ")
numheader = int(input())
print("Enter AU at which you'd like the output: ")
theAU = float(input())
print("Enter a descriptive tag for this spectrum (please use underscores): ")
descriptive_tag = input()

longest_euvm_wavelength = 189.51
print(f"Scaling to {theAU} above wavelength={longest_euvm_wavelength} nm. Shortwards of that, EUVM data are used.")
print(f"loading {solarfile}")
solarspec = np.loadtxt(solarfile, skiprows=numheader)

solarspec_df = pd.DataFrame(solarspec, columns=["wavelength (nm)", "irradiance (W/m^2/nm)"])
solarspec_df_tidy = interpolate_solar_spectrum(solarspec_df, theAU, show_plots=True, scale_above=longest_euvm_wavelength, desctag=descriptive_tag)

with open(f"marssolarphotonflux_{descriptive_tag}.dat", "a") as f:
    f.write("# wavelength (nm)\\tphoton flux (phot/s/cm^2/nm)"+'\\n')
    f.write(f'# Original input has been scaled to {theAU} above wavelength={longest_euvm_wavelength} nm\\n')
    solarspec_df_tidy.to_csv(f, sep='\t', float_format="%.2f", index=False)

solarspec_df_tidy.to_csv(f"marssolarphotonflux_{descriptive_tag}.dat", sep='\t', float_format="%.2f", columns=["wavelength (nm)", "photon flux (phot/s/cm^2/nm)"], index=False)
