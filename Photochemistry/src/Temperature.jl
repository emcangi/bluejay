# **************************************************************************** #
#                                                                              #
#                          Temperature functions                               #
#                                                                              #
# **************************************************************************** #

function T(z, Tsurf, Tmeso, Texo, sptype::String; lapserate=-1.4e-5)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    =#
    
    z_meso_bottom = alt[searchsortednearest(alt, (Tmeso-Tsurf)/(lapserate))]
    z_meso_top = 108e5  # height of the tropopause top
    
    # These are the altitudes at which we "stitch" together the profiles 
    # from fitting the tanh profile in Ergun+2015,2021 to Hanley+2021 DD8
    # ion and electron profiles, and the somewhat arbitary profiles defined for
    # the region roughly between z_meso_top and the bottom of the fitted profiles.
    stitch_alt_electrons = 142e5
    stitch_alt_ions = 165e5

    function T_upper_atmo_neutrals(zee)
        return Texo - (Texo - Tmeso)*exp(-((zee - z_meso_top)^2)/(8e10*Texo))
    end
    
    function bobs_profile(z, TH, TL, z0, H0)
        #=
        Functional form from Ergun+ 2015 and 2021, but fit to data for electrons and ions
        in Hanley+ 2021, DD8 data.
        =#
        return ((TH+TL)/2) + ((TH-TL)/2) * tanh(((z/1e5)-z0)/H0)
    end
    
    function T_upper_atmo_ions(z)
        #=
        This is a totally arbitary functional form for the region from [z_meso_top, 138]
        =#
        # Values used for Gwen's DD* profile, SZA 40, Hanley+2021:
        #M = 47/13
        #B = -3480/13
        # New values for fit to SZA 60,Hanley+2022:
        M = 3.39157034
        B = -286.48716122
        return M*(z/1e5) + B
    end
    
    function T_meso_ions_byeye(z)
        # This is completely made up! Not fit to any data!
        # It is only designed to make a smooth curve between the upper atmospheric temperatures,
        # which WERE fit to data, and the mesosphere, where we demand the ions thermalize.

        # Values used for Gwen's DD8 profile, SZA 40:  170/49 * (z/1e5) -11990/49
        return 136/49 * (z/1e5) -9000/49 # These values are for the new fit to SZA 60, Hanley+2022. (11/2/22)
    end
    
    # In the lower atmosphere, neutrals, ions, and electrons all 
    # have the same temperatures. 
    if z < z_meso_bottom
        return Tsurf + lapserate*z
    elseif z_meso_bottom <= z <= z_meso_top 
        return Tmeso
    
    # Near the top of the isothermal mesosphere, profiles diverge.        
    elseif z > z_meso_top
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="electron"
            # This region connects the upper atmosphere with the isothermal mesosphere
            if z_meso_top <= z < stitch_alt_electrons
                return bobs_profile(z, -1289.05806755, 469.31681082, 72.24740123, -50.84113252)
            
            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            elseif z >= stitch_alt_electrons
                return bobs_profile(z, 1409.23363494, 292.20319103, 191.39012079, 36.64138724)
            end
        elseif sptype=="ion"
            # This is similar to the electron handling, but for the ion profile in Hanley+2021, DD8.
            if z_meso_top < z <= stitch_alt_ions
                return T_meso_ions_byeye(z) < Tmeso ? Tmeso : T_meso_ions_byeye(z)
            elseif z > stitch_alt_ions
                return T_upper_atmo_ions(z)
            end
        end
    end
end

function T_updated_old(z, Tsurf, Tmeso, Texo, sptype::String)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    z_meso_bottom = alt[searchsortednearest(alt, (Tmeso-Tsurf)/(lapserate))]
    z_meso_top = 110e5  # height of the tropopause top
    
    # These are the altitudes at which we "stitch" together the profiles 
    # from fitting the tanh profile in Ergun+2015,2021 to Hanley+2021 DD8
    # ion and electron profiles, and the somewhat arbitary profiles defined for
    # the region roughly between z_meso_top and the bottom of the fitted profiles.
    stitch_alt_electrons = 142e5
    stitch_alt_ions = 135e5

    function T_upper_atmo_neutrals(zee)
        return Texo - (Texo - Tmeso)*exp(-((zee - z_meso_top)^2)/(8e10*Texo))
    end
    
    function bobs_profile(z, TH, TL, z0, H0)
        #=
        Functional form from Ergun+ 2015 and 2021, but fit to data for electrons and ions
        in Hanley+ 2021, DD8 data.
        =#
        return ((TH+TL)/2) + ((TH-TL)/2) * tanh(((z/1e5)-z0)/H0)
    end
    
    function T_ions_thermalize_region(z)
        #=
        This is a totally arbitary functional form for the region from [z_meso_top, 138]
        =#
        M = 6.06308414
        B = -538.97784139
        return M*(z/1e5) + B
    end
    
    # In the lower atmosphere, neutrals, ions, and electrons all 
    # have the same temperatures. 
    if z < z_meso_bottom
        return Tsurf + lapserate*z
    elseif z_meso_bottom <= z <= z_meso_top 
        return Tmeso
    
    # Near the top of the isothermal mesosphere, profiles diverge.        
    elseif z > z_meso_top
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="electron"
            # This region connects the upper atmosphere with the isothermal mesosphere
            if z_meso_top <= z < stitch_alt_electrons
                return bobs_profile(z, -1289.05806755, 469.31681082, 72.24740123, -50.84113252)
            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            elseif z >= stitch_alt_electrons
                return bobs_profile(z, 1409.23363494, 292.20319103, 191.39012079, 36.64138724)
            end
        elseif sptype=="ion"
            # This is similar to the electron handling, but for the ion profile in Hanley+2021, DD8.
            if z_meso_top < z <= stitch_alt_ions
                return T_ions_thermalize_region(z) < Tmeso ? Tmeso : T_ions_thermalize_region(z)
            elseif z > stitch_alt_ions
                return bobs_profile(z, 4.87796600e+06, 2.15643719e+02, 7.83610155e+02, 1.16129872e+02)
            end
        end
    end
end

function T_all(z, Tsurf, Tmeso, Texo, sptype::String)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in K
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of mesosphere. This code allows it to vary for surface or 
    # mesosphere experiments.
    ztropo_bot = (Tmeso-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    function find_offset(z1, T1, A, B, C)
        #=
        (z1, T1): the new point that satisfies T1 = T_upper_atmo_neutrals(z1), which we
                  also want to solve T_upper_atmo_ions.
        =#
        Q = A*(exp(-B*exp(-C*(z1/1e5))))
        
        D = T1 - Q
        return D
    end

    function T_upper_atmo_neutrals(zee)
        return Texo - (Texo - Tmeso)*exp(-((zee - ztropo)^2)/(8e10*Texo))
    end

    function T_upper_atmo_electrons(z1, T1)
        # for electrons, this is valid above 130 km in height.
        A = 2.24100983e+03
        B = 1.84024165e+01
        C = 1.44238590e-02
        D = find_offset(z1, T1, A, B, C)  
        return A*(exp(-B*exp(-C*z/1e5))) + D
    end

    function T_upper_atmo_ions(z1, T1)
        # valid above 160 km 
        A = 9.96741381e2
        B = 1.66317054e2
        C = 2.49434339e-2
        D = find_offset(z1, T1, A, B, C)  #
        return A*(exp(-B*exp(-C*(z/1e5 + 0)))) + D
    end

    # Used for correctly configuring the ion and electron profiles
    z1_ion = 160e5
    T1_ion = T_upper_atmo_neutrals(z1_ion)
    z1_e = 130e5
    T1_e = T_upper_atmo_neutrals(z1_e)

    # this is the WOOO--OOORRRST but someday I'll make it better MAYBE???
    if z >= 160e5 
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="ion"
            return T_upper_atmo_ions(z1_ion, T1_ion) # only valid up here
        elseif sptype=="electron"
            T_upper_atmo_electrons(z1_e, T1_e)
        end
    elseif 130e5 < z < 160e5
        if sptype=="neutral"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="ion"
            return T_upper_atmo_neutrals(z)
        elseif sptype=="electron"
            T_upper_atmo_electrons(z1_e, T1_e)
        end
    elseif ztropo < z <= 130e5   # upper atmosphere until neutrals, ions, electrons diverge
        return Texo - (Texo - Tmeso)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif (ztropo - ztropowidth) < z <= ztropo  # tropopause
        return Tmeso
    elseif z <= ztropo-ztropowidth  # lower atmosphere; <= makes it work for isothermal atm
        return Tsurf + lapserate*z
    end
end

function Tpiecewise(z, Tsurf, Tmeso, Texo)
    #= 
    FOR USE WITH FRACTIONATION FACTOR PROJECT ONLY.

    DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 "half-Gaussian" function for temperatures 
    altitudes above the tropopause, with a constant lapse rate (1.4K/km) 
    in the lower atmosphere. The tropopause width is allowed to vary
    in certain cases.

    z: altitude above surface in cm
    Tsurf: Surface temperature in K
    Tmeso: tropopause/mesosphere tempearture
    Texo: exobase temperature
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of mesosphere. This code allows it to vary for surface or 
    # mesosphere experiments.
    ztropo_bot = (Tmeso-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Tmeso)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Tmeso
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end