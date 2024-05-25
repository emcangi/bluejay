#===============================================================================#
#                        Temperature profile function                           #
#===============================================================================#

# TO DO: Ideally the temperature profile would be generic, but currently that's hard to do because
# we are loading stuff from a file for Venus. 
function T_Mars(Tsurf, Tmeso, Texo; lapserate=-1.4e-5, z_meso_top=108e5, weird_Tn_param=8, globvars...)
    #= 
    Input:
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
    Opt inputs:
        lapserate: adiabatic lapse rate for the lower atmosphere. 1.4e-5 from Zahnle+2008, accounts for dusty atmo.
        z_meso_top: height of the top of the mesosphere (sometimes "tropopause")
        weird_Tn_param: part of function defining neutral temp in upper atmo. See Krasnopolsky 2002, 2006, 2010.
    Output: 
        Arrays of temperatures in K for neutrals, ions, electrons
    =#
    GV = values(globvars)
    required = [:alt]
    check_requirements(keys(GV), required)

    # Altitudes at which various transitions occur -------------------------------
    z_meso_bottom = GV.alt[searchsortednearest(GV.alt, (Tmeso-Tsurf)/(lapserate))]
    
    # These are the altitudes at which we "stitch" together the profiles 
    # from fitting the tanh profile in Ergun+2015,2021 to Hanley+2021 DD8
    # ion and electron profiles, and the somewhat arbitary profiles defined for
    # the region roughly between z_meso_top and the bottom of the fitted profiles.
    z_stitch_electrons = 142e5
    z_stitch_ions = 164e5
    
    # Various indices to define lower atmo, mesosphere, and atmo -----------------
    i_lower = findall(z->z < z_meso_bottom, GV.alt)
    i_meso = findall(z->z_meso_bottom <= z <= z_meso_top, GV.alt)
    i_upper = findall(z->z > z_meso_top, GV.alt)
    i_meso_top = findfirst(z->z==z_meso_top, GV.alt)
    i_stitch_elec = findfirst(z->z==z_stitch_electrons, GV.alt)
    i_stitch_ions = findfirst(z->z==z_stitch_ions, GV.alt)

    function NEUTRALS()
        function upper_atmo_neutrals(z_arr)
            @. return Texo - (Texo - Tmeso)*exp(-((z_arr - z_meso_top)^2)/(weird_Tn_param*1e10*Texo))
        end
        Tn = zeros(size(GV.alt))

        if Tsurf==Tmeso==Texo # isothermal case
            Tn .= Tsurf
        else
            Tn[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
            Tn[i_meso] .= Tmeso
            Tn[i_upper] .= upper_atmo_neutrals(GV.alt[i_upper])
        end

        return Tn 
    end 

    function ELECTRONS() 
        function upper_atmo_electrons(z_arr, TH, TL, z0, H0)
            #=
            Functional form from Ergun+ 2015 and 2021, but fit to data for electrons and ions
            in Hanley+ 2021, DD8 data.
            =#

            @. return ((TH + TL) / 2) + ((TH - TL) / 2) * tanh(((z_arr / 1e5) - z0) / H0)
        end
        Te = zeros(size(GV.alt))

        if Tsurf==Tmeso==Texo # isothermal case
            Te .= Tsurf
        else
            Te[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
            Te[i_meso] .= Tmeso

            # This region connects the upper atmosphere with the isothermal mesosphere
            Te[i_meso_top+1:i_stitch_elec] .= upper_atmo_electrons(GV.alt[i_meso_top+1:i_stitch_elec], -1289.05806755, 469.31681082, 72.24740123, -50.84113252)

            # Don't let profile get lower than specified meso temperature
            Te[findall(t->t < Tmeso, Te)] .= Tmeso

            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            Te[i_stitch_elec:end] .= upper_atmo_electrons(GV.alt[i_stitch_elec:end], 1409.23363494, 292.20319103, 191.39012079, 36.64138724)
        end 

        return Te
    end

    function IONS()
        function upper_atmo_ions(z_arr)
            #=
            fit to Gwen's DD8 profile, SZA 40, Hanley+2021, with values M = 47/13, B = -3480/13
            =#
            # New values for fit to SZA 60,Hanley+2022:
            M = 3.40157034 
            B = -286.48716122
            @. return M*(z_arr/1e5) + B
        end
        
        function meso_ions_byeye(z_arr)
            # This is completely made up! Not fit to any data!
            # It is only designed to make a smooth curve between the upper atmospheric temperatures,
            # which WERE fit to data, and the mesosphere, where we demand the ions thermalize.

            # Values used for Gwen's DD8 profile, SZA 40:  170/49 * (z/1e5) -11990/49
            @. return 136/49 * (z_arr/1e5) -9000/49 # These values are for the new fit to SZA 60, Hanley+2022. (11/2/22)
        end

        Ti = zeros(size(GV.alt))

        if Tsurf==Tmeso==Texo # isothermal case
            Ti .= Tsurf
        else
            Ti[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
            Ti[i_meso] .= Tmeso

            # try as an average of neutrals and electrons. There is no real physical reason for this.
            Ti[i_meso_top:i_stitch_ions] = meso_ions_byeye(GV.alt[i_meso_top:i_stitch_ions])

            # Don't let profile get lower than specified meso temperature
            Ti[findall(t->t < Tmeso, Ti)] .= Tmeso

            # This next region is a fit of the tanh electron temperature expression in Ergun+2015 and 2021 
            # to the electron profile in Hanley+2021, DD8
            Ti[i_stitch_ions:end] .= upper_atmo_ions(GV.alt[i_stitch_ions:end])
        end 
        
        return Ti
    end 

    return Dict("neutrals"=>NEUTRALS(), "ions"=>IONS(), "electrons"=>ELECTRONS())
end 

function T_Venus(Tsurf::Float64, Tmeso::Float64, Texo::Float64, file_for_interp; z_meso_top=80e5, lapserate=-8e-5, weird_Tn_param=8, globvars...)
    #= 
    Input:
        z: altitude above surface in cm
        Tsurf: Surface temperature in KT
        Tmeso: tropopause/mesosphere tempearture
        Texo: exobase temperature
        sptype: "neutral", "ion" or "electron". NECESSARY!
    Output: 
        A single temperature value in K.
    
    Uses the typical temperature structure for neutrals, but interpolates temperatures
    for ions and electrons above 108 km according to the profiles in Fox & Sung 2001.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:alt])
    
    # Subroutines -------------------------------------------------------------------------------------

    function upperatmo_i_or_e(z, particle_type; Ti_array=Ti_interped, Te_array=Te_interped, select_alts=new_a)
        #=
        Finds the index for altitude z within new_a
        =#
        i = something(findfirst(isequal(z), select_alts), 0)
        returnme = Dict("electron"=>Te_array[i], "ion"=>Ti_array[i])
        return returnme[particle_type]
    end

    function NEUTRALS()
        Tn = zeros(size(GV.alt))

        Tn[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Tn[i_meso] .= Tmeso

        # upper atmo
        Tfile = readdlm(file_for_interp, ',')
        T_n = Tfile[2,:]
        alt_i = Tfile[1,:] .* 1e5
        interp_ion = LinearInterpolation(alt_i, T_n)
        Tn_interped = [interp_ion(a) for a in new_a];
        Tn[i_upper] .= Tn_interped # upper_atmo_neutrals(GV.alt[i_upper])

        return Tn 
    end 

    function ELECTRONS(;spc="electron") 
        Te = zeros(size(GV.alt))

        Te[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Te[i_meso] .= Tmeso

        # Upper atmo
        Tfile = readdlm(file_for_interp, ',')
        T_e = Tfile[4,:]
        alt_e = Tfile[1,:] .* 1e5
        interp_elec = LinearInterpolation(alt_e, T_e)
        Te_interped = [interp_elec(a) for a in new_a];
        Te[i_upper] .= Te_interped

        return Te
    end

    function IONS(;spc="ion") 
        Ti = zeros(size(GV.alt))

        Ti[i_lower] .= Tsurf .+ lapserate*GV.alt[i_lower]
        Ti[i_meso] .= Tmeso

        # Upper atmo
        Tfile = readdlm(file_for_interp, ',')
        T_i = Tfile[3,:]
        alt_i = Tfile[1,:] .* 1e5
        interp_ion = LinearInterpolation(alt_i, T_i)
        Ti_interped = [interp_ion(a) for a in new_a];
        Ti[i_upper] .= Ti_interped

        return Ti
    end

    # Define mesosphere scope.
    z_meso_bottom = GV.alt[searchsortednearest(GV.alt, (Tmeso-Tsurf)/(lapserate))]

    # Various indices to define lower atmo, mesosphere, and atmo -----------------
    i_lower = findall(z->z < z_meso_bottom, GV.alt)
    i_meso = findall(z->z_meso_bottom <= z <= z_meso_top, GV.alt)
    i_upper = findall(z->z > z_meso_top, GV.alt)
    # i_meso_top = findfirst(z->z==z_meso_top, GV.alt)

    # For interpolating upper atmo temperatures from Fox & Sung 2001
    new_a = collect(90e5:2e5:250e5) # TODO: Remove hard coded values

    return Dict("neutrals"=>NEUTRALS(), "ions"=>IONS(), "electrons"=>ELECTRONS())
end