# **************************************************************************** #
#                                                                              #
#                            Unit conversions                                  #
#                                                                              #
# **************************************************************************** #


function DH_conversion(; deltaD=0, VSMOW_multiplier=1, DH=1.56e-4)
    #=
    Convert between a D/H ratio, a multiplier of VSMOW (e.g. Mars atmosphere = 5.5 x VSMOW) 
    and delta D, which is quoted in milles and is either negative or positive.

    =#

    if (deltaD!=0) & ((VSMOW_multiplier!=1) | (DH != 1.56e-4))
        throw("InputError: Enter either deltaD, or (VSMOW_multiplier or DH), but not both.")
    end


    if (deltaD==0) & (VSMOW_multiplier==1.56)
        println(L"No conversion requested. VSMOW = 1.56e-4, which is $\delta$D = 0.") 
    elseif (deltaD != 0)
        # this is for when delta D is quoted as some number > 1 away from 0, not a fraction.
        multiplier = deltaD/1000 + 1
        return "D/H = $(round(multiplier*1.56, digits=2))e-4 which is $(multiplier)x VSMOW"
    elseif (VSMOW_multiplier!=1)
        return 1000*(VSMOW_multiplier - 1) # returns deltaD
    elseif (DH != 1.56e-4)
        println("Input: D/H ratio = $DH")
        return 1000*((DH/1.56e-4) - 1) # returns deltaD also
    end
end

function GEL_to_molecule(GEL, HorH2O)
    #=
    Converts a global equivalent layer in meters to molecules per cm^2. 
    GEL: layer of water in meters
    HorH2O: string, either "H" or "H2O" to specify which molecule to return
    =#
    
    molec_H2O_per_cm2 = (GEL * 100 * density_water_cgs) / (18 * mH)
    
    if HorH2O == "H"
        return 2 * molec_H2O_per_cm2
    elseif HorH2O == "H2O"
        return molec_H2O_per_cm2
    end
end

function molec_to_GEL(molec_per_cm2, HorH2O)
    #=
    Converts molecules per cm^2 to a global equivalent layer in meters. 
    molec_per_cm2: number of molec_per_cm2 of species HorH2O
    HorH2O: string, either "H" or "H2O" to specify which type the molecule is
    =#
    if HorH2O == "H"
        molec_per_cm2 = molec_per_cm2 / 2
    end
        
    return (molec_per_cm2 * 18 * mH) * (1 / density_water_cgs) * (1#=m=#/100 #=cm=#)
end

function total_escape_to_GEL(esc_per_s, t, HorH2O)
    #=
    esc_per_s: total escape of either H or H2O per second, for the whole planet
    t: years over which the escape occurs
    HorH2O: "H" or "H2O" to denote the escaping molecules
    =#

    esc_per_cm2 = esc_per_s * (t * s_per_yr) * (1/SA_Mars) 
    return molec_to_GEL(esc_per_cm2, HorH2O)
end

function total_escape_to_area_escape(esc_per_s)
    #=
    esc_per_s: total escape of some molecule over the whole planet 
    =#

    return esc_per_s / SA_Mars
end

function by_to_sec(y)
    return y*1e9*3.154e7
end

function prum_to_ppm(sp, prum, atmdict; globvars...)
    #=
    Converts pr μm to ppm, only works for water!
    Inputs
        sp: H2O or HDO
        prum: water pr μm
        atmdict: the atmosphere against which you will be comparing
    Output:
        water total in ppm. 
    =#


    GV = values(globvars)
    required =  [:molmass, :all_species]
    check_requirements(keys(GV), required)

    colabund = colabund_from_prum(sp, prum; globvars...)
    nt_sum = sum(n_tot(atmdict, ihoriz; GV.all_species) for ihoriz in 1:n_horiz) * dz
    
    return (colabund/nt_sum) / 1e-6 # return ppm
end

function rayleigh_fractionation(solve_for; water_lost=nothing, water_today=nothing, f=nothing, DH_original=nothing, DH_today=nothing)
    #=
    Solve the Rayleigh fractionation equation. Units don't matter as everything occurs in fractions. 
    
    Inputs
        solve_for: Any of the optional variable names
        water_lost: Water lost over the time period bookended by the two DH values.
        water_today: Current water inventory, or water inventory at the time designated by DH_today.
        f: fractionation factor, assumed constant
        DH_original: D/H value at some time in the past
        DH_today: D/H value "today" (or at some time you're interested in calculating)
    Outputs:
        Whichever variable wasn't provided.
    =#
    if solve_for=="water_lost"
        @assert water_lost == nothing
        @assert water_today != nothing
        @assert f != nothing
        @assert DH_original != nothing
        @assert DH_today != nothing
        return water_today * ((DH_today / DH_original)^(1/(1-f)) -1 )
    end

    if solve_for=="water_today"
        @assert water_lost != nothing
        @assert water_today == nothing
        @assert f != nothing
        @assert DH_original != nothing
        @assert DH_today != nothing
        return water_lost / ((DH_today / DH_original)^(1/(1-f)) -1 ) 
    end

    if solve_for=="f"
        @assert water_lost != nothing
        @assert water_today != nothing
        @assert f == nothing
        @assert DH_original != nothing
        @assert DH_today != nothing
        return 1 - ( (log(DH_today) - log(DH_original)) / (log(water_lost/water_today + 1)) )
    end

    if solve_for=="DH_today"
        @assert water_lost != nothing
        @assert water_today != nothing
        @assert f != nothing
        @assert DH_original != nothing
        @assert DH_today == nothing
        return DH_original * (water_lost / water_today + 1) ^ (1-f)
    end

    if solve_for=="DH_original"
        @assert water_lost != nothing
        @assert water_today != nothing
        @assert f != nothing
        @assert DH_original == nothing
        @assert DH_today != nothing
        return DH_today / ((water_lost / water_today + 1) ^ (1-f))
    end

end


#                       Converting joules and eletron volts                     #
#===============================================================================#

kJ_to_eV(kj) = kj * 6.242e21

ev_per_molecule(eV_per_mole) = eV_per_mole/6.022e23

