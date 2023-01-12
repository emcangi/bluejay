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
    
    molec_H2O = (GEL * 999.89) / (1e4 * 18 * 1.67e-27)
    
    if HorH2O == "H"
        return 2 * molec_H2O
    elseif HorH2O == "H2O"
        return molec_H2O
    end
end



function molec_to_GEL(molecules, HorH2O)
    #=
    Converts molecules of H2O per cm^2 to a global equivalent layer in meters. 
    molecules: number of molecules of species HorH2O
    HorH2O: string, either "H" or "H2O" to specify which type the molecule is
    =#
    if HorH2O == "H"
        molecules = molecules / 2
    end
        
    return molecules * 1e4 * 18 * 1.67e-27 / 999.89
end
