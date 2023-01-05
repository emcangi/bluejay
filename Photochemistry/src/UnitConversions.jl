# **************************************************************************** #
#                                                                              #
#                            Unit conversions                                  #
#                                                                              #
# **************************************************************************** #


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
