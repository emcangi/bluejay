# **************************************************************************** #
#                                                                              #
#           Functions to format and load the chemical reaction network         #
#                                                                              #
# **************************************************************************** #

function filter_network(sp::Symbol, rate_type::String, species_role::String; globvars...)
    #= 
    Given a reaction network, filter it down by:
    sp: A specific species, produced by
    rate_type reactions (Jrates, krates, or all), where the species is a
    species_role (reactant, product, or any).

    Returns the reaction network but without reactions that don't meet the criteria.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:reaction_network])

    # Select either photodissociation or bi-/tri-molecular reactions
    if rate_type=="Jrates"
        selected_rxns = filter(x->occursin("J", string(x[3])), GV.reaction_network)
    elseif rate_type=="krates" 
        selected_rxns = filter(x->!occursin("J", string(x[3])), GV.reaction_network)
    elseif rate_type=="all"
        selected_rxns = deepcopy(GV.reaction_network)
    end

    # now select only the reactions with the species in question (there is probably a less repetative way to write this)
    filtered_rxn_list = Any[]

    # need to make a regular expression containing the species that will ignore cases where
    # the species string appears as a subset of another species string, i.e. when O2pl is 
    # detected as being in the reaction CO2 -> CO2pl.
    species_re = r"\b"*string(sp)*r"\b"

    if species_role=="reactant"
        found_rxns = filter(x->(occursin(species_re, string(x[1]))), selected_rxns)
        filtered_rxn_list = vcat(filtered_rxn_list, found_rxns)
    elseif species_role=="product"
        found_rxns = filter(x->(occursin(species_re, string(x[2]))), selected_rxns)
        filtered_rxn_list = vcat(filtered_rxn_list, found_rxns)
    elseif species_role=="both"
        found_rxns = filter(x->(occursin(species_re, string(x[1])) || occursin(species_re, string(x[2]))), selected_rxns)
        filtered_rxn_list = vcat(filtered_rxn_list, found_rxns)
    else
        throw("No handling found for species_role=$(species_role)")
    end

    return unique(filtered_rxn_list)  # gets rid of duplicates since occursin() is greedy
end

function format_Jrates(spreadsheet, used_species, return_what; saveloc=nothing, write_rxns=false, hot_atoms=false, ions_on=true)
    #=
    This formats Jrate symbols from the reaction entrires in the spreadsheet.
    Input
        spreadsheet: .xlsx file with chemical and photochemical reaction entries
                    return_what: "Jratelist" to return just a list of Jrate symbols, or
                     "Jrate network" to return the network of photodissociation/
                     photoionization reactions. 
    Output: as described in return_what
    =#
    # Lists to store things
    Jrates = [] # The symbols for Jrates which are used like "species" in the model
    Jrate_network = [] # The reaction network for unimolecular reactions
    newJrates = [] # Keeps track of whether the rates are newly introduced, which is used to initialize them as 0.
                   # MUST be set as new in the spreadsheet.

    # Do photodissociation handling all the time, for neutral-only and all-inclusive simulations ==========
    # Handle the neutral photodissociation rates
    photodissociation = DataFrame(XLSX.readtable(spreadsheet, "Photodissociation"))
    if !ions_on
        photodissociation = filter(row -> row.Neutrals=="Yes", photodissociation)
    end 
    replace!(photodissociation."P2", missing=>"none");
    replace!(photodissociation."P3", missing=>"none");
    replace!(photodissociation."hotH", missing=>"No");
    replace!(photodissociation."hotD", missing=>"No");
    replace!(photodissociation."hotH2", missing=>"No");
    replace!(photodissociation."hotHD", missing=>"No");

    # Filter out reactions for species we aren't using
    photodissociation = filter(row->all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), photodissociation)
    unused_photodissociation = filter(row->!all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), photodissociation)

    # Set up arrays to hold the hot H and hot D producing reactions
    hot_H_J = []
    hot_D_J = []
    hot_H2_J = []
    hot_HD_J = []
    
    for r in eachrow(photodissociation)
        products = [Symbol(i) for i in [r.P1, r.P2, r.P3] if i!="none"]
        Jrate = get_Jrate_symb(r.R1, products)
        rxn = [[Symbol(r.R1)], products, Jrate]
        push!(Jrate_network, rxn)

        # Add either to the "new" list or the normal list
        if r.Status=="New"
            push!(newJrates, Jrate)
        elseif r.Status=="Conv"
            push!(Jrates, Jrate)
        end 

        # Keep lists of photodissociation that produces hot H, D
        if r.hotH=="Yes"
            push!(hot_H_J, rxn)
        end 
        if r.hotD=="Yes"
            push!(hot_D_J, rxn)
        end 
        if r.hotH2=="Yes"
            push!(hot_H2_J, rxn)
        end 
        if r.hotHD=="Yes"
            push!(hot_HD_J, rxn)
        end 
    end

    # Do photoionization if the simulation uses ions 
    if ions_on
        photoionization = DataFrame(XLSX.readtable(spreadsheet, "Photoionization"))
        replace!(photoionization."P2", missing=>"none");
        replace!(photoionization."P3", missing=>"none");
        replace!(photoionization."hotH", missing=>"No");
        replace!(photoionization."hotD", missing=>"No");
        replace!(photoionization."hotH2", missing=>"No");
        replace!(photoionization."hotHD", missing=>"No");

        # Filter out reactions for species we aren't using
        photoionization = filter(row->all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), photoionization)
        unused_photoionization = filter(row->!all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), photoionization)

        for r in eachrow(photoionization)
            products = [Symbol(i) for i in [r.P1, r.P2, r.P3] if i!="none"]
            Jrate = get_Jrate_symb(r.R1, products)
            rxn = [[Symbol(r.R1)], products, Jrate]
            push!(Jrate_network, rxn)
            if r.Status=="New"
                push!(newJrates, Jrate)
            elseif r.Status=="Conv"
                push!(Jrates, Jrate)
            end 
            if r.hotH=="Yes"
                push!(hot_H_J, rxn)
            end 
            if r.hotD=="Yes"
                push!(hot_D_J, rxn)
            end 
            if r.hotH2=="Yes"
                push!(hot_H2_J, rxn)
            end 
            if r.hotHD=="Yes"
                push!(hot_HD_J, rxn)
            end 
        end
    end 

    if return_what=="Jratelist"
        return Jrates, newJrates
    elseif return_what=="Jrate network"
        if write_rxns
            try 
                @assert saveloc != nothing 
            catch AssertionError
                throw("Please pass in the location to save the spreadsheet of active reactions to parameter 'saveloc'")
            end
            # Write the in-use reactions to a new file for the simulation 
            log_reactions(photodissociation, "Photodissociation", saveloc)
            log_reactions(unused_photodissociation, "Unused photodissociation", saveloc)

            if ions_on
                # Write the in-use reactions to a new file for the simulation 
                log_reactions(photoionization, "Photoionization", saveloc)
                log_reactions(unused_photoionization, "Unused photoionization", saveloc)
            end
        end
        
        returndict = Dict("all"=>Jrate_network)
    
        if hot_atoms
            returndict["hotH"] = hot_H_J
            returndict["hotD"] = hot_D_J
            returndict["hotH2"] = hot_H2_J
            returndict["hotHD"] = hot_HD_J
        end

        return returndict
    end
end

function format_neutral_network(reactions_spreadsheet, used_species; saveloc=nothing, write_rxns=false, verbose=false)
    #=
    Inputs:
        reactions_spreadsheet: Properly formatted .xlsx of chemical reactions with 
                               sheet names: Neutral reactions, Ion reactions.
        used_species: List of species in use, really just used once to find unneeded
                      reactions
    Outputs:
        lists of several reaction types, all neutrals, in form [[R1, R2], [P1, P2], :(k)].
        3 reactants or products are also possible.
    =#
    n_table_raw = DataFrame(XLSX.readtable(reactions_spreadsheet, "Neutral reactions"))
    
    # Replace missing values with proper stuff 
    replace!(n_table_raw."R2", missing=>"none");
    replace!(n_table_raw."R3", missing=>"none");
    replace!(n_table_raw."P2", missing=>"none");
    replace!(n_table_raw."P3", missing=>"none");
    replace!(n_table_raw."M2", missing=>1);
    replace!(n_table_raw."M1", missing=>1);
    replace!(n_table_raw."pow", missing=>0);
    replace!(n_table_raw."BR", missing=>1);

    # Get rid of reactions with a species we don't use
    n_table = filter(row -> all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2), Symbol(row.R3), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), n_table_raw)
    ununsed_reactions = filter(row->!all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2), Symbol(row.R3), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), n_table_raw) 
    if verbose
        println("Removed reactions: $(ununsed_reactions)")
    end

    # Write the in-use reactions to a new file for the simulation 
    if write_rxns
        try 
            @assert saveloc != nothing 
            log_reactions(n_table, "Neutral reactions", saveloc)
            log_reactions(ununsed_reactions, "Unused neutral reactions", saveloc)
        catch AssertionError
            throw("Please pass in the location to save the spreadsheet of active reactions to parameter 'saveloc'. If you did this, it may be that the file already exists.")
        end
    end
    
    # Set up the array to store neutral network 
    neutral_network = vec(Array{Any}(undef, nrow(n_table)))

    # Just some quick error handling
    @assert all(x->x>=0, n_table."kA")
    @assert all(x->x>=0, n_table."k0A")
    
    # Get all the counts of each type 
    counts = Dict([i=>count(x->x==i, n_table[:, "type"]) for i in 1:6])

    n = 1 # Keeps track of iteration, since some counts may be 0.

    for i in 1:1+counts[1]-1 # Type 1: Pressure independent unimolecular rxns, high pressure limit
        @assert n_table[i, "type"]==1
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]

        neutral_network[i] = [reactants, products, :($(make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn", 
                                                                   n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])))]
        n += 1
    end

    for i in n:n+counts[2]-1 # Type 2: P independent bimolecular; use k (same as *k inf), units cm^3 s^-1
        @assert n_table[i, "type"]==2
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        neutral_network[i] = [reactants, products, make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn",
                                                               n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])]
        n += 1
    end 

    for i in n:n+counts[3]-1 # Type 3: P dependent bimolecular.
        @assert n_table[i, "type"]==3
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        neutral_network[i] = [reactants, products, make_Troe([n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"]], 
                                                             [n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"]],
                                                              n_table[i, "F"], n_table[i, "M2"], n_table[i, "M1"], 
                                                              n_table[i, "pow"], n_table[i, "BR"])]
        n =+ 1
    end 

    for i in n:n+counts[4]-1 # Type 4: P dependent association rxns
        @assert n_table[i, "type"]==4
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]        
        neutral_network[i] = [reactants, products, make_modified_Troe([n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"]], 
                                                                      [n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"]],
                                                                      [n_table[i, "kradA"], n_table[i, "kradB"], n_table[i, "kradC"]],
                                                                       n_table[i, "F"], n_table[i, "M2"], n_table[i, "M1"], 
                                                                       n_table[i, "pow"], n_table[i, "BR"])]
        n += 1
    end 

    for i in n:n+counts[5]-1 # Three body reactions 
        @assert n_table[i, "type"]==5
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        
        k0 = make_k_expr(n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"], "Tn", 
                         n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])             
        kinf = make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn", 
                           n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])

        neutral_network[i] = [reactants, products, :($k0 .* M ./ (1 .+ $k0 .* M ./ $kinf).*0.6 .^ ((1 .+ (log10.($k0 .* M ./ $kinf)) .^2).^-1.0))]
        n += 1
    end 

    for i in n:n+counts[6]-1 # Three body reactions 
        @assert n_table[i, "type"]==6
        reactants = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "R1"], n_table[i, "R2"], n_table[i, "R3"]])]
        products = [Symbol(j) for j in filter!(j->j!="none", [n_table[i, "P1"], n_table[i, "P2"], n_table[i, "P3"]])]
        k0 = make_k_expr(n_table[i, "k0A"], n_table[i, "k0B"], n_table[i, "k0C"], "Tn", 
                         n_table[i, "M2"], n_table[i, "M2"], n_table[i, "pow"], n_table[i, "BR"])             
        kinf = make_k_expr(n_table[i, "kA"], n_table[i, "kB"], n_table[i, "kC"], "Tn", 
                           n_table[i, "M2"], n_table[i, "M1"], n_table[i, "pow"], n_table[i, "BR"])

        neutral_network[i] = [reactants, products, :($k0 ./ (1 .+ $k0 ./ ($kinf ./ M)).*0.6 .^ ((1 .+ (log10.($k0 ./ ($kinf .* M))) .^2).^-1.0))]
        n += 1
    end

    return neutral_network
end

function format_ion_network(reactions_spreadsheet, used_species; saveloc=nothing, write_rxns=false, hot_atoms=false, verbose=false, special_condition=nothing)
    #=
    Inputs:
        reactions_spreadsheet: Properly formatted .xlsx of chemical reactions with 
                               sheet names: Neutral reactions, Ion reactions.
        used_species: active species in the model
    Outputs:
        returndict, which has three keys: "all", "hotH", and "hotD", with each pointing to a
                    list of all the relevant ion reactions in format [[R1, R2], [P1, P2], :(k)].
                    3 reactants or products are also possible.
    =#
    
    ion_table_raw = DataFrame(XLSX.readtable(reactions_spreadsheet, "Ion reactions"))
    
    # Replace missing values with proper stuff 
    replace!(ion_table_raw."R2", missing=>"none");
    replace!(ion_table_raw."P2", missing=>"none");
    replace!(ion_table_raw."P3", missing=>"none");
    replace!(ion_table_raw."M2", missing=>1);
    replace!(ion_table_raw."M1", missing=>1);
    replace!(ion_table_raw."pow", missing=>0);
    replace!(ion_table_raw."BR", missing=>1);
    replace!(ion_table_raw."hotH", missing=>"No");
    replace!(ion_table_raw."hotD", missing=>"No");
    replace!(ion_table_raw."hotH2", missing=>"No");
    replace!(ion_table_raw."hotHD", missing=>"No");

    # Filter out reactions with reactants we don't use
    ion_table = filter(row->all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), ion_table_raw)
    unused_reactions = filter(row->!all(x->x in union(used_species, [:E, :M, :none]), [Symbol(row.R1), Symbol(row.R2), Symbol(row.P1), Symbol(row.P2), Symbol(row.P3)]), ion_table_raw)
    if special_condition != nothing 
        ion_table = filter(special_condition, ion_table)
    end

    if verbose
        println("Removed reactions: $(unused_reactions)")
    end 

    # Write the in-use reactions to a new file for the simulation 
    if write_rxns
        try 
            @assert saveloc != nothing 
            log_reactions(ion_table, "Ion reactions", saveloc)
            log_reactions(unused_reactions, "Unused ion reactions", saveloc)
        catch AssertionError
            throw("Please pass in the location to save the spreadsheet of active reactions to parameter 'saveloc'")
        end
    end
    
    # Set up arrays to hold the hot H and hot D producing reactions
    if hot_atoms
        num_hotH = size(filter(row->(row.hotH in ["Yes"]), ion_table))[1]
        hot_H_network = Array{Any}(undef, num_hotH, 1)
        Hi = 1

        num_hotD = size(filter(row->(row.hotD in ["Yes"]), ion_table))[1]
        hot_D_network = Array{Any}(undef, num_hotD, 1)
        Di = 1

        num_hotH2 = size(filter(row->(row.hotH2 in ["Yes"]), ion_table))[1]
        hot_H2_network = Array{Any}(undef, num_hotH2, 1)
        H2i = 1

        num_hotHD = size(filter(row->(row.hotHD in ["Yes"]), ion_table))[1]
        hot_HD_network = Array{Any}(undef, num_hotHD, 1)
        HDi = 1
    end

    numrows = size(ion_table)[1]
    ion_network = Array{Any}(undef, numrows, 1)
    
    # Just some quick error handling
    @assert all(x->x>=0, ion_table."kA")
    @assert all(x->x>=0, ion_table."k0A")
    @assert all(x->x>=0, ion_table."kradA")

    # Fill in the network
    for i in range(1, stop=numrows)

        # Set the temperature term
        Tstr = ion_table[i, "type"] == -4 ? "Te" : "Ti"
            
        kA = ion_table[i, "kA"]
        kB = ion_table[i, "kB"]
        kC = ion_table[i, "kC"]
        M2 = ion_table[i, "M2"]
        M1 = ion_table[i, "M1"]
        pow = ion_table[i, "pow"]
        BR = ion_table[i, "BR"]
       
        # Fill out the entry
        ion_network[i] = [[Symbol(ion_table[i, "R1"]), Symbol(ion_table[i, "R2"])],  # REACTANT LIST
                          [Symbol(j) for j in filter!(j->j!="none", [ion_table[i, "P1"], ion_table[i, "P2"], ion_table[i, "P3"]])], # PRODUCT LIST
                          make_k_expr(kA, kB, kC, Tstr, M2, M1, pow, BR)] # RATE COEFFICIENT EXPRESSION
        
        if hot_atoms
            if ion_table[i, "hotH"] in ["Yes"]
                hot_H_network[Hi] = ion_network[i]
                Hi += 1
            end
            if ion_table[i, "hotD"] in ["Yes"]
                hot_D_network[Di] = ion_network[i]
                Di += 1
            end
            if ion_table[i, "hotH2"] in ["Yes"]
                hot_H2_network[H2i] = ion_network[i]
                H2i += 1
            end
            if ion_table[i, "hotHD"] in ["Yes"]
                hot_HD_network[HDi] = ion_network[i]
                HDi += 1
            end
        end
        
    end
    
    returndict = Dict("all"=>ion_network)
    
    if hot_atoms
        returndict["hotH"] = hot_H_network
        returndict["hotD"] = hot_D_network
        returndict["hotH2"] = hot_H2_network
        returndict["hotHD"] = hot_HD_network
    end

    return returndict
end

function load_reaction_network(spreadsheet; saveloc=nothing, write_rxns=false, to_return="all", ions_on=true, get_hot_rxns=false, globvars...)
    #=
    Inputs:
        spreadsheet: path and filename of the reaction network spreadsheet
        saveloc: Path in which to store the re-written spreadsheet of active reactions. This should be the simulation result folder. 
    Output:
        reaction_network, a vector of vectors in the format
        [Symbol[reactants...], Symbol[products...], :(rate coefficient expression)
    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species])

    if get_hot_rxns
        try 
            @assert ions_on==true 
        catch
            throw("You can't calculate non-thermal escape in a neutral only simulation, at least for now, try again with ions_on=true")
        end
    end

    if to_return=="all"
        #load neutrals
        neutral_net = format_neutral_network(spreadsheet, GV.all_species; saveloc=saveloc, write_rxns=write_rxns)

        # Get Jrates
        Jdict = format_Jrates(spreadsheet, GV.all_species, "Jrate network"; saveloc=saveloc, write_rxns=write_rxns, hot_atoms=true, ions_on=ions_on)

        # Load all ions
        if ions_on == true
            ionsdict = format_ion_network(spreadsheet, GV.all_species; saveloc=saveloc, write_rxns=write_rxns, hot_atoms=get_hot_rxns)

            @assert haskey(ionsdict, "hotH")
            @assert haskey(Jdict, "hotH")
            @assert haskey(ionsdict, "hotD")
            @assert haskey(Jdict, "hotD")
            @assert haskey(ionsdict, "hotH2")
            @assert haskey(Jdict, "hotH2")
            @assert haskey(ionsdict, "hotHD")
            @assert haskey(Jdict, "hotHD")

            hot_H_network = [Jdict["hotH"]..., ionsdict["hotH"]...]
            hot_D_network = [Jdict["hotD"]..., ionsdict["hotD"]...]
            hot_H2_network = [Jdict["hotH2"]..., ionsdict["hotH2"]...]
            hot_HD_network = [Jdict["hotHD"]..., ionsdict["hotHD"]...]
            whole_network = [Jdict["all"]..., neutral_net..., ionsdict["all"]...]
        else 
            # For neutral-only simulations we stlil need the two CO2pl attack on molecular hydrogen reactions
            ionsdict = format_ion_network(spreadsheet, GV.all_species; saveloc=saveloc, write_rxns=write_rxns, hot_atoms=get_hot_rxns,
                       special_condition=row->(row.R1=="CO2pl" && (row.R2=="HD" || row.R2=="H2") && (row.P1=="CO2"))) 
            hot_H_network = []
            hot_D_network = []
            hot_H2_network = []
            hot_HD_network = []
            whole_network = [Jdict["all"]..., neutral_net...]
        end 

        return whole_network, hot_H_network, hot_D_network, hot_H2_network, hot_HD_network
    elseif to_return=="hot"
        # Load all ions
        ionsdict = format_ion_network(spreadsheet, GV.all_species; saveloc=saveloc, write_rxns=write_rxns, hot_atoms=get_hot_rxns)
        Jdict = format_Jrates(spreadsheet, GV.all_species, "Jrate network"; saveloc=saveloc, write_rxns=write_rxns, hot_atoms=true, ions_on=ions_on)

        # Construct the networks
        @assert haskey(ionsdict, "hotH")
        @assert haskey(Jdict, "hotH")
        @assert haskey(ionsdict, "hotD")
        @assert haskey(Jdict, "hotD")
        @assert haskey(ionsdict, "hotH2")
        @assert haskey(Jdict, "hotH2")
        @assert haskey(ionsdict, "hotHD")
        @assert haskey(Jdict, "hotHD")

        hot_H_network = [Jdict["hotH"]..., ionsdict["hotH"]...]
        hot_D_network = [Jdict["hotD"]..., ionsdict["hotD"]...]
        hot_H2_network = [Jdict["hotH2"]..., ionsdict["hotH2"]...]
        hot_HD_network = [Jdict["hotHD"]..., ionsdict["hotHD"]...]

        return hot_H_network, hot_D_network, hot_H2_network, hot_HD_network
    end
end

function get_Jrate_symb(reactant::String, products::Array)::Symbol
    #=
    All inputs are strings or arrays of strings.
    Formats a Jrate symbol.
    =#
    return Symbol("J$(reactant)to" * join(products, "p"))
end

function log_reactions(df, sheetname, spreadsheetname)
    #=
    Writes out the reactions as listed in df to the sheetname
    "$(sheetname)" in spreadsheetname.xlsx.
    =#

    themode = isfile(spreadsheetname) ? "rw" : "w"

    XLSX.openxlsx(spreadsheetname, mode=themode) do xf
        XLSX.addsheet!(xf,"$(sheetname)")
        sheet = xf["$(sheetname)"]
        sheet[1, :] = names(df)
        for r in 2:size(df)[1]+1, c in 1:size(df,2)
             sheet[XLSX.CellRef(r, c)] = df[r-1,c]
        end
    end
end

function make_k_expr(A, B, C, T::String, M2, M1, pow, BR)
    #=
    Constructs the modified Arrhenius equation, with leading multipliers for
    mass scaling and branching ratios.

    Inputs:
        A, B, C: Arrhenius equation parameters
        T: Temperature string to use; can be Tn, Ti, Te.
        M2, M1: Masses of the heavier and lighter isotope respectively. 1 if no mass scaling
        pow: -0.5. Passed in, in case it needs to change later. 
        BR: Branching ratio.
    Output:
        k = BR * (M2/M1)^pow * A * T^B * exp(C/T)
    =#
    if (B == 0) & (C == 0) # constant rate coefficient
        k = :($(((M2/M1)^pow)*BR*A))
    elseif (B != 0) & (C == 0)
        k = :($(((M2/M1)^pow)*BR*A) .* $(Meta.parse(T)) .^ $(B))
    elseif (B != 0) & (C != 0)
        k = :($(((M2/M1)^pow)*BR*A) .* $(Meta.parse(T)) .^ $(B) .* exp.($(C) ./ $(Meta.parse(T))))
    elseif (B == 0) & (C != 0)
        k = :($(((M2/M1)^pow)*BR*A) .* exp.($(C) ./ $(Meta.parse(T))))
    end
    return k
end

function make_Troe(k0_ABC, kinf_ABC, F, M2, M1, pow, BR)
    #=
    Make expression for type 3 reactions, of which there are currently none in the network (Feb 2022).

    Input:
        k0_ABC, kinf_ABC, kR_ABC: A, B, C arrhenius parameters for k0 (low pressure limit), kinf (high pressure limit), krad.
        F: Troe parameter
        M2, M1: Masses for heavier and lighter isotope for mass scaling
        pow: power (usually -0.5) for mass scaling
        BR: Branching ratio 
    Output: a symbolic expression for the rate coefficient
    =#
    k0 = make_k_expr(k0_ABC..., "Tn", M2, M1, pow, BR)
    kinf = make_k_expr(kinf_ABC..., "Tn", M2, M1, pow, BR)
    
    if F == 0
        return :(($(k0) .* $(kinf) .* M) ./ ($(k0) .* M .+ $(kinf)))  # Updated to match Roger's code 4 Feb 2021
    else        
        FF = troe_expr(k0, kinf, F)
        return :($(FF) .* ($(k0) .* $(kinf)) ./ ($(k0) .* M .+ $(kinf)))  # Updated to match Roger's code 5 Feb 2021
    end
end

function make_modified_Troe(k0_ABC, kinf_ABC, kR_ABC, F, M2, M1, pow, BR)
    #=
    "Type 4" pressure dependent association reactions. 
    Roger has the result for F=0 in his code as k_tot = k2 + k0*kinf*M/(k0*M+kinf), with k2 (kR in the appendix from Vuitton)
    being the general form and using the third set of A, B, C coefficients in his file. Those columns are 0 in the file he gave me, but
    I have coded it like this in case I ever get one where it's not 0.
    
    Input:
        k0_ABC, kinf_ABC, kR_ABC: A, B, C arrhenius parameters for k0 (low pressure limit), kinf (high pressure limit), krad.
        F: Troe parameter
        M2, M1: Masses for heavier and lighter isotope for mass scaling
        pow: power (usually -0.5) for mass scaling
        BR: Branching ratio 
    Output: a symbolic expression for the rate coefficient
    =#
    k0 = make_k_expr(k0_ABC..., "Tn", M2, M1, pow, BR)
    kinf = make_k_expr(kinf_ABC..., "Tn", M2, M1, pow, BR)
    kR = make_k_expr(kR_ABC..., "Tn", M2, M1, pow, BR)
    
    if F == 0
        return :($(kR) .+ ($(k0) .* $(kinf) .* M) ./ ($(k0) .* M .+ $(kinf)))
    else
        FF = troe_expr(k0, kinf, F)
        
        # for this one, we need the minimum of k_inf and a more complicated expression.
        # Format for minimum of A and B is: :(min.($:(A), $:(B)))
        return :(min.($:($kinf), $:($(kR) .+ ($(FF) .* $(k0) .* $(kinf) .* M) ./ ($(k0) .* M .+ $(kinf)))))
    end
end

function troe_expr(k0, kinf, F)
    #=
    Put together a really nasty expression
    =#
    outer_numer = :(log10.($(F)))
    inner_numer = :(log10.(($(k0) .* M) ./ ($(kinf))) .- 0.4 .- 0.67 .* log10.($(F)))
    inner_denom = :(0.75 .- 1.27 .* log10.($(F)) .- 0.14 .* (log10.(($(k0) .* M) ./ ($(kinf))) .- 0.4 .- 0.67 .* log10.($(F))))
    outer_denom = :(1 .+ (($(inner_numer)) ./ ($(inner_denom))) .^ 2)
    FF = :(10 .^ (($(outer_numer)) ./ ($(outer_denom))))
    return FF
end