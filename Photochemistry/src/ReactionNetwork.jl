# **************************************************************************** #
#                                                                              #
#                   Reaction network formatting and loading                    #
#                                                                              #
# **************************************************************************** #

#===============================================================================#
#                      Functions to load and manipulate the network             #
#===============================================================================#

function calculate_and_write_column_rates(rxn_filename, atm_state, n_horiz; globvars...)
    #=
    Calculates the column rates for all reactions in rxn_filename based on species densities in atm_state
    and writes them back out into rxn_filename. Creates separate columns for each horizontal position (ihoriz).
    =#

    GV = values(globvars)
    required = [:all_species, :dz, :ion_species, :num_layers, :reaction_network,
                :results_dir, :sim_folder_name, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)

    flush(stdout)
    println("Writing out column rates to the reaction log...")
    flush(stdout)

    # Jrates needed for photochemical reactions
    Jratedict = Dict([j => atm_state[j] for j in keys(atm_state) if occursin("J", string(j))])

    # Loop through sheets in the spreadsheet
    the_spreadsheet_file = GV.results_dir * GV.sim_folder_name * "/" * rxn_filename
    active_rxns_original = XLSX.readxlsx(the_spreadsheet_file)
    original_sheets = XLSX.sheetnames(active_rxns_original)

    for sheet in original_sheets
        if occursin("Unused", sheet) || occursin("Sheet1", sheet)
            continue
        end

        flush(stdout)
        println("Working on sheet $(sheet)")
        flush(stdout)

        df = DataFrame(XLSX.readtable(the_spreadsheet_file, sheet))

        # Add new columns for each horizontal column separately
        for ihoriz in 1:n_horiz
            insertcols!(df, size(df)[2]+1, Symbol("ColumnRate_$(ihoriz)") => zeros(Float64, size(df)[1]))
        end

        # Loop through reactions
        for (row_index, row) in enumerate(eachrow(df))
            rcols, pcols = get_product_and_reactant_cols(df)

            these_reactants = [Symbol(row[r]) for r in rcols if row[r] != "none"]
            these_products = [Symbol(row[p]) for p in pcols if row[p] != "none"]

            # Find matching reaction
            rxn_i = findfirst(s -> (s[1] == these_reactants && s[2] == these_products), GV.reaction_network)
            this_rxn = GV.reaction_network[rxn_i]

            this_rxn_func = mk_function(:((Tn, Ti, Te, M) -> $(this_rxn[3])))

            # Loop over horizontal columns to calculate rates
            for ihoriz in 1:n_horiz
                # explicitly get temperatures [horizontal, vertical]
                Tn_col = GV.Tn[ihoriz, :]
                Ti_col = GV.Ti[ihoriz, :]
                Te_col = GV.Te[ihoriz, :]

                # Total atmospheric density explicitly for this horizontal column
                Mtot = n_tot(atm_state, ihoriz; GV.all_species)

                vol_rate_by_alt = get_volume_rates(
                    these_reactants[1], this_rxn, this_rxn_func,
                    atm_state, Mtot, ihoriz;
                    Jratedict,
                    Tn = Tn_col[2:end-1],
                    Ti = Ti_col[2:end-1],
                    Te = Te_col[2:end-1],
                    globvars...
                )

                # Calculate and store column rates for each horizontal column separately
                col_rate = sum(vol_rate_by_alt .* GV.dz)

                # Insert column-specific rate into its corresponding new column
                df[row_index, Symbol("ColumnRate_$(ihoriz)")] = col_rate

                if col_rate == 0
                    println("$(this_rxn) in column $(ihoriz) has a 0 col rate. This may or may not be expected.")
                end
            end
        end

        # Write the updated DataFrame back to the spreadsheet (now multiple columns)
        for ihoriz in 1:n_horiz
            col_symbol = Symbol("ColumnRate_$(ihoriz)")
            add_column(df[!, col_symbol], string(col_symbol), sheet, the_spreadsheet_file)
        end

        println("Completed sheet $(sheet)")
    end

    println("Finished adding column rates to $(the_spreadsheet_file)")
end

function filter_network(sp::Symbol, rate_type::String, species_role::String; globvars...)
    #= 
    Given a reaction network, filter it down by:
    sp: A specific species, produced by
    rate_type reactions (Jrates, krates, or all), where the species is a
    species_role (reactant, product, or any).

    Returns the reaction network but without reactions that don't meet the criteria.
    =#

    GV = values(globvars)
    required = [:reaction_network]
    check_requirements(keys(GV), required)

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

function find_duplicates(rxnnet)
    #=
    Can be used periodically when reactions are added to check for duplicates. 
    
    Note that it will print the row numbers of the duplicates, but that doesn't necessarily correspond to the spreadsheet, but rather 
    the position within the Julia reaction_network array.
    
    Inputs:
        rxnnet: the Julia array of reactions. 
    
    Outputs:
        Print statements. your responsibility to remove duplicates from the spreadsheet!
    =#
    
    i = 1
    for rxn in rxnnet
        reactants = get_counts(rxn[1])
        products = get_counts(rxn[2])
        
        j = i + 1
        
        j_for_printing = j
        for other_rxn in rxnnet[j:end]
            other_reactants = get_counts(other_rxn[1])
            other_products = get_counts(other_rxn[2])
            
            if (Set(reactants) == Set(other_reactants)) & (Set(products) == Set(other_products))
                println("Found Duplicates:")
                println("Entry $(i) ", rxn)
                println("Entry $(j_for_printing) ", other_rxn)
                println()
            end
            j_for_printing += 1
        end
        i += 1
    end       
end

function load_network_and_make_functions(rxn_sheet; globvars...)
    #=
    Using the reaction spreadsheet, generates symbolic chemistry network and sub networks for hot atoms.

    Input: 
        rxn_sheet: path to a spreasheet containing chemical reactions.

    Outputs: 
        reaction_network: array of the symbolic network of reactions in the format [[reactants...], [products...], :(rate)]
        hHnet, hDnet, hH2net, hHDnet: smaller arrays containing the reactions which produce hot atoms of the specified type
        hot_H_rc_funcs, etc: special functions for each chemical reaction which can be evaluated on the fly in a fast manner.
                             critical for enabling the calculation of non-thermal escape of hot atoms.
    =#
    GV = values(globvars)
    required = [:all_species]
    check_requirements(keys(GV), required)

    reaction_network, hHnet, hDnet, hH2net, hHDnet = load_reaction_network(rxn_sheet; get_hot_rxns=true, GV.all_species);
    hot_H_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hHnet]);
    hot_D_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hDnet]);
    hot_H2_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hH2net]);
    hot_HD_rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in hHDnet]);
    
    return reaction_network, hHnet, hDnet, hH2net, hHDnet, hot_H_rc_funcs, hot_D_rc_funcs, hot_H2_rc_funcs, hot_HD_rc_funcs
end

function load_reaction_network(spreadsheet; saveloc=nothing, write_rxns=false, to_return="all", ions_on=true, get_hot_rxns=true, globvars...)
    #=
    Inputs:
        spreadsheet: path and filename of the reaction network spreadsheet
        saveloc: Path in which to store the re-written spreadsheet of active reactions. This should be the simulation result folder. 
    Output:
        reaction_network, a vector of vectors in the format
        [Symbol[reactants...], Symbol[products...], :(rate coefficient expression)
    =#
    GV = values(globvars)
    required = [:all_species]
    check_requirements(keys(GV), required)

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

function log_reactions(df, sheetname, spreadsheetname)
    #=
    Writes out the reactions as listed in df to the sheetname
    "$(sheetname)" in spreadsheetname.xlsx.
    =#

    themode = isfile(spreadsheetname) ? "rw" : "w"

    XLSX.openxlsx(spreadsheetname, mode=themode) do xf
        if "Sheet1" in XLSX.sheetnames(xf)
            sheet1 = xf[1]
            XLSX.rename!(sheet1, sheetname)
        else
            XLSX.addsheet!(xf,"$(sheetname)")
        end
        
        sheet = xf["$(sheetname)"]
        sheet = xf["$(sheetname)"]
        sheet[1, :] = names(df)
        for r in 2:size(df)[1]+1, c in 1:size(df,2)
             sheet[XLSX.CellRef(r, c)] = df[r-1,c]
        end
    end
end

function rxns_where_species_is_observer(sp, chemnet)
    #=
    Finds reactions where a given chemical species is a non-reacting third body ("observer")
    
    Input:
        sp: Species which may be an observer
        chemnet: chemistry reaction network
    Output: either:
        nothing: if the species was not found to be an observer in any reaction
        twosides: List of reactions where the species is found on both LHS and RHS
    =#

    krate_rxns = filter(x->!occursin("J", string(x[3])), chemnet)
    
    twosides = []

    flag = false
    for rxn in krate_rxns
        if in(sp, rxn[1]) && in(sp, rxn[2])
            flag = true
            twosides = push!(twosides, rxn)
        end
    end

    if flag == true 
        return twosides
    else
        return nothing
    end
end


#===============================================================================#
#              Functions to initially format the network object                 #
#===============================================================================#

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
        n += 1
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

function get_Jrate_symb(reactant::String, products::Array)::Symbol
    #=
    All inputs are strings or arrays of strings.
    Formats a Jrate symbol.
    =#
    return Symbol("J$(reactant)to" * join(products, "p"))
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

#                        Functions to calculate enthalpies                      #
#                                                                               #
# ALERT: YOU *MUST* RUN modify_rxn_spreadsheet once the first time you try to   #
# generate results for a new planet because it has a DIFFERENT ESCAPE V!        #
# Also run if you add new reactions to the spreadsheet or change rate coefs!    #
#===============================================================================#

function escape_velocity(; globvars...)
    #=
    Calculates escape velocity for the planet described by the supplied global parameters.
    this kinda doesn't make sense but it's the way to do it that is most consistent with 
    the rest of the module.
    =#
    GV = values(globvars)
    required = [:M_P, :R_P, :zmax]
    check_requirements(keys(GV), required)

    return sqrt(bigG * GV.M_P / (GV.R_P + GV.zmax)) / 100 # 100 converts from cm/s to m/s, bc escape energy function is written in MKS.
end 

function escape_energy(z; globvars...)
    #=
    Calculates escape energy from planet P for a molecule containing z protons.
    =#
    
    kJ_to_eV((0.5 * z * (1.67e-27 #=kg=#) * (escape_velocity(; globvars...) #=m/s=#)^2) / 1000 #=kJ/J=#)
end

function enthalpy_of_reaction(reactants, products, enthalpy_dict; globvars...)
    #=
    Calculates the total enthalpy of reaction, amount needed to escape hot atoms and molecules in the products,
    and the total excess to figure out if its exothermic or not. Note that it just automatically looks for H, D, H2, HD right now...

    Inputs:
        reactants, products: lists of exactly what it sounds like, with symbols for the species.
        enthalpy_dict: dictionary containing enthalpies of formation for many species found in the network. 
    =#
    
    dfH_products = 0
    dfH_reactants = 0

    for p in products
        dfH_products += enthalpy_dict[p]
    end
    for r in reactants
        dfH_reactants += enthalpy_dict[r]
    end

    # Calculate the total enthalpy of reaction:
    enthalpy_of_rxn = dfH_products - dfH_reactants
    
    # this block adds up the total mass that needs to escape, allows for reactions that produce, e.g., both H and D
    m = 0 
    flag = 0
    if :H in products
        m += 1
        flag += 1
    end
    if :H2 in products
        m += 2
        flag += 1
    end
    if :D in products
        m += 2
        flag += 1
    end
    if :HD in products
        m += 3
        flag += 1
    end

    # Determine if exothermic or endothermic, using the standard convention:
    # ΔfH0 < 0 ==> exothermic, ΔfH0 > 0 ==> endothermic. 
        
    # Convert the total enthalpy of reaction to eV.
    # Change the sign by multiplying by -1, since ΔfH0 < 0 (negative, exothermic) means there
    # excess energy (i.e. heat) is a product, but ΔfH0 > 0 (positive, endothermic means energy
    # is consumed in the reaction.
    total_reaction_enthalpy_ev = -1*ev_per_molecule(kJ_to_eV(enthalpy_of_rxn))
    
    # Calculate how much energy goes into escaping the hot atoms:
    energy_required_to_escape_all_hot = escape_energy(m; globvars...)
    
    # Calculate how much excess energy we have
    excess_energy = round(total_reaction_enthalpy_ev - energy_required_to_escape_all_hot, digits=2)

    # Calculate whether endothermic or exothermic
    can_drive_nonthermal_escape = excess_energy > 0 ? true : false
    
    if flag > 1
        println("Flag! Reaction produces two hot atoms")
    end
    
    return can_drive_nonthermal_escape, total_reaction_enthalpy_ev, energy_required_to_escape_all_hot, excess_energy
end

function get_product_and_reactant_cols(df)
    # Determine how many reactant and product columns there are 
    possible_Rcols = ["R1", "R2", "R3"]
    possible_Pcols = ["P1", "P2", "P3"]
    rcols = Symbol.(possible_Rcols[possible_Rcols .∈ Ref(names(df))])
    pcols = Symbol.(possible_Pcols[possible_Pcols .∈ Ref(names(df))])
    
    return rcols, pcols
end

function calculate_enthalpies(df; species=[:H, :D, :H2, :HD], new_cols=nothing, insert_i=nothing, globvars...)
    #=
    Inputs:
        df: A dataframe containing the contents of a reaction spreadsheet to which we must add excess energies
        species: species for which we want to identify reactions producing hot ones.
        new_cols: strings containing names for new columns to put into the df/spreadsheet.
        insert_i: first integer at which we will begin inserting columns.
    =#
    
    #Enthalpy of formation 
    enthalpy_df = DataFrame(XLSX.readtable("Resources/Enthalpies_of_Formation.xlsx", "enthalpy"))
    
    enthalpy = Dict([Symbol(k)=>df_lookup(enthalpy_df, "Species", k, "Enthalpy of formation (kJ/mol)")[1] for k in enthalpy_df."Species"])

    # Add new columns
    if new_cols != nothing
        j = 0
        for name in new_cols
            insertcols!(df, insert_i+j, "$(name)"=>["" for x in collect(1:size(df)[1])])
            j += 1
        end
    end
    
    rcols, pcols = get_product_and_reactant_cols(df)
    
    # Count endothermic and exothermic
    not_escape_driver_count = 0
    escape_driver_count = 0
    
    for row in eachrow(df)
        reactants = [Symbol(row.:($r)) for r in rcols]
        products = filter!(x->x!=:none, [Symbol(row.:($p)) for p in pcols])
        
        # Ignore reactions for which we don't have enthalpies of a species:
        if any(x->!(x in keys(enthalpy)), union(reactants, products))
            continue
        end

        if any(x->x in products, species)
            drives_NT_escape, total_enthalpy, loss_energy, excess_energy = enthalpy_of_reaction(reactants, products, enthalpy; globvars...)
            row.rxnEnthalpy = string(total_enthalpy)
            row.totalEscE = string(loss_energy)
            row.excessE = string(excess_energy)
            
            if drives_NT_escape == true
                escape_driver_count += 1
                row.NTEscape = "Yes"

                if :D in products
                    row.hotD = "Yes"
                end
                if :HD in products
                    row.hotHD = "Yes"
                end
                
                if :H in products
                    row.hotH = "Yes"
                end
                if :H2 in products
                    row.hotH2 = "Yes"
                end
            else # if drives_NT_escape == false
                not_escape_driver_count += 1

                # Handle special case where charge exchange should always be counted as producing hot atoms
                println("reactants $(Set(reactants)), products $(Set(products))" )
                println()
                

                if (Set(reactants)==Set([:Hpl, :H])) & (Set(products)==Set([:Hpl, :H])) # Resonant
                    row.NTEscape = "Yes"
                    row.hotH = "Yes"
                    not_escape_driver_count -= 1
                    escape_driver_count += 1
                elseif (Set(reactants)==Set([:Dpl, :H])) & (Set(products)==Set([:Hpl, :D])) # deuterated
                    row.NTEscape = "Yes"
                    row.hotD = "Yes"
                    not_escape_driver_count -= 1
                    escape_driver_count += 1
                elseif (Set(reactants)==Set([:Hpl, :O])) & (Set(products)==Set([:Opl, :H])) # with O
                    row.NTEscape = "Yes"
                    row.hotH = "Yes"
                    not_escape_driver_count -= 1
                    escape_driver_count += 1
                elseif (Set(reactants)==Set([:Dpl, :O])) & (Set(products)==Set([:Opl, :D])) # D with O
                    row.NTEscape = "Yes"
                    row.hotD = "Yes"
                    not_escape_driver_count -= 1
                    escape_driver_count += 1
                else
                    row.NTEscape = "No"
                    if :D in products
                        row.hotD = ""
                    end
                    if :HD in products
                        row.hotHD = ""
                    end
                    
                    if :H in products
                        row.hotH = ""
                    end
                    if :H2 in products
                        row.hotH2 = ""
                    end
                end
            end
        end
    end
    
    println("Drives non-thermal escape (exothermic): $(escape_driver_count)\nDoes not drive non-thermal escape (endothermic): $(not_escape_driver_count)")
    return df
end

function modify_rxn_spreadsheet(spreadsheet; new_file="REACTION_NETWORK_NEW.xlsx", spc=[:H, :D, :H2, :HD], new_cols=nothing, insert_i=[0,7,8,8], globvars...)
    #=
    Inputs:
        spreadsheet: A starting spreadsheet with reaction rate data.
        spc: spc for which we want to determine if hot ones are produced.
        new_cols: strings containing names for new columns to put into the df/spreadsheet.
        insert_i: a list of integers at which to insert columns. The order in this list corresponds to the order of sheets in the workbook.
                  so, if you have 4 sheets, and you want to insert the new columns at col index 2, 4, 6, and then 8 in sheets 1, 2, 3, and 4, 
                insert_i = [2, 4, 6, 8]. Note that it only works for sheets involving ions so you have to put in a zero to skip the neutral sheet, which 
                is usually first. Typical usage is insert_i=[0, 7, 8, 8].
    =#
    xf = XLSX.readxlsx(spreadsheet)
    original_sheets = XLSX.sheetnames(xf)
    
    
    if spc != [:H, :D, :H2, :HD]
        throw("Error: The code that calculates hot atom excess energies is not set up to handle extra species beyond H, D, H2, HD.")
    end
    
    # Go through each available sheet and open it as a dataframe for modification
    for (j, sheet) in enumerate(original_sheets)
        df = DataFrame(XLSX.readtable(spreadsheet, sheet));
               
        # Now process:        
        if sheet in ["Ion reactions", "Photodissociation", "Photoionization"]
            # Replace any missing values so we don't get rid of all our reactions before processing:
            rcols, pcols = get_product_and_reactant_cols(df)
            for r in rcols
                replace!(df.:($r), missing=>"none")
            end
            for p in pcols
                replace!(df.:($p), missing=>"none")
            end
            
            if insert_i != nothing
                df_to_write = calculate_enthalpies(df; species=spc, new_cols=new_cols, insert_i=insert_i[j], globvars...)
            else 
                df_to_write = calculate_enthalpies(df; species=spc, globvars...)
            end
        else
            df_to_write = df
        end
        println()
        log_reactions(df_to_write, sheet, new_file)
    end
end


