# **************************************************************************** #
#                                                                              #
#                       File input/output and logging                          #
#                                                                              #
# **************************************************************************** #
function create_folder(foldername::String, parentdir::String)
    #=
    Checks to see if foldername exists within parentdir. If not, creates it.
    =#
    # println("Checking for existence of $(foldername) folder in $(parentdir)")
    dircontents = readdir(parentdir)
    if foldername in dircontents
        println("Folder $(foldername) already exists")
    else
        mkdir(parentdir*foldername)
        println("Created folder ", foldername)
    end
end

function delete_old_h5file(filename::String)
    #= 
    the HDF5 module doesn't allow implicit deletion of existing names in the function
    h5write, so this function will delete previously xisting files with the name 'filename',
    i.e. in cases where a run did not complete successfully and left behind garbage. 
    =#
    current_folder = pwd()
    dircontents = readdir(current_folder)
    if filename in dircontents
        println("File $(filename) exists, deleting")
        rm(filename)
    end
end

function get_elapsed_time(filepath)
    #=
    Given an .h5 file located at filepath, this extracts the total elapsed time in seconds. 
    =# 
    return parse(Float64, split(h5read(filepath, "info")[end])[1])
end

function get_ncurrent(readfile::String)
    #=
    Input:
        readfile: HDF5 file containing a matrix with atmospheric density profiles
    Output:
        n_current: dictionary of atmospheric density profiles by species 
    =#

    # This accounts for the old format of some files. 
    try
        global n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
        global n_current_mat = h5read(readfile,"n_current/n_current_mat");
    catch
        global n_current_tag_list = map(Symbol, h5read(readfile,"ncur/species"))
        global n_current_mat = h5read(readfile,"ncur/ncur_mat");
    end
    
    n_current = Dict{Symbol, Array{ftype_ncur, 1}}()

    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies], length(n_current_mat[:, ispecies]))
    end
    return n_current
end

# searches path for key
searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function search_subfolders(path::String, key; type="folders")
    #=
    path: a folder containing subfolders and files.
    key: the text pattern in folder names that you wish to find.
    type: whether to look for "folders" or "files". 

    Searches the top level subfolders within path for all folders matching a 
    certain regex given by key. Does not search files or sub-subfolders.
    =#
    folderlist = []
    filelist = []
    if type=="folders"
        for (root, dirs, files) in walkdir(path)
            if root==path
                for dir in dirs
                    push!(folderlist, joinpath(root, dir)) # path to directories
                end
            end
        end
    elseif type=="files"
        for (root, dirs, files) in walkdir(path)
            if root==path
                for fil in files
                    push!(filelist, fil) # files
                end
            end
        end
    end

    if type=="folders"
        folderlist = filter(x->occursin(key, x), folderlist)
        return folderlist
    elseif type=="files"
        filelist = filter(x->occursin(key, x), filelist)
        return filelist
    end
end

function write_atmosphere(atmdict::Dict{Symbol, Vector{ftype_ncur}}, filename::String; t=0, globvars...) 
    #=
    Writes out the current atmospheric state to an .h5 file

    Input: 
        atmdict: atmospheric density profile dictionary
        filename: filename to write to
    =# 
    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:alt, :num_layers, :hrshortcode, :rshortcode])

    sorted_keys = sort(collect(keys(atmdict)))
    atm_mat = Array{Float64}(undef, GV.num_layers, length(sorted_keys));

    for ispecies in [1:length(sorted_keys);]
        for ialt in [1:GV.num_layers;]
            atm_mat[ialt, ispecies] = convert(Float64, atmdict[sorted_keys[ispecies]][ialt])
        end
    end
    delete_old_h5file(filename)
    h5open(filename, "w") do f # this syntax is ok because we never write multiple times to a file.
        write(f, "n_current/n_current_mat", atm_mat)
        write(f, "n_current/alt", GV.alt)
        write(f, "n_current/species", map(string, sorted_keys))
        write(f, "info", ["hrshortcode:" "$(GV.hrshortcode)"; "random shortcode:" "$(GV.rshortcode)"; "elapsed t:" " $(t) sec"])
    end
end

function write_final_state(atmdict, thedir, thefolder, fname; globvars...)
    #=
    Write out the final atmosphere to a file, first making sure the current Jrates are included.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:alt, :external_storage, :num_layers, :hrshortcode, :Jratedict, :rshortcode])

    # Make sure we have the updated Jrates
    println("Adding stored Jrates to the final output file")
    for j in keys(GV.Jratedict)
        atmdict[j] = GV.external_storage[j]
    end

    # Write out final atmosphere
    write_atmosphere(atmdict, thedir*thefolder*"/"*fname; globvars...)
    println("Saved final atmospheric state")
end

#                             Logging functions                                 #
#===============================================================================#


function get_param(param, df)
    #=
    Retrieve a particular parameter from a parameter log spreadsheet opened as dataframe df.
    This assumes that the parameter you want is listed in the first column. 
    =#
    key = names(df)[1]
    entry = names(df)[2]
    return filter(row -> row.:($key)==param, df).:($entry)[1]
end

function load_bcdict_from_paramdf(df)
    #=
    Special function to load the boundary condition variable from the log spreadsheet.
    Assumes that the column names are 'Species', 'Type', 'Lower' and 'Upper.'

    Output:
        speciesbclist_reconstructed: full speciesbclist dictionary object
    =#

    # Some of the older runs might have just "flux" in the log, which really means thermal flux alone
    typerevdict = Dict("thermal flux"=>"f", "flux"=>"f", "nonthermal flux"=>"ntf", "velocity"=>"v", "density"=>"n")
    speciesbclist_reconstructed = Dict()

    for s in unique(df.Species)
        speciesbclist_reconstructed[Symbol(s)] = Dict()
        # get just the rows with the species
        thisspecies = filter(row->row.Species==s, df)

        # The nonthermal flux isn't in here as a function because it would be hard to specify it there. This could be fixed later probably
        for row in eachrow(thisspecies)
            try 
                speciesbclist_reconstructed[Symbol(s)][typerevdict[row.Type]] = [parse(Float64, row.Lower), parse(Float64, row.Upper)]
            catch y
                speciesbclist_reconstructed[Symbol(s)][typerevdict[row.Type]] = [parse(Float64, row.Lower), row.Upper]
            end
        end
    end
    return speciesbclist_reconstructed
end

function load_from_paramlog(folder; quiet=true, globvars...)
    #=
    Given a folder containing simulation results, this will open the parameter log spreadsheet, 
    load as dataframe, and extract all the entries so as to return the global parameters that were used for
    that simulation. 
    =#

    # Load the workbook
    paramlog_wb = XLSX.readxlsx("$(folder)PARAMETERS.xlsx")

    # Basic variables
    df_gen = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "General"));
    ions_included = get_param("IONS", df_gen)
    hrshortcode = get_param("RSHORTCODE", df_gen)
    rshortcode = get_param("HRSHORTCODE", df_gen)
    rxn_spreadsheet = get_param("RXN_SOURCE", df_gen)

    # Species lists
    df_splists = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "SpeciesLists"));
    neutral_species = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.Neutrals)]
    ion_species = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.Ions)]
    all_species = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.AllSpecies)]
    no_transport_species = [Symbol(x) for x in [filter(x->typeof(x)==String, df_splists.NoTransport)]...];
    no_chem_species = [Symbol(x) for x in [filter(x->typeof(x)==String, df_splists.NoChem)]...];
    transport_species = setdiff(all_species, no_transport_species);
    chem_species = setdiff(all_species, no_chem_species);

    # Atmospheric conditions
    df_atmcond = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "AtmosphericConditions"));

    ## Temperatures first
    if "TemperatureArrays" in XLSX.sheetnames(paramlog_wb)
        df_temps = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "TemperatureArrays"));
        Tn_arr = df_temps.Neutrals
        Ti_arr = df_temps.Ions
        Te_arr = df_temps.Electrons
    else 
        GV = values(globvars)
        @assert all(x->x in keys(GV), [:alt])
        if quiet==false
            println("WARNING: Reconstructing temperature profiles with default options based on logged control temperatures. It is POSSIBLE the reconstruction could be wrong if the temp function changed.")
        end
        T_dict = T(get_param("TSURF", df_atmcond), get_param("TMESO", df_atmcond), get_param("TEXO", df_atmcond); alt)
        Tn_arr = T_dict["neutrals"]
        Ti_arr = T_dict["ions"]
        Te_arr = T_dict["electrons"]
    end

    Tplasma_arr = Ti_arr .+ Te_arr;
    Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr);
    Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
    Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)]; molmass) for sp in all_species]); 
    water_bdy = get_param("WATER_BDY", df_atmcond) * 1e5 # It's stored in km but we want it in cm

    # Boundary conditions
    df_bcs = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "BoundaryConditions"));
    speciesbclist = load_bcdict_from_paramdf(df_bcs);
    
    vardict = Dict("ions_included"=>ions_included,
                   "hrshortcode"=>hrshortcode,
                   "rshortcode"=>rshortcode,
                   "neutral_species"=>neutral_species,
                   "ion_species"=>ion_species,
                   "all_species"=>all_species,
                   "transport_species"=>transport_species,
                   "chem_species"=>chem_species,
                   "Tn_arr"=>Tn_arr,
                   "Ti_arr"=>Ti_arr,
                   "Te_arr"=>Te_arr,
                   "Tplasma_arr"=>Tplasma_arr,
                   "Tprof_for_Hs"=>Tprof_for_Hs,
                   "Tprof_for_diffusion"=>Tprof_for_diffusion,
                   "Hs_dict"=>Hs_dict,
                   "speciesbclist"=>speciesbclist,
                   "rxn_spreadsheet"=>rxn_spreadsheet,
                   "water_bdy"=>water_bdy)

    return vardict
end

function write_to_log(full_path::String, entries; mode="a")
    #=
    Inputs;
        full_path: full path to log file (folder structure and filename, with extension)
        entries: List of strings to write to the log file.
        Optional:
            mode: w or a for write or append.
    =#

    f = open(full_path, mode) 
    if isa(entries, Array)
        for e in entries 
            write(f, string(e)*"\n")
        end
    elseif isa(entries, String)
        write(f, entries*"\n")
    else
        throw("Wrong format for logging: $(typeof(entries))")
    end

    close(f)
end