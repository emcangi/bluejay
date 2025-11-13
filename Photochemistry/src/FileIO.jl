# **************************************************************************** #
#                                                                              #
#                       File input/output and logging                          #
#                                                                              #
# **************************************************************************** #
function add_column(data, colname, sheetname, filename)
    #=
    Adds a column containing the data within data to the sheet sheetname inside the excel spreadsheet filename.
    =#
    
    themode = isfile(filename) ? "rw" : "w"
    
    XLSX.openxlsx(filename, mode=themode) do xf
        sheet = xf["$(sheetname)"]
        
        df = DataFrame(XLSX.readtable(filename, sheetname))
        new_col_num = size(df)[2] + 1

        sheet[:, new_col_num] = [colname, data...]
    end
end

function create_folder(foldername::String, parentdir::String)
    #=
    Creates a folder within an enclosing directory.

    Inputs:
        foldername: folder to be created
        parentdir: directory in which to create foldername; will be created if it doesn't already exist.

    Outputs:
        None
    =#
    # Create the parent directory path if it does not exist
    if !isdir(parentdir)
        mkpath(parentdir)
    end

    target = joinpath(parentdir, foldername)

    if isdir(target)
        println("Folder $(foldername) already exists")
    else
        mkpath(target)
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

function get_ncurrent(readfile::String; globvars...)
    #=
    Input:
        readfile: HDF5 file containing a matrix with atmospheric density profiles
        globvars: keyword arguments including n_horiz and num_layers
    Output:
        n_current: dictionary of atmospheric density profiles by species 
    =#

    GV = values(globvars)
    
    # This accounts for the old format of some files. 
    try
        global species_list = map(Symbol, h5read(readfile,"n_current/species"))
        global n_current_mat = h5read(readfile,"n_current/n_current_mat");
    catch
        global species_list = map(Symbol, h5read(readfile,"ncur/species"))
        global n_current_mat = h5read(readfile,"ncur/ncur_mat");
    end
    
    n_current = Dict{Symbol, Vector{Array{ftype_ncur}}}()

    # Determine desired number of altitude layers from globvars if available,
    # otherwise fall back to Main.num_layers, or finally the file dimensions.
    if :num_layers in keys(GV)
        target_layers = GV.num_layers
    elseif isdefined(Main, :num_layers)
        target_layers = Main.num_layers
    else
        # For 2D arrays: (layers, species). For 3D arrays: (n_horiz, layers, species)
        target_layers = ndims(n_current_mat) == 3 ? size(n_current_mat, 2) : size(n_current_mat, 1)
    end
    
    # Determine desired number of horizontal columns from globvars if available,
    # otherwise fall back to Main.n_horiz, or finally default to 1.
    if :n_horiz in keys(GV)
        n_horiz = GV.n_horiz
    elseif isdefined(Main, :n_horiz)
        n_horiz = Main.n_horiz
    else
        n_horiz = 1
    end

    # Slice or pad the profiles so that their length matches the altitude grid
    for ispecies in eachindex(species_list)
        profile = convert(Vector{ftype_ncur}, n_current_mat[:, ispecies])

        if length(profile) > target_layers
            profile = profile[1:target_layers]
        elseif length(profile) < target_layers
            profile = vcat(profile, fill(last(profile), target_layers - length(profile)))
        end

        # Duplicate the profile for each column so that updates to one column do not modify the others by aliasing the same array.
        n_current[species_list[ispecies]] = [copy(profile) for _ in 1:n_horiz]
    end
    if isdefined(Main, :all_species)
        zero_prof = fill(ftype_ncur(0.0), target_layers)
        for sp in setdiff(Main.all_species, keys(n_current))
            n_current[sp] = fill(copy(zero_prof), n_horiz)
        end
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

function write_atmosphere(atmdict::Dict{Symbol, Vector{Array{ftype_ncur}}}, filename::String; t=0, globvars...) 
    #=
    Writes out the current atmospheric state to an .h5 file

    Input: 
        atmdict: atmospheric density profile dictionary
        filename: filename to write to
        globvars: keyword arguments including n_horiz, num_layers, alt, etc.
    =# 
    GV = values(globvars)
    required =  [:alt, :num_layers, :hrshortcode, :rshortcode]
    check_requirements(keys(GV), required)
    
    # Determine number of horizontal columns from globvars if available,
    # otherwise fall back to Main.n_horiz, or default to 1.
    if :n_horiz in keys(GV)
        n_horiz = GV.n_horiz
    elseif isdefined(Main, :n_horiz)
        n_horiz = Main.n_horiz
    else
        n_horiz = 1
    end
    
    sorted_keys = sort(collect(keys(atmdict)))
    atm_mat = Array{Float64}(undef, n_horiz, GV.num_layers, length(sorted_keys));

    for ihoriz in 1:n_horiz
    	for ispecies in [1:length(sorted_keys);]
            for ialt in [1:GV.num_layers;]
            	atm_mat[ihoriz, ialt, ispecies] = convert(Float64, atmdict[sorted_keys[ispecies]][ihoriz][ialt])
	    end
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
    required = [:alt, :external_storage, :num_layers, :hrshortcode, :Jratedict, :rshortcode]
    check_requirements(keys(GV), required)

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
        speciesbclist_reconstructed: full speciesbclist_vert dictionary object
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

    GV = values(globvars)
    required = [:molmass] 
    check_requirements(keys(GV), required)

    # Load the workbook
    paramlog_wb = XLSX.readxlsx("$(folder)PARAMETERS.xlsx")

    # Read as many of the parameter sheets as we can right now
    df_gen = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "General"));
    df_atmcond = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "AtmosphericConditions"));
    df_bcs = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "BoundaryConditions"));

    # Planet, M_P, R_P, and altitude, which were not logged in older runs.
    if ~(:M_P in keys(GV)) | ~(:R_P in keys(GV))
        try
            global planet = get_param("PLANET", df_gen)
            global M_P = get_param("M_P", df_gen)
            global R_P = get_param("R_P", df_gen)
        catch y
            println("WARNING: Exception: $(y) - you are trying to load parameters which aren't logged. File probably made before module updates.")
            println("Please pass in the following global variables: planet, M_P, R_P, then re-run this command.")
            println()
        end
    end

    if ~(:alt in keys(GV))
        try 
            global df_alt = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "AltGrid"));
            global df_altinfo = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "AltInfo"));
        catch y
            println("WARNING: Exception: $(y) - you tried to load the altitude grid but it's not logged. File probably made before module updates. Please pass in alt manually")
            println()
        end
    else
        global alt = GV.alt # This seems stupid but has to be done because apparently the if/else blocks can't access the keyword arg alt??
    end

    # AltInfo
    if (@isdefined df_altinfo) & (@isdefined df_alt)
        global alt = df_alt.Alt
        global non_bdy_layers = [parse(Float64, f) for f in collect(skipmissing(df_alt.non_bdy_layers))]
        global zmin = get_param("zmin", df_altinfo)
        global dz = get_param("dz", df_altinfo)
        global zmax = get_param("zmax", df_altinfo)
        global n_all_layers = get_param("n_all_layers", df_altinfo)
        global num_layers = get_param("num_layers", df_altinfo)
        global upper_lower_bdy = get_param("upper_lower_bdy", df_altinfo)
        global upper_lower_bdy_i = get_param("upper_lower_bdy_i", df_altinfo)
    else 
        # In this case, alt should have been passed in manually.
        global non_bdy_layers = alt[2:end-1]
        global zmin = alt[1]
        global dz = alt[2] - alt[1]
        global zmax = alt[end]
        global n_all_layers = length(alt)
        global num_layers = length(non_bdy_layers)
        # In older versions of the parameter logging, this variable was named "WATER_BDY".
        global upper_lower_bdy = get_param("WATER_BDY", df_atmcond) * 1e5
        # Ideally we would not repeat the code for n_alt_index here, but oh well.
        global upper_lower_bdy_i = Dict([z=>clamp((i-1),1, num_layers) for (i, z) in enumerate(alt)])[upper_lower_bdy]

    end

    # Codes
    hrshortcode = get_param("RSHORTCODE", df_gen)
    rshortcode = get_param("HRSHORTCODE", df_gen)

    # Basics - more
    rxn_spreadsheet = get_param("RXN_SOURCE", df_gen)
    DH = get_param("DH", df_gen)
    
    # Species lists
    ions_included = get_param("IONS", df_gen)
    df_splists = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "SpeciesLists"));
    neutral_species = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.Neutrals)]
    ion_species = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.Ions)]
    all_species = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.AllSpecies)]
    no_transport_species = [Symbol(x) for x in [filter(x->typeof(x)==String, df_splists.NoTransport)]...];
    no_chem_species = [Symbol(x) for x in [filter(x->typeof(x)==String, df_splists.NoChem)]...];
    transport_species = setdiff(all_species, no_transport_species);
    chem_species = setdiff(all_species, no_chem_species);
    Jratelist = [Symbol(x) for x in filter(x->typeof(x)==String, df_splists.Jratelist)]

    ## Temperatures
    if "TemperatureArrays" in XLSX.sheetnames(paramlog_wb)
        df_temps = DataFrame(XLSX.readtable("$(folder)PARAMETERS.xlsx", "TemperatureArrays"));
        Tn_arr = reshape(df_temps.Neutrals, GV.n_horiz, GV.num_layers+2)
        Ti_arr = reshape(df_temps.Ions, GV.n_horiz, GV.num_layers+2)
        Te_arr = reshape(df_temps.Electrons, GV.n_horiz, GV.num_layers+2)
    else 
        if quiet==false
            println("WARNING: Reconstructing temperature profiles with default options based on logged control temperatures. It is POSSIBLE the reconstruction could be wrong if the temp function changed.")
        end
        # Need to have some logic for constructing the temperature profile from the logged control parameters, and it's 
        # different by planet. 
        if planet=="Venus"
            T_dict = T_Venus(get_param("TSURF", df_atmcond), get_param("TMESO", df_atmcond), get_param("TEXO", df_atmcond), 
                             "Venus-Inputs/FoxandSung2001_temps_mike.txt"; GV.alt)
        elseif planet=="Mars"
            T_dict = T_Mars(get_param("TSURF", df_atmcond), get_param("TMESO", df_atmcond), get_param("TEXO", df_atmcond); GV.alt)
        end
        Tn_arr = T_dict["neutrals"]
        Ti_arr = T_dict["ions"]
        Te_arr = T_dict["electrons"]
    end

    # Note: Here, all_species is loaded above so it is NOT a global variable. Molmass needs to be passed in though.
    Tplasma_arr = Ti_arr .+ Te_arr;
    Tprof_for_Hs = Dict("neutral"=>Tn_arr, "ion"=>Ti_arr);
    Tprof_for_diffusion = Dict("neutral"=>Tn_arr, "ion"=>Tplasma_arr)
    try
        global Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(alt, sp, Tprof_for_Hs[charge_type(sp)]; M_P, R_P, globvars...) for sp in all_species]); 
    catch UndefVarError
        # hope the user has passed it in
        global Hs_dict = Dict{Symbol, Vector{Float64}}([sp=>scaleH(GV.alt, sp, Tprof_for_Hs[charge_type(sp)]; GV.M_P, GV.R_P, globvars...) for sp in all_species]); 
    end
   
    # Boundary conditions
    speciesbclist_vert = load_bcdict_from_paramdf(df_bcs);
    
    vardict = Dict("DH"=>DH, 
                   "ions_included"=>ions_included,
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
                   "speciesbclist_vert"=>speciesbclist_vert,
                   "rxn_spreadsheet"=>rxn_spreadsheet,
                   "Jratelist"=>Jratelist,
                   "non_bdy_layers"=>non_bdy_layers,
                   "zmin"=>zmin,
                   "zmax"=>zmax,
                   "dz"=>dz,
                   "n_all_layers"=>n_all_layers,
                   "num_layers"=>num_layers,
                   "upper_lower_bdy"=>upper_lower_bdy,
                   "upper_lower_bdy_i"=>upper_lower_bdy_i
                   )
    try
        vardict["alt"] = alt
        vardict["M_P"] = M_P
        vardict["R_P"] = R_P
        vardict["planet"] = planet
    catch UndefVarError
        vardict["alt"] = GV.alt
        vardict["M_P"] = GV.M_P
        vardict["R_P"] = GV.R_P
        vardict["planet"] = GV.planet
    end 
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