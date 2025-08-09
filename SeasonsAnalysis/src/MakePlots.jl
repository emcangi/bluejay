# **************************************************************************** #
#                                                                              #
#                             PLOTTING FUNCTIONS                               #
#                                                                              #
# **************************************************************************** #

#===============================================================================#
#                        Basic model mechanical functions                       #
#===============================================================================#

function make_seasonal_cycle_plots_inclusive(season_folders; mainalt=250, savepath=nothing, resolution="high", extrafn="",
                                                        make_main_cycle_plot=true, make_frac_plot=false, make_flux_plot=false, make_3panel=false,
                                                        make_6panel=false, make_reincorp_cycle_plot=false, make_flux_and_pct_plot=false,
                                                        make_DH_of_escaping=false, plot_pct=false,
                                                        flux_colormap="plasma", show_all_flux=false,
                                                        alt2=200, alt3=250, subplot_lbl_loc=[0, 0.95], subplot_lbl_sz=20, ihoriz=1, globvars...)
    #=
    Does all the work to load files and make plots for temp seasonal cycle. Made because we have two cases, one standard and one with 
    identical crosssections for HDO and H2O.
     Input:
        season_folders: a list of folder paths containing the simulation results, IN SEASON ORDER.
    Optional:
        mainalt: main altitude to plot for
        savepath: place to save the figures
        resolution: typically always high, can be done in low res and select fewer foldres. not needed since I figured out I can save the escape values
                    in a separate XLSX (I'm smart)
        plot_water_profile: if true, plots water vertical profile instead of pr um (typically true)
        spclbl_x, spclbl_loc, subplot_lbl_loc, subplot_lbl_sz: described elsewhere
        alt_legend_loc: location for legend on alttitude axis
        flloc: described elsewhere
        extrafn: more to put in the filename
        plot_pct: described elsewhere
        flux_colormap: colormap to use when plotting flux instead of other stuff
        plot_type: used to specify heatmap when plotting flux
        show_all_flux: if true, will print thermal and nonthermal flux lines.
        make_main_cycle_plot: make plot showing escape, densities, and D/H over time.
        make_frac_plot: same thing but shows fractionation factor.
        make_flux_and_pct: same thing as main_cycle but with an extra axis for percent of total column escape
        make_flux_plot: plot showing vertical flux in atmosphere
        make_3panel: make results plot with altitude as y axis and time as seprate lines.
        make_6panel: plot with D/H ratio in 6 different species
        make_DH_of_escaping: show D/H ratio of escaping atoms over time
    =#
    
    GV = values(globvars)
    required = [:alt, :all_species, :dz, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc, :speciesstyle, :M_P, :R_P, :planet]
    check_requirements(keys(GV), required)
    
    # Get necessary input for plots - temperature
    atm_states, state_files, Tn_all, Ti_all, Te_all = collect_atmospheres(season_folders, resolution, "temp"; globvars...)
    temp_vs_time, change_indices = generate_indvar_vs_time_array(state_files)
    
    # Get the water input also
    atm_states, atm_states_init, state_files, Tn_all, Ti_all, Te_all = collect_atmospheres(season_folders, resolution, "water"; globvars...)


    # println("State files for inclusive run: $(state_files)")

    tempcols_3 = get_colors(3, "plasma"; stp=0.8)

    # 3 panel plot
    if make_3panel
        # more colors
        tempcols_cycle = Dict("hot"=>tempcols_3[3, :], "cold"=>tempcols_3[1, :], "mean"=>tempcols_3[2, :])
        println("Making 3 panel")

        make_3panel_figure([atm_states["lowT"], atm_states["midT2"], atm_states["highT"]], tempcols_cycle, "temp", ["cold", "mean", "hot"]; 
                           savepath=savepath, fn="3panel_tempcycle.png", lloc=(0.02, 0.38), 
                           subplot_lbl_loc=[0, 0.92], subplot_lbl_sz=24, Tn_all, Ti_all, Te_all, GV.alt, GV.speciesstyle)
    end

    # 6 panel plot
    if make_6panel
        println("Making 6 panel")
        DH_6panel([atm_states["lowT"], atm_states["midT2"], atm_states["highT"]], savepath; 
                       fn="DH_profiles_vs_Texo_cycle", lines_mean = [L"\mathrm{T_{exo}=175}"*" K", L"\mathrm{T_{exo}=225}"*" K", L"\mathrm{T_{exo}=275}"*" K"], 
                       tempcols=tempcols_3, subplot_lbl_loc=subplot_lbl_loc, subplot_lbl_sz=subplot_lbl_sz)
    end
    
    # Temperature, D/H, density 
    if make_main_cycle_plot
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making cycle figure")
        flush(stdout)

        # Get the water array
        H2O_profile_array = Array{Float64}(undef, GV.num_layers, length(state_files))
        HDO_profile_array = deepcopy(H2O_profile_array)

        for i in 1:length(state_files)
            nc = get_ncurrent(state_files[i])
            ntot = n_tot(nc, ihoriz; GV.all_species, GV.non_bdy_layers)
            H2O_profile_array[:, i] = nc[:H2O] ./ ntot
            HDO_profile_array[:, i] = nc[:HDO] ./ ntot
        end
        ivar_lbl2 = "Water profile"
        water_to_plot = H2O_profile_array

        # Get the insolation and store some lyman alpha values
        Lyalpha_array = Array{Float64}(undef, length(state_files))
        Lya_values = [1.7616747457786e11, 2.0589378068516e11, 1.7616747457786e11, 1.501199487753e11, 1.7616747457786e11] # FOUND MANUALLY!
        Lya_ind = 0
        for i in 1:length(state_files)
            if i in change_indices
                Lya_ind += 1
            end
            Lyalpha_array[i] = Lya_values[Lya_ind]
        end
        ivar_lbl3 = L"Lyman $\alpha$" * "\n" * L"($\gamma$cm$^2$s$^{-1}$)"

        # Make plot skeleton
        DHax, ffax, nax, fax, fax2, bigfig = seasonal_cycling_figure_skeleton(state_files, temp_vs_time, mainalt, change_indices; 
                                                                          multi_ivar=true, IVAR_array2=water_to_plot, IVAR_array3=Lyalpha_array,
                                                                          ivar_label2=ivar_lbl2, ivar_label3=ivar_lbl3, 
                                                                          plot_pct=plot_pct, subplot_lbl_loc=[-0.3, subplot_lbl_loc[2]], subplot_lbl_sz=subplot_lbl_sz,                                              
                                                                          globvars...)

        # Make the plot
        seasonal_cycling_figure_original(state_files, mainalt, DHax, ffax, nax, fax, fax2; fn="tempcycle_$(resolution)res$(extrafn)", show_all_flux=show_all_flux,
                                         savepath=savepath, plot_pct=plot_pct, subplot_lbl_loc=[-0.3, subplot_lbl_loc[2]], subplot_lbl_sz=subplot_lbl_sz, 
                                         DH_ylims=[2.5e-4, 7e-3], ff_ylims=[0.01, 0.2], DHmult=[2, 5, 10, 20, 34],                                        
                                         globvars...)
        
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end
    
    if make_flux_plot
        println("Now doing flux plot")
        
        
        required = [:q, :molmass, :polarizability, :collision_xsect, # CONSTANTS.jl
                                      :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, # CUSTOMIZATIONS.jl
                                      :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc] # Simulation-unique stuff 
        check_requirements(keys(GV), required)
                                    
        flux_vs_time(state_files, temp_vs_time, change_indices, :H; fn="flux_vs_time_H$(resolution)res$(extrafn)", savepath=savepath, colormap=flux_colormap,
                                                                 num_layers=GV.num_layers, hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, 
                                                                 hot_HD_network=GV.hHDnet, hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, 
                                                                 hot_HD_rc_funcs=GV.hHDrc, globvars...)
        flux_vs_time(state_files, temp_vs_time, change_indices, :D; fn="flux_vs_time_D$(resolution)res$(extrafn)", savepath=savepath, colormap=flux_colormap,
                                                                 num_layers=GV.num_layers, hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, 
                                                                 hot_HD_network=GV.hHDnet, hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, 
                                                                 hot_HD_rc_funcs=GV.hHDrc, globvars...)
        println("finished with flux plot")
    end
    
    # Temperature, fractionation factor
    if make_frac_plot
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making fractionation figure")
        flush(stdout)
        f_vs_time(state_files, temp_vs_time, change_indices; fn="f_vs_time$(resolution)res$(extrafn)", savepath=savepath, globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end 

    # D/H of escaping atoms 
    if make_DH_of_escaping
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making D/H of escaping figure")
        flush(stdout)
        DH_of_escaping_vs_time(state_files, temp_vs_time, change_indices; fn="DH_escaping_vs_time$(resolution)res$(extrafn)", 
                                                                           show_all_flux=show_all_flux, savepath=savepath, globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end
    
    # Temperature, flux, reincorporation, and flux/reincorp
    if make_reincorp_cycle_plot
        required = [:ion_species, :num_layers, :reaction_network]
        check_requirements(keys(GV), required)
        flush(stdout)
        seasonal_cycling_flux_and_recomb(state_files, temp_vs_time, mainalt, change_indices;  fn="tempcycle_reincorp_$(resolution)res$(extrafn)", savepath=savepath,                                                       
            Tn=Tn_all, Ti=Ti_all[2:end-1, 1], Te=Te_all[2:end-1, 1], globvars...)
        flush(stdout)
    end
    
    if make_flux_and_pct_plot
        required = [:ion_species, :num_layers, :reaction_network]
        check_requirements(keys(GV), required)
        flush(stdout)
        seasonal_cycling_flux_and_pct(state_files, temp_vs_time, mainalt, change_indices; #=alt2=200, alt3=250,=# fn="tempcycle_flux_and_pct_$(resolution)res$(extrafn)", savepath=savepath,                                                       
                                      Tn=Tn_all, Ti=Ti_all[2:end-1, 1], Te=Te_all[2:end-1, 1], globvars...)
        flush(stdout)
    end
end

function make_equilibrium_plots_temp(case_folders; savepath=nothing, extrafn="", globvars...)
    #=
    Not used
    =#
    
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
               :q, :molmass, :polarizability, :collision_xsect,
               # From CUSTOMIZATIONS.jl
               :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, :speciesstyle, 
               # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc]
    check_requirements(keys(GV), required)
    
    Tn_all = Array{Float64}(undef, length(GV.alt), 3)
    Ti_all = Array{Float64}(undef, length(GV.alt), 3)
    Te_all = Array{Float64}(undef, length(GV.alt), 3)
    
    atm_keys = Dict(1=>"lowT", 2=>"midT", 3=>"highT")
    atm_states = Dict("lowT"=>Dict(), "midT"=>Dict(), "highT"=>Dict())
    vardicts = Dict("lowT"=>Dict(), "midT"=>Dict(), "highT"=>Dict()) 
    
    for i in 1:length(case_folders)
        thefolder = case_folders[i]
        atm_states[atm_keys[i]] = get_ncurrent(thefolder*"final_atmosphere.h5")

        vardicts[atm_keys[i]] = load_from_paramlog(thefolder; GV.alt)

        # Store the first column (ihoriz=1) for each temperature array for backward compatibility
        Tn_all[:, i] = vardicts[atm_keys[i]]["Tn_arr"][1, :];
        Ti_all[:, i] = vardicts[atm_keys[i]]["Ti_arr"][1, :];
        Te_all[:, i] = vardicts[atm_keys[i]]["Te_arr"][1, :];
    end
    
    esc_df_lowT = final_escape(case_folders[1], "final_atmosphere.h5"; globvars...)
    esc_df_midT = final_escape(case_folders[2], "final_atmosphere.h5"; globvars...)
    esc_df_highT = final_escape(case_folders[3], "final_atmosphere.h5"; globvars...);
    
    DH_of_escaping(esc_df_lowT, case_folders; t=L"T$_{exo}=175$ K" )
    DH_of_escaping(esc_df_midT, case_folders; t=L"T$_{exo}=225$ K" )
    DH_of_escaping(esc_df_midT, case_folders; t=L"T$_{exo}=400$ K" )
    
    make_3panel_figure([atm_states["lowT"], atm_states["midT"], atm_states["highT"]], tempcols_dict, "temp", ["cold", "mean", "hot"]; 
                       savepath=savepath, fn="3panel_equilibrium_temps$(extrafn).png", Tn_all, Ti_all, Te_all, GV.alt, GV.speciesstyle)
    
    DH_6panel([atm_states["lowT"], atm_states["midT"], atm_states["highT"]], savepath; 
                   fn="DH_profiles_vs_Texo$(extrafn)", lines_mean = [L"\mathrm{T_{exo}=175}"*" K", L"\mathrm{T_{exo}=225}"*" K", L"\mathrm{T_{exo}=400}"*" K"], tempcols=tempcols_3)
end

function make_seasonal_cycle_plots_temp(season_folders; mainalt=250, savepath=nothing, resolution="high", extrafn="", 
                                                        make_main_cycle_plot=true, make_frac_plot=false, make_flux_plot=false, make_3panel=false,
                                                        make_6panel=false, make_reincorp_cycle_plot=false, make_flux_and_pct_plot=false, 
                                                        make_DH_of_escaping=false, plot_pct=false, 
                                                        flux_colormap="plasma", show_all_flux=false,
                                                        alt2=200, alt3=250, subplot_lbl_loc=[0, 0.95], subplot_lbl_sz=20, globvars...)
    #=
    Does all the work to load files and make plots for temp seasonal cycle. Made because we have two cases, one standard and one with 
    identical crosssections for HDO and H2O.
     Input:
        season_folders: a list of folder paths containing the simulation results, IN SEASON ORDER.
    Optional:
        mainalt: main altitude to plot for
        savepath: place to save the figures
        resolution: typically always high, can be done in low res and select fewer foldres. not needed since I figured out I can save the escape values
                    in a separate XLSX (I'm smart)
        plot_water_profile: if true, plots water vertical profile instead of pr um (typically true)
        spclbl_x, spclbl_loc, subplot_lbl_loc, subplot_lbl_sz: described elsewhere
        alt_legend_loc: location for legend on alttitude axis
        flloc: described elsewhere
        extrafn: more to put in the filename
        plot_pct: described elsewhere
        flux_colormap: colormap to use when plotting flux instead of other stuff
        plot_type: used to specify heatmap when plotting flux
        show_all_flux: if true, will print thermal and nonthermal flux lines.
        make_main_cycle_plot: make plot showing escape, densities, and D/H over time.
        make_frac_plot: same thing but shows fractionation factor.
        make_flux_and_pct: same thing as main_cycle but with an extra axis for percent of total column escape
        make_flux_plot: plot showing vertical flux in atmosphere
        make_3panel: make results plot with altitude as y axis and time as seprate lines.
        make_6panel: plot with D/H ratio in 6 different species
        make_DH_of_escaping: show D/H ratio of escaping atoms over time
    =#
    
    GV = values(globvars)
    required = [:alt, :all_species, :dz, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc, :speciesstyle, :M_P, :R_P, :planet]
    check_requirements(keys(GV), required)
    
    # Get necessary input for plots
    atm_states, state_files, Tn_all, Ti_all, Te_all = collect_atmospheres(season_folders, resolution, "temp"; globvars...)
    temp_vs_time, change_indices = generate_indvar_vs_time_array(state_files)

    tempcols_3 = get_colors(3, "plasma"; stp=0.8)

    # 3 panel plot
    if make_3panel
        # more colors
        tempcols_cycle = Dict("hot"=>tempcols_3[3, :], "cold"=>tempcols_3[1, :], "mean"=>tempcols_3[2, :])
        println("Making 3 panel")

        make_3panel_figure([atm_states["lowT"], atm_states["midT2"], atm_states["highT"]], tempcols_cycle, "temp", ["cold", "mean", "hot"]; 
                           savepath=savepath, fn="3panel_tempcycle.png", lloc=(0.02, 0.38), 
                           subplot_lbl_loc=[0, 0.92], subplot_lbl_sz=24, Tn_all, Ti_all, Te_all, GV.alt, GV.speciesstyle)
    end

    # 6 panel plot
    if make_6panel
        println("Making 6 panel")
        DH_6panel([atm_states["lowT"], atm_states["midT2"], atm_states["highT"]], savepath; 
                       fn="DH_profiles_vs_Texo_cycle", lines_mean = [L"\mathrm{T_{exo}=175}"*" K", L"\mathrm{T_{exo}=225}"*" K", L"\mathrm{T_{exo}=275}"*" K"], 
                       tempcols=tempcols_3, subplot_lbl_loc=subplot_lbl_loc, subplot_lbl_sz=subplot_lbl_sz)
    end
    
    # Temperature, D/H, density 
    if make_main_cycle_plot
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making cycle figure")
        flush(stdout)
        DHax, ffax, nax, fax, fax2, bigfig = seasonal_cycling_figure_skeleton(state_files, temp_vs_time, mainalt, change_indices; fn="tempcycle_$(resolution)res$(extrafn)", 
                                                                              savepath=savepath, plot_pct=plot_pct, subplot_lbl_loc=[-0.3, subplot_lbl_loc[2]], subplot_lbl_sz=subplot_lbl_sz,                                              
                                                                              globvars...)
        seasonal_cycling_figure_original(state_files, mainalt, DHax, ffax, nax, fax, fax2; fn="tempcycle_$(resolution)res$(extrafn)", show_all_flux=show_all_flux,
                                         savepath=savepath, plot_pct=plot_pct, subplot_lbl_loc=[-0.3, subplot_lbl_loc[2]], subplot_lbl_sz=subplot_lbl_sz, 
                                         DH_ylims=[3e-4, 5e-3], DHmult=[2, 5, 10, 15, 23], ff_ylims=[0.01, 0.2],                                        
                                         globvars...)
        
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end
    
    if make_flux_plot
        println("Now doing flux plot")
        
        
        required = [:q, :molmass, :polarizability, :collision_xsect, # CONSTANTS.jl
                                      :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, # CUSTOMIZATIONS.jl
                                      :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc] # Simulation-unique stuff 
        check_requirements(keys(GV), required)
                                    
        flux_vs_time(state_files, temp_vs_time, change_indices, :H; fn="flux_vs_time_H$(resolution)res$(extrafn)", savepath=savepath, colormap=flux_colormap,
                                                                 num_layers=GV.num_layers, hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, 
                                                                 hot_HD_network=GV.hHDnet, hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, 
                                                                 hot_HD_rc_funcs=GV.hHDrc, globvars...)
        flux_vs_time(state_files, temp_vs_time, change_indices, :D; fn="flux_vs_time_D$(resolution)res$(extrafn)", savepath=savepath, colormap=flux_colormap,
                                                                 num_layers=GV.num_layers, hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, 
                                                                 hot_HD_network=GV.hHDnet, hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, 
                                                                 hot_HD_rc_funcs=GV.hHDrc, globvars...)
        println("finished with flux plot")
    end
    
    # Temperature, fractionation factor
    if make_frac_plot
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making fractionation figure")
        flush(stdout)
        f_vs_time(state_files, temp_vs_time, change_indices; fn="f_vs_time$(resolution)res$(extrafn)", savepath=savepath, globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end 

    # D/H of escaping atoms 
    if make_DH_of_escaping
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making D/H of escaping figure")
        flush(stdout)
        DH_of_escaping_vs_time(state_files, temp_vs_time, change_indices; DHcol="#6449E9", fn="DH_escaping_vs_time$(resolution)res$(extrafn)", 
                                                                           show_all_flux=show_all_flux, savepath=savepath, globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end
    
    # Temperature, flux, reincorporation, and flux/reincorp
    if make_reincorp_cycle_plot
        required = [:ion_species, :num_layers, :reaction_network]
        check_requirements(keys(GV), required)
        flush(stdout)
        seasonal_cycling_flux_and_recomb(state_files, temp_vs_time, mainalt, change_indices;  fn="tempcycle_reincorp_$(resolution)res$(extrafn)", savepath=savepath,                                                       
            Tn=Tn_all, Ti=Ti_all[2:end-1, 1], Te=Te_all[2:end-1, 1], globvars...)
        flush(stdout)
    end
    
    if make_flux_and_pct_plot
        required = [:ion_species, :num_layers, :reaction_network]
        check_requirements(keys(GV), required)
        flush(stdout)
        seasonal_cycling_flux_and_pct(state_files, temp_vs_time, mainalt, change_indices; #=alt2=200, alt3=250,=# fn="tempcycle_flux_and_pct_$(resolution)res$(extrafn)", savepath=savepath,                                                       
                                      Tn=Tn_all, Ti=Ti_all[2:end-1, 1], Te=Te_all[2:end-1, 1], globvars...)
        flush(stdout)
    end
end

function make_equilibrium_plots_water(case_folders; extrafn="", water_where="mesospheric", proflabels=[L"10.38 pr $\mu$m", L"10.5 pr $\mu$m", L"10.94 pr $\mu$m"], 
                                      spclbl_x=[0.8, 0.45], spclbl_loc=[0.75 0.25; 0.35 0.25], globvars...)
    #=,
    Not used.
    =#
    
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
               :q, :molmass, :polarizability, :collision_xsect,
               # From CUSTOMIZATIONS.jl
               :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, :speciesstyle, 
               # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc]
    check_requirements(keys(GV), required)
     
    
    atm_keys = Dict(1=>"low", 2=>"mean", 3=>"high")
    atm_states = Dict("low"=>Dict{Symbol, Vector{Float64}}(), "mean"=>Dict{Symbol, Vector{Float64}}(), "high"=>Dict{Symbol, Vector{Float64}}())
    atm_states_init = Dict("low"=>Dict{Symbol, Vector{Float64}}(), "mean"=>Dict{Symbol, Vector{Float64}}(), "high"=>Dict{Symbol, Vector{Float64}}())
    
    for i in 1:length(case_folders)
        thefolder = case_folders[i]
        atm_states[atm_keys[i]] = get_ncurrent(thefolder*"final_atmosphere.h5");
        atm_states_init[atm_keys[i]] = get_ncurrent(thefolder*"initial_atmosphere.h5");
    end
    
    # get escape
    esc_df_low = final_escape(case_folders[1], "final_atmosphere.h5"; globvars...)
    esc_df_mean = final_escape(case_folders[2], "final_atmosphere.h5"; globvars...)
    esc_df_high = final_escape(case_folders[3], "final_atmosphere.h5"; globvars...)
    
    DH_of_escaping(esc_df_low, case_folders; t="Dry" )
    DH_of_escaping(esc_df_mean, case_folders; t="Mean" )
    DH_of_escaping(esc_df_mean, case_folders; t="Wet" )
    
    # Define some colors
    watercolors = get_colors(3, "Blues"; strt=0.4, stp=1) # starts with lighter one and progresses to darker. 
    watercols_dict = Dict("low"=>watercolors[1, :], "mean"=>watercolors[2, :], "high"=>watercolors[3, :])

    
    # 3 panel plot
    make_3panel_figure([atm_states["low"], atm_states["mean"], atm_states["high"]], watercols_dict, "water", ["low", "mean", "high"]; 
                   initial_atms=[atm_states_init["low"], atm_states_init["mean"], atm_states_init["high"]], 
                   panel1labels=proflabels, spclbl_x=spclbl_x, spclbl_loc=spclbl_loc,
                   esclbl="Atomic", savepath=savepath, fn="3_panel_equilibrium_water$(extrafn).png",
                   titl="How $(water_where) water vapor content affects the density and D/H ratio of atomic H and D", globvars...)
    
    # 6 panel plot
    DH_6panel([atm_states["low"], atm_states["mean"], atm_states["high"]], savepath; fn="DH_profiles_vs_water_eq$(water_where)$(extrafn)", lines_mean=["Dry", "Mean", "Wet"], 
                    tempcols=watercolors)
end

function make_seasonal_cycle_plots_water(season_folders; mainalt=250, savepath=nothing, resolution="high", plot_water_profile=false, spclbl_x=[0.7, 0.45], 
                                         spclbl_loc=[0.75 0.25; 0.35 0.25], subplot_lbl_loc=[0,0.92], subplot_lbl_sz=24, 
                                         alt_legend_loc="upper right", flloc="upper right",  extrafn="", 
                                         plot_pct=false, flux_colormap="plasma", plot_type="heatmap", show_all_flux=false,
                                         # Which plots to make on this call
                                         make_main_cycle_plot=true, make_frac_plot=false, make_flux_and_pct=false, make_flux_plot=false, 
                                         make_3panel=true, make_6panel=true, make_DH_of_escaping=false, globvars...)
    #=
    Does all the work to load files and make plots for water seasonal cycle. Made because we have two cases, one standard and one with 
    identical crosssections for HDO and H2O.

    Input:
        season_folders: a list of folder paths containing the simulation results, IN SEASON ORDER.
    Optional:
        mainalt: main altitude to plot for
        savepath: place to save the figures
        resolution: typically always high, can be done in low res and select fewer foldres. not needed since I figured out I can save the escape values
                    in a separate XLSX (I'm smart)
        plot_water_profile: if true, plots water vertical profile instead of pr um (typically true)
        spclbl_x, spclbl_loc, subplot_lbl_loc, subplot_lbl_sz: described elsewhere
        alt_legend_loc: location for legend on alttitude axis
        flloc: described elsewhere
        extrafn: more to put in the filename
        plot_pct: described elsewhere
        flux_colormap: colormap to use when plotting flux instead of other stuff
        plot_type: used to specify heatmap when plotting flux
        show_all_flux: if true, will print thermal and nonthermal flux lines.
        make_main_cycle_plot: make plot showing escape, densities, and D/H over time.
        make_frac_plot: same thing but shows fractionation factor.
        make_flux_and_pct: same thing as main_cycle but with an extra axis for percent of total column escape
        make_flux_plot: plot showing vertical flux in atmosphere
        make_3panel: make results plot with altitude as y axis and time as seprate lines.
        make_6panel: plot with D/H ratio in 6 different species
        make_DH_of_escaping: show D/H ratio of escaping atoms over time
    =#
    
    GV = values(globvars)
    required = [:all_species, :dz, :num_layers, :non_bdy_layers, :n_alt_index, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc, 
                :molmass, :M_P, :R_P, :planet]
    check_requirements(keys(GV), required)
    
    # Get necessary inputs for plots --------------------------------------------------------------------- #
    atm_states, atm_states_init, state_files, Tn_all, Ti_all, Te_all = collect_atmospheres(season_folders, resolution, "water"; globvars...)
    water_prum, change_indices = generate_indvar_vs_time_array(state_files; val_order=[10.5, 12.73, 10.5, 10.38, 10.5])

    # Define some colors
    watercolors = get_colors(3, "Blues"; strt=0.4, stp=1) # starts with lighter one and progresses to darker. 
    watercols_dict = Dict("low"=>watercolors[1, :], "mean"=>watercolors[2, :], "high"=>watercolors[3, :])

    # 3 panel plot
    if make_3panel
        make_3panel_figure([atm_states["low"], atm_states["mean2"], atm_states["high"]], watercols_dict, "water", ["low", "mean", "high"]; 
                       initial_atms=[atm_states_init["low"], atm_states_init["mean2"], atm_states_init["high"]], subplot_lbl_loc=subplot_lbl_loc,
                       panel1labels=["Dry", "Mean", "Wet"], spclbl_x=spclbl_x, spclbl_loc=spclbl_loc, lloc=(0.05, 0.38), figsz=(20, 7), 
                       esclbl="Atomic", savepath=savepath, fn="3_panel_water_cycle$(extrafn)", globvars...)
    end
    
    # 6 panel plot -------------------------------------------------------------------------------- #
    if make_6panel
        DH_6panel([atm_states["low"], atm_states["mean"], atm_states["high"]], savepath; fn="DH_profiles_vs_water_cycle_mesosphere$(extrafn)", 
                       lloc=[0.35, 0.4],
                       lines_mean=["Dry", "Mean", "Wet"], subplot_lbl_loc=subplot_lbl_loc, tempcols=watercolors, globvars...)
    end

    if plot_water_profile==true
        # Fill in the water profiles into an array for plotting
        H2O_profile_array = Array{Float64}(undef, GV.num_layers, length(state_files))
        HDO_profile_array = deepcopy(H2O_profile_array)

        for i in 1:length(state_files)
            nc = get_ncurrent(state_files[i])
            ntot = n_tot(nc, ihoriz; GV.all_species, GV.non_bdy_layers)
            H2O_profile_array[:, i] = nc[:H2O] ./ ntot
            HDO_profile_array[:, i] = nc[:HDO] ./ ntot
        end
        ivar_lbl = "Water profile"
        water_to_plot = H2O_profile_array
    else
        ivar_lbl = "Water\n"*L"(H$_2$O pr. $\mu$m)"
        water_to_plot = water_prum
    end

    # Main cycle plot -------------------------------------------------------------------------------- #
    if make_main_cycle_plot
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making cycle figure")
        flush(stdout)

        DHax, ffax, nax, fax, fax2, bigfig = seasonal_cycling_figure_skeleton(state_files, water_to_plot, mainalt, change_indices;  
                                                                              fn="watercycle_$(resolution)res$(extrafn)", savepath=savepath, plot_pct=plot_pct, 
                                                                              subplot_lbl_loc=[-0.3, subplot_lbl_loc[2]], subplot_lbl_sz=subplot_lbl_sz,                                              
                                                                              globvars...)

        seasonal_cycling_figure_original(state_files, mainalt, DHax, ffax, nax, fax, fax2; thefig=bigfig, fn="watercycle_$(resolution)res$(extrafn)", show_all_flux=show_all_flux,
                                        ivar_label=ivar_lbl, savepath=savepath, flloc=flloc, alt_legend_loc=alt_legend_loc, subplot_lbl_loc=[-0.3, 0.92], 
                                        plot_pct, pct_ylims=[8e-7, 1e-4], DH_ylims=[3e-4, 5e-3], ff_ylims=[0.01, 0.2], DHmult=[2, 5, 10, 15, 23], globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end

    # Vertical transport plot -------------------------------------------------------------------------------- #
    if make_flux_plot
        println("Now doing flux plot")
        
        required = [:q, :molmass, :polarizability, :collision_xsect, # CONSTANTS.jl
                                      :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, # CUSTOMIZATIONS.jl
                                      :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc] # Simulation-unique stuff 
        check_requirements(keys(GV), required)
                                    
        flux_vs_time(state_files, water_to_plot, change_indices, :H; fn="flux_vs_time_H$(resolution)res$(extrafn)", savepath=savepath, 
                     colormap=flux_colormap, plot_type=plot_type, 
                     num_layers=GV.num_layers, ivar_label=ivar_lbl, hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, 
                     hot_HD_network=GV.hHDnet, hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, 
                     hot_HD_rc_funcs=GV.hHDrc, globvars...)
        flux_vs_time(state_files, water_to_plot, change_indices, :D; fn="flux_vs_time_D$(resolution)res$(extrafn)", savepath=savepath, 
                     colormap=flux_colormap, plot_type=plot_type, 
                     num_layers=GV.num_layers, ivar_label=ivar_lbl, hot_H_network=GV.hHnet, hot_D_network=GV.hDnet, hot_H2_network=GV.hH2net, 
                     hot_HD_network=GV.hHDnet, hot_H_rc_funcs=GV.hHrc, hot_D_rc_funcs=GV.hDrc, hot_H2_rc_funcs=GV.hH2rc, 
                     hot_HD_rc_funcs=GV.hHDrc, globvars...)
        println("finished with flux plot")
    end
    
    # Fractionation factor plot -------------------------------------------------------------------------------- #
    if make_frac_plot
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making fractionation figure")
        flush(stdout)
        f_vs_time(state_files, water_to_plot, change_indices; fn="f_vs_time_water_$(resolution)res$(extrafn)", subplot_lbl_loc=[-0.3, 0.92], 
                                                              ivar_label=ivar_lbl, savepath=savepath, subplot_lbl_sz=subplot_lbl_sz,
                                                              globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end

    # D/H of escaping atoms 
    if make_DH_of_escaping
        println("$(Dates.format(now(), "(HH:MM:SS)")) Start making D/H escaping figure")
        flush(stdout)
        DH_of_escaping_vs_time(state_files, water_to_plot, change_indices; fn="DH_escaping_vs_time$(resolution)res$(extrafn)", 
                                                                           show_all_flux=show_all_flux, savepath=savepath, globvars...)
        println("$(Dates.format(now(), "(HH:MM:SS)")) Done")
        flush(stdout)
    end
    
    # DO NOT USE
    if make_flux_and_pct
        required = [:ion_species, :num_layers, :reaction_network]
        check_requirements(keys(GV), required)
        flush(stdout)
        seasonal_cycling_flux_and_pct(state_files, water_to_plot, mainalt, change_indices; fn="watercycle_fluxes_$(resolution)res$(extrafn)", typeflag="water",
                                         savepath=savepath, ivar_label=ivar_lbl, Tn=Tn_all[2:end-1], Te=Te_all[2:end-1], Ti=Ti_all[2:end-1], globvars...)
        flush(stdout)
    end  
end

function make_insolation_plots(case_folders, fluxfiles; savepath=nothing, extrafn="", subplot_lbl_loc=[0.01, 0.92], subplot_lbl_sz=24, globvars...)
    #=
    Make a few plots for insolation, but not many since it didn't produce interesting results

    Inputs:
        case_folders: folders containing model results for 3 cases. must be in order of low, med, high insolation
    Optional:
        idk they're all defined elsewhere in this file and I'm tired.
    =#
    
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
               :q, :molmass, :polarizability, :collision_xsect,
               # From CUSTOMIZATIONS.jl
               :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, :speciesstyle, 
               # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc]
    check_requirements(keys(GV), required)
    
    atm_keys = Dict(1=>"equinox1", 2=>"perihelion", 3=>"equinox2", 4=>"aphelion", 5=>"equinox3")
    atm_states = Dict("equinox1"=>Dict{Symbol, Vector{Float64}}(), "equinox2"=>Dict{Symbol, Vector{Float64}}(), "equinox3"=>Dict{Symbol, Vector{Float64}}(),
                      "perihelion"=>Dict{Symbol, Vector{Float64}}(), "aphelion"=>Dict{Symbol, Vector{Float64}}())
    
    for i in 1:length(case_folders)
        thefolder = case_folders[i]
        atm_states[atm_keys[i]] = get_ncurrent(thefolder*"final_atmosphere.h5");
    end
    
    # 3 panel plot
    solcols = get_colors(3, "inferno"; stp=0.8)
    solcols_dict = Dict("ap"=>solcols[1, :], "eq"=>solcols[2, :], "peri"=>solcols[3, :])

    ap_flux_file, eq_flux_file, peri_flux_file = fluxfiles

    make_3panel_figure([atm_states["aphelion"], atm_states["equinox2"], atm_states["perihelion"]], solcols_dict, "insolation", 
                       ["ap", "eq", "peri"]; fluxfiles=[ap_flux_file, eq_flux_file, peri_flux_file], 
                        shy=false, fn="3_panel_insolation_cycle", savepath=savepath, panel1labels=["Aphelion", "Mean distance", "Perihelion"],
                        line_lbl_x = 0.64,
                        subplot_lbl_loc, subplot_lbl_sz, globvars...)
    
    # 6 panel plot
    DH_6panel([atm_states["aphelion"], atm_states["equinox2"], atm_states["perihelion"]], savepath; fn="DH_profiles_vs_insolation", 
                     lines_mean=["Aphelion", "Mean distance", "Perihelion"], lloc=[0.4, 0.7], 
                    tempcols=solcols, subplot_lbl_loc, subplot_lbl_sz, globvars...)
end

function plot_limiting_flux(atm_states, atm_keys, Tn_all; savepath=nothing, simple=false, all_carriers=true, fnextra="", globvars...)
    #=
    Plots the limiting flux 

    This is not a reasonable plot to use because limiting flux doesn't make sense to evalutae at all altitudes. Ask Roger Yelle why
    =#
    GV = values(globvars)
    required = [:all_species, :alt, :molmass, :n_alt_index, :non_bdy_layers, :plot_grid, :zmax]
    check_requirements(keys(GV), required)

    if simple
        therange = 2:2
    else
        therange = 1:3
    end
        
    # Limiting flux
    if all_carriers
        
        limflux_H = [limiting_flux(:H, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     limiting_flux(:H2, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     2 * limiting_flux(:H2O, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     limiting_flux(:HDO, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     limiting_flux(:OH, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) for i in therange]
        limflux_D = [limiting_flux(:D, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     limiting_flux(:HD, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     limiting_flux(:HDO, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) +
                     limiting_flux(:OD, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) for i in therange]
    else
        limflux_H = [limiting_flux(:H, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) for i in therange]
        limflux_D = [limiting_flux(:D, atm_states[atm_keys[i]], Tn_all[2:end-1, i]; ihoriz=i, globvars...) for i in therange]
    end

    # Typical escape fluxes for these files:
    max_H_flux = 1e9#effusion_velocity(275, 1; zmax) * atm_states["highT"][:H][end]
    min_H_flux = effusion_velocity(175, 1; GV.zmax) * atm_states["lowT"][:H][end]
    max_D_flux = 3e4#effusion_velocity(275, 2; zmax) * atm_states["highT"][:D][end]
    min_D_flux = 1e4#effusion_velocity(175, 2; zmax) * atm_states["lowT"][:D][end]

    set_rc_params(fs=18, axlab=20, xtls=18, ytls=18, sansserif=sansserif_choice, monospace=monospace_choice)

    fig, ax = subplots(figsize=(8,6))
    plot_bg(ax)

    descs = ["cold", "mean", "hot"]
    lbls = [L"\mathrm{T_{exo}=175 K}", L"\mathrm{T_{exo}=225 K}", L"\mathrm{T_{exo}=400 K}"]
    if simple
        ax.plot(limflux_H[1], GV.plot_grid, color=tempcols_dict[descs[2]], label=lbls[2], zorder=3)
        ax.plot(limflux_D[1], GV.plot_grid, color=tempcols_dict[descs[2]], linestyle="--", zorder=3)
    else
        for i in therange
            ax.plot(limflux_H[i], GV.plot_grid, color=tempcols_dict[descs[i]], label=lbls[i], zorder=3)
            ax.plot(limflux_D[i], GV.plot_grid, color=tempcols_dict[descs[i]], linestyle="--", zorder=3)
        end
    end

    ax.text(0.35, 0.5, L"$\phi_{\ell,D}$", transform=ax.transAxes, color="#262626", fontsize=18)
    ax.text(0.75, 0.5, L"$\phi_{\ell,H}$", transform=ax.transAxes, color="#262626", fontsize=18)


    ax.fill_betweenx([0, 250], min_H_flux, x2=max_H_flux, zorder=2, color="xkcd:dark gray", alpha=0.5)
    ax.fill_betweenx([0, 250], min_D_flux, x2=max_D_flux, zorder=2, color="xkcd:dark gray", alpha=0.5)

    ax.text(0.65, 0.85, "Typical range\nof H escape", transform=ax.transAxes)#, fontsize=14)
    ax.text(0.1, 0.85, "Typical range\nof D escape", transform=ax.transAxes)#, fontsize=14)

    ax.set_xlabel(L"Limiting upward flux (cm$^{-2}$s$^{-1}$)")
    ax.set_ylabel("Altitude (km)")

    ax.set_xscale("log")
    ax.set_xlim(6e3, 1e12)
    ax.set_xticks([i for i in logrange(1e4, 1e12, 5)])
    majorticks = ax.get_xticks()
    minor_locator = matplotlib.ticker.FixedLocator([i for i in logrange(1e5, 1e11, 4)])
    ax.xaxis.set_minor_locator(minor_locator)
    ax.tick_params(axis="x", which="minor",  labelbottom=false) # 

    ax.legend(fontsize=14, loc=(0.67, 0.04))

    savefig(savepath*"limiting_flux$(fnextra)_400Khigh.png", dpi=300, bbox_inches="tight")
    show()
end


# Make figures ==============================================================================


function make_3panel_figure(atms, colors, exptype, atm_state_order; fn="3panel", savepath=nothing, shy=true, initial_atms=nothing, 
                            figsz=(20, 5),  specialH2=false, panel1labels=["Low", "Mean", "High"], line_lbl_x=0.01,
                            spclbl_x=[0.8, 0.45], spclbl_loc=[0.75 0.25; 0.35 0.25], lloc=(0, 0.4),
                            fluxfiles=nothing,
                            subplot_lbl_loc=[0.05, 0.9], subplot_lbl_sz=20, plotfmt="pdf",  globvars...)
    #=
    This makes the 3 panel figure which shows the inputs, the D/H vs. altitude, and the densities of H and D.
    
    Inputs:
        atms: list of atmospheric state dictionaries
        colors: dictionary where keys describe exobase state and entry is a 3-item array with RGB values
        exptype: string, experiment type; "temp" "water" or "insolation"
        atm_state_order: ordering of the entries in atms, helps label things correctly.

    Optional inputs:
        fn: filename
        savepath: path to save to 
        shy: share y. i guess at one point I wanted it set to false???
        initial_atms: optional initial atmospheric dictionaries, useful when plotting the water  modeling results.
        figsz: figure size
        specialH2: whether to show H2 in the density panel instead of H (yes if true)
        panel1labels: line labels to print in the first panel showing the independent variable
        line_lbl_x: x-coordinate to print the line labels in panel 1 
        spclbl_x: x-coordinate to print the species names H2O and HDO in the water panel
        spclbl_loc: location to print the labels for H and D in the density panel
        lloc: location for legend
        subplot_lbl_loc: where to print the subplot panel labels
        subplot_lbl_sz: font size for said labels
    =#
    
    GV = values(globvars)
    # required = [:zmax]
    # check_requirements(keys(GV), required)
    
    set_rc_params(fs=18, axlab=20, xtls=18, ytls=18, sansserif=sansserif_choice, monospace=monospace_choice)
    
    fig, ax = subplots(1, 3, sharey=shy, figsize=figsz)
    subplots_adjust(wspace=0.1)
    
    for i in 1:length(ax)
        ax[i].text(subplot_lbl_loc..., "$(Unicode.julia_chartransform('a') + (i-1))"*")", fontsize=subplot_lbl_sz, transform=ax[i].transAxes)
    end
    
    # PANEL 1
    colors_3 = Array{Float64}(undef, 3, 3) 
    if exptype=="temp"
        keys = ["cold", "mean", "hot"]
    elseif exptype=="water"
        keys = ["low", "mean", "high"]
    elseif exptype=="insolation"
        keys = ["ap", "eq", "peri"]
    end
    colors_3[1, :] .= colors[keys[1]]
    colors_3[2, :] .= colors[keys[2]]
    colors_3[3, :] .= colors[keys[3]]
    if exptype == "temp"
        make_temperature_panel(ax[1], GV.Tn_all, GV.Ti_all, GV.Te_all; GV.alt, tempcols=colors_3)
    elseif exptype=="water"
        make_water_panel(ax[1], initial_atms, colors_3, panel1labels; spclbl_x=spclbl_x, line_lbl_x=line_lbl_x, globvars...)
    elseif exptype=="insolation"
        make_insolation_panel(ax[1], colors_3, panel1labels, fluxfiles; line_lbl_x, GV.speciesstyle)
    end
    
    # PANEL 2
    if shy==false
        ax[2].set_ylabel("Altitude (km)")
        subplots_adjust(wspace=0.2)
    end

    DH_alt_profile_singlepanels(ax[2], atms, colors, atm_state_order;
                            heavysp=:D, lightsp=:H, wl=[0.29, 0.9], ol=[0.6, 0.5], lloc=lloc,
                            xlims=[1e-4, 1e-2], legendon=true, globvars...)
    # PANEL 3
    make_density_panel(ax[3], atms, colors, atm_state_order; specialH2, spclbl_loc=spclbl_loc, globvars...)

    savefig("$(savepath)$(fn).$(plotfmt)", format=plotfmt, dpi=300, bbox_inches="tight")
    show()
end

function DH_6panel(atmdict_list, savepath; lines_mean=["Solar minimum", "Solar mean", "Solar maximum"],  
                        lloc=[0.3, 0.7], fn="DH_profiles", cutoff=[1 1 1; 1 1 25], 
                        subplot_lbl_loc=[0, 0.95], subplot_lbl_sz=20, globvars...)
    #= 
    Plots the D/H ratio in 6 different species vs. altitude for several atmospheres and saves it.
    
    Inputs:
        atmdict_list: a list of several dictionary objects containing atmospheric states.
        savepath: Where to save the plot

    Optional:
        lines_mean: line labels 
        lloc: location to draw the legend for the lines 
        fn: Filename 
        cutoff: altitude index below which not to plot. Useful because the HCO+ D/H ratios 
                go crazy low down because of really small numbers dividing other small numbers. 
                The shape is the shape of the final grid of plots.
        subplot_lbl_loc: where to draw the subplot panel labels, in transform coords
        subplot_lbl_sz: font size for these labels

    Output:
        A 6-panel plot of D/H ratio in different species.
    =#
    
    GV = values(globvars)
    required =  [:tempcols, :zmax, :plot_grid]
    check_requirements(keys(GV), required)
    
    set_rc_params(fs=14, axlab=16, xtls=14, ytls=14, sansserif=sansserif_choice, monospace=monospace_choice)

    # get D/H in water
    DH_toplot = Array{Array}(undef, 2, 3)
    
    # Fill in the stuff
    DH_toplot[1, 1] = hcat([n[:HDO] ./ (2 .* n[:H2O]) for n in atmdict_list]...)
    DH_toplot[1, 2] = hcat([n[:D] ./ (n[:H]) for n in atmdict_list]...)
    DH_toplot[1, 3] = hcat([n[:HD] ./ (2 .* n[:H2]) for n in atmdict_list]...)
    DH_toplot[2, 1] = hcat([n[:OD] ./ (n[:OH]) for n in atmdict_list]...)
    DH_toplot[2, 2] = hcat([n[:H2DOpl] ./ (3 .* n[:H3Opl]) for n in atmdict_list]...)
    DH_toplot[2, 3] = hcat([n[:DCOpl] ./ (n[:HCOpl]) for n in atmdict_list]...)
    
    #heavysp
    heavysp = [:HDO :D :HD; :OD :H2DOpl :DCOpl]
    lightsp = [:H2O :H :H2; :OH :H3Opl :HCOpl]
    sppair = ["Water" "Atomics" "Molecular"; "Hydroxl" "Hydronium" "Formyl cation"]
    
    # Make a scale to show D/H in terms of earth value
    smow = 1.6e-4
    DHtoSMOW(x) = x ./ smow
    SMOWtoDH(x) =  x .* smow
    
    # plot setup
    r = 2
    c = 3
    fig, ax = subplots(2, 3, figsize=(10, 10))
    subplots_adjust(wspace=0.08, hspace=0.23)
    ax[2, 2].set_xlabel("log(D/H ratio)")
    ax[1, 1].set_ylabel("Altitude (km)")
    ax[2, 1].set_ylabel("Altitude (km)")
    
    lbl_i = 1
    for a in ax
        plot_bg(a)
        a.text(subplot_lbl_loc..., "$(Unicode.julia_chartransform('a') + (lbl_i-1))"*")", fontsize=subplot_lbl_sz, transform=a.transAxes)
        lbl_i+=1
        a.set_xscale("log")
        a.set_xlim(1e-4, 1e-2)
        a.set_xticklabels(["", "-4", "-3", "-2"])
        a.set_yticks(ticks=collect(0:50:GV.zmax/1e5))
        a.tick_params(axis="x", which="both", bottom=true, top=true, labelbottom=true, labeltop=false)
        multiplier = [1, 2, 5, 10, 25, 50]
        draw_DH_lines(a, multiplier; DHfs=12, which="vertical")

    end

    # plot actual content
    for i in 1:r # rows 
        for j in 1:c # columns
            if j > 1
                ax[i, j].tick_params(axis="y", labelleft=false)
            end
            for sc in 1:length(atmdict_list) 
                plotme = DH_toplot[i, j][:, sc]
                if (i==2) & (j==3)
                    plotme[1:cutoff[i, j]] .= 0
                end
                
                ax[i, j].plot(plotme, GV.plot_grid, zorder=10, color=GV.tempcols[sc, :])
                
                # Labels 
                heavysp_str = string_to_latexstr(string(heavysp[i, j]); dollarsigns=false)
                lightsp_str = string_to_latexstr(string(lightsp[i, j]); dollarsigns=false)
                
                basic_lightspstr = string(lightsp[i, j]) # Used for determining the factor to attach to the H species

                if occursin("H2", basic_lightspstr)
                    Hfactor_light = 2
                    lbl = L"%$(sppair[i, j]): $\frac{[%$heavysp_str]}{%$Hfactor_light[%$lightsp_str]}$"
                elseif occursin("H3", basic_lightspstr)
                    Hfactor_light = 3
                    lbl = L"%$(sppair[i, j]): $\frac{[%$heavysp_str]}{%$Hfactor_light[%$lightsp_str]}$"
                else
                    Hfactor_light = 1
                    lbl = L"%$(sppair[i, j]): $\frac{[%$heavysp_str]}{[%$lightsp_str]}$"
                end
                
                ax[i, j].set_title(lbl, size=18)#, va="top", size=14)
            end
        end
    end
    
    # LEgend
    L2D = PyPlot.matplotlib.lines.Line2D
    lines = [L2D([0], [0], color=GV.tempcols[i, :]) for i in 1:length(lines_mean)]  # linewidth=lw[i]
    legend(lines, [lines_mean...], fontsize=12, bbox_to_anchor=lloc)
    
    if savepath != nothing
        savefig(savepath*"/$(fn).png", bbox_inches="tight", dpi=300)
    end
end

function DH_of_escaping_vs_time(thefiles, IVAR_array, change_indices; DHcol="#6449E9", show_all_flux=false, lloc=nothing, fn="DH_escaping_vs_time", ivar_label="Temp. (K)", 
                     savepath=nothing, colorbar_tix=collect(logrange(1e-10,1e-2,5)), subplot_lbl_loc=[-0.3,0.92], subplot_lbl_sz=24,
                     ivar_cmap="RdBu", globvars...)

    #=
    Plots D/H of escaping atoms as a function of time

    thefiles: list of final atmosphere files to use for each datapoint in the plot.
    IVAR_array: array of values to plot for temperature cycling
    change_indices: indices of IVAR_array at which the value changes (i.e., start of seasons)

    Optional:
    fn: Name of image file to save
    ivar_label: y-axis label for the top panel of the plot, showing the independent variable. Temp. (K) or Water
    savepath: path in which to save the image
    colorbar_tix: ticks to assign to the colorbar for the independent variable plot (if water).
    subplot_lbl_loc: axis-based coordinates in which to place subfigure labels like "a)", etc
    subplot_lbl_sz: fontsize to use for said labels.
    ivar_cmap: colormap to use for the colorbar, must be valid in matplotlib.
    =#
    
    GV = values(globvars)

    total_num_files = length(thefiles)
    t = 1:total_num_files

    # Collect escape info...
    thefolders = parent_folders_from_full_path(thefiles)
    all_esc_df = get_escape_for_all_atmospheres(thefolders, thefiles; globvars...)

    
    # calculate D/H (total)
    DH_t = all_esc_df."Dt" ./ all_esc_df."Ht"
    DH_n = all_esc_df."Dn" ./ all_esc_df."Hn"
    DH_tn = all_esc_df."Dtn" ./ all_esc_df."Htn" 

    # plot
    set_rc_params(fs=18, axlab=20, xtls=18, ytls=18, sansserif=sansserif_choice, monospace=monospace_choice)
    
    gridspec = matplotlib.gridspec
    fig = figure(figsize=(8, 6), constrained_layout=true)
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 2])
    
    IVax = subplot(gs.new_subplotspec((0, 0)))
    DHax = subplot(gs.new_subplotspec((1, 0)))
    setp(IVax.get_xticklabels(), visible=false)
    
    subplots_adjust(hspace=0.08)

    # Independent variable axis
    stages = set_up_independent_var_ax(IVax, thefiles, t, IVAR_array, ivar_label, ivar_cmap, colorbar_tix; thefig=fig, GV.alt)

    # D/H axis
    plot_bg(DHax)
    DHax.plot(t, DH_tn, color=DHcol)

    if show_all_flux
        DHax.plot(t, DH_t, color=DHcol, linestyle="--")
        DHax.plot(t, DH_n, color=DHcol, linestyle=":")

        L2D = PyPlot.matplotlib.lines.Line2D
        lines = [L2D([0], [0], color=DHcol), 
                 L2D([0], [0], color=DHcol, linestyle="--"),
                 L2D([0], [0], color=DHcol, linestyle=":")]  
        DHax.legend(lines, ["total", "thermal", "nonthermal"], fontsize=14, loc=lloc)
    end
    
    DHax.set_yscale("log")
    
    # D/H lines
    multiplier = [0.01, 0.1, 1, 5, 10]
    draw_DH_lines(DHax, multiplier)

    axlist = [IVax, DHax]
    for i in 1:length(axlist)
        axlist[i].xaxis.set_major_locator(matplotlib.ticker.FixedLocator(change_indices))
        axlist[i].set_xlim(change_indices[1]-2, total_num_files+2)
        axlist[i].text(subplot_lbl_loc..., "$(Unicode.julia_chartransform('a') + (i-1))"*")", fontsize=subplot_lbl_sz, transform=axlist[i].transAxes)
    end
    format_time_axis(DHax, change_indices, length(thefiles), stages)
    
    DHax.set_ylabel("D/H of escaping atoms")

    # SAVE IT
    if savepath != nothing
        savefig("$(savepath)$(fn).png", bbox_inches="tight", dpi=300)
    end
    
    show()
end

function f_vs_time(thefiles, IVAR_array, change_indices; fn="f_vs_time", ivar_label="Temp. (K)", 
                     savepath=nothing, colorbar_tix=collect(logrange(1e-10,1e-2,5)), subplot_lbl_loc=[0,0.92], subplot_lbl_sz=24,
                     ivar_cmap="RdBu", globvars...)
    #=
    Plots the fractionation factor as a function of time. 

    thefiles: list of final atmosphere files to use for each datapoint in the plot.
    IVAR_array: array of values to plot for temperature cycling
    change_indices: indices of IVAR_array at which the value changes (i.e., start of seasons)

    Optional:
    fn: Name of image file to save
    ivar_label: y-axis label for the top panel of the plot, showing the independent variable. Temp. (K) or Water
    savepath: path in which to save the image
    colorbar_tix: ticks to assign to the colorbar for the independent variable plot (if water).
    subplot_lbl_loc: axis-based coordinates in which to place subfigure labels like "a)", etc
    subplot_lbl_sz: fontsize to use for said labels.
    ivar_cmap: colormap to use for the colorbar, must be valid in matplotlib.
    =#
      
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
               :q, :molmass, :polarizability, :collision_xsect, :DH,
               # From CUSTOMIZATIONS.jl
               :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, :speciesstyle, 
               # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc]
    check_requirements(keys(GV), required)
    
    # Setup figure
    set_rc_params(fs=18, axlab=20, xtls=18, ytls=18, sansserif=sansserif_choice, monospace=monospace_choice)
    
    gridspec = matplotlib.gridspec
    fig = figure(figsize=(8, 6), constrained_layout=true)
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 2])
    
    IVax = subplot(gs.new_subplotspec((0, 0)))
    fracax = subplot(gs.new_subplotspec((1, 0)))
    setp(IVax.get_xticklabels(), visible=false)    
    subplots_adjust(hspace=0.08)

    # Ind Var ax --------------------------------------------------------------------------------------------------------
    total_num_files = length(thefiles)
    t = 1:total_num_files
    stages = set_up_independent_var_ax(IVax, thefiles, t, IVAR_array, ivar_label, ivar_cmap, colorbar_tix; GV.alt)
    
    # Frac axis -----------------------------------------------------------------------------------------------------
    println("$(Dates.format(now(), "(HH:MM:SS)")) Working on frac axis")
    
    # Tease apart folder and filenames because it's needed for getting flux
    folder_pattern = r".+v\d\/"
    file_pattern = r"/[a-z]+_.+\.h5"
    thefolders = [String(match(folder_pattern, f).match) for f in thefiles]
    filenames = [String(match(file_pattern, f).match) for f in thefiles]

    # Get escape and calculate f
    all_esc_df = get_escape_for_all_atmospheres(thefolders, thefiles; globvars...)
    
    f_all = (all_esc_df."Dtn" ./ all_esc_df."Htn") ./ DH
    f_t = (all_esc_df."Dt" ./ all_esc_df."Ht") ./ DH
    f_n = (all_esc_df."Dn" ./ all_esc_df."Hn") ./ DH

    # Plot it
    fcol = "xkcd:merlot"
    plot_bg(fracax)
    fracax.plot(t, f_all, color=fcol)
    fracax.plot(t, f_t, color=fcol, linestyle="--")
    fracax.plot(t, f_n, color=fcol, linestyle=":")
    fracax.set_yscale("log")
    
    L2D = PyPlot.matplotlib.lines.Line2D
    lines = [L2D([0], [0], color=fcol), 
             L2D([0], [0], color=fcol, linestyle="--"),
             L2D([0], [0], color=fcol, linestyle=":")]  
    fracax.legend(lines, ["total", "thermal", "nonthermal"], fontsize=16, loc=(0.55,0.2))
    
    println("$(Dates.format(now(), "(HH:MM:SS)")) Fracfac axis completed")

    # -------------------------------------------------------------------------------------------------------------

    # Set up the x ticks
    axlist = [IVax, fracax]
    for i in 1:length(axlist)
        axlist[i].xaxis.set_major_locator(matplotlib.ticker.FixedLocator(change_indices))
        axlist[i].set_xlim(change_indices[1]-2, total_num_files+2)
        axlist[i].text(subplot_lbl_loc..., "$(Unicode.julia_chartransform('a') + (i-1))"*")", fontsize=subplot_lbl_sz, transform=axlist[i].transAxes)
    end
    format_time_axis(fracax, change_indices, length(thefiles), stages)
    
    fracax.set_ylabel("Fractionation factor")

    # SAVE IT
    if savepath != nothing
        savefig("$(savepath)$(fn).png", bbox_inches="tight", dpi=300)
    end
        
    show()
end

function flux_vs_time(thefiles, IVAR_array, change_indices, sp; fn="flux_vs_time", ivar_label="Temp. (K)", plot_type="line", flux_alt=80,
                     savepath=nothing,  subplot_lbl_loc=[0,0.92], subplot_lbl_sz=24, logthresh=[5, 2], logstep=1, linscale=1, 
                     colormap="plasma", colorbar_tix=collect(logrange(1e-10,1e-2,5)), globvars...)
    #=
    Plots the net production/loss by altitude due to transport as a function of time. 

    thefiles: list of final atmosphere files to use for each datapoint in the plot.
    IVAR_array: array of values to plot for temperature cycling
    change_indices: indices of IVAR_array at which the value changes (i.e., start of seasons)
    sp: species for which to calculate the net production/loss due to transport.

    Optional:
    fn: Name of image file to save
    ivar_label: y-axis label for the top panel of the plot, showing the independent variable. Temp. (K) or Water
    savepath: path in which to save the image
    colorbar_tix: ticks to assign to the colorbar for the independent variable plot (if water).
    subplot_lbl_loc: axis-based coordinates in which to place subfigure labels like "a)", etc
    subplot_lbl_sz: fontsize to use for said labels.
    colormap: colormap to use for the colorbar, must be valid in matplotlib.
    =#
      
    GV = values(globvars)
    required = [:alt, :all_species, :dz, :num_layers]
    check_requirements(keys(GV), required)
    
    stages=["Spring", "Summer", "Autumn", "Winter", "Spring"]

    # Setup figure
    set_rc_params(fs=18, axlab=20, xtls=18, ytls=18, sansserif=sansserif_choice, monospace=monospace_choice)

    gridspec = matplotlib.gridspec
    fig = figure(figsize=(8, 6), constrained_layout=true)
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 1.5])
    
    IVax = subplot(gs.new_subplotspec((0, 0)))
    fluxax = subplot(gs.new_subplotspec((1, 0)))
    setp(IVax.get_xticklabels(), visible=false)
    
    subplots_adjust(hspace=0.08)

    
    # Ind Var ax --------------------------------------------------------------------------------------------------------
    stages = set_up_independent_var_ax(IVax, thefiles, t, IVAR_array, ivar_label, ivar_cmap, colorbar_tix; GV.alt)

    # Flux axis -----------------------------------------------------------------------------------------------------
    println("$(Dates.format(now(), "(HH:MM:SS)")) Working on flux axis")
    flush(stdout)
    
    # Tease apart folder and filenames because it's needed for getting flux
    folder_pattern = r".+v\d\/"
    file_pattern = r"/[a-z]+_.+\.h5"
    thefolders = [String(match(folder_pattern, f).match) for f in thefiles]
    filenames = [String(match(file_pattern, f).match) for f in thefiles]
    
    # Make a zero array
    flux_heatmap = zeros(GV.num_layers, length(filenames))
    
    c = 1 # column counter
    for (fold, fl) in zip(thefolders, filenames)
        # Open atmosphere
        thisatm = get_ncurrent(fold*fl)
        
        # Get things we need to know
        vardict = load_from_paramlog(fold; GV.alt)
        
        # Get jratdict
        Jratelist = format_Jrates(fold*"active_rxns.xlsx", GV.all_species, "Jratelist"; hot_atoms=true, ions_on=true)[1];
        Jratedict = Dict([j=>thisatm[j] for j in Jratelist])
        
        # Get flux
        n_horiz = length(thisatm[GV.all_species[1]])
        flux_heatmap[:, c] = get_directional_fluxes(sp, thisatm, n_horiz; nonthermal=true, Jratedict, speciesbclist=vardict["speciesbclist"],
                                                         Tprof_for_Hs=vardict["Tprof_for_Hs"], Tprof_for_diffusion=vardict["Tprof_for_diffusion"],
                                                         Hs_dict = vardict["Hs_dict"], Tn = vardict["Tn_arr"], Ti = vardict["Ti_arr"],
                                                         Te = vardict["Te_arr"], Tp = vardict["Tplasma_arr"],
                                                         transport_species=vardict["transport_species"], globvars...)[2:end-1]
        c += 1
    end
   
    # Get the largest absolute value so we know what range ot use for cbar
    max_abs = maximum(abs.(flux_heatmap))

    maxlog=Int(ceil(log10(maximum(flux_heatmap))))
    minlog=-Int(ceil(log10(abs(minimum(flux_heatmap)))))
    
    # PLOT
    plot_bg(fluxax)

    if plot_type=="heatmap"
        if sp == :H
            logthresh=logthresh[1]
        elseif sp==:D
            logthresh=logthresh[2]
        end
        img_flux = fluxax.imshow(flux_heatmap, cmap=colormap, extent=(t[1], t[end], alt[1]/1e5, alt[end]/1e5),
                                  interpolation="nearest", origin="lower", aspect="auto", zorder=10, 
                                  norm=matplotlib.colors.SymLogNorm(10.0^(logthresh), linscale=linscale))# matplotlib.colors.CenteredNorm())# 
        fluxax.set_ylabel("Altitude (km)")

        # Colorbar
        fig.subplots_adjust(right=0.9)
        
        fluxax_bbox = fluxax.get_position().get_points()
        xmin = fluxax_bbox[1, 1]
        xmax = fluxax_bbox[2, 1]
        ymin = fluxax_bbox[1, 2]
        ymax = fluxax_bbox[2, 2]
        cbar_ax = fig.add_axes([xmax + 0.02, ymin, (xmax-xmin)/10, ymax-ymin]) # left, bottom, width, height

        cbar = fig.colorbar(img_flux,  cax=cbar_ax, label="Directional flux "*L"(cm$^{-2}$s$^{-1}$)")
        cbar.ax.tick_params(labelsize=12)
    elseif plot_type=="line"
        selected_fluxes = flux_heatmap[GV.n_alt_index[flux_alt*1e5], :]
        fluxax.plot(t, selected_fluxes)
        fluxax.set_ylabel("Flux at $(flux_alt) km (cm^-2s^-1)")
    end
    
    println("$(Dates.format(now(), "(HH:MM:SS)")) Flux axis completed")

    # -------------------------------------------------------------------------------------------------------------

    # Set up the x ticks
    axlist = [IVax, fluxax]
    for i in 1:length(axlist)
        axlist[i].xaxis.set_major_locator(matplotlib.ticker.FixedLocator(change_indices))
        axlist[i].set_xlim(change_indices[1]-2, length(filenames)+2)
        axlist[i].text(subplot_lbl_loc..., "$(Unicode.julia_chartransform('a') + (i-1))"*")", fontsize=subplot_lbl_sz, transform=axlist[i].transAxes)
    end

    format_time_axis(fluxax, change_indices, length(thefiles), stages)

    # SAVE IT
    if savepath != nothing
        savefig("$(savepath)$(fn).png", bbox_inches="tight", dpi=300)
    end
        
    show()
end

function seasonal_cycling_figure_original(thefiles, A, DHax, ffax, nax, fax, fax2;  fn="vstime_multialt", #ivar_label=L"T$_{exo}$ (K)", 
                                 savepath=nothing, W_array=nothing, alt2=nothing, alt3=nothing, show_all_flux=false, plot_pct=false,
                                 # ylimits 
                                 Hflux_ylims=[1e7, 1.2e9],  Dflux_ylims=[1e1, 5e4], pct_ylims=[1e-6, 1e-4], DH_ylims=[1e-4, 1e-2], ff_ylims=[0.01, 0.1], 
                                 # ticks
                                 Hflux_ticks=[1e7, 1e8, 1e9], Dflux_ticks=[1e1, 1e2, 1e3, 1e4, 1e5], density_ticks=[1e3, 1e4, 1e5, 1e6, 1e7], colorbar_tix=collect(logrange(1e-10,1e-2,5)),
                                 # DH guidelines
                                 DHmult=[2, 5, 10, 20, 25], 
                                 # Figure stuff
                                 alt_legend_loc=(0.72, 0.4), flloc="lower left", subplot_lbl_loc=[0,0.95], subplot_lbl_sz=20,
                                 ivar_cmap="RdBu", plotfmt="pdf", globvars...)
    #=
    
    Inputs:
        thefiles: files containing final atmospheric state to use for each datapoint in the plot
        A: an altitude to choose
        DHax, nax, fax, fax2: Axes objects returned by seasonal_cycling_figure_skeleton 
    =#
      
    GV = values(globvars)
    required = [:alt, :all_species, :dz]
    check_requirements(keys(GV), required)

    # Define colors
    faxcol = "#FF6F59" # "#F26430"# #e66101"
    densitycol = "#9C0D38" # 
    DHcol = "#6449E9"# 

    altcols = get_colors(3, "Purples"; strt=0.4, stp=1)

    # twin axis colors
    pctcol = "#357DED" # "xkcd:cornflower"
    fraccol = "#009B72"
    
    # Needed for flux axis
    total_num_files = length(thefiles)
    t = 1:total_num_files
    
    # Flux axes (H and D) --------------------------------------------------------------------------------------------
    # println("$(Dates.format(now(), "(HH:MM:SS)")) Working on flux axes")
    flush(stdout)
    
    total_H_flux, total_D_flux = esc_flux_vs_time(fax, fax2, t, thefiles, [faxcol, faxcol]; show_all_flux=show_all_flux, lloc=flloc, globvars...)

    if Hflux_ylims!=nothing
        fax.set_ylim(Hflux_ylims[1], Hflux_ylims[2])
        fax.tick_params(axis="y", which="minor", labelleft=false, labelcolor="white")
    end
    if Dflux_ylims!=nothing
        fax2.set_ylim(Dflux_ylims[1], Dflux_ylims[2])
        fax2.tick_params(axis="y", which="minor", labelleft=false, labelcolor="white")
    end

    # Set ticks
    fax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(Hflux_ticks))
    fax2.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(Dflux_ticks))

    # hopefully log minor ticks?:
    y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], numticks = 10)
    fax2.yaxis.set_minor_locator(y_minor)
    fax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    
    # Add percent which escapes for H ------------------------------------------------------------
    if plot_pct==true
        total_H_col = [sum(get_ncurrent(f)[:H]) .* GV.dz for f in thefiles]
        # twin ax for ratio of total col
        fax_pct = fax.twinx()
        turn_off_borders(fax_pct)
        
        fax_pct.set_ylabel("% of total column",  color=pctcol)
        fax_pct.tick_params(axis="y", labelcolor=pctcol)        
        Hpct = esc_pct_of_total_col(total_H_flux, total_H_col)
        fax_pct.plot(t, Hpct, color=pctcol) # H
        fax_pct.set_yscale("log")
        fax_pct.set_ylim(pct_ylims...)
    end
    
    # Add percent which escapes for H -------------------------------------------------------------
    if plot_pct==true
        total_D_col = [sum(get_ncurrent(f)[:D]) .* GV.dz for f in thefiles]
        # TWIN AX for ratio
        fax2_pct = fax2.twinx()
        turn_off_borders(fax2_pct)
        fax2_pct.set_ylabel("% of total column", color=pctcol)
        fax2_pct.tick_params(axis="y", labelcolor=pctcol)
        Dpct = esc_pct_of_total_col(total_D_flux, total_D_col)
        fax2_pct.plot(t, Dpct, color=pctcol) # D
        fax2_pct.set_yscale("log")
        fax2_pct.set_ylim(pct_ylims...) # same ylimits so they can be compared
    end

    # println("$(Dates.format(now(), "(HH:MM:SS)")) Flux axis completed")

    # Density axis --------------------------------------------------------------------------------------------------
    # println("$(Dates.format(now(), "(HH:MM:SS)")) Working on density axis")
    density_vs_time(nax, thefiles, densitycol, A; alt2=alt2, alt3=alt3, altcols=altcols, density_ticks=density_ticks, globvars...)
    # println("$(Dates.format(now(), "(HH:MM:SS)")) Density axis completed")

    # LEGEND for density and D/H
    L2D = PyPlot.matplotlib.lines.Line2D
    if (alt2!=nothing) & (alt3!=nothing)
        L=3
        lbls = ["Alt. $(Int64(A)) km", "Alt. $(Int64(alt2)) km", "Alt. $(Int64(alt3)) km"]
        lines = [L2D([0], [0], color=altcols[i, :]) for i in 1:L]  
        nax.legend(lines, lbls, fontsize=14, loc=alt_legend_loc)
    end
    
    # Fractionation factor - reviewer request  ----------------------------------------------------------
    # println("$(Dates.format(now(), "(HH:MM:SS)")) Working on ff axis")
    f_vs_time_panel(ffax, thefiles; fraccol, globvars...)
    ffax.set_ylim(ff_ylims...) 
    # println("$(Dates.format(now(), "(HH:MM:SS)")) ff axis completed")

    # D/H axis -----------------------------------------------------------------------------------------------------
    # println("$(Dates.format(now(), "(HH:MM:SS)")) Working on D/H axis")
    thiscol = alt2 == nothing ? DHcol : altcols[1, :]
    DH_vs_time(DHax, thefiles, A, [thiscol]; lw=2, DHmult, globvars...)
    
    if alt2 != nothing
        DH_vs_time(DHax, thefiles, alt2, [altcols[2, :]]; lw=2, DHmult)
    end
    if alt3 != nothing
        DH_vs_time(DHax, thefiles, alt3, [altcols[3, :]];  lw=2, DHmult)
    end

    DHax.set_ylim(DH_ylims...)
    
    # println("$(Dates.format(now(), "(HH:MM:SS)")) D/H axis completed")

    # ------------------------------------------------------------------------------------------------------------

    # SAVE IT
    if savepath != nothing
        savefig("$(savepath)$(fn).$(plotfmt)", format=plotfmt, bbox_inches="tight", dpi=300)
    end
        
    show()
end

function seasonal_cycling_figure_skeleton(thefiles, IVAR_array, A, change_indices; 
                                          # Additional independent variables 
                                          multi_ivar=false, ivar_cmap="RdBu", IVAR_array2=nothing, IVAR_array3=nothing,
                                          ivar_label=L"T$_{exo}$ (K)", ivar_label2="", ivar_label3="", 
                                          # Additional altitudes at which to plot
                                          alt2=nothing, alt3=nothing, 
                                          # Y limits
                                          pct_ylims=[9e-7, 1e-4], DH_ylims=[4e-4,5e-3],  ff_ylims=[0.01, 0.1], 
                                          # Ticks
                                          Hflux_ticks=[1e7, 1e8, 1e9], Dflux_ticks=[1e1, 1e2, 1e3, 1e4], density_ticks=[1e3, 1e4, 1e5, 1e6, 1e7], 
                                          colorbar_tix=collect(logrange(1e-10,1e-2,5)), # ticks for water abundance colorbar
                                          # Plot stylings
                                          flloc="upper center", plot_pct=false, subplot_lbl_loc=[0,0.95], subplot_lbl_sz=20, 
                                          globvars...)
    #=
    Inputs:
        thefiles: list of files to plot results for
        IVAR_array: array of values to plot for independent variable (temp or water)
        A: an altitude to choose
        change_indices: indices within thefiles that signify a change to a new simulation (change in independent variable)

    Optional:
        fn: filename extension to sae
        ivar_label: label for the first panel Y-axis
        alt2, alt3: other optional alts for which to plot densities and D/H
        pct_ylims: y-limits for the percent of total column escape axis
        DH_ylims: y-limits or D/H axis
        Hflux_ticks: ticks to use for the H flux y axis
        Dflux_ticks: likewise for D
        density_ticks: likewise for density axis
        colorbar_tix: ticks for the independent variable colorbar if plotting water
        flloc: flux legend location
        plot_pct: whether to make the axis showing percent of total colkumn that escapes
        subplot_lbl_loc: location of labels for subplot panels
        subplot_lbl_sz: font for previous
        ivar_cmap: Colormap to use for water colorbar
    =#
      
    # Setup figure -----------------------------
    set_rc_params(fs=16, axlab=18, xtls=16, ytls=16, sansserif=sansserif_choice, monospace=monospace_choice)

    gridspec = matplotlib.gridspec

    if multi_ivar == false
        gs_rows = 6
        subplot_rel_sizes = [1.25, 2.5, 2.5, 2.5, 1.5, 3]
        num_ivars = 1
    else  # 3 independent variables
        gs_rows = 8
        subplot_rel_sizes = [1, 1, 1, 2.5, 2.5, 2.5, 1.5, 3]
        num_ivars = 3
    end

    figheight = gs_rows*3 + 1

    fig = figure(figsize=(8, figheight), constrained_layout=true)
    gs = gridspec.GridSpec(gs_rows, 1, figure=fig, height_ratios=subplot_rel_sizes)
    
    for row in 1:gs_rows  # has to be done this way to get the indexing right, more julia/python confusion
        axi = fig.add_subplot(gs[row, 0])
    end 

    ax = gcf().axes

    if multi_ivar == false
        IVax = ax[1]
        fax =  ax[2]
        fax2 = ax[3]
        nax =  ax[4]
        ffax = ax[5]
        DHax = ax[6]
    else  # 3 independent variables
        IVax = ax[1]
        IVax2 = ax[2]
        IVax3 = ax[3]
        fax =  ax[4]
        fax2 = ax[5]
        nax =  ax[6]
        ffax = ax[7]
        DHax = ax[8]
    end
    
    for a in ax[1:end-1]
        setp(a.get_xticklabels(), visible=false)
    end
    
    # setp(IVax.get_xticklabels(), visible=false)
    subplots_adjust(hspace=0.1)

    # Needed to set xlim bounds correctly and create an independent variable array against which to plot 
    total_num_files = length(thefiles)
    t = 1:total_num_files
    
    # Ind Var ax --------------------------------------------------------------------------------------------------------
    stages = set_up_independent_var_ax(IVax, thefiles, t, IVAR_array, ivar_label, ivar_cmap, colorbar_tix; thefig=fig, globvars...)

    if multi_ivar==true 
        stages = set_up_independent_var_ax(IVax2, thefiles, t, IVAR_array2, ivar_label2, ivar_cmap, colorbar_tix; thefig=fig, globvars...)
        stages = set_up_independent_var_ax(IVax3, thefiles, t, IVAR_array3, ivar_label3, ivar_cmap, colorbar_tix; thefig=fig, globvars...)
    end 
    
    # Flux axes (H and D) -------------------------------------------------------------------------------------------- #
    
    # H axis configure ------------------------------------------------------------ #
    plot_bg(fax)
    fax.set_ylabel("H escape\n"*L"(cm$^{-2}$s$^{-1}$)")
    fax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(Hflux_ticks))
    fax.set_yscale("log")


    # D axis configure ------------------------------------------------------------- #
    plot_bg(fax2)
    fax2.set_ylabel("D escape\n"*L"(cm$^{-2}$s$^{-1}$)")#, rotation=0, labelpad=50)
    fax2.set_yscale("log")
    
    # Density axis --------------------------------------------------------------------------------------------------# 
   
    plot_bg(nax)
    nax.set_yscale("log")
    nax.set_ylabel(L"Density (cm$^{-3}$)")

    # Frac fac axis -----------------------------------------------------------------------------------------------------#
    plot_bg(ffax)
    ffax.set_ylabel("Fractionation\nfactor")
    ffax.set_yscale("log")
    ffax.set_ylim(ff_ylims...)

    # D/H axis -----------------------------------------------------------------------------------------------------#
    plot_bg(DHax)
    DHax.set_ylabel("Atomic D/H ratio")
    DHax.set_yscale("log")
    DHax.set_ylim(DH_ylims...)


    # ------------------------------------------------------------------------------------------------------------- #

    # Set up the x ticks
    for i in 1:length(ax)
        ax[i].xaxis.set_major_locator(matplotlib.ticker.FixedLocator(change_indices))
        ax[i].set_xlim(change_indices[1]-2, total_num_files+2)
        ax[i].text(subplot_lbl_loc..., "$(Unicode.julia_chartransform('a') + (i-1))"*")", fontsize=subplot_lbl_sz, transform=ax[i].transAxes)
    end

    format_time_axis(DHax, change_indices, length(thefiles), stages)
    
    show()
        
    return DHax, ffax, nax, fax, fax2, fig
end

# Individual panels ========================================================================

# By altitude ------------------------------------------------------------------------------

function make_temperature_panel(ax, Tn_all, Ti_all, Te_all; globvars...)
    #=
    Fills in ax with a plot showing the temperature profiles:

    Inputs:
        ax: existing axis object to plot onto.
        Tn_all, Ti_all, Te_all: neutrals, ions, electrons
    =#
    GV = values(globvars)
    
    nicegray = "#262626"
    
    # AX 1: TEMPERATURE ===================================================
    ioncol = "xkcd:tea"
    ecol = "xkcd:greyish blue"
    plot_bg(ax)
    
    # Find last isothermal altitude
    i_isotherm = findlast(Tn_all[:, 1] .== Ti_all .== Te_all)[1]

    # Find which column to use for each
    cold_col_i = argmin(Tn_all[end, :])
    hot_col_i = argmax(Tn_all[end, :])
    col_opts_for_mid = setdiff(1:size(Tn_all)[1], [cold_col_i, hot_col_i])
    col_mid_i = col_opts_for_mid[1] # just choose the first since Tn is the same for all moderate T profs

    # plot
    ax.plot(Tn_all[:, cold_col_i], GV.alt./1e5, label="Neutrals, solar min", color=GV.tempcols[1, :])
    ax.plot(Tn_all[:, col_mid_i], GV.alt./1e5, label="Neutrals, solar mean", color=GV.tempcols[2, :])
    ax.plot(Tn_all[:, hot_col_i], GV.alt./1e5, label="Neutrals, solar max", color=GV.tempcols[3, :])

    ax.plot(Ti_all[:, 1], GV.alt./1e5, label="Ions", color=ioncol)
    ax.plot(Te_all[:, 1], GV.alt./1e5, label="Electrons", color=ecol)
    ax.plot(Tn_all[1:i_isotherm, 1], GV.alt[1:i_isotherm] ./ 1e5, color=nicegray)

    ax.set_xscale("log")
    

    # Label the control temps
    # surface
    ax.scatter(Tn_all[:, 1][1], 0, marker="o", color=nicegray, zorder=10)
    ax.text(Tn_all[:, 1][1]+10, 0,"$(Int64(round(Tn_all[:, 1][1], digits=0))) K ", color=nicegray)

    # mesosphere
    meso_ind = findfirst(x->x==minimum(Tn_all[:, 1]), Tn_all[:, 1])
    ax.scatter(Tn_all[:, 1][meso_ind], GV.alt[meso_ind+5]/1e5, marker="o", color=nicegray, zorder=10)
    ax.text(Tn_all[:, 1][meso_ind]+5, GV.alt[meso_ind+5]/1e5, "$(Int64(round(Tn_all[:, 1][meso_ind], digits=0))) K ", color=nicegray)

    # exosphere
    Tn_lbl_y = 255
    ax.scatter(Tn_all[end, cold_col_i], 250, marker="o", zorder=10, color=GV.tempcols[1, :])
    ax.text(Tn_all[end, cold_col_i], Tn_lbl_y, "$(Int64(round(Tn_all[:, cold_col_i][end], digits=0))) K ", ha="right", fontsize=16, color=GV.tempcols[1, :])
    ax.scatter(Tn_all[end, col_mid_i], 250, marker="o", zorder=10, color=GV.tempcols[2, :])
    ax.text(Tn_all[end, col_mid_i], Tn_lbl_y, "$(Int64(round(Tn_all[:, col_mid_i][end], digits=0))) K", ha="center", fontsize=16, color=GV.tempcols[2, :])
    ax.scatter(Tn_all[end, hot_col_i], 250, marker="o", zorder=10, color=GV.tempcols[3, :])
    ax.text(Tn_all[end, hot_col_i], Tn_lbl_y, "$(Int64(round(Tn_all[:, hot_col_i][end], digits=0))) K ", fontsize=16, color=GV.tempcols[3, :])

    # final labels
    ax.set_ylabel("Altitude (km)")
    ax.set_yticks(collect(0:50:Int64(GV.alt[end]/1e5)))
    ax.set_xlabel("Temperature (K)")
    ax.set_xlim(95, 2e3)
    ax.tick_params(which="both", axis="x", top=true)
    
    ax.text(Tn_all[end, 1]-5, 270, "Neutrals")
    ax.text(Ti_all[end]+10, 230, "Ions", color=ioncol)
    ax.text(Ti_all[end]+10, 150, "Electrons", color=ecol)
end

function make_water_panel(ax, atms, colors, txtlbls; spclbl_x = [0.8, 0.45], line_lbl_x=0.01, globvars...)
    #=
    Plots water mixing ratio vs. altitude. Includes a conversion axis for ppm.

    Inputs
        ax: axis object on which to plot
        atms: atmospheric state dictionaries to use in plotting
        colors: colors to use for the three profiles
        txtlbls: strings to print and label the lines
    Optional inputs:
        spclbl_x: optional x-coordinate for printing the species labels (h2O and HDO)
        line_lbl_x: optional x-coordinate to start printing the labels for the lines.

    Outputs: n/a
    =#

    GV = values(globvars)

    plot_bg(ax)
    for w in 1:length(atms)
        ax.plot(atms[w][:H2O][1:GV.upper_lower_bdy_i]./n_tot(atms[w], ihoriz; GV.all_species, GV.alt)[1:GV.upper_lower_bdy_i], GV.plot_grid[1:GV.upper_lower_bdy_i],
                color=colors[w, :])
        ax.plot(atms[w][:HDO][1:GV.upper_lower_bdy_i]./n_tot(atms[w], ihoriz; GV.all_species, GV.alt)[1:GV.upper_lower_bdy_i], GV.plot_grid[1:GV.upper_lower_bdy_i],
                color=colors[w, :], linestyle=GV.speciesstyle[:HDO])
    end
    ax.set_xscale("log")
    ax.set_xlabel(L"Water mixing ratio (cm$^{-3}$/cm$^{-3}$)")
    ax.set_xlim(1e-13, 1e-3)

    # make a conversion ax for PPM ------------------ 
    ppmax = ax.twiny()
    subplots_adjust(bottom=0.4)
    ppmax.xaxis.set_ticks_position("bottom")
    ppmax.xaxis.set_label_position("bottom")
    ppmax.spines["bottom"].set_position(("axes", -0.25))
    turn_off_borders(ppmax)
    ppmax.spines["bottom"].set_visible(true)
    ppmax.set_xticks(ax.get_xticks() ./ 1e-6)
    ppmax.set_xlim(1e-13/1e-6, 1e-3/1e-6)
    ppmax.set_xscale("log")
    ppmax.set_xlabel("Water mixing ratio (ppm)")

    # Print labels for lines
    for t in 1:length(txtlbls)
        ax.text(line_lbl_x, 0.75-0.08*t, txtlbls[t], color=colors[t, :], transform=ax.transAxes)
    end
    
    ax.text(1, 0.95, "Water abundances above $(GV.plot_grid[GV.upper_lower_bdy_i]) km\nare solved for, not prescribed.", ha="right", va="top", color="black",
            transform=ax.transAxes)
    
    ax.text(spclbl_x[1], 0.1, L"H$_2$O", color="black", transform=ax.transAxes)
    ax.text(spclbl_x[2], 0.1, "HDO", color="black", transform=ax.transAxes)
    ax.set_ylabel("Altitude (km)")
end

function make_insolation_panel(ax, colors, txtlbls, photonfluxfiles; line_lbl_x=0.7, globvars...)
    #=
    Plots solar flux by wavelength, similar to make_temperature_panel

    Inputs:
        ax: axis object on which to plot
        colors: colors to use for the three lines
        txtlbls: labels to apply to each line.
    =#

    GV = values(globvars)
    plot_bg(ax)

    ap_flux_file, eq_flux_file, peri_flux_file = photonfluxfiles
    
    ap_flux = readdlm(ap_flux_file, '\t', Float64, comments=true, comment_char='#')[1:2000,:]
    eq_flux = readdlm(eq_flux_file, '\t', Float64, comments=true, comment_char='#')[1:2000,:]
    peri_flux = readdlm(peri_flux_file, '\t', Float64, comments=true, comment_char='#')[1:2000,:]
    
    ax.plot(ap_flux[1:2000, 1], ap_flux[1:2000, 2], color=colors[1, :])
    ax.plot(eq_flux[1:2000, 1], eq_flux[1:2000, 2], color=colors[2, :])
    ax.plot(peri_flux[1:2000, 1], peri_flux[1:2000, 2], color=colors[3, :])
    
    ax.set_yscale("log")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Insolation\n"*L"($\gamma$ s$^{-1}$ cm$^{-2}$ nm$^{-1}$)")
    ax.set_xlim(-1, 500)
    
    for t in 1:length(txtlbls)
        ax.text(line_lbl_x, 0.7-0.1*t, txtlbls[t], color=colors[t, :], transform=ax.transAxes)
    end
end

function DH_alt_profile_singlepanels(ax, atmdict_list, linecols, atm_state_order; lines_mean=nothing,
                        heavysp=:D, lightsp=:H, species_pair="atomics", wl=[0.1, 0.9], ol=[0.6, 0.8],
                        titl="", xlims=[1e-4, 1e-2], cutoff=nothing,fn="DH_profiles", lloc=(0.05,0.4),
                        opttext=[], opttext_loc=[()], opttext_col=[], legendon=false, globvars...)
    #=
    Also plots D/H by altitude, but does so on the ax object for part of a larger plot.
    
    =#

    GV = values(globvars)
    required = [:zmax, :plot_grid]
    check_requirements(keys(GV), required)

    
    # get D/H in water
    DH_water = [n[:HDO] ./ (2 .* n[:H2O]) for n in atmdict_list]
    
    if occursin("H2", string(lightsp))
        Hfactor_light = 2
    elseif occursin("H3", string(lightsp))
        Hfactor_light = 3
    else
        Hfactor_light = 1
    end
    
    DH_other = [n[heavysp] ./ (Hfactor_light * n[lightsp]) for n in atmdict_list]

    # plot setup
    plot_bg(ax)
    ax.set_xlabel("D/H ratio")
    
    ax.set_xscale("log")
    ax.set_xlim(xlims[1], xlims[2])
    ax.set_yticks(ticks=collect(0:50:GV.zmax/1e5))
    ax.tick_params(axis="x", which="both", bottom=true, top=true, labelbottom=true, labeltop=false)
    ax.set_title(titl)
   
    # plot actual content   
    for i in 1:length(atmdict_list)
        ax.plot(DH_water[i], GV.plot_grid, zorder=10, linewidth=2, color=linecols[atm_state_order[i]], linestyle=":")
        
        if cutoff != nothing
            DH_other[i][1:cutoff] .= 0
        end
        ax.plot(DH_other[i], GV.plot_grid, zorder=10, color=linecols[atm_state_order[i]], linewidth=2, linestyle="-")#ls[i])
    end
    
    # Show SMOW lines for orientation
    multiplier = [1, 2, 5, 10, 15, 25, 50]
    draw_DH_lines(ax, multiplier; which="vertical")

    # Optional text
    for t in 1:length(opttext)
        ax.text(opttext_loc[t]..., opttext[t], color=opttext_col[t], transform=ax.transAxes, fontsize=14)
    end
    # Text about temperatures
    heavysp_str = "[$(string_to_latexstr(string(heavysp); dollarsigns=false))]"
    lightsp_str = "[$(string_to_latexstr(string(lightsp); dollarsigns=false))]"
    if Hfactor_light > 1
        lightsp_str = "$(Hfactor_light)"*lightsp_str
    end
    # Legend
    if legendon==true
        L2D = PyPlot.matplotlib.lines.Line2D
        lines = [L2D([0], [0], color="black", linewidth=2, linestyle="-"),
                 L2D([0], [0], color="black", linewidth=2, linestyle=":")]
        if species_pair=="atomics"
            printme = L"R$_{atomic}$"
        end
        ax.legend(lines, ["$(printme)", L"R$_{water}$"], fontsize=16, loc=lloc)
    end
end

function make_density_panel(ax, atmdict_list, colors, atm_state_order; spclbl_loc=[0.75 0.25; 0.35 0.25], specialH2=false, globvars...)
    #=
    Plots densities of H and D vs. altitude from the atmospheres in atmdict_list.

    Inputs:
        ax: axis object on which to plot. Should be part of an externally-created figure.
        atmdict_list: atmospheric state dictionaries to plot for (should be 3)
        colors: colors to use for each line (3)
        atm_state_order: tells what order the dictionaries appear in so things plot correctly

    Optional:
        spclbl_loc: optional location for species labels (H and D) on density panel
        specialH2: whethe to plot H2 instead of H (troubleshooting reasons)
    =#
    GV = values(globvars)
    required = [:plot_grid]
    check_requirements(keys(GV), required)
    plot_bg(ax)
    
    for a in 1:length(atmdict_list)
        if specialH2
            ax.plot(atmdict_list[a][:H2], GV.plot_grid, color=colors[atm_state_order[a]])
        else
            ax.plot(atmdict_list[a][:H], GV.plot_grid, color=colors[atm_state_order[a]])
        end
        ax.plot(atmdict_list[a][:D], GV.plot_grid, color=colors[atm_state_order[a]], linestyle=GV.speciesstyle[:D])
    end

    ax.set_xscale("log")
    ax.set_xlabel(L"Final density (cm$^{-3}$)")
    
    ax.text(spclbl_loc[2, :]..., "D", transform=ax.transAxes)
    if specialH2
        ax.text(spclbl_loc[1, :]..., "H2", transform=ax.transAxes)
        ax.set_xlim(1e1, 1e14)
    else 
        ax.text(spclbl_loc[1, :]..., "H", transform=ax.transAxes)
        ax.set_xlim(1e1, 1e10)
    end
end

# Vs time ------------------------------------------------------------------------------------

function set_up_independent_var_ax(the_ax, thefiles, t, the_arr, ivar_label, ivar_cmap, colorbar_tix; thefig=nothing, globvars...)
    #=
    Make the independent variable axis for the time figures. 
    =#
    GV = values(globvars)
    tempcol = "#2A2D34"#
    watercol = "xkcd:true blue"
    lya_col = "#2A2D34" # 
    plot_bg(the_ax)
    stages = ["Spring", "Summer", "Autumn", "Winter", "Spring"]

    # Setup figure
    set_rc_params(fs=16, axlab=18, xtls=16, ytls=16, sansserif=sansserif_choice, monospace=monospace_choice)

    # Needed to set up x axis correctly.
    total_num_files = length(thefiles)
    t = 1:total_num_files

    # Set the y label, will be overwritten if needed
    the_ax.set_ylabel(ivar_label)
    
    if occursin("(K)", ivar_label)
        the_ax.plot(t, the_arr, color=tempcol)
        the_ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([minimum(the_arr), mode(the_arr), maximum(the_arr)]))
        # the_ax.set_ylabel(ivar_label)
    elseif occursin("Water", ivar_label)
        if !isa(the_arr, Vector)
            img = the_ax.imshow(the_arr, cmap=ivar_cmap, extent=(t[1], t[end], GV.alt[1]/1e5, GV.alt[end]/1e5),
                              interpolation="nearest", origin="lower", aspect="auto", norm=matplotlib.colors.LogNorm(vmin=1e-10, vmax=1e-2), zorder=10)
            the_ax.set_yticks([0, 50, 100, 150, 200, 250])
            the_ax.set_ylabel("Altitude (km)")
            # the_ax.tick_params(which="both", axis="both", labelsize=16)

            thefig.subplots_adjust(right=0.9)
            cbar_ax = thefig.add_axes(get_cbar_size(the_ax)) 
            cbar = thefig.colorbar(img, cax=cbar_ax, label="Water MR", ticks=colorbar_tix)
            cbar.ax.tick_params(labelsize=12)
        else 
            the_ax.plot(t, the_arr, color=watercol)
            the_ax.set_ylim(9, 13)
            the_ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([9, 10, 11, 12, 13]))
            # the_ax.set_ylabel(ivar_label)
        end
        stages = [" ", "Perihelion", " ", "Aphelion", " "]
    elseif occursin("Ly", ivar_label)
        the_ax.plot(t, the_arr, color=lya_col)
        the_ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([minimum(the_arr), mode(the_arr), maximum(the_arr)]))
        stages = [" ", "Perihelion", " ", "Aphelion", " "]
    end

    return stages
end

function esc_flux_vs_time(fax, fax2, t, thefiles, cols; lloc="lower left", show_all_flux=false, globvars...)
    #=
    Plots escape fluxes over time.

    Inputs:
        fax: First escape axis object
        fax2: Second escape axis object
        t: x-axis array of times
        thefiles: files to plot from
        cols: colors to use, should be the same color currently

    Optional:
        lloc: legend location, string or tuple
        show_all_flux: if set to false, will not print the thermal and non-thermal lines, just total
    =#
    
    GV = values(globvars)
    required = [# Things from CONSTANTS.jl
                :q, :molmass, :polarizability, :collision_xsect,
                # From CUSTOMIZATIONS.jl
                :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, :speciesstyle, 
                # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc, :ion_species, :Tn, :Ti, :Te]
    check_requirements(keys(GV), required)
    
    # Collect escape information
    thefolders = parent_folders_from_full_path(thefiles)
    all_esc_df = get_escape_for_all_atmospheres(thefolders, thefiles; globvars...)

    # Plot it
    fax.plot(t, all_esc_df."Htn", color=cols[1])
    fax2.plot(t, all_esc_df."Dtn", color=cols[2])

    if show_all_flux
        fax.plot(t, all_esc_df."Ht", color=cols[1], linestyle="--")
        fax.plot(t, all_esc_df."Hn", color=cols[1], linestyle=":")
        fax2.plot(t, all_esc_df."Dt", color=cols[2], linestyle="--")
        fax2.plot(t, all_esc_df."Dn", color=cols[2], linestyle=":")

        L2D = PyPlot.matplotlib.lines.Line2D
        lines = [L2D([0], [0], color=cols[1]), 
                 L2D([0], [0], color=cols[1], linestyle="--"),
                 L2D([0], [0], color=cols[1], linestyle=":")]  
        fax2.legend(lines, ["total", "thermal", "nonthermal"], fontsize=14, loc=lloc)
    end
 
    return all_esc_df."Htn", all_esc_df."Dtn"    
end

function density_vs_time(n_ax, thefiles, sole_color, A; alt2=nothing, alt3=nothing, altcols=nothing, density_ticks=[1e3, 1e4, 1e5, 1e6, 1e7], globvars...)
    #=
    Makes panel of H and D densities over time.
    Inputs:
        n_ax: axis object on which to plot densities
        thefiles: all the filenames for which to plot results
        sole_color: color to use when plotting only one altitude
        A: main altitude for which to plot results.
        

    Optional:
        alt2, alt3: second/third altitude option
        altcols: an array of shape (3, 3) which can be used to plot values for up to 3 altitudes.
        density_ticks: ticks for the y axis 

    Output: n/a
    =#
    GV = values(globvars)
    required = [:n_alt_index]
    check_requirements(keys(GV), required)
    
    
    nH = Array{Float64}(undef, length(thefiles), 3)
    nD = Array{Float64}(undef, length(thefiles), 3)

    # Collect densities into new arrays
    for (i, a) in enumerate(thefiles)
        n = get_ncurrent(a)
        nH[i, 1] = n[:H][GV.n_alt_index[A*1e5]]
        nD[i, 1] = n[:D][GV.n_alt_index[A*1e5]]

        if alt2 != nothing
            nH[i, 2] = n[:H][GV.n_alt_index[alt2*1e5]]
            nD[i, 2] = n[:D][GV.n_alt_index[alt2*1e5]]
        end
        if alt3 != nothing
            nH[i, 3] = n[:H][GV.n_alt_index[alt3*1e5]]
            nD[i, 3] = n[:D][GV.n_alt_index[alt3*1e5]]
        end
    end

    # PLOT
    x = 1:length(thefiles)
    thiscol = alt2 == nothing ? sole_color : altcols[1, :]
    n_ax.plot(x, nH[:, 1], color=thiscol)
    n_ax.plot(x, nD[:, 1], color=thiscol)

    if alt2 != nothing
        n_ax.plot(x, nH[:, 2], color=altcols[2, :])
        n_ax.plot(x, nD[:, 2], color=altcols[2, :])
    end
    if alt3 != nothing
        n_ax.plot(x, nH[:, 3], color=altcols[3, :])
        n_ax.plot(x, nD[:, 3], color=altcols[3, :])
    end

    # Text labels
    printtext=["H", "D"]
    text_loc=[[0.1, 0.95], [0.1, 0.25]]
    for (ti, l) in zip(printtext, text_loc)
        n_ax.text(l..., ti, va="top", size=20, color="black", transform=n_ax.transAxes)
    end

    n_ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(density_ticks))
end

function f_vs_time_panel(fracax, thefiles; fraccol="#009B72", subplot_lbl_loc=[0,0.92], subplot_lbl_sz=24, globvars...)
    #=
    Plots the fractionation factor as a function of time. 

    thefiles: list of files to use for plotting
    ffax: axis on which to plot

    Optional:
    colorbar_tix: ticks to assign to the colorbar for the independent variable plot (if water).
    subplot_lbl_loc: axis-based coordinates in which to place subfigure labels like "a)", etc
    subplot_lbl_sz: fontsize to use for said labels.
    =#
      
    GV = values(globvars)
    required = [ # Things from CONSTANTS.jl
               :q, :molmass, :polarizability, :collision_xsect, :DH,
               # From CUSTOMIZATIONS.jl
               :alt, :dz, :num_layers, :n_alt_index, :non_bdy_layers, :speciesstyle, 
               # Simulation-unique stuff 
                :all_species, :hHnet, :hDnet, :hH2net, :hHDnet, :hHrc, :hDrc, :hH2rc, :hHDrc]
    check_requirements(keys(GV), required)
    
    # Setup figure
    set_rc_params(fs=16, axlab=18, xtls=16, ytls=16, sansserif=sansserif_choice, monospace=monospace_choice)
    
    # What to plot against -----------------------------------------------------------------------------------------
    total_num_files = length(thefiles)
    t = 1:total_num_files
    
    # Frac axis -----------------------------------------------------------------------------------------------------
    # Tease apart folder and filenames because it's needed for getting flux
    folder_pattern = r".+v\d\/"
    file_pattern = r"/[a-z]+_.+\.h5"
    thefolders = [String(match(folder_pattern, f).match) for f in thefiles]
    filenames = [String(match(file_pattern, f).match) for f in thefiles]

    # Get escape and calculate f
    all_esc_df = get_escape_for_all_atmospheres(thefolders, thefiles; globvars...)

    # Plot it
    fcol = "xkcd:merlot"
    # plot_bg(fracax)
    fracax.plot(t, (all_esc_df."Dtn" ./ all_esc_df."Htn") ./ GV.DH, color=fraccol)
    fracax.set_yscale("log")
end

function DH_vs_time(ax, atmfile_list, whichalt, col; DHmult = [2, 5, 10, 20, 25],
                        heavysp=:D, lightsp=:H, species_pair="atomics", #ylims=[1e-4, 1e-2], cutoff=nothing, fn="DH_profiles", 
                        lw=1, ls="-", plot_waterdh=false, globvars...)
    #=
    Plots D/H over the anual cycle in a single panel as part of a larger figure.
    
    Inputs:
        ax: axis object on which to plot
        atmdict_list: List of atmospheric state dictionaries to use, in order.
        whichalt: altitude at which to plot
        col: color to use for the main plot line
    Optional:
        heavysp: heavier of the two species, usually :D
        lightsp: usually :H
        species_pair: isotope pairs; "atomics", "water", etc. 
        wmark, omark: markers to use in plotting the lines for water and "other", usually unused.
        lw: linewidth
        ls: linestyle
        plot_waterdh: whether to plot D/H in water, almost never used

        
    
    wl, ol: where to place labels for the lines
    
    =#

    GV = values(globvars)
    required = [:n_alt_index]
    check_requirements(keys(GV), required)
        
    # Set a factor used to identify number of H atoms
    if occursin("H2", string(lightsp))
        Hfactor_light = 2
    elseif occursin("H3", string(lightsp))
        Hfactor_light = 3
    else
        Hfactor_light = 1
    end
    
    # Populate the plottable D/H arrays
    DH_water = Array{Float64}(undef, length(atmfile_list), 1)
    DH_other = Array{Float64}(undef, length(atmfile_list), 1)
    
    for (i, a) in enumerate(atmfile_list)
        n = get_ncurrent(a)
        # get D/H in water
        DH_water[i] = n[:HDO][GV.n_alt_index[whichalt*1e5]] ./ (2 .* n[:H2O][GV.n_alt_index[whichalt*1e5]])
        DH_other[i] = n[heavysp][GV.n_alt_index[whichalt*1e5]] ./ (Hfactor_light * n[lightsp][GV.n_alt_index[whichalt*1e5]])
    end
    
    # plot actual content
    ax.plot(1:length(atmfile_list), DH_other, zorder=10, linewidth=lw, linestyle=ls, color=col[1], label=species_pair)
    if plot_waterdh
        ax.plot(1:length(atmfile_list), DH_water, zorder=10, linewidth=lw, linestyle=ls, color=col[2], label="water")
    end

    # Show SMOW lines for orientation
    # multiplier = DHmult
    draw_DH_lines(ax, DHmult)
    
    # Text showing which line is which
    heavysp_str = "[$(string_to_latexstr(string(heavysp); dollarsigns=false))]"
    lightsp_str = "[$(string_to_latexstr(string(lightsp); dollarsigns=false))]"
    if Hfactor_light > 1
        lightsp_str = "$(Hfactor_light)"*lightsp_str
    end
end

function reincorp_vs_time(fax, fax2, t, thefiles, cols, change_indices; typeflag="temp", globvars...)
    #=
    Not currently using this because whether it returns a believable value is suspect

    t: time array to plot against
    thefiles: list of filenames
    A: an altitude
    col: color array
    =#
    
    GV = values(globvars)
    
    # Tease apart folder and filenames because it's needed for getting flux
    folder_pattern = r".+v\d\/"
    file_pattern = r"/[a-z]+_.+\.h5"
    thefolders = [String(match(folder_pattern, f).match) for f in thefiles]
    filenames = [String(match(file_pattern, f).match) for f in thefiles]

    # Dataframes
    reincorp_array = Array{Float64}(undef, length(thefiles), 2)
    println(size(reincorp_array))
    
    i = 1
    ti = 1# index of temp profile to use; there are 5. 
    ci = 1 # change_index tracker
    for (foldy, f) in zip(thefolders, filenames)
        this_n = get_ncurrent(foldy*f)
        
        if (i == change_indices[ci]) & (i != 1)
            ti += 1
            ci += 1           
        end
        
        # Need to get the temperature arrays, but we have stored them all for every file
        if typeflag=="temp"
            # Extract the first column for backward compatibility with 1-D analysis
            this_Tn = GV.Tn[1, 2:end-1]
        else
            # Extract the first column for backward compatibility with 1-D analysis
            this_Tn = GV.Tn[1, :]
        end
            
        println(length(this_Tn))

        Hreinc, Dreinc = reincorporate_H_D(this_n; GV.all_species, GV.ion_species, GV.num_layers, GV.reaction_network, Tn=this_Tn, Ti=GV.Ti, Te=GV.Te)

        reincorp_array[i, 1] = Hreinc
        reincorp_array[i, 2] = Dreinc
        i += 1
    end

    # Plot it
    fax.plot(t, reincorp_array[:, 1], color=cols[1])
    fax2.plot(t, reincorp_array[:, 2], color=cols[2])
    
    return reincorp_array 
end

# Helpers
function get_cbar_size(the_ax)
    #=
    Creates and returns a list of parameters to properly size a colorbar so it has the same height as the_ax. 
    =#
    the_ax_bbox = the_ax.get_position().get_points() # left, bottom, width, height 
    xmin = the_ax_bbox[1, 1]
    xmax = the_ax_bbox[2, 1]
    ymin = the_ax_bbox[1, 2]
    ymax = the_ax_bbox[2, 2]
    return [xmax + 0.02, ymin, (xmax-xmin)/10, ymax-ymin]
end

function format_time_axis(the_ax, delta_inds, n, lbls)
    #=
    Formatting wizardry for the_ax, which will have the x ticks for a plot of stuff over time.
        delta_inds: wherever seasons change
        n: length(thefiles)
    =#
    the_ax.minorticks_on()
    major_tix = [delta_inds..., n]
    minor_tick_locs = [(major_tix[i] + major_tix[i+1])/2 for i in 1:length(major_tix)-1]
    the_ax.tick_params(axis="x", which="major", bottom=true, labelbottom=false)
    the_ax.tick_params(axis="x", which="minor", bottom=false, labelbottom=true)
    the_ax.set_xticks(minor_tick_locs, labels=lbls, minor=true)
end

function collect_atmospheres(season_folders, resolution, set; globvars...)
    #=
    Helper function to retrieve variables necessary to make complicated plots.
    Inputs:
        season_folders: ordered list of folders containing simulations for each Mars season
        resolution: "high" or "low"
        set: "temp" or "water"

    Outputs:
        atm_states: A dictionary containing sub-dictionaries which each contain the final atmospheric state for the season
        state_files: list of files, in order that they will be analyzed and added to plots
        Tn_all, Ti_all, Te_all: Arrays of shape (num_layers, len(state_files)) providing all the temperature profiles
    =#
    GV = values(globvars)
    required = [:alt, :M_P, :R_P, :planet]
    check_requirements(keys(GV), required)
    

    if set=="temp"
        atm_keys = Dict(1=>"midT1", 2=>"highT", 3=>"midT2", 4=>"lowT", 5=>"midT3")
        atm_states = Dict("lowT"=>Dict(), "midT1"=>Dict(), "highT"=>Dict(), "midT2"=>Dict(), "midT3"=>Dict()) # Final atmospheres
        vardicts = Dict("lowT"=>Dict(), "midT1"=>Dict(), "highT"=>Dict(), "midT2"=>Dict(), "midT3"=>Dict())
    elseif set=="water"
        atm_keys = Dict(1=>"mean", 2=>"high", 3=>"mean2", 4=>"low", 5=>"mean3")
        atm_states = Dict("low"=>Dict{Symbol, Vector{Float64}}(), "mean"=>Dict{Symbol, Vector{Float64}}(), "high"=>Dict{Symbol, Vector{Float64}}(), 
                          "mean2"=>Dict{Symbol, Vector{Float64}}(), "mean3"=>Dict{Symbol, Vector{Float64}}()) # final atmospheres 
        atm_states_init = Dict("low"=>Dict{Symbol, Vector{Float64}}(), "mean1"=>Dict{Symbol, Vector{Float64}}(), "high"=>Dict{Symbol, Vector{Float64}}(), 
                               "mean2"=>Dict{Symbol, Vector{Float64}}(), "mean3"=>Dict{Symbol, Vector{Float64}}())
        vardicts = Dict("low"=>Dict(), "mean1"=>Dict(), "high"=>Dict(), "mean2"=>Dict(), "mean3"=>Dict())
    end
    
    Tn_all = Array{Float64}(undef, length(GV.alt), length(season_folders))
    Ti_all = Array{Float64}(undef, length(GV.alt), length(season_folders))
    Te_all = Array{Float64}(undef, length(GV.alt), length(season_folders))
    
    for i in 1:length(season_folders)
        thefolder = season_folders[i]

        atm_states[atm_keys[i]] = get_ncurrent(thefolder*"final_atmosphere.h5")
        if set=="water"
            atm_states_init[atm_keys[i]] = get_ncurrent(thefolder*"initial_atmosphere.h5")
        end

        vardicts[atm_keys[i]] = load_from_paramlog(thefolder; globvars...)

        Tn_all[:, i] = vardicts[atm_keys[i]]["Tn_arr"];
        Ti_all[:, i] = vardicts[atm_keys[i]]["Ti_arr"];
        Te_all[:, i] = vardicts[atm_keys[i]]["Te_arr"];
    end

    if resolution=="high"
        # Get a list of all atmospehric state files and collect the water data for the independent variable axis
        state_files = collect(Iterators.flatten([[thisfol*a for a in search_subfolders(thisfol, ".h5"; type="files")[1:end-2]] for thisfol in season_folders]));
    elseif resolution=="low"
        # Same but for select atmospheres
        state_files = [collect(Iterators.flatten([[thisfol*a for a in search_subfolders(thisfol, r"(0[1-9]|19|29)\.h5"; type="files")[1:end-2]] for thisfol in season_folders]))...];
    end

    if set=="temp"
        return atm_states, state_files, Tn_all, Ti_all, Te_all
    else
        return atm_states, atm_states_init, state_files, Tn_all, Ti_all, Te_all
    end
end

function get_escape_for_all_atmospheres(thefolders, thefiles; globvars...)
    #=
    Helper function that wil either load the XLSX file containing escape of H and D for all files in thefiles
    or make and save a new XLSX file. 

    Inputs:
        thefolders: list of parent folders for each entry in thefiles. 
        thefiles: Ordered list of files in thefolders. 

    Outputs: a dataframe containing the thermal, non-thermal, and total escape of H and D for each particular file in thefiles.

    Called in esc_flux_vs_time and f_vs_time
    =#
    # Check for a saved file with the escape values in it
    parent_folder = String(match(r".+(?=I-)", thefiles[1]).match)
    println("Searching $(parent_folder) for an escape spreadsheet")
    flush(stdout)

    if isfile(parent_folder*"escape.xlsx")
        println("Found the spreadsheet")
        flush(stdout)
        all_esc_df = DataFrame(XLSX.readtable(parent_folder*"escape.xlsx", "Escape"))
    else # If no saved file, generate the escapes
        println("No saved escape info found, generating escape for all folders and files in $(parent_folder) and creating spreadsheet")
        flush(stdout)

        file_pattern = r"/[a-z]+_.+\.h5"

        # This is what takes the most time and makes it take forever to do seasonal_cycling_figure.
        filenames = [String(match(file_pattern, f).match) for f in thefiles]
        escdfs = [final_escape(foldy, f; globvars...) for (foldy, f) in zip(thefolders, filenames)]

        H_tn = collect(Iterators.flatten([df_lookup(i, "EscapeType", "Total", "H") for i in escdfs]))
        H_t = collect(Iterators.flatten([df_lookup(i, "EscapeType", "Thermal", "H") for i in escdfs]))
        H_n = collect(Iterators.flatten([df_lookup(i, "EscapeType", "Nonthermal", "H") for i in escdfs]))
        
        D_tn = collect(Iterators.flatten([df_lookup(i, "EscapeType", "Total", "D") for i in escdfs]))
        D_t = collect(Iterators.flatten([df_lookup(i, "EscapeType", "Thermal", "D") for i in escdfs]))
        D_n = collect(Iterators.flatten([df_lookup(i, "EscapeType", "Nonthermal", "D") for i in escdfs]))

        all_esc_df = DataFrame("Ht"=>H_t, "Hn"=>H_n, "Htn"=>H_tn, "Dt"=>D_t, "Dn"=>D_n, "Dtn"=>D_tn)

        # Write out the escape to the file
        XLSX.writetable(parent_folder*"escape.xlsx", "Escape"=>all_esc_df)
        println("Finished generating and saving escape info")
        flush(stdout)
    end
    return all_esc_df
end
