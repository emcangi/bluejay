# **************************************************************************** #
#                                                                              #
#                             Plotting Functions                               #
#                                                                              #
# **************************************************************************** #

function get_colors(L::Int64, cmap; strt=0, stp=1)
    #=
    Generates some colors based on a non-gradient color map for use in plotting a 
    bunch of lines all at once.
    Input:
        L: number of colors to generate.
        cmap: color map name
    Output:
        An array of RGB values: [[R1, G1, B1], [R2, G2, B2]...]
    =#

    Cmap = get_cmap(Symbol(cmap))
    colors = [Cmap(x) for x in range(strt, stop=stp, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors))
        c[i, 1] = colors[i][1]
        c[i, 2] = colors[i][2]
        c[i, 3] = colors[i][3]
    end
    return c
end

function get_grad_colors(L::Int64, cmap; strt=0, stp=1)
    #=
    Generates some colors based on a GRADIENT color map for use in plotting a 
    bunch of lines all at once.

    Input:
        L: number of colors to generate.
        cmap: color map name
        strt and stp: By setting these to other values between 0 and 1 you can restrict 
                      the range of colors drawn from.
    Output:
        An array of RGB values: [[R1, G1, B1], [R2, G2, B2]...]

    AVAILABLE MAPS: blues, viridis, pu_or, magma, plasma, inferno
    =#

    colors = [cgrad(Symbol(cmap))[x] for x in range(strt, stop=stp, length=L)]
    c = Array{Float64}(undef, L, 3)

    for i in range(1, length=length(colors))
        c[i, 1] = red(colors[i])
        c[i, 2] = green(colors[i])
        c[i, 3] = blue(colors[i])
    end
    return c
end

function plot_atm(atmdict::Dict{Symbol, Vector{ftype_ncur}}, savepath::String, atol, E_prof; print_shortcodes=true, 
                  t="", showonly=false, xlab=L"Species concentration (cm$^{-3}$)", xlim_1=(1e-12, 1e18), xlim_2=(1e-5, 1e5), 
                  legloc=[0.8,1], globvars...)
    #=
    Makes a "spaghetti plot" of the species concentrations by altitude in the
    atmosphere. 

    Input:
        atmdict: dictionary of vertical profiles of a plottable quantity, keys are species names.
        savepath: path and name for saving resulting .png file
        atol: absolute tolerance to plot
        E_prof: E densities for plotting the electron line
        print_shortcodes: whether to print unique simulation IDs on plots
        t: title text for whole plot
        showonly: whether to just show() the plot instead of saving. If setting to true, send in junk string for savepath.
        xlab: override for x axis label
        xlim_1: override for x axis limits in column 1 or entire plot of length(species_lists)==1
        xlim_2: override for x axis limits in column 2
        legloc: legend location
    Output:
        Big atmospheric density plot
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:neutral_species, :plot_grid, :speciescolor, :speciesstyle, :zmax])

    if print_shortcodes
        @assert all(x->x in keys(GV),  [:hrshortcode, :rshortcode])
    end 

    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    
    # Establish logical groups for axes =======================================================

    # This dictionary was created by hand. Species may be missing. If so, an error will automatically be thrown by the code
    species_groups = Dict( # PRIMARY NEUTRALS + IONS
                            1=>[:CO2,:CO2pl,
                                :CO,:COpl,
                                :N2,:Nup2D,:N2pl,
                                :Ar,:Arpl,
                                :ArHpl,:ArDpl,
                                :HCO,:HCOpl,:DCO,:DCOpl,
                                :HOCpl,:DOCpl,
                                :HOCO,:DOCO, :HCO2pl,:DCO2pl,
                                :O1D,
                                :O,:Opl,
                                :O2,:O2pl,
                                :O3],
                            # ODD HYDROGEN + IONS
                            2=>[:H,:Hpl,:D,:Dpl,
                                :H2,:H2pl, :HD, :HDpl,
                                :H3pl,:H2Dpl, :HD2pl, 
                                :H2O,:H2Opl,:HDO,:HDOpl,
                                :H3Opl,:H2DOpl,
                                :H2O2,:HDO2,
                                :HO2,:HO2pl,:DO2,
                                :OH,:OHpl,:OD,:ODpl],
                            # NITROGEN NEUTRALS + IONS and RARE C MOLECULES
                            3=>[:C,:Cpl,:CH,:CHpl,
                                :CN,:CNpl,:HCN,:HCNpl,:HCNHpl,
                                :HNO,:HNOpl,:HN2Opl,
                                :N,:Npl,:N2O,:N2Opl,
                                :NH,:NHpl,:NH2,:NH2pl,:NH3pl,
                                :NO,:NOpl,:NO2,:NO2pl,
                                :N2Hpl,:N2Dpl]
                        );

    axes_by_sp = Dict()

    for k in keys(species_groups)
        for sp in species_groups[k]
            axes_by_sp[sp] = k
        end
    end

    # Plot neutrals and ions together =========================================================
    if haskey(GV, :ion_species) # neutrals and ions 
        
        # set up the overall plot -------------------------------------------------------------
        atm_fig, atm_ax = subplots(3, 2, sharex=false, sharey=true, figsize=(14, 16))
        subplots_adjust(wspace=0, hspace=0)
        tight_layout()
                
        # only the neutral-col axes
        atm_ax[1, 1].set_title("Neutrals")
        for i in 1:3
            plot_bg(atm_ax[i, 1])
            atm_ax[i, 1].set_xlim(xlim_1[1], xlim_1[2])
            atm_ax[i, 1].fill_betweenx(GV.plot_grid, xlim_1[1] .* ones(size(GV.plot_grid,)), x2=atol, alpha=0.1, color=medgray, zorder=10)
            atm_ax[i, 1].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
            atm_ax[i, 1].set_ylabel("Altitude (km)")
        end
        atm_ax[3, 1].set_xlabel(xlab)
        
        # only the ion-col axes
        atm_ax[1, 2].set_title("Ions")
        for i in 1:3
            plot_bg(atm_ax[i, 2])
            atm_ax[i, 2].set_xlim(xlim_2[1], xlim_2[2])
            atm_ax[i, 2].fill_betweenx(GV.plot_grid, xlim_2[1] .* ones(size(GV.plot_grid)), x2=atol, alpha=0.1, color=medgray, zorder=10)
            atm_ax[i, 2].tick_params(which="both", labeltop=false, top=true, labelbottom=true, bottom=true)
        end
        atm_ax[3, 2].set_xlabel(xlab)
         
        # plot the neutrals according to logical groups -------------------------------------------------------
        for sp in GV.neutral_species
            atm_ax[axes_by_sp[sp], 1].plot(convert(Array{Float64}, atmdict[sp]), GV.plot_grid, color=get(GV.speciescolor, sp, "black"),
                                           linewidth=2, label=string_to_latexstr(string(sp)), linestyle=get(GV.speciesstyle, sp, "-"), zorder=2)
        end
        
        # plot the ions according to logical groups ------------------------------------------------------------
        for sp in GV.ion_species
            atm_ax[axes_by_sp[sp], 2].plot(convert(Array{Float64}, atmdict[sp]), GV.plot_grid, color=get(GV.speciescolor, sp, "black"),
                                           linewidth=2, label=string_to_latexstr(string(sp)), linestyle=get(GV.speciesstyle, sp, "-"), zorder=2)
        end

        # plot electron profile --------------------------------------------------------------------------------
        atm_ax[1, 2].plot(convert(Array{Float64}, E_prof), GV.plot_grid, color="black", linewidth=2, linestyle=":", zorder=10, label=L"e$^-$")

        # stuff that applies to all axes
        for r in 1:size(atm_ax)[1]
            for c in 1:size(atm_ax)[2]
                atm_ax[r,c].set_ylim(0, GV.zmax/1e5)
                atm_ax[r,c].set_xscale("log")
                handles, labels = atm_ax[r,c].get_legend_handles_labels()
                if isempty(handles) == false
                    x, y = legloc
                    if (r==1) && (c==2)
                        x, y = 1.01, 1
                    end
                    atm_ax[r,c].legend(handles, labels, fontsize=12, bbox_to_anchor=[x,y], loc=2, borderaxespad=0)
                end
            end
        end

    # Plot only neutrals - to support the fractionation factor project ==========================================
    else # ion species is not defined 
        atm_fig, atm_ax = subplots(figsize=(16,6))
        tight_layout()
        for sp in GV.neutral_species
            atm_ax.plot(convert(Array{Float64}, atmdict[sp]), GV.plot_grid, color=get(GV.speciescolor, sp, "black"),
                        linewidth=2, label=sp, linestyle=get(GV.speciesstyle, sp, "-"), zorder=1)
            atm_ax.set_xlim(xlim_1[1], xlim_1[2])
            atm_ax.set_ylabel("Altitude [km]")
        end
        atm_ax.tick_params(which="both", labeltop=true, top=true)
        plot_bg(atm_ax)
        atm_ax.set_ylim(0, GV.zmax/1e5)
        atm_ax.set_xscale("log")
        atm_ax.set_xlabel(xlab)
        atm_ax.legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0, fontsize=16)
    end

    suptitle(t, y=1.03)

    # Shortcodes as watermarks
    if print_shortcodes
        text(0.9, 1.05, GV.hrshortcode, transform=gcf().transFigure, color="dimgrey", ha="right")
        text(0.9, 1.02, GV.rshortcode, transform=gcf().transFigure, color="dimgrey", ha="right")
    end

    if showonly==false  
        atm_fig.savefig(savepath, format="pdf", bbox_inches="tight", dpi=300)
        close(atm_fig)
    else
        show()
    end
end

function plot_bg(axob; bg="#ededed")
    #=
    Make plots not look ugly. 

    intended to immitate seaborn
    =#
    axob.set_facecolor(bg)
    axob.grid(zorder=-5, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function plot_extinction(solabs; fnextr="", path=nothing, tauonly=false, xsect_info=nothing, solflux=nothing, extra_t="", linth=1e-8, vm=1e-6, globvars...)
    #=

    TODO: This is out of date and needs updating

    Used to plot the extinction in the atmosphere by wavelength and altitude

    solabs: a ROW VECTOR (sadly) of solar absorption in the atmosphere. rows = altitudes, "columns" = wavelengths,
            but being a row vector it supposedly has only one column, but each element is a row.
    path: a path at which to save the plot
    tauonly: plot just tau rather than e^-tau
    xsect_info: a list containing two things:
                xsect: list of crosssections for some Jrate, shape (124, 2000)
                xsect_sp: string representation of a chemical species which is absorbing using xsect
    solflux: solar flux values of shape (1, 2000), optional
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:zmax, :plot_grid])

    
    fig, ax = subplots()

    plot_bg(ax)

    num_x = length(solabs[1])  # extinction is in an awkward format
    num_y = size(solabs)[1]
    println("num_x = $(num_x), num_y = $(num_y)")
    X = 1:num_x  
    Y = GV.plot_grid

    # Don't mess with this line. It is witchcraft that translates the row vector that is extinction
    # into an actual 2D matrix so that we can do normal math on it. 
    solabs = transpose(reshape(collect(Iterators.flatten(solabs[:,:])), (num_x,num_y)))

    # now convert it to actual extinction:
    if tauonly # plot the opacity (I think)
        # println(solabs)
        solabs[solabs .< 0] .= 0
        # z = log10.(solabs)
        z = solabs
        titlestr = L"\tau"
        savestr = "tau_$(fnextr).png"
        z_min = 0
        z_max = 1
        heatmap = ax.pcolor(X, Y, z, cmap="bone", vmin=z_min, vmax=z_max)
        cbarlbl = L"\tau"
    else  # this is if we want to plot the extinction.
        if xsect_info==nothing
            throw("Error! Do you want to plot tau or the extinction? Please pass either tauonly=true or xsect_info")
        end
        xsect, jr = xsect_info
        z = exp.(-solabs)
        
        titlestr = L"e^{-\tau}\sigma_{\lambda}"
        if xsect != nothing  # multiply by the crosssection if we want to plot the Jrates
            xsect = transpose(reshape(collect(Iterators.flatten(xsect[:,:])), (num_x,num_y)))
            z = z .* xsect
            titlestr = L"e^{-\tau}\sigma_{\lambda}" * ", $(string_to_latexstr(string(jr)))"
            savestr = "extinction_$(fnextr)_$(jr).png"
        end

        if solflux != nothing
            # The following multiplies the solar flux at every wavelength by z, which currently = exp(-solabs). 
            @time z = solflux[:, 2]' .* z # solflux is an an array (actually a vector, tbh) of shape (1, 2000);
                                    # exp.(-solabs) is an array of shape (124, 2000). This multiplies the solar flux
                                    # values by the extinction across the 2000 wavelengths. Unbelievably, this 
                                    # operation works in Julia in either direction, whether you do solflux .* z or reverse.

            # Actinic flux
            # @time for ialt in 1:num_layers
            #     z[ialt, :] = solflux[:,2] .* z[ialt, :]#exp.(-solarabs[ialt])
            # end


            titlestr = L"J_{\lambda} e^{-\tau}\sigma_{\lambda}" * ", $(string_to_latexstr(string(jr)))" * " $(extra_t)"
        end
        heatmap = ax.pcolor(X, Y, z, cmap="bone", norm=PyPlot.matplotlib.colors.SymLogNorm(vmin=0, vmax=vm, linthresh=linth))# 
        cbarlbl = "photons/sec"
    end

    cbar = fig.colorbar(heatmap, ax=ax)
    cbar.set_label(cbarlbl, rotation=270, fontsize=14, labelpad=15)
    xlabel("Wavelength (nm)")
    ylabel("Altitude (km)")
    yticks(0:50:(GV.zmax/1e5))

    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(25))

    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(25))

    title(titlestr)
    ax.set_xlim(0, 300)  # turn off to see full range to 2000 nm.
    show()
    if path != nothing
        savefig(path*savestr, bbox_inches="tight")
        close(fig)
    end
end

function plot_Jrates(sp, atmdict::Dict{Symbol, Vector{ftype_ncur}}; savedir=nothing, opt="", globvars...)                
    #=
    Plots the Jrates for each photodissociation or photoionizaiton reaction. Override for small groups of species.
    Input:
        sp: species for which to plot the Jrates
        atmdict: Present atmospheric state dictionary
        savedir: directory in which to save the plots
        Optional:
            opt: Extra string for the filename
    Output:
        Plot of Jrates by altitude
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:all_species, :ion_species, :num_layers, :plot_grid, :reaction_network, :speciesbclist, :Tn, :Ti, :Te])

    # Plot setup
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

    # --------------------------------------------------------------------------------
    # make plot
    rxd_prod, prod_rc = get_volume_rates(sp, atmdict; which="Jrates", globvars..., Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1]) 
    rxd_loss, loss_rc = get_volume_rates(sp, atmdict; which="Jrates", globvars..., Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1]) 

    fig, ax = subplots(figsize=(8,6))
    plot_bg(ax)

    minx = 1e10
    maxx = 0
    
    # Collect chem production equations and total 
    if !isempty(keys(rxd_prod))
        for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
            lbl = "$(kv[1])"

            # if source2 != nothing 
            #     if !all(x->x<=1e-10, abs.(kv[2] - rxd_prod2[kv[1]]))
            #         ax.semilogx(kv[2] - rxd_prod2[kv[1]], GV.plot_grid, linestyle="-", linewidth=1, label=lbl)
            #     else
            #         text(0.5, 0.5, "Everything is basically 0, nothing to plot", transform=ax.transAxes)
            #     end
            # else
            #     ax.semilogx(kv[2], GV.plot_grid, linestyle="-", linewidth=1, label=lbl)
            # end
            ax.semilogx(kv[2], GV.plot_grid, linestyle="-", linewidth=1, label=lbl)
            
            # set the xlimits
            if minimum(kv[2]) <= minx
                minx = minimum(kv[2])
            end
            if maximum(kv[2]) >= maxx
                maxx = maximum(kv[2])
            end
        end

        suptitle("J rates for $(sp), $(opt)", fontsize=20)
        ax.set_ylabel("Altitude (km)")
        ax.legend(bbox_to_anchor=(1.01, 1))
        ax.set_xlabel(L"Reaction rate (cm$^{-3}$s$^{-1}$)")
        if minx < 1e-14
            minx = 1e-14
        end
        maxx = 10^(ceil(log10(maxx)))
        ax.set_xlim([minx, maxx])

        if savedir!=nothing
            savefig(savedir*"J_rates_$(sp)_$(opt).png", bbox_inches="tight", dpi=300)
            close(fig)
        else
            show()
        end
    end
end

function plot_production_and_loss(final_atm, results_dir, thefolder; globvars...)
    #=

    =#
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :chem_species, :collision_xsect, :dz, :hot_D_rc_funcs, :hot_H_rc_funcs, 
                                  :hot_H2_rc_funcs, :hot_HD_rc_funcs, :Hs_dict, :hot_H_network, :hot_D_network, :hot_H2_network, :hot_HD_network,
                                  :hrshortcode, :ion_species, :Jratedict, :molmass, :neutral_species, :non_bdy_layers, :nonthermal,
                                  :num_layers, :n_all_layers, :n_alt_index, :polarizability, :plot_grid, :q, :rshortcode, :reaction_network,
                                  :speciesbclist, :Tn, :Ti, :Te, :Tp, :Tprof_for_Hs, :Tprof_for_diffusion, 
                                  :transport_species, :upper_lower_bdy_i, :upper_lower_bdy, :zmax])

    println("Creating production and loss plots to show convergence of species...")
    create_folder("chemeq_plots", results_dir*thefolder*"/")
    for sp in GV.all_species
        plot_rxns(sp, final_atm, results_dir; subfolder=thefolder,num="final_atmosphere", globvars...)
    end
    println("Finished convergence plots")
end

function plot_rxns(sp::Symbol, atmdict::Dict{Symbol, Vector{ftype_ncur}}, results_dir::String; 
                   nonthermal=true, shown_rxns=nothing, subfolder="", plotsfolder="chemeq_plots", dt=nothing, num="", extra_title="", 
                   plot_timescales=false, plot_total_rate_coefs=false, showonly=false, globvars...)
    #=
    Plots the production and loss rates from chemistry, transport, and both together by altitude for a given species, sp, 
    at a given snapshot in time of the atmosphere atmdict. 

    Inputs:
        sp: species to make plot for
        atmdict: the atmospheric state dictionary
        Tn, Ti, Te: atmospheric temperature arrays
        bcdict: Boundary conditions dictionary specified in parameters file
        reaction_network: Chemical reaction network
        plot_indiv_rxns: whether to plot lines for the individual chemical reactions. if false, only the total will be plotted.
        thresh: a threshhold of reaction rate. will only plot the reactions that have values above this threshhold
        subfolder: an optional subfolder within results_dir
        plotsfolder: specific folder within which the plots will go.
        dt: timestep, optional, to be included in the plot title
        num: optional string to add to filename
        extra_title: extra text for the title, e.g. "converged"
        plot_timescales: plots the timescale over which a species population changes. 
        plot_total_rate_coefs: plots chemical produciton and loss rates only (not multiplied by density)
        showonly: if true, will show the plots instead of saving them
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:all_species, :alt, :chem_species, :dz, :hrshortcode, :Hs_dict, :ion_species, 
                                    :molmass, :n_all_layers, :n_alt_index, :neutral_species, :num_layers, 
                                    :plot_grid, :polarizability, :q, :rshortcode, :reaction_network, :speciesbclist, 
                                    :Te, :Ti, :Tn, :Tp, :Tprof_for_Hs, :Tprof_for_diffusion, :transport_species, 
                                    :upper_lower_bdy, :upper_lower_bdy_i])

    # ================================================================================
    # Plot setup stuff
    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"]= 16
    rcParams["xtick.labelsize"] = 16
    rcParams["ytick.labelsize"] = 16

    if plot_timescales==true
        total_ax = 4
    else
        total_ax = 3
    end

    fig, ax = subplots(1, total_ax, sharey=true, figsize=(20,6))
    # subplot order will be chemistry, transport, total.
    subplots_adjust(wspace=0.1, bottom=0.25)
    for a in ax
        plot_bg(a)
    end

    # this will be used to determine the x limits. the three entries are the chemistry panel, transport panel, and the total panel
    # these are just starting values. 
    minx = [1e-8, 1e-8, 1e-8]  
    maxx = [1e2, 1e2, 1e2]

    # ================================================================================
    # Plot the chemical production and loss rates: Panel #1

    if shown_rxns != nothing
        shown_rxns = [format_chemistry_string(r[1], r[2]) for r in shown_rxns]
    end

    # Arrays to store the total reactions per second for this species of interest
    total_prod_rate = zeros(GV.num_layers)
    total_loss_rate = zeros(GV.num_layers)

    # Arrays to hold the total chemical production and loss 
    total_chem_prod = zeros(GV.num_layers)
    total_chem_loss = zeros(GV.num_layers)
    total_chem_prod_ratecoef = zeros(GV.num_layers)
    total_chem_loss_ratecoef = zeros(GV.num_layers)

    if sp in GV.chem_species
        # --------------------------------------------------------------------------------
        # calculate reaction rates x density for all reactions of the species at each level of the atmosphere.
        # Temperature arrays include boundary layers so that diffusion can be calculated later, so we have to 
        # index these in this way to make the evaluation of chemistry reaction rate coefficients work. 
        # Entering them separately from globvars allows us to keep passing globvars as a "packed" variable but 
        # use the most recent Tn, Ti, Te according to rightmost taking precedence.
        rxd_prod, rate_coefs_prod = get_volume_rates(sp, atmdict; species_role="product", globvars..., Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1])
        rxd_loss, rate_coefs_loss = get_volume_rates(sp, atmdict; species_role="reactant", globvars..., Tn=GV.Tn[2:end-1], Ti=GV.Ti[2:end-1], Te=GV.Te[2:end-1])

        # Water is turned off in the lower atmosphere, so we should represent that.
        if sp in [:H2O, :HDO] 
            for (prod_k, loss_k) in zip(keys(rxd_prod), keys(rxd_loss))
                rxd_prod[prod_k][1:GV.upper_lower_bdy_i] .= NaN
                rxd_loss[loss_k][1:GV.upper_lower_bdy_i] .= NaN
            end
        end

        # set up a selection of colors and linestyles so that individual reaction rate lines are distinguishable
        cols = ["#cc156d", "#7e5fe7", "#b26145", "#257950", "#265582", "#fc5931"]  # colorgorical
        ls = ["-", "--", (0, (3, 1, 1, 1)), ":"]

        # ---------------------------------------------------------------------------------
        col_i = 1
        ls_i = 1

        # Pattern to help identify whether there is more than one product of species sp in the reaction
        pat = r"((?<=--> ).+)"

        # Chemical production - add up total, and plot individual reactions if needed
        for kv in rxd_prod  # loop through the dict of format reaction => [rates by altitude]
            if shown_rxns != nothing
                if kv[1] in shown_rxns 
                    ax[1].semilogx(kv[2], GV.plot_grid, linestyle=ls[ls_i], marker=9, markevery=20, color=cols[col_i], linewidth=1, label=kv[1])
                    col_i = next_in_loop(col_i, length(cols))
                    ls_i = next_in_loop(ls_i, length(ls))
                end
            end
            # Add up the total chemical production. Accounts for cases where species is produced more than once. 
            prods = split(match(pat, kv[1])[1], " + ")
            num_created = count(i->(i==string(sp)), prods)
            total_chem_prod += num_created .* kv[2]
        end
        for kv in rate_coefs_prod
            total_chem_prod_ratecoef += kv[2]
        end

        # Chemical loss 
        for kv in rxd_loss
            if shown_rxns != nothing
                if kv[1] in shown_rxns
                    ax[1].semilogx(kv[2], GV.plot_grid, linestyle=ls[ls_i], marker=8, markevery=20, color=cols[col_i], linewidth=1, label=kv[1])
                    col_i = next_in_loop(col_i, length(cols))
                    ls_i = next_in_loop(ls_i, length(ls))
                end
            end
            total_chem_loss += kv[2]
        end
        for kv in rate_coefs_loss
            total_chem_loss_ratecoef += kv[2]
        end

        # Plot the totals 
        ax[1].semilogx(total_chem_prod, GV.plot_grid, color="xkcd:forest green", linestyle=(0, (4,2)), marker=9, markevery=20, linewidth=2, label="Total chemical production", zorder=5)
        ax[1].semilogx(total_chem_loss, GV.plot_grid, color="xkcd:shamrock", linestyle=(0, (4,2)), marker=8, markevery=20, linewidth=2, label="Total chemical loss", zorder=5)

        # set the x lims for chem axis
        prod_without_nans = filter(x->!isnan(x), total_chem_prod)
        loss_without_nans = filter(x->!isnan(x), total_chem_loss)
        minx[1] = minimum( [minimum(prod_without_nans), minimum(loss_without_nans)] )
        maxx[1] = maximum( [maximum(prod_without_nans), maximum(loss_without_nans)] )
        minx[1] = 10^(floor(log10(minx[1])))
        maxx[1] = 10^(ceil(log10(maxx[1])))

        ax[1].legend(fontsize=10)

        if plot_total_rate_coefs
            ax_1_2 = ax[1].twiny()
            ax_1_2.tick_params(axis="x", labelcolor="xkcd:royal purple")
            ax_1_2.set_ylabel("Rate coefficient (cm^3/s)", color="xkcd:royal purple")
            ax_1_2.semilogx(total_chem_prod_ratecoef, GV.plot_grid, color="xkcd:royal purple", linestyle=":", linewidth=3, label="Total chemical production rate coef", zorder=6)
            ax_1_2.semilogx(total_chem_loss_ratecoef, GV.plot_grid, color="xkcd:lavender", linestyle=":", linewidth=3, label="Total chemical loss rate coef", zorder=6)
            ax_1_2.legend()
        end

    else
        println()
        minx[1] = 0
        maxx[1] = 1
        ax[1].text(0.1, 225, "Chemical production/loss is off for $(sp).")
    end

    # ================================================================================
    # Plot the transport production and loss rates: panel #2

    # Calculate the transport fluxes for the species
    plottitle_ext = "" # no extra info in the plot title if flux==false 
    if sp in GV.transport_species
        transportPL = get_transport_PandL_rate(sp, atmdict; nonthermal=nonthermal, globvars...)
        # now separate into two different arrays for ease of addition.
        production_i = transportPL .>= 0  # boolean array for where transport entries > 0 (production),
        loss_i = transportPL .< 0 # and for where transport entries < 0 (loss).
        total_transport_prod = production_i .* transportPL
        total_transport_loss = loss_i .* abs.(transportPL)

        if sp in [:H2O, :HDO] # Water is turned off in the lower atmosphere, so we should represent that.
            total_transport_prod[1:GV.upper_lower_bdy_i] .= NaN
            total_transport_loss[1:GV.upper_lower_bdy_i] .= NaN
        end

        # set the x lims for transport axis. Special because total_transport_prod, and etc are incomplete arrays.
        minx[2] = 10.0^(floor(log10(minimum(abs.(transportPL)))))
        maxx[2] = 10.0^(ceil(log10(maximum(abs.(transportPL)))))

        # Plot the transport production and loss without the boundary layers
        ax[2].scatter(total_transport_prod, GV.plot_grid, color="red", marker=9, label="Total gain this layer", zorder=4)
        ax[2].scatter(total_transport_loss, GV.plot_grid, color="blue", marker=8, label="Total loss this layer", zorder=4)
        ax[2].set_xscale("log")

        ax[2].legend(fontsize=12)

        plottitle_ext = " by chemistry & transport"
    else
        total_transport_prod = zeros(GV.num_layers)
        total_transport_loss = zeros(GV.num_layers)
        minx[2] = 0
        maxx[2] = 1
        ax[2].text(0.1, 225, "Vertical transport is off for $(sp).")
    end

    # =================================================================================
    # Total production and loss from chemistry and transport: panel #3
    total_prod_rate = total_transport_prod .+ total_chem_prod
    total_loss_rate = total_transport_loss .+ total_chem_loss

    prod_without_nans = filter(x->!isnan(x), total_prod_rate)
    loss_without_nans = filter(x->!isnan(x), total_loss_rate)

    ax[3].semilogx(total_prod_rate, GV.plot_grid, color="black", marker=9, markevery=15, linewidth=2, label="Total production", zorder=3) 
    ax[3].semilogx(total_loss_rate, GV.plot_grid, color="gray", marker=8, markevery=15, linewidth=2, label="Total loss", zorder=3) 

    minx[3] = minimum([minimum(prod_without_nans), minimum(loss_without_nans)])
    maxx[3] = maximum([maximum(prod_without_nans), maximum(loss_without_nans)])
    minx[3] = 10^(floor(log10(minx[3])))
    maxx[3] = 10^(ceil(log10(maxx[3])))

    ax[3].legend(fontsize=12)

    # ---------------------------------------------------------------------------------
    # Plot the timescale of dn/dt if requested
    if plot_timescales==true
        timescale_of_change = atmdict[sp] ./ abs.(total_prod_rate .- total_loss_rate)
        c = "teal"
        ax[total_ax].set_xscale("log")
        ax[total_ax].plot(timescale_of_change, GV.plot_grid, color=c, zorder=5)
    end

    # Final plotting tasks ============================================================
    # check for and correct any ridiculously low limits
    for i in 1:length(minx)
        if minx[i] < 1e-12
            minx[i] = 1e-12
        end
    end

    titles = ["Chemistry", "Transport", "Total chem+trans"]
    notes = ["Chemical production/loss is off for $(sp) below $(Int64(GV.upper_lower_bdy /1e5)) km.", 
             "Vertical transport is off for $(sp) below $(Int64(GV.upper_lower_bdy /1e5)) km.", 
             "Abundance of $(sp) below $(Int64(GV.upper_lower_bdy /1e5)) km is fixed."] 
    if plot_timescales==true
        titles = ["Chemistry", "Transport", "Total chem+trans", "Timescale of dn/dt (s)"]
    end
    dtstr = dt == nothing ? "" : ", $(dt)"
    titlestr = extra_title == "" ? "" : ", "*extra_title
    for i in 1:length(ax)
        ax[i].set_title(titles[i])
        ax[i].set_xlim(minx[i], maxx[i])
        ax[i].set_ylim(0, zmax/1e5)
        if sp in [:H2O, :HDO]
            ax[i].text(2*minx[i], 5, notes[i])
        end
        
    end
    suptitle("Production & loss" * plottitle_ext * ", $(string(sp))" * dtstr * titlestr, fontsize=20)
    ax[1].set_ylabel("Altitude (km)")
    ax[1].set_xlabel("Rate ("*L"cm^{-3}s^{-1})")

    if num==""
        num=extra_title
    end

    
    path_folders = [results_dir[1:end-1], subfolder, plotsfolder, "chem_rates_$(sp)_$(num).png"]
    filter!(e->e≠"", path_folders)  # gets rid of empty names, in case subfolder or plotsfolder hasn't been passed in
    savepathname = join(path_folders, "/")
    
    # Shortcodes as watermarks
    text(0.9, 0.9, GV.hrshortcode, transform=gcf().transFigure, color="dimgrey")
    text(0.9, 0.85, GV.rshortcode, transform=gcf().transFigure, color="dimgrey")

    if showonly
        show()
    else
        savefig(savepathname, bbox_inches="tight", dpi=300)
        close(fig)
    end
end

function plot_temp_prof(Tprof_1; opt="", lbls=["Neutrals", "Ions", "Electrons"], Tprof_2=nothing, Tprof_3=nothing, savepath=nothing, showonly=false, globvars...)
    #=
    Creates a .png image of the tepmeratures plotted by altitude in the atmosphere

    Inputs:
        Tprof_1: an array of neutral temperature by altitude
        Tprof_2: same, but for the ions
        Tprof_3: same but for the electrons
        opt: optional extension to the filename, must start with _ to not look stupid
        savepath: where to save the resulting .png image
        showonly: whether to just show() the figure. If false, must be accompanied by a value for savepath.
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV),  [:alt])

    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18

    fig, ax = subplots(figsize=(4,6))
    plot_bg(ax)

    plot(Tprof_1, GV.alt./1e5, label="Neutrals", color=medgray)

    if Tprof_2 != nothing
        ax.plot(Tprof_2, GV.alt./1e5, label="Ions", color="xkcd:bright orange")
        ax.legend(fontsize=16)
        ax.set_xscale("log")
    end
    if Tprof_3 != nothing
        ax.plot(Tprof_3, GV.alt./1e5, label="Electrons", color="cornflowerblue")
        ax.legend(fontsize=16, loc=(.45,.3))#"center right")
        ax.set_xscale("log")
    end

    # plot the control temps

    # surface
    ax.scatter(Tprof_1[1], 0, marker="o", color=medgray, zorder=10)
    ax.text(Tprof_1[1]+10, 0, #=L"\mathrm{T}_{\mathrm{surface}}"*=#"$(Int64(round(Tprof_1[1], digits=0))) K ")

    # mesosphere
    meso_ind = findfirst(x->x==minimum(Tprof_1), Tprof_1)
    ax.scatter(Tprof_1[meso_ind], GV.alt[meso_ind+5]/1e5, marker="o", color=medgray, zorder=10)
    ax.text(Tprof_1[meso_ind]+5, GV.alt[meso_ind+5]/1e5, #=L"\mathrm{T}_{\mathrm{meso}}"*=#"$(Int64(round(Tprof_1[meso_ind], digits=0))) K ")

    # exosphere
    ax.scatter(Tprof_1[end], 250, marker="o", color=medgray, zorder=10)
    ax.text(Tprof_1[end]*1.05, 240, #=L"\mathrm{T}_{\mathrm{exo}}"*=#"$(Int64(round(Tprof_1[end], digits=0))) K ")
    
    # final labels
    ax.set_ylabel("Altitude [km]")
    ax.set_yticks(collect(0:50:Int64(GV.alt[end]/1e5)))
    ax.set_xlabel("Temperature [K]")
    ax.set_xlim(95, 2e3)
    ax.tick_params(which="both", axis="x", top=true, labeltop=true)

    if showonly==true
        show()
    else
        try
            savefig(savepath*"/temp_profile$(opt).png", bbox_inches="tight", dpi=300) 
        catch
            println("Error: You asked to save the figure but didn't provide a savepath")
        end
    end
end

function plot_water_profile(atmdict, savepath::String; showonly=false, watersat=nothing, H2Oinitf=nothing, prev_profs=nothing, globvars...)  # H2Oinitf, HDOinitf, nH2O, nHDO, 
    #=
    Plots the water profile in mixing ratio and number densities, in two panels.

    Inputs:
        H2Oinitf: initial fraction of H2O in the atmosphere
        HDOinitf: same for HDO
        nH2O: number density of H2O
        nHDO: same for HDO
        savepath: where to place the image
        showonly: show rather than save image
        watersat: optional. must be a list of the saturation fractions with HDO second, 
                  i.e. [H2Osatfrac, HDOsatfrac]
    =#

    GV = values(globvars)
    @assert all(x->x in keys(GV), [:plot_grid, :all_species, :non_bdy_layers, :speciescolor, :speciesstyle])

    rcParams = PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 16
    rcParams["axes.labelsize"]= 18
    rcParams["xtick.labelsize"] = 14
    rcParams["ytick.labelsize"] = 14

    fig, ax = subplots(1, 3, sharey=true, figsize=(14,6))
    for a in ax
        plot_bg(a)
        a.tick_params(axis="x", which="minor", bottom=true, top=true)
    end
    subplots_adjust(wspace=0.12)

    prevcol = "#666"
    
    # mixing ratio axis
    # to get in ppmv, divide the mixing ratio by 1e-6. 
    if prev_profs != nothing
        ax[1].semilogx(prev_profs[1] ./ n_tot(atmdict; globvars...), GV.plot_grid, color=prevcol)
        ax[1].semilogx(prev_profs[2] ./ n_tot(atmdict; globvars...), GV.plot_grid, color=prevcol)
    end
    ax[1].semilogx(atmdict[:H2O] ./ n_tot(atmdict; globvars...), GV.plot_grid, color=GV.speciescolor[:H2O], linewidth=2)
    ax[1].semilogx(atmdict[:HDO] ./ n_tot(atmdict; globvars...), GV.plot_grid, color=GV.speciescolor[:HDO], linestyle=GV.speciesstyle[:HDO], linewidth=2)
    ax[1].set_xlabel("Mixing Ratio")
    ax[1].set_ylabel("Altitude (km)")
    ax[1].set_xticks(collect(logrange(1e-12, 1e-2, 6)))
    
    # Number density 
    if prev_profs != nothing
        ax[2].semilogx(prev_profs[1], GV.plot_grid, color=prevcol)
        ax[2].semilogx(prev_profs[2], GV.plot_grid, color=prevcol, linestyle=GV.speciesstyle[:HDO])
    end
    ax[2].semilogx(atmdict[:H2O], GV.plot_grid, color=GV.speciescolor[:H2O], linewidth=2, label=string_to_latexstr("H2O"))
    ax[2].semilogx(atmdict[:HDO], GV.plot_grid, color=GV.speciescolor[:HDO], linestyle=GV.speciesstyle[:HDO], linewidth=2, label="HDO")
    ax[2].set_xlabel(L"Number density (cm$^{-3}$)")
    ax[2].set_xticks(collect(logrange(1e-4, 1e16, 6)))

    # ppm
    if prev_profs != nothing
        ax[3].semilogx((prev_profs[1] ./ n_tot(atmdict; globvars...)) ./ 1e-6, GV.plot_grid, color=prevcol)
        ax[3].semilogx((prev_profs[2] ./ n_tot(atmdict; globvars...)) ./ 1e-6, GV.plot_grid, color=prevcol, linestyle=GV.speciesstyle[:HDO])
    end
    ax[3].semilogx((atmdict[:H2O] ./ n_tot(atmdict; globvars...)) ./ 1e-6, GV.plot_grid, color=GV.speciescolor[:H2O], linewidth=2)
    ax[3].semilogx((atmdict[:HDO] ./ n_tot(atmdict; globvars...)) ./ 1e-6, GV.plot_grid, color=GV.speciescolor[:HDO], linestyle=GV.speciesstyle[:HDO], linewidth=2)
    ax[3].set_xlabel("ppmv")

    # Title and legend
    suptitle(L"Initial H$_2$O and HDO vertical profiles", y=1.05)

    L2D = PyPlot.matplotlib.lines.Line2D
    if prev_profs != nothing
        lines = [L2D([0], [0], color=GV.speciescolor[:H2O]),
                 L2D([0], [0], color=GV.speciescolor[:HDO], linestyle=GV.speciesstyle[:HDO]),
                 L2D([0], [0], color=prevcol)]
        lbls = ["$(string_to_latexstr("H2O")) (new)", "HDO (new)", "Initial"]
    else  
        lines = [L2D([0], [0], color=GV.speciescolor[:H2O]),
                 L2D([0], [0], color=GV.speciescolor[:HDO], linestyle=GV.speciesstyle[:HDO])]
        lbls = ["$(string_to_latexstr("H2O"))", "HDO"]
    end 
    ax[1].legend(lines, lbls, fontsize=12)
    
    # save it
    if showonly==true
        show()
    else
        savefig(savepath*"/water_profiles.png", dpi=300, bbox_inches="tight")
        close(fig)
    end

    # Will make a plot of the water profile and the saturation vapor pressure curve on the same axis for comparison
    if (watersat != nothing) & (H2Oinitf!=nothing)
        throw("Need to update the code to plot water saturation")
        fig, ax = subplots(figsize=(4,6))
        plot_bg(ax)
        semilogx(convert(Array{Float64}, H2Oinitf), GV.plot_grid, color=col, linewidth=3, label=L"H$_2$O initial fraction")
        semilogx(convert(Array{Float64}, watersat[2:end-1]), GV.plot_grid, color="black", alpha=0.5, linewidth=3, label=L"H$_2$O saturation")
        xlabel("Mixing ratio", fontsize=18)
        ylabel("Altitude [km]", fontsize=18)
        xlim(1e-8, 1)
        ax.set_xticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-0, 1])
        title(L"H$_2$O saturation fraction", fontsize=20)
        ax.tick_params("both",labelsize=16)
        legend()
        if showonly==true 
            show()
        else
            savefig(savepath*"/water_MR_and_saturation.png", dpi=300, bbox_inches="tight")
        end
    end
end

function top_mechanisms(x, sp, atmdict, p_or_r; savepath=nothing, filename_extra="", y0=100, lowerlim=nothing, upperlim=nothing, globvars...) 
    #=
    Reports the top x dominant mechanisms for production or loss of species sp, and shows a plot.

    Input:
        x: number of top mechanisms to print in the plots and tables.
        sp: Species for which to calculate most important mechanisms (symbol)
        atmdict: atmosphere dictionary
        p_or_r: product or reactant
        savepath: somewhere to put the figure
    Output: 
        Plot of the production profiles of the 5 reactions that produce the most flux
        of hot H and D.
    =#
    
    # Collect global variables
    GV = values(globvars)
    @assert all(x->x in keys(GV), [:all_species, :alt, :collision_xsect, :ion_species, :Jratedict, :molmass, :non_bdy_layers, :num_layers,  
                           :n_alt_index, :reaction_network, :Tn, :Ti, :Te, :dz, :zmax])
    
    # String used for various labels
    rxntype = p_or_r == "product" ? "production" : "loss"
    
    # Build an evalutable network
    relevant_reactions = filter_network(sp, "all", p_or_r; GV.reaction_network)
    rc_funcs = Dict([rxn => mk_function(:((Tn, Ti, Te, M) -> $(rxn[3]))) for rxn in relevant_reactions]);
    
    # Atmospheric density
    Mtot = n_tot(atmdict; GV.all_species)
    
    # Reaction strings used for labeling dataframes
    rxn_strings = vec([format_chemistry_string(r[1], r[2]) for r in relevant_reactions])

    
    # Get volume rates by altitude 
    by_alt = volume_rate_wrapper(sp, relevant_reactions, rc_funcs, atmdict, Mtot; globvars...) # array format--by alt
    by_alt_df = DataFrame(by_alt, rxn_strings)
    
    # Get the column value and its sorted equivalent
    column_val = sum(by_alt .* GV.dz, dims=1)
    column_val_df = DataFrame("Rxn"=>rxn_strings, "Value"=>vec(column_val))
    sorted_column_val = sort(column_val_df, [:Value], rev=true)
    
    # Top number of reactions, limit x
    if nrow(sorted_column_val) < x
        L = nrow(sorted_column_val)
    else
        L = x
    end

    println("Top $(L) $(rxntype) reactions sorted by highest column value: $(sorted_column_val[1:L, :])")

    rcParams = PyCall.PyDict(matplotlib."rcParams")
    rcParams["font.sans-serif"] = ["Louis George Caf?"]
    rcParams["font.monospace"] = ["FreeMono"]
    rcParams["font.size"] = 18
    rcParams["axes.labelsize"]= 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    
    
    # PLOT -----------------------------------------------------
    fig, ax = subplots(figsize=(7.5, 5))
    plot_bg(ax)
    ax.set_ylabel("Altitude (km)")
    if lowerlim!=nothing
        ax.set_xlim(left=lowerlim)
    end
    if upperlim!=nothing
        ax.set_xlim(right=upperlim) # H plot limits
    end
    ax.tick_params(which="both", labeltop=true, labelbottom=true, top=true)
    ax.set_ylim(y0, 250)
    ax.set_xscale("log")
    spstr = string_to_latexstr(string(sp))
    ax.set_title(L"Top %$(L) %$(rxntype) mechanisms, %$(spstr) %$(filename_extra)", size=16)
    ax.set_xlabel(L"Rate (cm$^{-3}$ s$^{-1}$)")
    
    # get the reaction strings for the top L reactions of each panel
    top5_rxn_strs = sorted_column_val.Rxn[1:L]
    
    for row in eachrow(sorted_column_val)[1:L]
        ax.plot(by_alt_df[!, row.Rxn], plot_grid, label=string_to_latexstr(row.Rxn), linewidth=2)#, color=thiscol)
    end
    ax.legend(loc=(1.01, 0.5))
    
    if savepath==nothing
        show()
    else
        savefig(savepath*"top$(L)_$(rxntype)_$(sp)$(filename_extra).png", bbox_inches="tight", dpi=300)
    end
end