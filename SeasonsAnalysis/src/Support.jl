# **************************************************************************** #
#                                                                              #
#                               SUPPORT FUNCTIONS                              #
#                                                                              #
# **************************************************************************** #

function esc_pct_of_total_col(flux, col_without_flux)
    #=
    Given an escape flux and column total (not including the flux), 
    this calculates the percent of the total column that the escape makes up. 
    =#
    return 100 .* (flux ./ (col_without_flux + flux))
end

function cyclical_color_array(;N=5, map="coolwarm")
    #=
    Used to generate N colors from map, such that they form a sort of intuitive cycle for plotting purposes.
    =#
    colors = get_colors(N, map; strt=0, stp=1)

    colors[1:2, :] = get_colors(2, map; strt=0.5, stp=1)
    colors[3:4, :] = get_colors(2, map; strt=0.5, stp=0)
    colors[5, :] = get_colors(1, map; strt=0.5, stp=0.5)

    for j in 1:2:5
        colors[j, :] .-= 0.05*j
    end
    return colors
end

function parent_folders_from_full_path(list_of_full_paths)
    #=
    Given a list of full file pathnames, return their parent directories.
    Using `dirname` avoids reliance on fragile regular expressions.
    =#
    return [dirname(f) * "/" for f in list_of_full_paths]
end

function DH_of_escaping(df, thefiles; t="")
    
    SMOW = 1.6e-4

    thefolders = parent_folders_from_full_path(thefiles)
    all_esc_df = get_escape_for_all_atmospheres(thefolders, thefiles)

    
    DH_t = (df_lookup(df, "EscapeType", "Thermal", "D") / df_lookup(df, "EscapeType", "Thermal", "H"))[1]
    DH_n = (df_lookup(df, "EscapeType", "Nonthermal", "D") / df_lookup(df, "EscapeType", "Nonthermal", "H"))[1]
    DH_tn = (df_lookup(df, "EscapeType", "Total", "D") / df_lookup(df, "EscapeType", "Total", "H"))[1]
    
    fig, ax = subplots()
    plot_bg(ax)
    bar([1, 2, 3], [DH_t, DH_n, DH_tn], zorder=3, width=0.5)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["Thermal escape","Nonthermal escape","Total escape"])
    ax.set_ylabel("D/H of escaping atoms")
    yscale("log")
    
    axhline(0.01*SMOW, color="gray")
    axhline(0.1*SMOW, color="gray")
    axhline(SMOW, color="gray")
    axhline(5*SMOW, color="gray")
    axhline(10*SMOW, color="gray")
    
    ax.text(0.9, 0.01*SMOW*1.1, "0.01xSMOW")
    ax.text(0.9, 0.1*SMOW*1.1, "0.1xSMOW")
    ax.text(0.9, SMOW*1.1, "SMOW")
    ax.text(0.9, 5*SMOW*1.1, "5xSMOW")
    ax.text(0.9, 10*SMOW*1.1, "10xSMOW")
    
    ax.set_title(t)
    
    show()
end

function draw_DH_lines(ax, multipliers; which="horizontal", smowcol="gray", DHfs=16)
    smow = 1.6e-4
    for m in multipliers
        if which=="horizontal"
            ax.axhline(m*smow, color=smowcol, linestyle=":", linewidth=1)
            ax.text(1, m*smow, "$(m)x VSMOW", color=smowcol, fontsize=14)
        else
            ax.axvline(m*smow, color=smowcol, linestyle=":", linewidth=1)
            ax.text(m*smow, 1, "$(m)x VSMOW", color=smowcol, fontsize=DHfs, rotation=90, ha="right")
        end
    end
end

function generate_indvar_vs_time_array(filelist; val_order=[225 275 225 175 225])
    # This identifies which files point to a temperature change 
    change = zeros(size(filelist))
    for (i,f) in enumerate(filelist)
        if f[end-4:end-3] == "01"
            change[i] = 1
        end
    end
    
    change_indices = []
    for (i,t) in enumerate(change)
        if t==1
            append!(change_indices, i)
        end
    end
    
    # make array of temps for plotting
    the_array = fill(val_order[1], length(filelist))

    c = 0
    for i in 1:length(filelist)
        if i in change_indices
            c += 1
        end
        if i > 1
            the_array[i] = val_order[c]
        end
    end
    
    return the_array, change_indices
end
