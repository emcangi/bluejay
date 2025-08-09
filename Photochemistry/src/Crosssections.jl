# **************************************************************************** #
#                                                                              #
#                     PHOTOCHEMISTRY CROSS SECTIONS                            #
#                                                                              #
# This file contains all code which does the work of reading from the input    #
# files which contain cross section info and are provided to the model. Also,  # 
# it includes populate_xsect_dict, which makes the grand dictionary needed to  # 
# run the model. There are probably many inefficiencies here, although it has  #
# been improved over time. They're probably not worth fixing because they only #
# add a small amount of time at the beginning of the model run.                #
#                                                                              # 
# **************************************************************************** #


function populate_xsect_dict(pd_dataf, xsecfolder; ion_xsects=true, globvars...)
    #=
    Creates a dictionary of the 1-nm photodissociation or photoionization
    cross-sections important in the atmosphere. keys are symbols found in
    Jratelist. each entry is an array of arrays, yielding the wavelengths
    and cross-sections for each altitude in the atmosphere.

    NOTE: jspecies refers to the photodissociation or photoionization
    cross section for a particular species which produces a UNIQUE SET OF
    PRODUCTS. In this sense, xsect_dict has already folded in quantum
    efficiency considerations (branching ratios).

    Input
        pd_dataf: Dictionary of filenames associated with a given species for its photodissociation reactions.
        Tn_array: Neutral temperatures for certain photolysis processes.
        ion_xsects: Whether to fill in crosssection information for ion species
        Jrates: the list of Jrates as defined in the PARAMETER file.
    Output
        crosssections: dictionary of cross sections by wavelength for each species. 
    =#

    GV = values(globvars)
    # required = [:Tn, :n_all_layers]
    required = [:Tn, :n_all_layers, :n_horiz]
    check_requirements(keys(GV), required)

    # Set up =======================================================================
    # xsect_dict = Dict{Symbol, Array{Array{Float64}}}()
    xsect_dict = Dict{Symbol, Vector{Vector{Array{Float64}}}}()
    compute_by_T(f) = [ [f(GV.Tn[ihoriz, ialt]) for ialt  in 1:GV.n_all_layers] for ihoriz in 1:GV.n_horiz ]

    # Loading Data =================================================================
    # CO2 photodissociation --------------------------------------------------------
    # temperature-dependent between 195-295K
    co2xdata = readdlm(xsecfolder*pd_dataf[:CO2]["main"],'\t', Float64, comments=true, comment_char='#')

    # H2O & HDO --------------------------------------------------------------------
    h2oxdata = readdlm(xsecfolder*pd_dataf[:H2O]["main"],'\t', Float64, comments=true, comment_char='#')

    # These xsect_dicts for HDO are for 298K.
    hdoxdata = readdlm(xsecfolder*pd_dataf[:HDO]["main"],'\t', Float64, comments=true, comment_char='#')

    # H2O2 + HDO2 ------------------------------------------------------------------
    # the data in the following table cover the range 190-260nm
    h2o2xdata = readdlm(xsecfolder*pd_dataf[:H2O2]["main"],'\t', Float64, comments=true, comment_char='#')
    hdo2xdata = readdlm(xsecfolder*pd_dataf[:HDO2]["main"],'\t', Float64, comments=true, comment_char='#')

    # O3 ---------------------------------------------------------------------------
    # including IR bands which must be resampled from wavenumber
    o3xdata = readdlm(xsecfolder*pd_dataf[:O3]["main"],'\t', Float64, comments=true, comment_char='#')
    global o3ls = o3xdata[:,1]
    global o3xs = o3xdata[:,2]
    o3chapxdata = readdlm(xsecfolder*pd_dataf[:O3]["chapman"],'\t', Float64, comments=true, comment_char='#')
    o3chapxdata[:,1] = map(p->1e7/p, o3chapxdata[:,1])
    for i in [round(Int, floor(minimum(o3chapxdata[:,1]))):round(Int, ceil(maximum(o3chapxdata))-1);]
        posss = getpos(o3chapxdata, x->i<x<i+1)
        dl = diff([map(x->o3chapxdata[x[1],1], posss); i])
        x = map(x->o3chapxdata[x[1],2],posss)
        ax = reduce(+,map(*,x, dl))/reduce(+,dl)
        global o3ls = [o3ls; i+0.5]
        global o3xs = [o3xs; ax]
    end
    o3xdata = reshape([o3ls; o3xs],length(o3ls),2)
    

    # O2 ---------------------------------------------------------------------------
    o2xdata = readdlm(xsecfolder*pd_dataf[:O2]["main"],'\t', Float64, comments=true, comment_char='#')
    o2schr130K = readdlm(xsecfolder*pd_dataf[:O2]["schr_short"],'\t', Float64, comments=true, comment_char='#')
    o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
    o2schr130K = binupO2(o2schr130K)
    o2schr190K = readdlm(xsecfolder*pd_dataf[:O2]["schr_mid"],'\t', Float64, comments=true, comment_char='#')
    o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
    o2schr190K = binupO2(o2schr190K)
    o2schr280K = readdlm(xsecfolder*pd_dataf[:O2]["schr_long"],'\t', Float64, comments=true, comment_char='#')
    o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
    o2schr280K = binupO2(o2schr280K)

    # HO2 & DO2 --------------------------------------------------------------------
    ho2xsect = [190.5:249.5;]
    ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)
    do2xsect = deepcopy(ho2xsect)

    # H2 & HD ----------------------------------------------------------------------
    h2xdata = readdlm(xsecfolder*pd_dataf[:H2]["main"],',',Float64, comments=true, comment_char='#')
    hdxdata = readdlm(xsecfolder*pd_dataf[:HD]["main"],',',Float64, comments=true, comment_char='#')

    # OH & OD ----------------------------------------------------------------------
    ohxdata = readdlm(xsecfolder*pd_dataf[:OH]["main"],',',Float64, comments=true, comment_char='#')
    ohO1Dxdata = readdlm(xsecfolder*pd_dataf[:OH]["O1D+H"],',',Float64, comments=true, comment_char='#')
    odxdata = readdlm(xsecfolder*pd_dataf[:OD]["main"],',',Float64, comments=true, comment_char='#')


    # Populating the dictionary ======================================================
    # CO2 photodissociation -------------------------------------------------------------
    # CO2+hv->CO+O
    # xsect_dict[get_Jrate_symb("CO2", ["CO", "O"])] = map(xs->quantumyield(xs,((l->l>167, 1), (l->95>l, 0.5))), map(t->co2xsect(co2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("CO2", ["CO", "O"])] = compute_by_T(t->quantumyield(co2xsect(co2xdata, t), ((l->l>167, 1), (l->95>l, 0.5))))

    # CO2+hv->CO+O1D
    # xsect_dict[get_Jrate_symb("CO2", ["CO", "O1D"])] = map(xs->quantumyield(xs,((l->95<l<167, 1), (l->l<95, 0.5))), map(t->co2xsect(co2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("CO2", ["CO", "O1D"])] = compute_by_T(t->quantumyield(co2xsect(co2xdata, t), ((l->95<l<167, 1), (l->l<95, 0.5))))

    # CO₂ + hν -> C + O + O; JCO2toCpOpO  
    thisjr = get_Jrate_symb("CO2", ["C", "O", "O"])
    # xsect_dict[thisjr] = fill(readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#'), GV.n_all_layers)
    xsect_dict[thisjr] = [fill(readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#'), GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # CO2 + hν -> C + O₂; JCO2toCpO2
    thisjr = get_Jrate_symb("CO2", ["C", "O2"])
    # xsect_dict[thisjr] = fill(readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#'), GV.n_all_layers)
    xsect_dict[thisjr] = [fill(readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#'), GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # O2 photodissociation ---------------------------------------------------------
    # O2+hv->O+O
    # xsect_dict[get_Jrate_symb("O2", ["O", "O"])] = map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), GV.Tn))
    xsect_dict[get_Jrate_symb("O2", ["O", "O"])] = compute_by_T(t->quantumyield(o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), ((x->x>175, 1),)))

    # O2+hv->O+O1D
    # xsect_dict[get_Jrate_symb("O2", ["O", "O1D"])] = map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), GV.Tn))
    xsect_dict[get_Jrate_symb("O2", ["O", "O1D"])] = compute_by_T(t->quantumyield(o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, t), ((x->x<175, 1),)))

    # O3 photodissociation ---------------------------------------------------------
    # O3+hv->O2+O
    # xsect_dict[get_Jrate_symb("O3", ["O2", "O"])] = map(t->quantumyield(o3xdata,
    #                                                                     (
    #                                                                      (l->l<193, 1-(1.37e-2*193-2.16)),
    #                                                                      (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
    #                                                                      (l->225<=l<306, 0.1),
    #                                                                      (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
    #                                                                      (l->328<=l<340, 0.92),
    #                                                                      (l->340<=l, 1.0)
    #                                                                     )), GV.Tn)
    xsect_dict[get_Jrate_symb("O3", ["O2", "O"])] = compute_by_T(t->quantumyield(o3xdata,
                                                                        (
                                                                         (l->l<193, 1-(1.37e-2*193-2.16)),
                                                                         (l->193<=l<225, l->(1 .- (1.37e-2*l-2.16))),
                                                                         (l->225<=l<306, 0.1),
                                                                         (l->306<=l<328, l->(1 .- O3O1Dquantumyield(l, t))),
                                                                         (l->328<=l<340, 0.92),
                                                                         (l->340<=l, 1.0)
                                                                        )))

    # xsect_dict[get_Jrate_symb("O3", ["O2", "O1D"])] = map(t->quantumyield(o3xdata,
    #                                    (
    #                                     (l->l<193, 1.37e-2*193-2.16),
    #                                     (l->193<=l<225, l->(1.37e-2*l-2.16)),
    #                                     (l->225<=l<306, 0.9),
    #                                     (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
    #                                     (l->328<=l<340, 0.08),
    #                                     (l->340<=l, 0.0)
    #                                    )), GV.Tn)
    xsect_dict[get_Jrate_symb("O3", ["O2", "O1D"])] = compute_by_T(t->quantumyield(o3xdata,
                                       (
                                        (l->l<193, 1.37e-2*193-2.16),
                                        (l->193<=l<225, l->(1.37e-2*l-2.16)),
                                        (l->225<=l<306, 0.9),
                                        (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                                        (l->328<=l<340, 0.08),
                                        (l->340<=l, 0.0)
                                       )))

    # xsect_dict[get_Jrate_symb("O3", ["O", "O", "O"])] =  fill(quantumyield(o3xdata,((x->true, 0.),)),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("O3", ["O", "O", "O"])] =  [fill(quantumyield(o3xdata,((x->true, 0.),)),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # TODO: These are sort of redundant with the loop on line 254. This could be tidied up. 
    # H2 and HD photodissociation --------------------------------------------------
    # H2+hv->H+H
    # xsect_dict[get_Jrate_symb("H2", ["H", "H"])] = fill(h2xdata, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("H2", ["H", "H"])] = [fill(h2xdata, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # HD+hν -> H+D 
    # xsect_dict[get_Jrate_symb("HD", ["H", "D"])] = fill(hdxdata, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("HD", ["H", "D"])] = [fill(hdxdata, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # OH and OD photodissociation --------------------------------------------------
    # OH+hv->O+H
    # xsect_dict[get_Jrate_symb("OH", ["O", "H"])] = fill(ohxdata, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("OH", ["O", "H"])] = [fill(ohxdata, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # OH + hv -> O(¹D) + H
    # xsect_dict[get_Jrate_symb("OH", ["O1D", "H"])] = fill(ohO1Dxdata, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("OH", ["O1D", "H"])] = [fill(ohO1Dxdata, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # OD + hv -> O+D  
    # xsect_dict[get_Jrate_symb("OD", ["O", "D"])] = fill(odxdata, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("OD", ["O", "D"])] = [fill(odxdata, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # OD + hν -> O(¹D) + D 
    # xsect_dict[get_Jrate_symb("OD", ["O1D", "D"])] = fill(ohO1Dxdata, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("OD", ["O1D", "D"])] = [fill(ohO1Dxdata, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # HO2 and DO2 photodissociation ------------------------------------------------
    # HO2 + hν -> OH + O
    # xsect_dict[get_Jrate_symb("HO2", ["OH", "O"])] = fill(ho2xsect, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("HO2", ["OH", "O"])] = [fill(ho2xsect, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # DO2 + hν -> OD + O
    # xsect_dict[get_Jrate_symb("DO2", ["OD", "O"])] = fill(do2xsect, GV.n_all_layers)
    xsect_dict[get_Jrate_symb("DO2", ["OD", "O"])] = [fill(do2xsect, GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # H2O and HDO photodissociation ------------------------------------------------
    # H2O+hv->H+OH
    # xsect_dict[get_Jrate_symb("H2O", ["H", "OH"])] = fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("H2O", ["H", "OH"])] = [fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # H2O+hv->H2+O1D
    # xsect_dict[ get_Jrate_symb("H2O", ["H2", "O1D"])] = fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),GV.n_all_layers)
    xsect_dict[ get_Jrate_symb("H2O", ["H2", "O1D"])] = [fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # H2O+hv->H+H+O
    # xsect_dict[get_Jrate_symb("H2O", ["H", "H", "O"])] = fill(quantumyield(h2oxdata,((x->true, 0),)),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("H2O", ["H", "H", "O"])] = [fill(quantumyield(h2oxdata,((x->true, 0),)),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # HDO + hν -> H + OD
    # xsect_dict[get_Jrate_symb("HDO", ["H", "OD"])] = fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("HDO", ["H", "OD"])] = [fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # HDO + hν -> D + OH
    # xsect_dict[get_Jrate_symb("HDO", ["D", "OH"])] = fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("HDO", ["D", "OH"])] = [fill(quantumyield(hdoxdata,((x->x<145, 0.5*0.89),(x->x>145, 0.5*1))),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # HDO + hν -> HD + O1D
    # xsect_dict[get_Jrate_symb("HDO", ["HD", "O1D"])] = fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("HDO", ["HD", "O1D"])] = [fill(quantumyield(hdoxdata,((x->x<145, 0.11),(x->x>145, 0))),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # HDO + hν -> H + D + O
    # xsect_dict[get_Jrate_symb("HDO", ["H", "D", "O"])] = fill(quantumyield(hdoxdata,((x->true, 0),)),GV.n_all_layers)
    xsect_dict[get_Jrate_symb("HDO", ["H", "D", "O"])] = [fill(quantumyield(hdoxdata,((x->true, 0),)),GV.n_all_layers) for ihoriz in 1:GV.n_horiz]

    # H2O2 and HDO2 photodissociation ----------------------------------------------
    # H2O2+hν->OH+OH
    # xsect_dict[get_Jrate_symb("H2O2", ["OH", "OH"])] = map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))), map(t->h2o2xsect(h2o2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("H2O2", ["OH", "OH"])] = compute_by_T(t->quantumyield(h2o2xsect(h2o2xdata, t), ((x->x<230, 0.85),(x->x>230, 1))))

    # H2O2+hv->HO2+H
    # xsect_dict[get_Jrate_symb("H2O2", ["HO2", "H"])] = map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))), map(t->h2o2xsect(h2o2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("H2O2", ["HO2", "H"])] = compute_by_T(t->quantumyield(h2o2xsect(h2o2xdata, t), ((x->x<230, 0.15),(x->x>230, 0))))

    # H2O2+hv->H2O+O1D
    # xsect_dict[get_Jrate_symb("H2O2", ["H2O", "O1D"])] = map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(h2o2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("H2O2", ["H2O", "O1D"])] = compute_by_T(t->quantumyield(h2o2xsect(h2o2xdata, t), ((x->true, 0),)))

    # HDO2 + hν -> OH + OD
    # xsect_dict[get_Jrate_symb("HDO2", ["OH", "OD"])] = map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))), map(t->hdo2xsect(hdo2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("HDO2", ["OH", "OD"])] = compute_by_T(t->quantumyield(hdo2xsect(hdo2xdata, t), ((x->x<230, 0.85),(x->x>230, 1))))

    # HDO2 + hν-> DO2 + H
    # xsect_dict[get_Jrate_symb("HDO2", ["DO2", "H"])] = map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))), map(t->hdo2xsect(hdo2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("HDO2", ["DO2", "H"])] = compute_by_T(t->quantumyield(hdo2xsect(hdo2xdata, t), ((x->x<230, 0.5*0.15),(x->x>230, 0))))

    # HDO2 + hν-> HO2 + D
    # xsect_dict[get_Jrate_symb("HDO2", ["HO2", "D"])] = map(xs->quantumyield(xs,((x->x<230, 0.5*0.15),(x->x>230, 0))), map(t->hdo2xsect(hdo2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("HDO2", ["HO2", "D"])] = compute_by_T(t->quantumyield(hdo2xsect(hdo2xdata, t), ((x->x<230, 0.5*0.15),(x->x>230, 0))))

    # HDO2 + hν -> HDO + O1D
    # xsect_dict[get_Jrate_symb("HDO2", ["HDO", "O1D"])] = map(xs->quantumyield(xs,((x->true, 0),)), map(t->hdo2xsect(hdo2xdata, t), GV.Tn))
    xsect_dict[get_Jrate_symb("HDO2", ["HDO", "O1D"])] = compute_by_T(t->quantumyield(hdo2xsect(hdo2xdata, t), ((x->true, 0),)))

    # The following reactions have associated files listing cross sections.
    # TODO: Ideally this should not be hard coded in; ought to be passed in based on the reaction spreadsheet.
    reactant_product_sets = Dict("CO2"=>[["CO2pl"], ["CO2plpl"], ["Cplpl", "O2"], ["Cpl", "O2"], ["COpl", "Opl"], ["COpl", "O"], ["Opl", "CO"], ["Opl", "Cpl", "O"], ["C", "O", "O"], ["C", "O2"]],
                                 "CO"=>[["COpl"],  ["C", "Opl"],  ["O", "Cpl"], ["C", "O"]],
                                 "H2O"=>[["H2Opl"], ["Opl", "H2"], ["Hpl", "OH"], ["OHpl", "H"]],
                                 "HDO"=>[["HDOpl"], ["Opl", "HD"], ["Hpl", "OD"],  ["Dpl", "OH"], ["OHpl", "D"],  ["ODpl", "H"]],
                                 "N2"=>[["N2pl"], ["Npl", "N"]],
                                 "NO2"=>[["NO2pl"], ["NO", "O"]],
                                 "NO"=>[["NOpl"], ["N", "O"]],
                                 "N2O"=>[["N2Opl"], ["N2", "O1D"]],
                                 "H"=>[["Hpl"]],
                                 "D"=>[["Dpl"]],
                                 "H2"=>[["H2pl"], ["Hpl", "H"]],
                                 "HD"=>[["HDpl"], ["Dpl", "H"], ["Hpl", "D"]],
                                 "H2O2"=>[["H2O2pl"]],
                                 "HDO2"=>[["HDO2pl"]],
                                 "O"=>[["Opl"]],
                                 "O2"=>[["O2pl"]],
                                 "O3"=>[["O3pl"]],

                                # Cl containing species
                                "HCl"=>[["H","Cl"]],
                                "DCl"=>[["D","Cl"]],
                                #S containing species
                                "SO2"=>[["SO","O"], ["S","O2"]],
                                "SO3"=>[["SO2","O"]], #there may be more possible products than what I have listed, but JPL didn't give a recommended quantum yield
                                "H2SO4"=>[["SO3","H2O"]], # same issue as SO3, JPL did not give a recommended quantum yield
                                "HDSO4"=>[["SO3","HDO"]]) # uses mass scaling of H2SO4
    
        

    for r in keys(reactant_product_sets)
        for ps in reactant_product_sets[r]
            thisjr = get_Jrate_symb(r, ps)
            # xsect_dict[thisjr] = fill(readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#'), GV.n_all_layers)
            xsect_dict[thisjr] = [fill(readdlm(xsecfolder*"$(thisjr).csv",',',Float64, comments=true, comment_char='#'), GV.n_all_layers) for ihoriz in 1:GV.n_horiz]
        end
    end

    # Validate that cross sections have shape (n_horiz, n_all_layers)
    for (jr, cols) in xsect_dict
        @assert length(cols) == GV.n_horiz "cross section column count mismatch for $(jr)"
        @assert all(length(profile) == GV.n_all_layers for profile in cols) "cross section altitude count mismatch for $(jr)"
    end

    return xsect_dict
end

function binupO2(list)
    #=
    Mike originally wrote this
    =#
    ret = Float64[];
    for i in [176:203;]
        posss = getpos(list[:,1],x->i<x<i+1)
        dl = diff([map(x->list[x[1],1],posss); i])
        x0 = map(x->list[x[1],2], posss)
        x1 = map(x->list[x[1],3], posss)
        x2 = map(x->list[x[1],4], posss)
        ax0 = reduce(+,map(*,x0, dl))/reduce(+, dl)
        ax1 = reduce(+,map(*,x1, dl))/reduce(+, dl)
        ax2 = reduce(+,map(*,x2, dl))/reduce(+, dl)
        append!(ret,[i+0.5, ax0, ax1, ax2])
    end
    return transpose(reshape(ret, 4, 203-176+1))
end

function co2xsect(co2xdata, T)
    #=
    Makes an array of CO2 cross sections at temperature T, of format 
    [wavelength in nm, xsect]. 

    Data are available at 195K and 295K, so crosssections at temperatures between 
    these are created by first calculating the fraction that T is along the 
    line from 195 to 295, and then calculating a sort of weighted average of 
    the crosssections at 195 and 295K based on that information. 
    =#
    clamp(T, 195, 295)
    Tfrac = (T-195)/(295-195)

    arr = [co2xdata[:,1]; (1-Tfrac)*co2xdata[:,2] + Tfrac*co2xdata[:,3]]
    reshape(arr, length(co2xdata[:,1]),2)
end

function h2o2xsect_l(l, T) 
    #=
    from 260-350 the following analytic calculation fitting the
    temperature dependence is recommended by Sander 2011.

    Analytic calculation of H2O2 cross section using temperature dependencies
    l: wavelength in nm
    T: temperature in K
    =#
    l = clamp(l, 260, 350)
    T = clamp(T, 200, 400)

    A = [64761., -921.70972, 4.535649,
         -0.0044589016, -0.00004035101,
         1.6878206e-7, -2.652014e-10, 1.5534675e-13]
    B = [6812.3, -51.351, 0.11522, -0.000030493, -1.0924e-7]

    lpowA = map(n->l^n,[0:7;])
    lpowB = map(n->l^n,[0:4;])

    expfac = 1.0/(1+exp(-1265/T))

    return 1e-21*(expfac*reduce(+, map(*, A, lpowA))+(1-expfac)*reduce(+, map(*, B, lpowB)))
end

function h2o2xsect(h2o2xdata, T) 
    #=
    stitches together H2O2 cross sections, some from Sander 2011 table and some
    from the analytical calculation recommended for 260-350nm recommended by the
    same.
    T: temperature in K
    =#
    retl = h2o2xdata[:,1]
    retx = 1e4*h2o2xdata[:,2] # factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    return reshape([retl; retx], length(retl), 2)
end

function hdo2xsect(hdo2xdata, T) 
    #=
    Currently this is a direct copy of h2o2xsect because no HDO2 crosssections
    exist. In the future this could be expanded if anyone ever makes those
    measurements.
    =#
    retl = hdo2xdata[:,1]
    retx = 1e4*hdo2xdata[:,2] # factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    reshape([retl; retx], length(retl), 2)
end

function ho2xsect_l(l) 
    #= 
    compute HO2 cross-section as a function of wavelength l in nm, as given by 
    Sander 2011 JPL Compilation 
    =#
    a = 4.91
    b = 30612.0
    sigmamed = 1.64e-18
    vmed = 50260.0
    v = 1e7/l;
    if 190<=l<=250
        return HO2absx = sigmamed / ( 1 - b/v ) * exp( -a * log( (v-b)/(vmed-b) )^2 )
    else
        return 0.0
    end
end

function o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, T)
    #=
    For O2, used in quantum yield calculations,
    including temperature-dependent Schumann-Runge bands.

    ...details missing
    =# 
    
    o2x = deepcopy(o2xdata);
    # fill in the schumann-runge bands according to Minschwaner 1992
    T = clamp(T, 130, 500)
    if 130<=T<190
        o2schr = o2schr130K
    elseif 190<=T<280
        o2schr = o2schr190K
    else
        o2schr = o2schr280K
    end

    del = ((T-100)/10)^2

    for i in [176.5:203.5;]
        posO2 = something(findfirst(isequal(i), o2x[:, 1]), 0)
        posschr = something(findfirst(isequal(i), o2schr[:, 1]), 0)
        o2x[posO2, 2] += 1e-20*(o2schr[posschr, 2]*del^2
                                + o2schr[posschr, 3]*del
                                + o2schr[posschr, 4])
    end

    # add in the herzberg continuum (though tiny)
    # measured by yoshino 1992
    for l in [192.5:244.5;]
        posO2 = something(findfirst(isequal(l), o2x[:, 1]), 0)
        o2x[posO2, 2] += 1e-24*(-2.3837947e4
                            +4.1973085e2*l
                            -2.7640139e0*l^2
                            +8.0723193e-3*l^3
                            -8.8255447e-6*l^4)
    end

    return o2x
end

function O3O1Dquantumyield(lambda, temp) 
    #=
    Quantum yield for O(1D) from O3 photodissociation. 
    "The quantum yield of O1D from ozone photolysis is actually well-studied! 
    This adds some complications for processing." - Mike
    =#
    if lambda < 306. || lambda > 328.
        return 0.
    end
    temp=clamp(temp, 200, 320)#expression is only valid in this T range

    X = [304.225, 314.957, 310.737];
    w = [5.576, 6.601, 2.187];
    A = [0.8036, 8.9061, 0.1192];
    v = [0.,825.518];
    c = 0.0765;
    R = 0.695;
    q = exp.(-v/(R*temp))
    qrat = q[1]/(q[1]+q[2])

    (q[1]/sum(q)*A[1]*exp.(-((X[1]-lambda)/w[1])^4.)
     +q[2]/sum(q)*A[2]*(temp/300.)^2 .* exp.(-((X[2]-lambda)/w[2])^2.)
     +A[3]*(temp/300.)^1.5*exp.(-((X[3]-lambda)/w[3])^2.)
     +c)
end

function padtosolar(solarflux, crosssection::Array{Float64, 2}) 
    #=
    a function to take an Nx2 array crosssection and pad it with zeroes until it's the
    same length as the solarflux array. Returns the cross sections only, as
    the wavelengths are shared by solarflux 
    =#
    positions = map(x->something(findfirst(isequal(x), solarflux[:,1]), 0), crosssection[:,1])
    retxsec = fill(0.,length(solarflux[:,1]))
    retxsec[positions] = crosssection[:,2]
    return retxsec
end

function quantumyield(xsect, arr)
    #= 
    function to assemble cross-sections for a given pathway. 

    xsect: an Nx2 array, format [wavelength, crosssection], where N is 
           the number of wavelengths available.

    arr: a tuple of tuples with a condition and a quantum yield multiplicative 
         factor, either constant or a function of wavelength in the given regime. 
         Return is an array with all of the matching wavelengths and the scaled cross-sections.
    =#
    lambdas = Float64[];
    rxs = Float64[];
    for (cond, qeff) in arr
        places = findall(cond, xsect[:,1])  # locate wavelengths that match the cond condition
        append!(lambdas, xsect[places, 1])  # put those wavelengths into a new list
        # if we have a number then map to a function
        isa(qeff, Function) ? (qefffn = qeff) : (qefffn = x->qeff)
        append!(rxs, map(*,map(qefffn, xsect[places, 1]),xsect[places, 2]))
    end

    return reshape([lambdas; rxs],length(lambdas),2)
end