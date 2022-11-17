#Eq. (22)
function constraint_cc_filter_voltage_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vk_s  = [_PM.var(pm, n, :vk_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vk_s, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vk_s, mop) <= vmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(vk_s, T2)
                                <=
                               ((_PCE.mean(vk_s, mop) - vmin^2) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vk_s, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vk_s, mop)) / λmax)^2
                    )
end

#Eq. (23)
function constraint_cc_converter_voltage_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vc_s  = [_PM.var(pm, n, :vc_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vc_s, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vc_s, mop) <= vmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(vc_s, T2)
                                <=
                               ((_PCE.mean(vc_s, mop) - vmin^2) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vc_s, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vc_s, mop)) / λmax)^2
                    )
end

#Eq. (33)
function constraint_cc_transformer_current_from_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    iik_s  = [_PM.var(pm, n, :iik_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iik_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iik_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iik_s, mop)) / λmax)^2
                    )
end

#Eq. (34)
function constraint_cc_transformer_current_to_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    iki_s  = [_PM.var(pm, n, :iki_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iki_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iki_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iki_s, mop)) / λmax)^2
                    )
end

#Eq. (35)
function constraint_cc_reactor_current_from_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    ikc_s  = [_PM.var(pm, n, :ikc_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(ikc_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(ikc_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(ikc_s, mop)) / λmax)^2
                    )
end

#Eq. (36)
function constraint_cc_reactor_current_to_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    ick_s  = [_PM.var(pm, n, :ick_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(ick_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(ick_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(ick_s, mop)) / λmax)^2
                    )
end

#=
#Eq. (40)
function constraint_cc_converter_current_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    ic_s  = [_PM.var(pm, n, :ic_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(ic_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(ic_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(ic_s, mop)) / λmax)^2
                    )
end
=#

#Eq. (40)
function constraint_cc_converter_current_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    iconv_lin_s  = [_PM.var(pm, n, :iconv_lin_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_lin_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iconv_lin_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iconv_lin_s, mop)) / λmax)^2
                    )
end


function constraint_cc_dc_branch_current(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop)
    i_dc_fr = [_PM.var(pm, n, :igrid_dc, f_idx) for n in sorted_nw_ids(pm)]
    i_dc_to = [_PM.var(pm, n, :igrid_dc, t_idx) for n in sorted_nw_ids(pm)]

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(i_dc_fr, mop) <= Imax)

    JuMP.@constraint(pm.model, Imin <= _PCE.mean(i_dc_fr, mop))

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr, T2)
                               <=
                                ((Imax - _PCE.mean(i_dc_fr, mop)) / λmax)^2
                    )

    
    JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr, T2)
                                <=
                               ((_PCE.mean(i_dc_fr, mop) - Imin) / λmin)^2
    )
    

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(i_dc_to, mop) <= Imax)
    JuMP.@constraint(pm.model, Imin <= _PCE.mean(i_dc_to, mop))

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(i_dc_to, T2)
                               <=
                                ((Imax - _PCE.mean(i_dc_to, mop)) / λmax)^2
                    )
    
    JuMP.@constraint(pm.model,  _PCE.var(i_dc_to, T2)
                                <=
                               ((_PCE.mean(i_dc_to, mop) - Imin) / λmin)^2
    )
    

end

function constraint_cc_iconv_lin_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    iconv_lin_s  = [_PM.var(pm, n, :iconv_lin_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_lin_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iconv_lin_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iconv_lin_s, mop)) / λmax)^2
                    )
end

function constraint_cc_iconv_lin(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, T2, mop)
    iconv_lin  = [_PM.var(pm, n, :iconv_lin, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_lin, mop) <= Imax)
    JuMP.@constraint(pm.model, Imin <= _PCE.mean(iconv_lin, mop))
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iconv_lin, T2)
                               <=
                                ((Imax - _PCE.mean(iconv_lin, mop)) / λmax)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(iconv_lin, T2)
                                <=
                                ((_PCE.mean(iconv_lin, mop) - Imin) / λmin)^2
                    )
end


function constraint_cc_conv_ac_power_real(pm, i, Pacmin, Pacmax, λmin, λmax, T2, mop)
    pconv_ac  = [_PM.var(pm, n, :pconv_ac, i) for n in sorted_nw_ids(pm)]


    # bounds on the expectation
    JuMP.@constraint(pm.model, Pacmin <= _PCE.mean(pconv_ac, mop))
    JuMP.@constraint(pm.model, _PCE.mean(pconv_ac, mop) <= Pacmax)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(pconv_ac, T2)
                            <=
                            ((_PCE.mean(pconv_ac, mop) - Pacmin) / λmin)^2
                    )

    JuMP.@constraint(pm.model,  _PCE.var(pconv_ac, T2)
                            <=
                            ((Pacmax - _PCE.mean(pconv_ac, mop)) / λmax)^2
                    )

end

function constraint_cc_conv_ac_power_imaginary(pm, i, Qacmin, Qacmax, λmin, λmax, T2, mop)
    qconv_ac  = [_PM.var(pm, n, :qconv_ac, i) for n in sorted_nw_ids(pm)]



    # bounds on the expectation
    JuMP.@constraint(pm.model, Qacmin <= _PCE.mean(qconv_ac, mop))
    JuMP.@constraint(pm.model, _PCE.mean(qconv_ac, mop) <= Qacmax)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(qconv_ac, T2)
                            <=
                            ((_PCE.mean(qconv_ac, mop) - Qacmin) / λmin)^2
                    )

    JuMP.@constraint(pm.model,  _PCE.var(qconv_ac, T2)
                            <=
                            ((Qacmax - _PCE.mean(qconv_ac, mop)) / λmax)^2
                    )

end

function constraint_cc_conv_dc_power(pm, i, Pdcmin, Pdcmax, λmin, λmax, T2, mop)
    pconv_dc  = [_PM.var(pm, n, :pconv_dc, i) for n in sorted_nw_ids(pm)]


    # bounds on the expectation
    JuMP.@constraint(pm.model, Pdcmin <= _PCE.mean(pconv_dc, mop))
    JuMP.@constraint(pm.model, _PCE.mean(pconv_dc, mop) <= Pdcmax)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(pconv_dc, T2)
                            <=
                            ((_PCE.mean(pconv_dc, mop) - Pdcmin) / λmin)^2
                    )

    JuMP.@constraint(pm.model,  _PCE.var(pconv_dc, T2)
                            <=
                            ((Pdcmax - _PCE.mean(pconv_dc, mop)) / λmax)^2
                    )

end




function constraint_cc_converter_dc_current(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, T2, mop)
    iconv_dc = [_PM.var(pm, n, :iconv_dc, i) for n in sorted_nw_ids(pm)]

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_dc, mop) <= Imax)

    JuMP.@constraint(pm.model, Imin <= _PCE.mean(iconv_dc, mop))

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iconv_dc, T2)
                               <=
                                ((Imax - _PCE.mean(iconv_dc, mop)) / λmax)^2
                    )

    
    JuMP.@constraint(pm.model,  _PCE.var(iconv_dc, T2)
                                <=
                               ((_PCE.mean(iconv_dc, mop) - Imin) / λmin)^2
    )
    
   

end

function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vms  = [_PM.var(pm, n, :vms, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vms, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vms, mop) <= vmax^2)
    # chance constraint bounds
    
    JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                                <=
                               ((_PCE.mean(vms, mop) - vmin^2) / λmin)^2
                    )
    
    JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vms, mop)) / λmax)^2
                    )
end

function constraint_cc_conv_voltage_magnitude(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    
    vdcm  = [_PM.var(pm, n, :vdcm, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin <= _PCE.mean(vdcm, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vdcm, mop) <= vmax)
    # chance constraint bounds
    
    JuMP.@constraint(pm.model,  _PCE.var(vdcm, T2)
                                <=
                               ((_PCE.mean(vdcm, mop) - vmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vdcm, T2)
                               <=
                                ((vmax - _PCE.mean(vdcm, mop)) / λmax)^2
                    )

    JuMP.@constraint(pm.model,  _PCE.var(vdcm, T2)
                                <=
                                (0.00001)
                    )
    

end


function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractACRModel, b, cmax, λcmax, T2, mop)
    cmss = [_PM.var(pm, nw, :cmss, b) for nw in sorted_nw_ids(pm)]

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= cmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                <=
                                ((cmax^2 - _PCE.mean(cmss,mop)) / λcmax)^2
                    )
end


function constraint_cc_gen_power_real(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop)
    pg  = [_PM.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

     # bounds on the expectation 
     JuMP.@constraint(pm.model,  pmin <= _PCE.mean(pg, mop))
     JuMP.@constraint(pm.model,  _PCE.mean(pg, mop) <= pmax)
     # chance constraint bounds
     
     
     JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                 <=
                                ((_PCE.mean(pg, mop) - pmin) / λmin)^2
                   )

    
     JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                 <=
                                 ((pmax - _PCE.mean(pg, mop)) / λmax)^2
                   )
end

function constraint_cc_gen_power_imaginary(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PM.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]

    # bounds on the expectation 
    JuMP.@constraint(pm.model,  qmin <= _PCE.mean(qg, mop))
    JuMP.@constraint(pm.model,  _PCE.mean(qg, mop) <= qmax)
    # chance constraint bounds
    
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qmin) / λmin)^2
                    )
    
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                               <=
                                ((qmax - _PCE.mean(qg,mop)) / λmax)^2
                    )
end