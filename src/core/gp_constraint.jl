################################################################################
# Copyright 2023, Kaan Yurtseven                                               #
################################################################################
# StochasticPowerModelsACDC.jl                                                 #
# An extention package of PowerModels.jl and StochasticPowerModels.jl for      #
#                                 Stochastic Optimal Power Flow in AC/DC grids #
# See https://github.com/kaanyurtseven/StochasticPowerModelsACDC               #
################################################################################

function constraint_gp_ohms_dc_branch(pm::AbstractACRModel, n::Int, i, T2, T3, f_bus, t_bus, f_idx, t_idx, r, p)

    p_fr  = _PM.var(pm, n,  :p_dcgrid, f_idx)
    p_to  = _PM.var(pm, n,  :p_dcgrid, t_idx)

    #Dict(nw => _PM.var(pm, nw, :vk_r, i) for nw in _PM.nw_ids(pm))
    vmdc_fr = Dict(nw => _PM.var(pm, nw, :vdcm, f_bus) for nw in _PM.nw_ids(pm))
    vmdc_to = Dict(nw => _PM.var(pm, nw, :vdcm, t_bus) for nw in _PM.nw_ids(pm))
    i_dc_fr = Dict(nw => _PM.var(pm, nw, :igrid_dc, f_idx) for nw in _PM.nw_ids(pm))
    i_dc_to = Dict(nw => _PM.var(pm, nw, :igrid_dc, t_idx) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * p_fr ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vmdc_fr[n1] * i_dc_fr[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )


    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * p_to ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) *
                                                    (vmdc_to[n1] * i_dc_to[n2])
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )


end


function constraint_gp_filter_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    vk_s = _PM.var(pm, n, :vk_s, i)
    vk_r  = Dict(nw => _PM.var(pm, nw, :vk_r, i) for nw in _PM.nw_ids(pm))
    vk_i  = Dict(nw => _PM.var(pm, nw, :vk_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vk_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vk_r[n1] * vk_r[n2] + vk_i[n1] * vk_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_converter_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    vc_s = _PM.var(pm, n, :vc_s, i)
    vc_r  = Dict(nw => _PM.var(pm, nw, :vc_r, i) for nw in _PM.nw_ids(pm))
    vc_i  = Dict(nw => _PM.var(pm, nw, :vc_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vc_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vc_r[n1] * vc_r[n2] + vc_i[n1] * vc_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_transformer_current_from_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    iik_s = _PM.var(pm, n, :iik_s, i)
    iik_r  = Dict(nw => _PM.var(pm, nw, :iik_r, i) for nw in _PM.nw_ids(pm))
    iik_i  = Dict(nw => _PM.var(pm, nw, :iik_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * iik_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (iik_r[n1] * iik_r[n2] + iik_i[n1] * iik_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_transformer_current_to_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    iki_s = _PM.var(pm, n, :iki_s, i)
    iki_r  = Dict(nw => _PM.var(pm, nw, :iki_r, i) for nw in _PM.nw_ids(pm))
    iki_i  = Dict(nw => _PM.var(pm, nw, :iki_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * iki_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (iki_r[n1] * iki_r[n2] + iki_i[n1] * iki_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_reactor_current_from_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    ikc_s = _PM.var(pm, n, :ikc_s, i)
    ikc_r  = Dict(nw => _PM.var(pm, nw, :ikc_r, i) for nw in _PM.nw_ids(pm))
    ikc_i  = Dict(nw => _PM.var(pm, nw, :ikc_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * ikc_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (ikc_r[n1] * ikc_r[n2] + ikc_i[n1] * ikc_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_reactor_current_to_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    ick_s = _PM.var(pm, n, :ick_s, i)
    ick_r  = Dict(nw => _PM.var(pm, nw, :ick_r, i) for nw in _PM.nw_ids(pm))
    ick_i  = Dict(nw => _PM.var(pm, nw, :ick_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * ick_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (ick_r[n1] * ick_r[n2] + ick_i[n1] * ick_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_converter_current_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    ic_s = _PM.var(pm, n, :ic_s, i)
    ic_r  = Dict(nw => _PM.var(pm, nw, :ic_r, i) for nw in _PM.nw_ids(pm))
    ic_i  = Dict(nw => _PM.var(pm, nw, :ic_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * ic_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (ic_r[n1] * ic_r[n2] + ic_i[n1] * ic_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_iconv_lin_squared_1(pm::AbstractACRModel, n::Int, i, T2, T3)

    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    ic_r  = Dict(nw => _PM.var(pm, nw, :ic_r, i) for nw in _PM.nw_ids(pm))
    ic_i  = Dict(nw => _PM.var(pm, nw, :ic_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * iconv_lin_s
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (ic_r[n1] * ic_r[n2] + ic_i[n1] * ic_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end



function constraint_gp_converter_dc_power(pm::AbstractACRModel, n::Int, i, T2, T3, b_idx)

    pconv_dc = _PM.var(pm, n, :pconv_dc, i)
    iconv_dc  = Dict(nw => _PM.var(pm, nw, :iconv_dc, i) for nw in _PM.nw_ids(pm))
    vdcm  = Dict(nw => _PM.var(pm, nw, :vdcm, b_idx) for nw in _PM.nw_ids(pm))

    #Eq. (43)
    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * pconv_dc ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vdcm[n1] * iconv_dc[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

end

function constraint_gp_converter_ac_power(pm::_PM.AbstractIVRModel, n::Int, i::Int, T2, T3)
    vc_r = Dict(nw => _PM.var(pm, nw, :vc_r, i) for nw in _PM.nw_ids(pm))
    vc_i = Dict(nw => _PM.var(pm, nw, :vc_i, i) for nw in _PM.nw_ids(pm))
    ic_r = Dict(nw => _PM.var(pm, nw, :ic_r, i) for nw in _PM.nw_ids(pm))
    ic_i = Dict(nw => _PM.var(pm, nw, :ic_i, i) for nw in _PM.nw_ids(pm))
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)
    qconv_ac = _PM.var(pm, n, :qconv_ac, i)


    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * pconv_ac ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vc_r[n1] * ic_r[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                                                    +
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vc_i[n1] * ic_i[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * qconv_ac ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vc_i[n1] * ic_r[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                                                    -
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vc_r[n1] * ic_i[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

end



function constraint_gp_converter_losses(pm::_PM.AbstractIVRModel, n::Int, i, T2, T3, a, b, c, plmax)
    iconv_lin = _PM.var(pm, n, :iconv_lin, i)
    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)
    pconv_dc = _PM.var(pm, n, :pconv_dc, i)

    JuMP.@constraint(pm.model, pconv_ac + pconv_dc == a + b * iconv_lin + c * iconv_lin_s)

end




function constraint_gp_iconv_lin_squared_2(pm::AbstractACRModel, n::Int, i, T2, T3)

    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    iconv_lin  = Dict(nw => _PM.var(pm, nw, :iconv_lin, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * iconv_lin_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (iconv_lin[n1] * iconv_lin[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vms = _PM.var(pm, n, :vms, i)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vms 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_RES_power_real(pm::AbstractIVRModel, n::Int, i, p, pd, T2, T3, p_size; curt=0.0)
        
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    crd_RES = Dict(nw => _PM.var(pm, nw, :crd_RES, p) for nw in _PM.nw_ids(pm))
    cid_RES = Dict(nw => _PM.var(pm, nw, :cid_RES, p) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd * p_size *(1-curt)
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd_RES[n2] + vi[n1] * cid_RES[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
        
end

function constraint_gp_RES_power_imaginary(pm::AbstractIVRModel, n::Int, i, p, qd, T2, T3, q_size; curt=0.0)
    
    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in _PM.nw_ids(pm))
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in _PM.nw_ids(pm))

    crd_RES = Dict(n => _PM.var(pm, n, :crd_RES, p) for n in _PM.nw_ids(pm))
    cid_RES = Dict(n => _PM.var(pm, n, :cid_RES, p) for n in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd * q_size * (1-curt)
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd_RES[n2] - vr[n1] * cid_RES[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    cmss  = _PM.var(pm, n, :cmss, i)
    csr = Dict(nw => _PM.var(pm, nw, :csr, i) for nw in _PM.nw_ids(pm))
    csi = Dict(nw => _PM.var(pm, nw, :csi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmss
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

## generator
""
function constraint_gp_gen_power_real(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))

    pg  = _PM.var(pm, n, :pg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * crg[n2] + vi[n1] * cig[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end
""
function constraint_gp_gen_power_imaginary(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))

    qg  = _PM.var(pm, n, :qg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crg[n2] - vr[n1] * cig[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

## load
""
function constraint_gp_load_power_real(pm::AbstractIVRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    crd = Dict(nw => _PM.var(pm, nw, :crd, l) for nw in _PM.nw_ids(pm))
    cid = Dict(nw => _PM.var(pm, nw, :cid, l) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end
""
function constraint_gp_load_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in _PM.nw_ids(pm))
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in _PM.nw_ids(pm))

    crd = Dict(n => _PM.var(pm, n, :crd, l) for n in _PM.nw_ids(pm))
    cid = Dict(n => _PM.var(pm, n, :cid, l) for n in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_branch_series_current_on_off(pm::AbstractIVRModel, n::Int, i, T2, T3)

    csr_on_off  = _PM.var(pm, n, :csr_on_off, i)
    csi_on_off  = _PM.var(pm, n, :csi_on_off, i)

    csr  = _PM.var(pm, n, :csr, i)
    csi = _PM.var(pm, n, :csi, i)
    z_branch  = _PM.var(pm, 1, :z_branch, i)

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * csr_on_off
                                ==
                                T2.get([n-1,n-1]) * csr * z_branch
                    )

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * csi_on_off
                                ==
                                T2.get([n-1,n-1]) * csi * z_branch
                    )

end

function constraint_gp_branch_series_current_magnitude_squared_on_off(pm::AbstractIVRModel, n::Int, i, T2, T3)
    cmss  = _PM.var(pm, n, :cmss, i)
    
    csr_on_off = Dict(nw => _PM.var(pm, nw, :csr_on_off, i) for nw in _PM.nw_ids(pm))
    csi_on_off = Dict(nw => _PM.var(pm, nw, :csi_on_off, i) for nw in _PM.nw_ids(pm))


    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmss
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (csr_on_off[n1] * csr_on_off[n2] + csi_on_off[n1] * csi_on_off[n2])  
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )


end

function constraint_gp_ac_branch_indicator(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    err = 1e-8

    z_branch = _PM.var(pm, 1, :z_branch, i)

    JuMP.@constraint(pm.model, z_branch * (1 - z_branch) <= err)

end

function constraint_gp_dc_branch_indicator(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    err = 1e-8

    z_branch_dc = _PM.var(pm, 1, :z_branch_dc, i)

    JuMP.@constraint(pm.model, z_branch_dc * (1 - z_branch_dc) <= err)

    # JuMP.fix(z_branch_dc[1], 1; force=true)

end

function constraint_gp_ohms_dc_branch_on_off_part1(pm::AbstractACRModel, n::Int, b, T2, T3, f_bus, t_bus, f_idx, t_idx, r, p)

    p_fr  = _PM.var(pm, n,  :p_dcgrid, f_idx)
    p_to  = _PM.var(pm, n,  :p_dcgrid, t_idx)

    i_dc_fr_on_off = _PM.var(pm, n,  :igrid_dc_on_off, f_idx)
    i_dc_to_on_off = _PM.var(pm, n,  :igrid_dc_on_off, t_idx)
    
    i_dc_fr = _PM.var(pm, n,  :igrid_dc, f_idx)
    i_dc_to = _PM.var(pm, n,  :igrid_dc, t_idx)

    z_branch_dc  = _PM.var(pm, 1, :z_branch_dc, b)

    vmdc_fr = Dict(nw => _PM.var(pm, nw, :vdcm, f_bus) for nw in _PM.nw_ids(pm))
    vmdc_to = Dict(nw => _PM.var(pm, nw, :vdcm, t_bus) for nw in _PM.nw_ids(pm))


    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * p_fr ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vmdc_fr[n1] * i_dc_fr_on_off[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )


    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * p_to  ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) *
                                                    (vmdc_to[n1] * i_dc_to_on_off[n2])
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )


    vmdc_fr = _PM.var(pm, n,  :vdcm, f_bus)
    vmdc_to = _PM.var(pm, n,  :vdcm, t_bus)


    if r == 0
        JuMP.@constraint(pm.model, (i_dc_fr_on_off + i_dc_to_on_off) == 0)
    else
        JuMP.@constraint(pm.model, 0 <= vmdc_to - (vmdc_fr - (1/p * r * i_dc_fr)) + (1-z_branch_dc)*(1e-8))
        JuMP.@constraint(pm.model, 0 >= vmdc_to - (vmdc_fr - (1/p * r * i_dc_fr)) + (1-z_branch_dc)*(0))
        JuMP.@constraint(pm.model, 0 <= vmdc_fr - (vmdc_to - (1/p * r * i_dc_to)) + (1-z_branch_dc)*(1e-8))
        JuMP.@constraint(pm.model, 0 >= vmdc_fr - (vmdc_to - (1/p * r * i_dc_to)) + (1-z_branch_dc)*(0))

    end


end


function constraint_gp_ohms_dc_branch_on_off_part2(pm::AbstractACRModel, n::Int, b, T2, T3, f_bus, t_bus, f_idx, t_idx, r, p)

    i_dc_fr = _PM.var(pm, n,  :igrid_dc, f_idx)
    i_dc_to = _PM.var(pm, n,  :igrid_dc, t_idx)
    i_dc_fr_on_off = _PM.var(pm, n,  :igrid_dc_on_off, f_idx)
    i_dc_to_on_off = _PM.var(pm, n,  :igrid_dc_on_off, t_idx)

    z_branch_dc  = _PM.var(pm, 1, :z_branch_dc, b)


    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * i_dc_fr_on_off ==  
                                                        T2.get([n-1,n-1]) * i_dc_fr * z_branch_dc
                    )

    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * i_dc_to_on_off ==  
                                                        T2.get([n-1,n-1]) * i_dc_to * z_branch_dc
                    )

end

