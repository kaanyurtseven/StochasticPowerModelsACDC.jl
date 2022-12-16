################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# general constraints
## bus
""
function constraint_bus_voltage_ref(pm::AbstractACRModel, n::Int, i::Int)

    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    vn = ifelse(n == 1, 1.0, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end

#####################################

function constraint_conv_transformer(pm::_PM.AbstractIVRModel, n::Int, i::Int, rtf, xtf, acbus, tm, transformer)
    vi_r = _PM.var(pm, n, :vr, acbus)
    vi_i = _PM.var(pm, n, :vi, acbus)
    vk_r = _PM.var(pm, n, :vk_r, i)
    vk_i = _PM.var(pm, n, :vk_i, i)

    iik_r = _PM.var(pm, n, :iik_r, i)
    iik_i = _PM.var(pm, n, :iik_i, i)
    iki_r = _PM.var(pm, n, :iki_r, i)
    iki_i = _PM.var(pm, n, :iki_i, i)

    #TODO add transformation ratio.....
    if transformer
        JuMP.@constraint(pm.model, vk_r == vi_r - rtf * iik_r + xtf * iik_i) #(24)
        JuMP.@constraint(pm.model, vk_i == vi_i - rtf * iik_i - xtf * iik_r) #(25)
        JuMP.@constraint(pm.model, vi_r == vk_r - rtf * iki_r + xtf * iki_i) #reverse
        JuMP.@constraint(pm.model, vi_i == vk_i - rtf * iki_i - xtf * iki_r) #reverse
        
    else
        JuMP.@constraint(pm.model, vk_r == vi_r)
        JuMP.@constraint(pm.model, vk_i == vi_i)
        JuMP.@constraint(pm.model, iik_r + iki_r == 0)
        JuMP.@constraint(pm.model, iik_i + iki_i == 0)
    end

end

function constraint_conv_reactor(pm::_PM.AbstractIVRModel, n::Int, i::Int, rc, xc, reactor)
    vk_r = _PM.var(pm, n, :vk_r, i)
    vk_i = _PM.var(pm, n, :vk_i, i)
    vc_r = _PM.var(pm, n, :vc_r, i)
    vc_i = _PM.var(pm, n, :vc_i, i)

    ikc_r = _PM.var(pm, n, :ikc_r, i)
    ikc_i = _PM.var(pm, n, :ikc_i, i)
    ick_r = _PM.var(pm, n, :ick_r, i)
    ick_i = _PM.var(pm, n, :ick_i, i)
    ic_r = _PM.var(pm, n, :ic_r, i)
    ic_i = _PM.var(pm, n, :ic_i, i)

    JuMP.@constraint(pm.model, ick_r + ic_r == 0) #(20)
    JuMP.@constraint(pm.model, ick_i + ic_i == 0) #(21)

    if reactor
        JuMP.@constraint(pm.model, vc_r == vk_r - rc * ikc_r + xc * ikc_i) #(28)
        JuMP.@constraint(pm.model, vc_i == vk_i - rc * ikc_i - xc * ikc_r) #(29)
        JuMP.@constraint(pm.model, vk_r == vc_r - rc * ick_r + xc * ick_i) #reverse
        JuMP.@constraint(pm.model, vk_i == vc_i - rc * ick_i - xc * ick_r) #reverse
    else
        JuMP.@constraint(pm.model, vk_r == vc_r)
        JuMP.@constraint(pm.model, vk_i == vc_i)
        JuMP.@constraint(pm.model, ikc_r + ick_r == 0)
        JuMP.@constraint(pm.model, ikc_i + ick_i == 0)
    end
end

function constraint_conv_filter(pm::_PM.AbstractIVRModel, n::Int, i::Int, bv, filter)
    iki_r = _PM.var(pm, n, :iki_r, i)
    iki_i = _PM.var(pm, n, :iki_i, i)
    ikc_r = _PM.var(pm, n, :ikc_r, i)
    ikc_i = _PM.var(pm, n, :ikc_i, i)

    vk_r = _PM.var(pm, n, :vk_r, i)
    vk_i = _PM.var(pm, n, :vk_i, i)

    JuMP.@constraint(pm.model,   iki_r + ikc_r + bv * filter * vk_i == 0)
    JuMP.@constraint(pm.model,   iki_i + ikc_i - bv * filter * vk_r == 0)
end

# Kirchhoff's current law for DC nodes
function constraint_current_balance_dc(pm::_PM.AbstractIVRModel, n::Int, bus_arcs_dcgrid, bus_convs_dc, pd)
    
    igrid_dc = _PM.var(pm, n, :igrid_dc)
    iconv_dc = _PM.var(pm, n, :iconv_dc)


    JuMP.@constraint(pm.model, sum(igrid_dc[a] for a in bus_arcs_dcgrid) + sum(iconv_dc[c] for c in bus_convs_dc) == 0) # deal with pd

end
   

function constraint_gp_pv_power_real(pm::AbstractIVRModel, n::Int, i, p, pd, T2, T3, p_size; curt=0.0)
        
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    crd_pv = Dict(nw => _PM.var(pm, nw, :crd_pv, p) for nw in _PM.nw_ids(pm))
    cid_pv = Dict(nw => _PM.var(pm, nw, :cid_pv, p) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd * p_size *(1-curt)
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd_pv[n2] + vi[n1] * cid_pv[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
        
end

function constraint_gp_pv_power_imaginary(pm::AbstractIVRModel, n::Int, i, p, qd, T2, T3, q_size; curt=0.0)
    
    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in _PM.nw_ids(pm))
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in _PM.nw_ids(pm))

    crd_pv = Dict(n => _PM.var(pm, n, :crd_pv, p) for n in _PM.nw_ids(pm))
    cid_pv = Dict(n => _PM.var(pm, n, :cid_pv, p) for n in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd * q_size * (1-curt)
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd_pv[n2] - vr[n1] * cid_pv[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end


