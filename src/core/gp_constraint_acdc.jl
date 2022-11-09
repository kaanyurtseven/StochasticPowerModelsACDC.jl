
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

function constraint_ohms_dc_branch(pm::AbstractACRModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, p)

    vmdc_fr = _PM.var(pm, n,  :vdcm, f_bus)
    vmdc_to = _PM.var(pm, n,  :vdcm, t_bus)
    i_dc_fr = _PM.var(pm, n,  :igrid_dc, f_idx)
    i_dc_to = _PM.var(pm, n,  :igrid_dc, t_idx)


    if r == 0
        JuMP.@constraint(pm.model, i_dc_fr + i_dc_to == 0)
    else
        JuMP.@constraint(pm.model, vmdc_to ==  vmdc_fr - 1/p * r * i_dc_fr)
        JuMP.@constraint(pm.model, vmdc_fr ==  vmdc_to - 1/p * r * i_dc_to)
    end

   # JuMP.@constraint(pm.model, i_dc_fr + i_dc_to == 0)


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

function constraint_gp_converter_current_iconv_lin_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

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



function constraint_gp_converter_limits(pm::AbstractACRModel, n::Int, i, T2, T3, b_idx)

#=

    ic_s = _PM.var(pm, n, :ic_s, i)
    iconv_lin  = Dict(nw => _PM.var(pm, nw, :iconv_lin, i) for nw in _PM.nw_ids(pm))

    #iconv_lin = _PM.var(pm, n, :iconv_lin)[i]



    #Eq. (48)
    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * ic_s ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (iconv_lin[n1] * iconv_lin[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

=#

    pconv_dc = _PM.var(pm, n, :pconv_dc, i)
    iconv_dc  = Dict(nw => _PM.var(pm, nw, :iconv_dc, i) for nw in _PM.nw_ids(pm))
    vc  = Dict(nw => _PM.var(pm, nw, :vdcm, b_idx) for nw in _PM.nw_ids(pm))
    #=
    vc = _PM.var(pm, n, :vdcm)[b_idx]
    pconv_dc = _PM.var(pm, n, :pconv_dc)[i]
    iconv_dc = _PM.var(pm, n, :iconv_dc)[i]
    =#
    #Eq. (43)
    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * pconv_dc ==  
                                                    sum(T3.get([n1-1,n2-1,n-1]) * 
                                                    (vc[n1] * iconv_dc[n2]) 
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

end

#=
function constraint_gp_converter_losses(pm::_PM.AbstractIVRModel, n::Int, i::Int, T2, T3, a, b, c, plmax)
    iconv_lin = _PM.var(pm, n, :iconv_lin, i)
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)
    pconv_dc = _PM.var(pm, n, :pconv_dc, i)
    #Eq. (47)
    JuMP.@constraint(pm.model, T2.get([n-1,n-1]) * pconv_ac +
                                 T2.get([n-1,n-1]) * pconv_dc == 
                                                            a + b*T2.get([n-1,n-1])*iconv_lin + 
                                                            c * sum(T3.get([n1-1,n2-1,n-1]) * 
                                                            (iconv_lin[n1] * iconv_lin[n2]) 
                                                            for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
    
                    )
end

=#
function constraint_gp_converter_current(pm::_PM.AbstractIVRModel, n::Int, i::Int, T2, T3, Umax, Imax)
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

    JuMP.@NLconstraint(pm.model, pconv_ac + pconv_dc == a + c*iconv_lin_s)

end




function constraint_gp_converter_iconv_lin_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    iconv_lin  = Dict(nw => _PM.var(pm, nw, :iconv_lin, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * iconv_lin_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (iconv_lin[n1] * iconv_lin[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end
