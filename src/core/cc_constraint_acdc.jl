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


