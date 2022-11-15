function variable_dc_converter(pm::_PM.AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true,report::Bool=false, kwargs...)
    
    
    _PMACDC.variable_filter_voltage_real(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_filter_voltage_imaginary(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_voltage_real(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_voltage_imaginary(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_real_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_real_to(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_imaginary_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_imaginary_to(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_real_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_real_to(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_imaginary_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_imaginary_to(pm, nw=nw, bounded=bounded,kwargs...)
    _PMACDC.variable_converter_current_real(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_converter_current_imaginary(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_current_dc(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_current_lin(pm, nw=nw, bounded=true, kwargs...)
    _PMACDC.variable_converter_active_power(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_converter_reactive_power(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_dcside_power(pm, nw=nw, bounded=bounded,  kwargs...)


    ################
    _PMACDC.variable_conv_tranformer_flow(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_conv_reactor_flow(pm, nw=nw, bounded=bounded, kwargs...)

    _PMACDC.variable_converter_active_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMACDC.variable_converter_reactive_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMACDC.variable_acside_current(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_dcside_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMACDC.variable_converter_firing_angle(pm, nw=nw, bounded=bounded, kwargs...)

    _PMACDC.variable_converter_filter_voltage(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_internal_voltage(pm, nw=nw, bounded=bounded, kwargs...)

    _PMACDC.variable_converter_to_grid_active_power(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_to_grid_reactive_power(pm, nw=nw, bounded=bounded, kwargs...)



end

function variable_dc_converter_squared(pm::_PM.AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true,report::Bool=true, kwargs...)
    

    variable_filter_voltage_squared(pm, nw=nw, bounded=bounded, kwargs...)                  # v ct gp cc
    variable_converter_voltage_squared(pm, nw=nw, bounded=bounded, kwargs...)               # v ct gp cc
    variable_transformer_current_from_squared(pm, nw=nw, bounded=bounded,  kwargs...)       # v ct gp cc
    variable_transformer_current_to_squared(pm, nw=nw, bounded=bounded,  kwargs...)         # v ct gp cc
    variable_reactor_current_from_squared(pm, nw=nw, bounded=bounded,  kwargs...)           # v ct gp cc
    variable_reactor_current_to_squared(pm, nw=nw, bounded=bounded,  kwargs...)             # v ct gp cc
    variable_converter_current_lin_squared(pm, nw=nw, bounded=bounded,  kwargs...)          # v
    
end

function variable_filter_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)



    vk_s = _PM.var(pm, nw)[:vk_s] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_vk_s",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "vr_start", 1.0)
    )

    if bounded
        for (i, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(vk_s[i], convdc["Vmmin"]^2)
            JuMP.set_upper_bound(vk_s[i], convdc["Vmmax"]^2)
        end
    end
    report && _PM.sol_component_value(pm, nw, :convdc, :vk_s, _PM.ids(pm, nw, :convdc), vk_s)
    #report && _IM.sol_component_value(pm, _PM.pm_it_sym, nw, :convdc, :vk_s, _PM.ids(pm, nw, :convdc), vk_s)
end

function variable_converter_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    vc_s = _PM.var(pm, nw)[:vc_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_vc_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "v_start", 1.0)
    )

    if bounded
        for (i, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(vc_s[i], convdc["Vmmin"]^2)
            JuMP.set_upper_bound(vc_s[i], convdc["Vmmax"]^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :vc_s, _PM.ids(pm, nw, :convdc), vc_s)

end

function variable_transformer_current_from_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 1;
    vpu = 1;
    iik_s = _PM.var(pm, nw)[:iik_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_iik_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(iik_s[c],  0)
            JuMP.set_upper_bound(iik_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :iik_s, _PM.ids(pm, nw, :convdc), iik_s)
end


function variable_transformer_current_to_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 1;
    vpu = 1;
    iki_s = _PM.var(pm, nw)[:iki_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_iki_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(iki_s[c],  0)
            JuMP.set_upper_bound(iki_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :iki_s, _PM.ids(pm, nw, :convdc), iki_s)
end

function variable_reactor_current_from_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 1;
    vpu = 1;
    ikc_s = _PM.var(pm, nw)[:ikc_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_ikc_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(ikc_s[c],  0)
            JuMP.set_upper_bound(ikc_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :ikc_s, _PM.ids(pm, nw, :convdc), ikc_s)
end

function variable_reactor_current_to_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 1;
    vpu = 1;
    ick_s = _PM.var(pm, nw)[:ick_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_ick_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(ick_s[c],  0)
            JuMP.set_upper_bound(ick_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :ick_s, _PM.ids(pm, nw, :convdc), ick_s)
end

function variable_converter_current_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 1;
    vpu = 1;
    ic_s = _PM.var(pm, nw)[:ic_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_ic_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(ic_s[c],  0)
            JuMP.set_upper_bound(ic_s[c],  (convdc["Imax"]^2) * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :ic_s, _PM.ids(pm, nw, :convdc), ic_s)
end


function variable_converter_current_lin_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 1;
    vpu = 1;
    iconv_lin_s = _PM.var(pm, nw)[:iconv_lin_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_iconv_lin_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(iconv_lin_s[c],  0)
            JuMP.set_upper_bound(iconv_lin_s[c],  (convdc["Imax"]^2) * vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :iconv_lin_s, _PM.ids(pm, nw, :convdc), iconv_lin_s)
end










