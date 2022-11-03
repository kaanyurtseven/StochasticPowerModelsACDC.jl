function constraint_gp_filter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_filter_voltage_squared(pm, nw, i, T2, T3)
end

function constraint_gp_converter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_voltage_squared(pm, nw, i, T2, T3)
end

function constraint_gp_transformer_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_transformer_current_from_squared(pm, nw, i, T2, T3)
end

function constraint_gp_transformer_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_transformer_current_to_squared(pm, nw, i, T2, T3)
end

function constraint_gp_reactor_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_reactor_current_from_squared(pm, nw, i, T2, T3)
end

function constraint_gp_reactor_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_reactor_current_to_squared(pm, nw, i, T2, T3)
end

function constraint_gp_converter_current_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_current_squared(pm, nw, i, T2, T3)
end

function constraint_gp_converter_limits(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_limits(pm, nw, i, T2, T3)
end

function constraint_cc_filter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :bus, i, "vmin")
    vmax = _PM.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PM.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PM.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_filter_voltage_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_converter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :bus, i, "vmin")
    vmax = _PM.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PM.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PM.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_voltage_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_transformer_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Imax = _PM.ref(pm, nw, :convdc, i, "Imax")
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_transformer_current_from_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_transformer_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Imax = _PM.ref(pm, nw, :convdc, i, "Imax")
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_transformer_current_to_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_reactor_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Imax = _PM.ref(pm, nw, :convdc, i, "Imax")
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_reactor_current_from_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_reactor_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Imax = _PM.ref(pm, nw, :convdc, i, "Imax")
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_reactor_current_to_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_converter_current_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Imax = _PM.ref(pm, nw, :convdc, i, "Imax")
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_current_squared(pm, i, Imax, λmax, T2, mop)
end


function constraint_gp_converter_limits(pm::_PM.AbstractIVRModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    b_idx =   conv["busdc_i"]

    constraint_gp_converter_limits(pm, nw, i, T2, T3, b_idx)
end


function constraint_gp_converter_losses(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    
    conv = _PM.ref(pm, nw, :convdc, i)
    a = conv["LossA"]
    b = conv["LossB"]
    c = conv["LossCinv"]
    plmax = conv["LossA"] + conv["LossB"] * conv["Pacrated"] + conv["LossCinv"] * (conv["Pacrated"])^2
    constraint_gp_converter_losses(pm, nw, i, T2, T3, a, b, c, plmax)
end

function constraint_gp_converter_current(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    
    conv = _PM.ref(pm, nw, :convdc, i)
    Vmax = conv["Vmmax"]
    Imax = conv["Imax"]
    constraint_gp_converter_current(pm, nw, i, T2, T3, Vmax, Imax)
end