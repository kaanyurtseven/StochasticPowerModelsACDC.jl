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

function constraint_gp_converter_current_iconv_lin_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_current_iconv_lin_squared(pm, nw, i, T2, T3)
end

function constraint_gp_converter_iconv_lin_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_iconv_lin_squared(pm, nw, i, T2, T3)
end




function constraint_cc_filter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :convdc, i, "Vmmin")
    vmax = _PM.ref(pm, nw, :convdc, i, "Vmmax")
    
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_filter_voltage_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_converter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :convdc, i, "Vmmin")
    vmax = _PM.ref(pm, nw, :convdc, i, "Vmmax")
    
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_voltage_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_transformer_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_transformer_current_from_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_transformer_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_transformer_current_to_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_reactor_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_reactor_current_from_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_reactor_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_reactor_current_to_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_converter_current_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
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

function constraint_cc_dc_branch_current(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    vpu = 1;
    branch = _PM.ref(pm, nw, :branchdc, i)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    
    Imax = branch["rateA"]/vpu
    Imin = - branch["rateA"]/vpu
   
       
    λmax = _PM.ref(pm, nw, :branchdc, i, "λcmax")
    λmin = _PM.ref(pm, nw, :branchdc, i, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_dc_branch_current(pm, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop)
end

function constraint_cc_iconv_lin_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_iconv_lin_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_iconv_lin(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_iconv_lin(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_conv_ac_power(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Pacmin = _PM.ref(pm, nw, :convdc, i, "Pacmin")
    Pacmax = _PM.ref(pm, nw, :convdc, i, "Pacmax")
    Qacmin = _PM.ref(pm, nw, :convdc, i, "Qacmin")
    Qacmax = _PM.ref(pm, nw, :convdc, i, "Qacmax")

    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")

    #λqmin = _PM.ref(pm, nw, :convdc, g, "λqmin")
    #λqmax = _PM.ref(pm, nw, :convdc, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_conv_ac_power_real(pm, i, Pacmin, Pacmax, λmin, λmax, T2, mop)
    constraint_cc_conv_ac_power_imaginary(pm, i, Qacmin, Qacmax, λmin, λmax, T2, mop)
end

function constraint_cc_conv_dc_power(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Pdcmin = - _PM.ref(pm, nw, :convdc, i, "Pacrated")
    Pdcmax = _PM.ref(pm, nw, :convdc, i, "Pacrated")
    
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")


    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_conv_dc_power(pm, i, Pdcmin, Pdcmax, λmin, λmax, T2, mop)
end


function constraint_cc_converter_dc_current(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    Imin = - conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_dc_current(pm, i, Imax, Imin, λmax, λmin, T2, mop)
end