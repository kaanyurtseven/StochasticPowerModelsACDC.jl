


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
    Imin = 0
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmax") # All λ values are equal.
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_iconv_lin(pm, i, Imax, Imin, λmax, λmin, T2, mop)
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


function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :bus, i, "vmin")
    vmax = _PM.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PM.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PM.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_bus_voltage_magnitude_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_conv_voltage_magnitude(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
 
    vmin = _PM.ref(pm, nw, :busdc, i, "Vdcmin")
    vmax = _PM.ref(pm, nw, :busdc, i, "Vdcmax")
    
    λmin = _PM.ref(pm, nw, :busdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :busdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_conv_voltage_magnitude(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    cmax = _PM.ref(pm, nw, :branch, b, "cmax")
    λmax = _PM.ref(pm, nw, :branch, b, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_branch_series_current_magnitude_squared(pm, b, cmax, λmax, T2, mop)
end

function constraint_cc_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    pmin = _PM.ref(pm, nw, :gen, g, "pmin")
    pmax = _PM.ref(pm, nw, :gen, g, "pmax")
    qmin = _PM.ref(pm, nw, :gen, g, "qmin")
    qmax = _PM.ref(pm, nw, :gen, g, "qmax")

    λpmin = _PM.ref(pm, nw, :gen, g, "λpmin")
    λpmax = _PM.ref(pm, nw, :gen, g, "λpmax")
    λqmin = _PM.ref(pm, nw, :gen, g, "λqmin")
    λqmax = _PM.ref(pm, nw, :gen, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_gen_power_real(pm, g, pmin, pmax, λpmin, λpmax, T2, mop)
    constraint_cc_gen_power_imaginary(pm, g, qmin, qmax, λqmin, λqmax, T2, mop)
end

