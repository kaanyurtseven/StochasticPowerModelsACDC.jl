################################################################################
# Copyright 2023, Kaan Yurtseven                                               #
################################################################################
# StochasticPowerModelsACDC.jl                                                 #
# An extention package of PowerModels.jl and StochasticPowerModels.jl for      #
#                                 Stochastic Optimal Power Flow in AC/DC grids #
# See https://github.com/kaanyurtseven/StochasticPowerModelsACDC               #
################################################################################

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

function constraint_gp_iconv_lin_squared_1(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_iconv_lin_squared_1(pm, nw, i, T2, T3)
end

function constraint_gp_iconv_lin_squared_2(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_iconv_lin_squared_2(pm, nw, i, T2, T3)
end

function constraint_gp_converter_dc_power(pm::_PM.AbstractIVRModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    b_idx =   conv["busdc_i"]

    constraint_gp_converter_dc_power(pm, nw, i, T2, T3, b_idx)
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

function constraint_gp_converter_ac_power(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_ac_power(pm, nw, i, T2, T3)
end


function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_bus_voltage_magnitude_squared(pm, nw, i, T2, T3)
end


function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_magnitude_squared(pm, nw, b, T2, T3)
end

function constraint_gp_ohms_dc_branch(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    branch = _PM.ref(pm, nw, :branchdc, b)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)

    constraint_gp_ohms_dc_branch(pm, nw, b, T2, T3, f_bus, t_bus, f_idx, t_idx, branch["r"], p)
end


function constraint_gp_power_branch_to(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
end

function constraint_gp_power_branch_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
end

function constraint_gp_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :gen, g, "gen_bus")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_gen_power_real(pm, nw, i, g, T2, T3)
    constraint_gp_gen_power_imaginary(pm, nw, i, g, T2, T3)
end

function constraint_gp_load_power(pm::AbstractPowerModel, l::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :load, l, "load_bus") 

    pd  = _PM.ref(pm, nw, :load, l, "pd")
    qd  = _PM.ref(pm, nw, :load, l, "qd")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_load_power_real(pm, nw, i, l, pd, T2, T3)
    constraint_gp_load_power_imaginary(pm, nw, i, l, qd, T2, T3)
end

function constraint_gp_RES_power(pm::AbstractPowerModel, p::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :RES, p, "load_bus") 

    pd  = _PM.ref(pm, nw, :RES, p, "pd")
    qd  = _PM.ref(pm, nw, :RES, p, "qd")

    # p_size= _PM.var(pm, 1, :p_size, p)
    # q_size= _PM.var(pm, 1, :q_size, p)

    p_size = _PM.ref(pm, nw, :RES, p, "p_size")
    q_size = _PM.ref(pm, nw, :RES, p, "q_size")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    # c = pm.data["curt"]

    constraint_gp_RES_power_real(pm, nw, i, p, pd, T2, T3, p_size)
    constraint_gp_RES_power_imaginary(pm, nw, i, p, qd, T2, T3, q_size)
end

function constraint_gp_branch_series_current_on_off(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_on_off(pm, nw, b, T2, T3)
end

function constraint_gp_branch_series_current_magnitude_squared_on_off(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_magnitude_squared_on_off(pm, nw, b, T2, T3)
end

function constraint_gp_ohms_dc_branch_on_off_part1(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    branch = _PM.ref(pm, nw, :branchdc, b)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)

    constraint_gp_ohms_dc_branch_on_off_part1(pm, nw, b, T2, T3, f_bus, t_bus, f_idx, t_idx, branch["r"], p)
end

function constraint_gp_ohms_dc_branch_on_off_part2(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    branch = _PM.ref(pm, nw, :branchdc, b)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)

    constraint_gp_ohms_dc_branch_on_off_part2(pm, nw, b, T2, T3, f_bus, t_bus, f_idx, t_idx, branch["r"], p)
end