################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# general constraints
## reference
""
function constraint_bus_voltage_ref(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_bus_voltage_ref(pm, nw, i)
end

## bus
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(_PM.con(pm, nw), :kcl_cr)
        _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PM.con(pm, nw), :kcl_ci)
        _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
end


function constraint_current_balance_ac(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(_PM.con(pm, nw), :kcl_cr)
        _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PM.con(pm, nw), :kcl_ci)
        _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance_ac(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs)
end


""
function constraint_power_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus_arcs   = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens   = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads  = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance(pm, nw, i, bus_arcs, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end

function constraint_ohms_dc_branch(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)

    branch = _PM.ref(pm, nw, :branchdc, b)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)

    constraint_ohms_dc_branch(pm, nw, b, f_bus, t_bus, f_idx, t_idx, branch["r"], p)
end

function constraint_conv_reactor(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_reactor(pm, nw, i, conv["rc"], conv["xc"], Bool(conv["reactor"]))
end

#
function constraint_conv_filter(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_filter(pm, nw, i, conv["bf"], Bool(conv["filter"]) )
end

#
function constraint_conv_transformer(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_transformer(pm, nw, i, conv["rtf"], conv["xtf"], conv["busac_i"], conv["tm"], Bool(conv["transformer"]))
end

function constraint_current_balance_dc(pm::_PM.AbstractIVRModel, i::Int; nw::Int=_PM.nw_id_default)
    bus_arcs_dcgrid = _PM.ref(pm, nw, :bus_arcs_dcgrid, i)
    bus_convs_dc = _PM.ref(pm, nw, :bus_convs_dc, i)
    pd = _PM.ref(pm, nw, :busdc, i)["Pdc"]
    constraint_current_balance_dc(pm, nw, bus_arcs_dcgrid, bus_convs_dc, pd)
end

# current balance with PV
""
function constraint_current_balance_with_PV(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
   if !haskey(_PM.con(pm, nw), :kcl_cr)
        _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PM.con(pm, nw), :kcl_ci)
        _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_PV = _PM.ref(pm, nw, :bus_loads, i)  #to ask Tom to get :PV it works now as all load has PV

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance_with_PV(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs, bus_PV)

end

function constraint_gp_pv_power(pm::AbstractPowerModel, p::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :PV, p, "load_bus") 

    pd  = _PM.ref(pm, nw, :PV, p, "pd")
    qd  = _PM.ref(pm, nw, :PV, p, "qd")

    # p_size= _PM.var(pm, 1, :p_size, p)
    # q_size= _PM.var(pm, 1, :q_size, p)

    p_size = _PM.ref(pm, nw, :PV, p, "p_size")
    q_size = _PM.ref(pm, nw, :PV, p, "q_size")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    # c = pm.data["curt"]

    constraint_gp_pv_power_real(pm, nw, i, p, pd, T2, T3, p_size)
    constraint_gp_pv_power_imaginary(pm, nw, i, p, qd, T2, T3, q_size)
end