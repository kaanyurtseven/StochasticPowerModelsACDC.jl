################################################################################
# Copyright 2023, Kaan Yurtseven                                               #
################################################################################
# StochasticPowerModelsACDC.jl                                                 #
# An extention package of PowerModels.jl and StochasticPowerModels.jl for      #
#                                 Stochastic Optimal Power Flow in AC/DC grids #
# See https://github.com/kaanyurtseven/StochasticPowerModelsACDC               #
################################################################################


function variable_bus_voltage(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)
    _PM.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_bus_voltage_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)

end

""
function variable_dcgrid_voltage_magnitude(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)
    
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)

end

function variable_active_dcbranch_flow(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)

    #DC grid variables
    _PMACDC.variable_active_dcbranch_flow(pm, nw=nw, bounded=bounded, report=report; kwargs...)

end

function variable_dcbranch_current(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)

    _PMACDC.variable_dcbranch_current(pm, nw=nw, bounded=bounded, report=report; kwargs...)

end


"variable: `vms[i]` for `i` in `bus`"
function variable_bus_voltage_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vms = _PM.var(pm, nw)[:vms] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vms",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vms_start", 1.0)
    )

    if bounded
        for (i, bus) in _PM.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vms[i], 0)
            JuMP.set_upper_bound(vms[i], 2.0 * bus["vmax"]^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :bus, :vms, _PM.ids(pm, nw, :bus), vms)
end

# branch
"variable: `cmss[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cmss = _PM.var(pm, nw)[:cmss] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_cmss",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "cmss_start", 0.0)
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        branch = _PM.ref(pm, nw, :branch)

        for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"] * b["tap"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                series_current = max(rate / bus[i]["vmin"], rate / bus[j]["vmin"])
                ub = series_current + shunt_current
            end
            if haskey(b, "c_rating_a")
                total_current = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                ub = total_current + shunt_current
            end

            if !isinf(ub)
                JuMP.set_lower_bound(cmss[l], 0)
                JuMP.set_upper_bound(cmss[l],  2.0 * ub^2)
            end
        end
    end

    report && _PM.sol_component_value(pm, nw, :branch, :cmss, _PM.ids(pm, nw, :branch), cmss)
end

# generator
""
function variable_gen_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)
    _PM.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_RES_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)

    variable_RES_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_RES_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# load current
"variable: `crd[j]` for `j` in `load`"
function variable_RES_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd_RES = _PM.var(pm, nw)[:crd_RES] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :RES)], base_name="$(nw)_crd_RES",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :RES, i), "crd_RES_start")
    )
    # if bounded
    #     for (i, RES) in _PM.ref(pm, nw, :RES)
    #         JuMP.set_lower_bound(crd_RES[i], 0)
    #     end
    # end
    
    report && _PM.sol_component_value(pm, nw, :RES, :crd_RES, _PM.ids(pm, nw, :RES), crd_RES)
end


"variable: `cid[j]` for `j` in `load`"
function variable_RES_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid_RES = _PM.var(pm, nw)[:cid_RES] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :RES)], base_name="$(nw)_cid_RES",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :RES, i), "cid_RES_start")
    )
    # if bounded
    #     for (i, RES) in _PM.ref(pm, nw, :RES)
    #         JuMP.set_lower_bound(cid_RES[i], 0)
    #     end
    # end

    report && _PM.sol_component_value(pm, nw, :RES, :cid_RES, _PM.ids(pm, nw, :RES), cid_RES)
end

# RES size
"variable: `p_size` for `j` in `load`"
function variable_RES_size(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    p_size = _PM.var(pm, nw)[:p_size] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :RES)], base_name="$(nw)_p_size",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :RES, i), "p_size_start", 1)
    )
    
    if bounded
        for (i, RES) in _PM.ref(pm, nw, :RES)
            if haskey(RES, "p_max") & haskey(RES,"p_min")
                # print(RES)
                JuMP.set_lower_bound(p_size[i], RES["p_min"])
                JuMP.set_upper_bound(p_size[i], RES["p_max"]) #2*RES["conn_cap_kW"])
            else
                JuMP.set_lower_bound(p_size[i], 0.1)
                JuMP.set_upper_bound(p_size[i], 10) #2*RES["conn_cap_kW"])
            end
        end
    end
    report && _PM.sol_component_value(pm, nw, :RES, :p_size, _PM.ids(pm, nw, :RES), p_size)
end

function variable_RES_size_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    q_size = _PM.var(pm, nw)[:q_size] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :RES)], base_name="$(nw)_q_size",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :RES, i), "q_size_start")
    )
    
    if bounded
        for (i, RES) in _PM.ref(pm, nw, :RES)
            if haskey(RES, "p_max") & haskey(RES,"p_min")
                JuMP.set_lower_bound(q_size[i], -0.2*RES["p_max"])
                JuMP.set_upper_bound(q_size[i], 0.2*RES["p_max"]) #2*RES["conn_cap_kW"])
            else

                #### Negatif bound verince sifira esitliyor kendini.

                JuMP.set_lower_bound(q_size[i], -1) 
                JuMP.set_upper_bound(q_size[i], 1) #2*RES["conn_cap_kW"])
            end
        end
    end
    report && _PM.sol_component_value(pm, nw, :RES, :q_size, _PM.ids(pm, nw, :RES), q_size)
end

function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    expression_variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cr = _PM.var(pm, nw)[:cr] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        cr[(l,i,j)] = (tr * csr_fr - ti * csi_fr + g_sh_fr * vr_fr - b_sh_fr * vi_fr) / tm^2
        cr[(l,j,i)] = -csr_fr + g_sh_to * vr_to - b_sh_to * vi_to

    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), cr)
end
"variable: `ci[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ci = _PM.var(pm, nw)[:ci] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        ci[(l,i,j)] = (tr * csi_fr + ti * csr_fr + g_sh_fr * vi_fr + b_sh_fr * vr_fr) / tm^2
        ci[(l,j,i)] = -csi_fr + g_sh_to * vi_to + b_sh_to * vr_to

    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), ci)
end

function variable_load_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `crd[j]` for `j` in `load`"
function variable_load_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd = _PM.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_crd",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "crd_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :crd, _PM.ids(pm, nw, :load), crd)
end
"variable: `cid[j]` for `j` in `load`"
function variable_load_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid = _PM.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cid",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cid_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :cid, _PM.ids(pm, nw, :load), cid)
end

function variable_gen_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end


function variable_ac_branch_indicator(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, relax::Bool=false, report::Bool=true)
    br = Dict()
    for l in _PM.ids(pm, nw, :branch)
        br[l] = _PM.ref(pm,nw,:branch,l)
        # display(br[l]["br_status_initial"])
    end 


    z_branch = _PM.var(pm, nw)[:z_branch] = JuMP.@variable(pm.model,
    [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_z_branch",
    binary = false,
    lower_bound = 0,
    upper_bound = 1,
    start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "z_branch_start", br[l]["br_status_initial"])
    # start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "z_branch_start", 0.5)
    )
    
    report && _PM.sol_component_value(pm, nw, :branch, :br_status, _PM.ids(pm, nw, :branch), z_branch)
end

function variable_dc_branch_indicator(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, relax::Bool=false, report::Bool=true)
 
    br = Dict()
    for l in _PM.ids(pm, nw, :branchdc)
        br[l] = _PM.ref(pm,nw,:branchdc,l)
    end 

    z_branch_dc = _PM.var(pm, nw)[:z_branch_dc] = JuMP.@variable(pm.model,
    [l in _PM.ids(pm, nw, :branchdc)], base_name="$(nw)_z_branch_dc",
    binary = false,
    lower_bound = 0,
    upper_bound = 1,
    start = _PM.comp_start_value(_PM.ref(pm, nw, :branchdc, l), "z_branch_dc_start", br[l]["br_status_initial"])
    )

    report && _PM.sol_component_value(pm, nw, :branchdc, :br_status_dc, _PM.ids(pm, nw, :branchdc), z_branch_dc)
end


function variable_branch_current_on_off_part1(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_real_on_off(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_imaginary_on_off(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_branch_series_current_real_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    csr_on_off = _PM.var(pm, nw)[:csr_on_off] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_csr_on_off",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "csr_on_off_start", 0.0)
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        branch = _PM.ref(pm, nw, :branch)

        for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"]*b["tap"]
                y_fr = abs(b["g_fr"] + im*b["b_fr"])
                y_to = abs(b["g_to"] + im*b["b_to"])
                shunt_current = max(y_fr*bus[i]["vmax"]^2, y_to*bus[j]["vmax"]^2)
                series_current = max(rate/bus[i]["vmin"], rate/bus[j]["vmin"])
                ub = series_current + shunt_current
            end
            if haskey(b, "c_rating_a")
                total_current = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im*b["b_fr"])
                y_to = abs(b["g_to"] + im*b["b_to"])
                shunt_current = max(y_fr*bus[i]["vmax"]^2, y_to*bus[j]["vmax"]^2)
                ub = total_current + shunt_current
            end

            if !isinf(ub)
                JuMP.set_lower_bound(csr_on_off[l], -ub)
                JuMP.set_upper_bound(csr_on_off[l],  ub)
            end
        end
    end


    report && _PM.sol_component_value(pm, nw, :branch, :csr_fr_on_off, _PM.ids(pm, nw, :branch), csr_on_off)
end

"variable: `csi[l,i,j] ` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_imaginary_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    csi_on_off = _PM.var(pm, nw)[:csi_on_off] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_csi_on_off",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "csi_on_off_start", 0.0)
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        branch = _PM.ref(pm, nw, :branch)

        for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"]*b["tap"]
                y_fr = abs(b["g_fr"] + im*b["b_fr"])
                y_to = abs(b["g_to"] + im*b["b_to"])
                shuntcurrent = max(y_fr*bus[i]["vmax"]^2, y_to*bus[j]["vmax"]^2)
                seriescurrent = max(rate/bus[i]["vmin"], rate/bus[j]["vmin"])
                ub = seriescurrent + shuntcurrent
            end
            if haskey(b, "c_rating_a")
                totalcurrent = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im*b["b_fr"])
                y_to = abs(b["g_to"] + im*b["b_to"])
                shuntcurrent = max(y_fr*bus[i]["vmax"]^2, y_to*bus[j]["vmax"]^2)
                ub = totalcurrent + shuntcurrent
            end

            if !isinf(ub)
                JuMP.set_lower_bound(csi_on_off[l], -ub)
                JuMP.set_upper_bound(csi_on_off[l],  ub)
            end

        end
    end


    report && _PM.sol_component_value(pm, nw, :branch, :csi_fr_on_off, _PM.ids(pm, nw, :branch), csi_on_off)
end


# function variable_active_dcbranch_flow_on_off(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
#     p_on_off = _PM.var(pm, nw)[:p_dcgrid_on_off] = JuMP.@variable(pm.model,
#     [(l,i,j) in _PM.ref(pm, nw, :arcs_dcgrid)], base_name="$(nw)_pdcgrid_on_off",
#     start = _PM.comp_start_value(_PM.ref(pm, nw, :branchdc, l), "p_on_off_start", 1.0)
#     )

#     if bounded
#         for arc in _PM.ref(pm, nw, :arcs_dcgrid)
#             l,i,j = arc
#             JuMP.set_lower_bound(p_on_off[arc], -_PM.ref(pm, nw, :branchdc, l)["rateA"])
#             JuMP.set_upper_bound(p_on_off[arc],  _PM.ref(pm, nw, :branchdc, l)["rateA"])
#         end
#     end

#     report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branchdc, :pf_on_off, :pt_on_off, _PM.ref(pm, nw, :arcs_dcgrid_from), _PM.ref(pm, nw, :arcs_dcgrid_to), p_on_off)
# end



function variable_dcbranch_current_on_off(pm::_PM.AbstractIVRModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    vpu = 1;
    igrid_dc_on_off = _PM.var(pm, nw)[:igrid_dc_on_off] = JuMP.@variable(pm.model,
    [(l,i,j) in _PM.ref(pm, nw, :arcs_dcgrid)], base_name="$(nw)_igrid_dc_on_off",
    start = (_PM.comp_start_value(_PM.ref(pm, nw, :branchdc, l), "p_start", 0.0) / vpu)
    )
    if bounded
        for arc in _PM.ref(pm, nw, :arcs_dcgrid)
            l,i,j = arc
            JuMP.set_lower_bound(igrid_dc_on_off[arc], -_PM.ref(pm, nw, :branchdc, l)["rateA"] / vpu)
            JuMP.set_upper_bound(igrid_dc_on_off[arc],  _PM.ref(pm, nw, :branchdc, l)["rateA"] / vpu)
        end
    end
    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branchdc, :if_on_off, :it_on_off, _PM.ref(pm, nw, :arcs_dcgrid_from), _PM.ref(pm, nw, :arcs_dcgrid_to), igrid_dc_on_off)
end

function variable_branch_current_on_off_part2(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    
    expression_variable_branch_current_real_on_off(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary_on_off(pm, nw=nw, bounded=bounded, report=report; kwargs...)

end

function expression_variable_branch_current_real_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cr = _PM.var(pm, nw)[:cr] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        z_branch = _PM.var(pm, 1, :z_branch, l)

        cr[(l,i,j)] = z_branch * ((tr * csr_fr - ti * csi_fr + g_sh_fr * vr_fr - b_sh_fr * vi_fr) / tm^2)
        cr[(l,j,i)] = z_branch * (-csr_fr + g_sh_to * vr_to - b_sh_to * vi_to)

    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), cr)
end


function expression_variable_branch_current_imaginary_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ci = _PM.var(pm, nw)[:ci] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)

        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)
    
        z_branch = _PM.var(pm, 1, :z_branch, l)

        ci[(l,i,j)] = z_branch * ((tr * csi_fr + ti * csr_fr + g_sh_fr * vi_fr + b_sh_fr * vr_fr) / tm^2)
        ci[(l,j,i)] = z_branch * (-csi_fr + g_sh_to * vi_to + b_sh_to * vr_to)

    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), ci)
end
