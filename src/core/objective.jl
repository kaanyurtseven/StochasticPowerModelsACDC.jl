################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

"expected cost of active power generation"
function objective_min_expected_generation_cost(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()

    T2 = pm.data["T2"]

    for (g, gen) in _PM.ref(pm, :gen, nw=1)
        pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _PM.nw_ids(pm))

        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*pg[1] + 
                          gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PM.nw_ids(pm)) + 
                          gen["cost"][2]*pg[1] + 
                          gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PM.ids(pm, :gen, nw=1))
    )
end

"expected cost of active power generation"
function objective_min_expected_generation_cost_PV(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()
    pv_cost = Dict()

    T2 = pm.data["T2"]

    for (g, gen) in _PM.ref(pm, :gen, nw=1)
        pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _PM.nw_ids(pm))

        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*pg[1] + 
                          gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PM.nw_ids(pm)) + 
                          gen["cost"][2]*pg[1] + 
                          gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    for (p, PV) in _PM.ref(pm, :PV, nw=1)

        p_size = _PM.ref(pm, 1, :PV, p, "p_size") #Dict(nw => _PM.ref(pm, nw, :PV, p))# for nw in _PM.nw_ids(pm))
        pv_cost[p] = p_size * 0.1

    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PM.ids(pm, :gen, nw=1)) + sum(pv_cost[p] for p in _PM.ids(pm, :PV, nw=1))
    )
end


"expected max PV generation"
function objective_max_PV(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()


    for (p, PV) in _PM.ref(pm, :PV, nw=1)
        p_size[p] = Dict(nw => _PM.var(pm, nw, :p_size, p) for nw in [1])
    end

    return JuMP.@objective(pm.model, Max,
            sum(p_size[p][1] for p in _PM.ids(pm, :PV, nw=1))
    )
end


function objective_max_PV_imaginary(pm::AbstractPowerModel; kwargs...)
    p_size = Dict()
    q_size = Dict()

    for (p, PV) in _PM.ref(pm, :PV, nw=1)
        p_size[p] = Dict(nw => _PM.var(pm, nw, :p_size, p) for nw in [1])
        q_size[p] = Dict(nw => _PM.var(pm, nw, :q_size, p) for nw in [1])
    end

    return JuMP.@objective(pm.model, Max,
            sum(p_size[p][1]^2 + q_size[p][1]^2 for p in _PM.ids(pm, :PV, nw=1))
    )
end
