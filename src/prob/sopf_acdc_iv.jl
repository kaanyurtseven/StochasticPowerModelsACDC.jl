################################################################################
#  Copyright 2021, Tom Van Acker, Frederik Geth                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function solve_sopf_acdc_iv(file::String, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)
    process_additional_data!(data)
    
    return solve_sopf_acdc_iv(data, model_constructor, optimizer; deg=deg,ref_extensions = [_PMACDC.add_ref_dcgrid!], solution_processors=solution_processors, kwargs...)
end

""
function solve_sopf_acdc_iv(data::Dict, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractIVRModel "This problem type only supports the IVRModel"
    
    sdata = build_stochastic_acdc_data(data, deg)
    result = _PM.solve_model(sdata, model_constructor, optimizer, build_sopf_acdc_iv; multinetwork=true, ref_extensions = [_PMACDC.add_ref_dcgrid!], solution_processors=solution_processors, kwargs...)
    result["mop"] = sdata["mop"]
    
    return result
end



""
function build_sopf_acdc_iv(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 

        _PM.variable_dcline_current(pm, nw=n)

        variable_bus_voltage(pm, nw=n)

        variable_branch_current(pm, nw=n)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)

        #DC grid variables
        variable_active_dcbranch_flow(pm, nw=n)
        variable_dcbranch_current(pm, nw=n)
        variable_dcgrid_voltage_magnitude(pm, nw=n)
        #DC converter variables
        variable_dc_converter(pm, nw=n, bounded=true)
        variable_dc_converter_squared(pm, nw=n, bounded=true)

    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance_ac(pm, i, nw=n)
            #constraint_current_balance(pm, i, nw=n)########################################################
            constraint_gp_bus_voltage_magnitude_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_branch_series_current_magnitude_squared(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end


        #DC Constraints        
        for i in _PM.ids(pm, :busdc, nw=n)
            constraint_current_balance_dc(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :branchdc, nw=n)
            constraint_gp_ohms_dc_branch(pm, i, nw=n) #######################################################
            constraint_ohms_dc_branch(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :convdc, nw=n)
            println("$i", i)
            constraint_gp_filter_voltage_squared(pm, i, nw=n)
            constraint_gp_converter_voltage_squared(pm, i, nw=n)
            constraint_gp_transformer_current_from_squared(pm, i, nw=n)
            constraint_gp_transformer_current_to_squared(pm, i, nw=n)
            constraint_gp_reactor_current_from_squared(pm, i, nw=n)
            constraint_gp_reactor_current_to_squared(pm, i, nw=n)
            #constraint_gp_converter_current_squared(pm, i, nw=n)
            constraint_gp_converter_current_iconv_lin_squared(pm, i, nw=n)

            constraint_gp_converter_limits(pm, i, nw=n)
            constraint_gp_converter_losses(pm, i, nw=n)
            constraint_gp_converter_current(pm, i, nw=n)

            constraint_conv_transformer(pm, i, nw=n)
            constraint_conv_reactor(pm, i, nw=n)
            constraint_conv_filter(pm, i, nw=n)

        end

    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_cc_branch_series_current_magnitude_squared(pm, b, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_cc_gen_power(pm, g, nw=1)
    end

    for i in _PM.ids(pm, :busdc, nw=1)
        constraint_cc_conv_voltage_magnitude(pm, i, nw=1)
        constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for i in _PM.ids(pm, :convdc, nw=1)

        constraint_cc_filter_voltage_squared(pm, i, nw=1)
        constraint_cc_converter_voltage_squared(pm, i, nw=1)
        constraint_cc_transformer_current_from_squared(pm, i, nw=1)
        constraint_cc_transformer_current_to_squared(pm, i, nw=1)
        constraint_cc_reactor_current_from_squared(pm, i, nw=1)
        constraint_cc_reactor_current_to_squared(pm, i, nw=1)
        constraint_cc_converter_current_squared(pm, i, nw=1)

    end

    objective_min_expected_generation_cost(pm)
end