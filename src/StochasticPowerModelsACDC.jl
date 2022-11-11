################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels(Distribution).jl for Stochastic (Optimal)#
# Power Flow                                                                   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

module StochasticPowerModelsACDC

    # import pkgs
    import InfrastructureModels
    import Ipopt
    import JuMP
    import KernelDensity
    import Memento
    import PolyChaos
    import PowerModels
    import PowerModelsACDC

    # import types
    import PowerModels: AbstractPowerModel, AbstractACRModel, AbstractIVRModel

    # pkgs const
    const _IM = InfrastructureModels
    const _KDE = KernelDensity
    const _PCE = PolyChaos
    const _PM = PowerModels
    const _SPMACDC = StochasticPowerModelsACDC
    const _PMACDC = PowerModelsACDC

    # memento logger
    function __init__()
        global _LOGGER = Memento.getlogger(PowerModels)
    end

    # const
    const nw_id_default = 1

    # funct
    sorted_nw_ids(pm) = sort(collect(_PM.nw_ids(pm)))

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("core/constraint.jl")
    #include("core/gp_constraint_acdc.jl")
    include("core/gp_constraint.jl")
    #include("core/cc_constraint_acdc.jl")
    include("core/cc_constraint.jl")

    include("core/constraint_template.jl")
    include("core/cc_constraint_template.jl")
    include("core/gp_constraint_template.jl")
    #include("core/cc_gp_constraint_template_acdc.jl")

    include("core/objective.jl")

    include("core/variable.jl")
    include("core/variableconv.jl")

    include("form/iv.jl")

    include("prob/sopf_iv.jl")
    include("prob/sopf_acdc_iv.jl")

    include("util/data.jl")
    include("util/util.jl")

    # export
    export BASE_DIR

    export solve_sopf_iv
    export solve_sopf_acdc_iv

    export build_stochastic_data
    export build_stochastic_acdc_data
    export extend_matlab_file
    export pce_coeff, sample, density, print_summary
end 
