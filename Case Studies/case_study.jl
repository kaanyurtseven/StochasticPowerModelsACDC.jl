using Pkg

Pkg.activate(".")


using JuMP
using Ipopt
using PowerModels
using StochasticPowerModelsACDC
using PowerModelsACDC
using Plots

# constants 
const _PM = PowerModels
const _SPMACDC = StochasticPowerModelsACDC
const _PMACDC = PowerModelsACDC

# solvers
#ipopt_solver = Ipopt.Optimizer
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_iter"=>3000, "fixed_variable_treatment" => "relax_bounds")
#ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_iter"=>3000)

# input
deg  = 1

case1 = "case67_AC_PV_SPMACDC.m"
file1  = joinpath(BASE_DIR, "test/data/matpower/Case Studies", case1)

case2 = "case67_ACDC_PV_SPMACDC.m"
file2  = joinpath(BASE_DIR, "test/data/matpower/Case Studies", case2)
 
#result_opf = _PM.solve_opf(file1, _PM.ACPPowerModel, ipopt_solver)

#result_spm = solve_sopf_iv(file, _PM.IVRPowerModel, ipopt_solver, deg=deg);

# s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)
# result_acdc = _PMACDC.run_acdcopf_iv(file, _PM.IVRPowerModel, ipopt_solver; setting = s)

result_spmacdc_case1 = solve_sopf_acdc_iv(file1, _PM.IVRPowerModel, ipopt_solver, deg=deg);

result_spmacdc_case2 = solve_sopf_acdc_iv(file2, _PM.IVRPowerModel, ipopt_solver, deg=deg);

# println("\n\n>>> OPF Results >>>")
# println(result_opf["primal_status"])
# print("Objective: ")
# print(result_opf["objective"])


# println("\n\n>>> SPM Results >>>")
# println(result_spm["primal_status"])
# print("Objective: ")
# print(result_spm["objective"])


# println("\n\n>>> ACDC Results >>>")
# println(result_acdc["primal_status"])
# print("Objective: ")
# print(result_acdc["objective"])


println("\n\n>>> SPMACDC Results - Case 1 >>>")
println(result_spmacdc_case1["primal_status"])
print("Objective: ")
print(result_spmacdc_case1["objective"])

println("\n\n>>> SPMACDC Results - Case 2 >>>")
println(result_spmacdc_case2["primal_status"])
print("Objective: ")
print(result_spmacdc_case2["objective"])




#=

idx = 1;
nw_idx = 1;

result_acdc["solution"]["gen"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["gen"]["$idx"]

result_acdc["solution"]["bus"]["$idx"]["vr"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["bus"]["$idx"]["vr"]

result_acdc["solution"]["branch"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["branch"]["$idx"]

result_acdc["solution"]["bus"]["$idx"]["vi"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["bus"]["$idx"]["vi"]

result_acdc["solution"]["busdc"]["$idx"]["vm"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["busdc"]["$idx"]["vm"]

result_acdc["solution"]["branchdc"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["branchdc"]["$idx"]



vm_sample1 = sample(result_spmacdc, "busdc", 1, "vm"; sample_size=1000)

q_sample1 = sample(result_spmacdc, "gen", 1, "qg"; sample_size=1000)

histogram(vm_sample1)

=#