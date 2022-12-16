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
#ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_iter"=>3000, "fixed_variable_treatment" => "relax_bounds")
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_iter"=>3000)

# input
deg  = 2
obj = Dict()
status = Dict()
#for p_size = 1:5
p_size = 50

case = "case5_acdc_SPMACDC.m"

#case = "deneme5case.m"

#case = "case39_acdc_SPMACDC.m"

#case = "case67_acdc_scopf_SPMACDC.m"

#case = "case67_ACDC_SPMACDC.m"


file  = joinpath(BASE_DIR, "test/data/matpower", case)

 
#result_opf = _PM.solve_opf(file, _PM.ACPPowerModel, ipopt_solver)

#result_spm = solve_sopf_iv(file, _PM.IVRPowerModel, ipopt_solver, deg=deg);

# s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)
# result_acdc = _PMACDC.run_acdcopf_iv(file, _PM.IVRPowerModel, ipopt_solver; setting = s)

result_spmacdc = solve_sopf_acdc_PV(file, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);

# obj[p_size] = result_spmacdc["objective"] 
# status[p_size] = result_spmacdc["primal_status"]
# end

# plot(obj)

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


println("\n\n>>> SPMACDC Results >>>")
println(result_spmacdc["primal_status"])
print("Objective: ")
print(result_spmacdc["objective"])


idx = 1;
nw_idx = 1;

# result_acdc["solution"]["gen"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["gen"]["$idx"]

# result_acdc["solution"]["bus"]["$idx"]["vr"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["bus"]["$idx"]["vr"]

# result_acdc["solution"]["branch"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["branch"]["$idx"]

# result_acdc["solution"]["bus"]["$idx"]["vi"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["bus"]["$idx"]["vi"]

# result_acdc["solution"]["busdc"]["$idx"]["vm"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["busdc"]["$idx"]["vm"]

# result_acdc["solution"]["branchdc"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["branchdc"]["$idx"]

result_spmacdc["solution"]["nw"]["$nw_idx"]["PV"]["$idx"]



q_sample1 = sample(result_spmacdc, "PV", 1, "cid_pv"; sample_size=1000)
q_sample2 = sample(result_spmacdc, "PV", 1, "crd_pv"; sample_size=1000)


dummy["solution"]["nw"]["1"]["PV"]["1"] = sdata["nw"]["1"]["PV"]["1"]
q_sample = sample(sdata, "PV", 1, "cid_pv"; sample_size=1000)


sample3 = _PCE.samplePCE(1000, [sdata["nw"]["$i"]["PV"]["4"]["pd"] for i=1:6], sdata["mop"])

histogram(sample3)
histogram(q_sample1)
histogram!(q_sample2)

vm_sample1 = sample(result_spmacdc, "bus", 2, "vr"; sample_size=1000)
pg_sample1 = sample(result_spmacdc, "gen", 1, "pg"; sample_size=1000)
crd_sample1 = sample(result_spmacdc, "branch", 2, "cr_fr"; sample_size=1000)
crd_sample2 = sample(result_spmacdc, "PV", 4, "crd_pv"; sample_size=1000)


histogram(vm_sample1)
histogram(crd_sample1)
histogram(crd_sample2)
histogram(pg_sample1)

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


pg_coeff = pce_coeff(result_spmacdc, "gen", 1, "pg") 

pg_coeff = pce_coeff(result_spmacdc, "gen", 1, "pg") 

pg_coeff = pce_coeff(result_spmacdc, "gen", 1, "crg") 

pg_coeff = pce_coeff(result_spmacdc, "gen", 1, "qg") 

pg_coeff = pce_coeff(result_spmacdc, "gen", 1, "cig") 



vr_sample1 = sample(result_spmacdc, "bus", idx, "vr"; sample_size=100)
histogram(vr_sample1)

pg_sample = sample(result_spmacdc, "gen", 2, "pg"; sample_size=100)

cur_sample1 = sample(result_spmacdc, "branch", 1, "cr_to"; sample_size=100)
cur_sample2 = sample(result_spmacdc, "branch", 1, "cr_fr"; sample_size=100)

dccur_sample1 = sample(result_spmacdc, "branchdc", 1, "cr_to"; sample_size=100)
dccur_sample2 = sample(result_spmacdc, "branchdc", 1, "cr_fr"; sample_size=100)

vm_sample1 = sample(result_spmacdc, "busdc", 1, "vm"; sample_size=100)
vm_sample2 = sample(result_spmacdc, "busdc", 2, "vm"; sample_size=100)
vm_sample3 = sample(result_spmacdc, "busdc", 3, "vm"; sample_size=100)

histogram(vm_sample1)


histogram(pg_sample)

histogram(cur_sample1)

histogram!(-cur_sample2)


=#


#=


@assert result_ivr["termination_status"] == PM.LOCALLY_SOLVED
#@assert result_acr["termination_status"] == PM.LOCALLY_SOLVED
#the optimal objective values (expectation) are 
obj_ivr = result_ivr["objective"] 
#obj_acr = result_acr["objective"] 
#@assert obj_ivr ≈ obj_acr

# print variables for all polynomial indices k
SPM.print_summary(result_ivr["solution"])

# print variables for a specific index k
k=1
SPM.print_summary(result_ivr["solution"]["nw"]["$k"])

# get polynomial chaos coefficients for specific component
pg_coeff = pce_coeff(result_ivr, "gen", 1, "pg") 

# obtain 10 samples of the generator active power output variable
pg_sample = sample(result_ivr, "gen", 1, "pg"; sample_size=10) 

# obtain an kernel density estimate of the generator active power output variable
pg_density = density(result_ivr, "gen", 1, "pg"; sample_size=10) 

#-----------------------------------
# alternatively, you can first read in PowerModels dict, 
# from a file with stochastic data extensions:
data  = PM.parse_file(file)

result_ivr2 = solve_sopf_iv(data, PM.IVRPowerModel, ipopt_solver; deg=deg)
@assert result_ivr2["termination_status"] == PM.LOCALLY_SOLVED
obj_ivr2 = result_ivr2["objective"]

#-----------------------------------
# finally, we can also build the multinetwork dictionary here.
# we run the stochastic replicate function, introduces polynomial index k
sdata = SPM.build_stochastic_data(data, deg)

# run the reduced IVR with auxiliary variables
result_ivr3 = PM.solve_model(sdata, PM.IVRPowerModel, ipopt_solver, SPM.build_sopf_iv; multinetwork=true, solution_processors=[PM.sol_data_model!])
@assert result_ivr3["termination_status"] == PM.LOCALLY_SOLVED
obj_ivr3 = result_ivr3["objective"]

@assert obj_ivr ≈ obj_ivr2 ≈ obj_ivr3




=#