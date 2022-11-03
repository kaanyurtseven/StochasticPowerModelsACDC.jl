using Pkg

Pkg.activate(".")


using JuMP
using Ipopt
using PowerModels
using StochasticPowerModelsACDC
using PowerModelsACDC

# constants 
const _PM = PowerModels
const _SPMACDC = StochasticPowerModelsACDC
const _PMACDC = PowerModelsACDC

# solvers
ipopt_solver = Ipopt.Optimizer

# input
deg  = 2

case = "case5_SPMACDC_spm.m"


file  = joinpath(BASE_DIR, "test/data/matpower", case)

result_ivr = solve_sopf_iv(file, _PM.IVRPowerModel, ipopt_solver, deg=deg);

result_ivracdc = solve_sopf_acdc_iv(file, _PM.IVRPowerModel, ipopt_solver, deg=deg);


#_SPMACDC.print_summary(result_ivr["solution"])
#_SPMACDC.print_summary(result_ivracdc["solution"])


println("\n\n>>> SPM Results >>>")
println(result_ivr["primal_status"])
print("Objective: ")
print(result_ivr["objective"])

println("\n\n>>> SPMACDC Results >>>")
println(result_ivracdc["primal_status"])
print("Objective: ")
print(result_ivracdc["objective"])









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