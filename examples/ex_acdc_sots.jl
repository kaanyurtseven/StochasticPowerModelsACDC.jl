using Pkg
Pkg.activate(".")

using JuMP
using Ipopt
# using Juniper
using PowerModels
using StochasticPowerModelsACDC
using PowerModelsACDC
using Memento
using InfrastructureModels
using Statistics
# using Plots
using PolyChaos

using XLSX

#Constants
const _PM = PowerModels
const _SPMACDC = StochasticPowerModelsACDC
const _PMACDC = PowerModelsACDC

Memento.setlevel!(Memento.getlogger(StochasticPowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PolyChaos), "error")

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "max_cpu_time" => 300.0, "max_iter"=>3000, "sb"=>"yes")


#gPC degree input
deg  = 2

#PV inputs
pen_level_start = 0.3
pen_level_step = 0.01
pen_level_end = 0.3

#Necessary initializations
obj_AC = Dict()
stat_AC = Dict()
time_AC = Dict()
obj_ACDC = Dict()
stat_ACDC = Dict()
time_ACDC = Dict()
obj_sopf = Dict()
stat_sopf = Dict()
time_sopf = Dict()
p_size_dict = Dict()

global result_sopf = Dict()
global result_sots_AC = Dict()
global result_sots_ACDC = Dict()

global feas_ctr_sopf = 0
global feas_ctr_AC_sots = 0
global feas_ctr_ACDC_sots = 0

solve_sopf = true
solve_case_AC = true
solve_case_ACDC = true

case = "case5_ACDC_mod_SPMACDC_95cc_RES.m"

file  = joinpath(BASE_DIR, "test/data/matpower/", case)

set = Dict("output" => Dict("duals" => false, "branch_flows" => true), "conv_losses_mp" => true)
data = _PM.parse_file(file)
data_AC = _PM.parse_file(file)
data_ACDC = _PM.parse_file(file)
data_OPF = _PM.parse_file(file)


for (b, branch) in data["branch"]
    data_AC["branch"][b]["br_status_initial"] = 1
    data_ACDC["branch"][b]["br_status_initial"] = 1
    data_OPF["branch"][b]["br_status"] = 1
end

for (b, branchdc) in data["branchdc"]
    data_ACDC["branchdc"][b]["br_status_initial"] = 1
end
    
_SPMACDC.process_additional_data!(data_AC)
_SPMACDC.process_additional_data!(data_ACDC)
_SPMACDC.process_additional_data!(data_OPF)

#Penetration level iteration
for pen_level = pen_level_start:pen_level_step:pen_level_end

    println("Penetration Level = $pen_level")

    if solve_case_AC
        
        total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
        p_size = pen_level * total_load / length(data["RES"]) #Calculate PV size for each load bus
        p_size_dict[pen_level] = p_size
    

        println("   AC SOTS: Solution progress: Solving...")
        global result_sots_AC[pen_level] = _SPMACDC.solve_sots_acdc_AC(data_AC, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
        
       
        #Store necessary values for reporting
        obj_AC[pen_level] = result_sots_AC[pen_level]["objective"] 
        stat_AC[pen_level] = string(result_sots_AC[pen_level]["primal_status"])
        time_AC[pen_level] = string(result_sots_AC[pen_level]["solve_time"])

        if string(result_sots_AC[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_sots_AC[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
            global feas_ctr_AC_sots += 1
        else
            global feas_ctr_AC_sots = 0
        end

        println("   AC SOTS: Solution progress: Solved! (", string(result_sots_AC[pen_level]["primal_status"]), ")")

        if feas_ctr_AC_sots >= 2
            println("   Case reached infeasibility on penetration level of $pen_level.")
            global solve_case_AC = false
        end

        # global sots_AC_vm_2 = sample(result_sots_AC[pen_level], "bus", 2, "vm"; sample_size=10000);
        # global sots_AC_vm_3 = sample(result_sots_AC[pen_level], "bus", 3, "vm"; sample_size=10000);
        # global sots_AC_vm_4 = sample(result_sots_AC[pen_level], "bus", 4, "vm"; sample_size=10000);
        # global sots_AC_vm_5 = sample(result_sots_AC[pen_level], "bus", 5, "vm"; sample_size=10000);
        # global sots_AC_pg_10 = sample(result_sots_AC[pen_level], "gen", 10, "pg"; sample_size=10000);
        # global sots_AC_csr_8 = sample(result_sots_AC[pen_level], "branch", 8, "csr_fr_on_off"; sample_size=10000);

    end

    if solve_case_ACDC
         
        println("   AC/DC SOTS: Solution progress: Solving...")

        total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
        p_size = pen_level * total_load / length(data["RES"]) #Calculate PV size for each load bus
        p_size_dict[pen_level] = p_size

        global result_sots_ACDC[pen_level] = _SPMACDC.solve_sots_acdc(data_ACDC, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
        
       
        #Store necessary values for reporting
        obj_ACDC[pen_level] = result_sots_ACDC[pen_level]["objective"] 
        stat_ACDC[pen_level] = string(result_sots_ACDC[pen_level]["primal_status"])
        time_ACDC[pen_level] = string(result_sots_ACDC[pen_level]["solve_time"])

        if string(result_sots_ACDC[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_sots_ACDC[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
            global feas_ctr_ACDC_sots += 1
        else
            global feas_ctr_ACDC_sots = 0
        end

        println("   AC/DC SOTS: Solution progress: Solved! (", string(result_sots_ACDC[pen_level]["primal_status"]), ")")

        if feas_ctr_ACDC_sots >= 2
            println("   Case reached infeasibility on penetration level of $pen_level.")
            global solve_case_ACDC = false
        end

    end

    # for (b, branch) in result_sots_ACDC[pen_level]["solution"]["nw"]["1"]["branch"]

    #     # data["branch"][b]["br_status_initial"] = branch["br_status"]

    #     if branch["br_status"] < 0.5
    #         data_OPF["branch"][b]["br_status"] = 0
    #     else
    #         data_OPF["branch"][b]["br_status"] = 1
    #     end
    # end

    
    # _SPMACDC.process_additional_data!(data_OPF)

    if solve_sopf
        
        total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
        p_size = pen_level * total_load / length(data["RES"]) #Calculate PV size for each load bus
        p_size_dict[pen_level] = p_size
    
    
        println("   SOPF: Solution progress: Solving...")
        global result_sopf[pen_level] = _SPMACDC.solve_sopf_acdc(data_OPF, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
        
       
        #Store necessary values for reporting
        obj_sopf[pen_level] = result_sopf[pen_level]["objective"] 
        stat_sopf[pen_level] = string(result_sopf[pen_level]["primal_status"])
        time_sopf[pen_level] = string(result_sopf[pen_level]["solve_time"])


        if string(result_sopf[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_sopf[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
            global feas_ctr_sopf += 1
        else
            global feas_ctr_sopf = 0
        end

        println("   SOPF: Solution progress: Solved! (", string(result_sopf[pen_level]["primal_status"]), ")")

        if feas_ctr_sopf >= 50
            println("   Case reached infeasibility on penetration level of $pen_level.")
            global solve_sopf = false
        end

    end

end

    
#Show results on the terminal
pen_level = pen_level_start

if solve_sopf || solve_case_AC || solve_case_ACDC
    println("\n\n>>> Penetration Level: $(pen_level*100)% >>>")
end

if solve_sopf
    println("\n>>> Results - SOPF >>>")
    println(result_sopf[pen_level]["primal_status"])
    print("Objective: ")
    print(result_sopf[pen_level]["objective"])
    print("\nSolve Time: ")
    print(result_sopf[pen_level]["solve_time"])
end

if solve_case_AC
    println("\n\n>>> Results - AC SOTS >>>")
    println(result_sots_AC[pen_level]["primal_status"])
    print("Objective: ")
    print(result_sots_AC[pen_level]["objective"])
    print("\nSolve Time: ")
    print(result_sots_AC[pen_level]["solve_time"])
end

if solve_case_ACDC
    println("\n\n>>> Results - ACDC SOTS >>>")
    println(result_sots_ACDC[pen_level]["primal_status"])
    print("Objective: ")
    print(result_sots_ACDC[pen_level]["objective"])
    print("\nSolve Time: ")
    print(result_sots_ACDC[pen_level]["solve_time"])
end
