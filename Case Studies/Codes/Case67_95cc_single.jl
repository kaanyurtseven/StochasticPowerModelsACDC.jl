using Pkg
Pkg.activate(".")

using JuMP
using Ipopt
using PowerModels
using StochasticPowerModelsACDC
using PowerModelsACDC
using Memento
using InfrastructureModels

using XLSX

#Constants
const _PM = PowerModels
const _SPMACDC = StochasticPowerModelsACDC
const _PMACDC = PowerModelsACDC

Memento.setlevel!(Memento.getlogger(StochasticPowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "max_iter"=>5000, "sb"=>"yes", "fixed_variable_treatment" => "relax_bounds")

#gPC degree input
deg  = 2

#PV inputs
pen_level_start = 0.00
pen_level_step = 0.01
pen_level_end = 2.00

#Report decision
Report = true

#Necessary initializations
obj_case1 = Dict()
obj_case2 = Dict()
stat_case1 = Dict()
stat_case2 = Dict()
p_size_dict = Dict()

global feas_ctr1 = 0
global feas_ctr2 = 0

global solve_case1 = true
global solve_case2 = true

#Case file and data reading
case1 = "case67_AC_SPMACDC_95cc_single.m"
file1  = joinpath(BASE_DIR, "Case Studies/Data", case1)

case2 = "case67_ACDC_SPMACDC_95cc_single.m"
# case2 = "case5_ACDC_SPMACDC.m"
file2  = joinpath(BASE_DIR, "Case Studies/Data", case2)

set = Dict("output" => Dict("duals" => false, "branch_flows" => true), "conv_losses_mp" => true)
data = _PM.parse_file(file2)

#Penetration level iteration
for pen_level = pen_level_start:pen_level_step:pen_level_end
    
    if solve_case1 || solve_case2
        println("Penetration Level = $pen_level")
        
        total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
        p_size = pen_level * total_load / length(data["load"]) #Calculate PV size for each load bus
        p_size_dict[pen_level] = p_size
        
        if solve_case1
            println("   Case 1 solution progress: Solving...")
            global result_spmacdc_case1 = solve_sopf_acdc_PV(file1, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
            
            #Store necessary values for reporting
            obj_case1[pen_level] = result_spmacdc_case1["objective"]
            stat_case1[pen_level] = string(result_spmacdc_case1["primal_status"])

            if string(result_spmacdc_case1["primal_status"]) != "FEASIBLE_POINT" && string(result_spmacdc_case1["primal_status"]) !="NEARLY_FEASIBLE_POINT"
                global feas_ctr1 += 1
            else
                global feas_ctr1 = 0
            end

            println("   Case 1 solution progress: Solved!")

            if feas_ctr1 >= 2
                global solve_case1 = false
                println("   Case 1 reached infeasibility on penetration level of $pen_level.")
            end

        end

        if solve_case2
            println("   Case 2 solution progress: Solving...")
            global result_spmacdc_case2 = solve_sopf_acdc_PV(file2, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
            
            #Store necessary values for reporting
            obj_case2[pen_level] = result_spmacdc_case2["objective"] 
            stat_case2[pen_level] = string(result_spmacdc_case2["primal_status"])

            if string(result_spmacdc_case2["primal_status"]) != "FEASIBLE_POINT" && string(result_spmacdc_case2["primal_status"]) !="NEARLY_FEASIBLE_POINT"
                global feas_ctr2 += 1
            else
                global feas_ctr2 = 0
            end

            println("   Case 2 solution progress: Solved!")

            if feas_ctr2 >= 2
                println("   Case 2 reached infeasibility on penetration level of $pen_level.")
                global solve_case2 = false
            end

        end

    end

end


#Show results on the terminal
# println("\n\n>>> SPMACDC Results - Case 1 >>>")
# println(result_spmacdc_case1["primal_status"])
# print("Objective: ")
# print(result_spmacdc_case1["objective"])

# println("\n\n>>> SPMACDC Results - Case 2 >>>")
# println(result_spmacdc_case2["primal_status"])
# print("Objective: ")
# print(result_spmacdc_case2["objective"])

# Reporting
if Report
    file_name1 = "Results\\gPC Results - $case1.xlsx"
    fid    = XLSX.openxlsx(file_name1, mode="w")
    header = ["Penetration Level";"Objective Value";"Status"]
    data   = [[collect(keys(obj_case1))];[collect(values(obj_case1))];[collect(values(stat_case1))] ]
    XLSX.writetable(file_name1, data,header)

    file_name2 = "Results\\gPC Results - $case2.xlsx"
    fid    = XLSX.openxlsx(file_name2, mode="w")
    header = ["Penetration Level";"Objective Value";"Status"]
    data   = [[collect(keys(obj_case2))];[collect(values(obj_case2))];[collect(values(stat_case2))] ]
    XLSX.writetable(file_name2, data,header)
end



