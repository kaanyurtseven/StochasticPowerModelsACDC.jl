using Pkg
Pkg.activate(".")

using JuMP
using Ipopt
using PowerModels
using StochasticPowerModelsACDC
using PowerModelsACDC
using Plots
using Memento
using InfrastructureModels
using Random, Distributions

using XLSX

#Constants 
const _PM = PowerModels
const _SPMACDC = StochasticPowerModelsACDC
const _PMACDC = PowerModelsACDC

Memento.setlevel!(Memento.getlogger(StochasticPowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")

Random.seed!(1234)

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "max_iter"=>3000, "sb"=>"yes")

#Monte-Carlo Simulation inputs
MC_size = 5;

#PV inputs
pen_level_start = 0.00
pen_level_step = 0.05
pen_level_end = 0.10

#Uncertain input parameters
load_std1 = 0.1;
load_std2 = 0.15;
α_g = 5;
β_g = 2;
PV_pf = 0.9;
global load_bus_gauss_1 = 3
global load_bus_gauss_2 = 4

#Desired Chance-Constraint level
global CC_level = 0.90

#Case file and data reading
case = "case5_ACDC.m"
file  = joinpath(BASE_DIR, "test/data/matpower/Case Studies", case)
data = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)
data_org = data #Store the original data

#Necessary initializations
load_means = Dict()
load_dist_gauss1 = Dict()
load_dist_gauss2 = Dict()
load_sample_gauss1 = Dict()
load_sample_gauss2 = Dict()
result_acdc = Dict()
feas_cost = Dict()
PV_sample = Dict()
global solve_case = true
global total_load = sum([data_org["load"]["$i"]["pd"] for i=1:length(data["load"])])

#Generate load samples
for i in keys(data["load"])

    load_means[i] = deepcopy(data["load"][i]["pd"])
    
    load_dist_gauss1[i] = Normal(load_means[i], load_means[i]*load_std1)
    load_sample_gauss1[i] = rand(load_dist_gauss1[i], MC_size)

    load_dist_gauss2[i] = Normal(load_means[i], load_means[i]*load_std2)
    load_sample_gauss2[i] = rand(load_dist_gauss2[i], MC_size)

end

#Setting for SPMACDC solver
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)

#Penetration level iteration
for pen_level_dummy = pen_level_start:pen_level_step:pen_level_end

    global pen_level = pen_level_dummy; #Necessary step for storing pen_level dictionaries
    
    #Initiliaze dictionaries
    result_acdc["$pen_level"] = Dict()
    result_acdc["$pen_level"]["MC Results"] = Dict()
    feas_cost["$pen_level"] = Dict()
    PV_sample["$pen_level"] = Dict()

    #Print the penetration level
    println("Penetration Level = $pen_level")
    
    PV_size = pen_level * total_load / length(data["load"]) #Calculate PV size for each load bus
    
    if pen_level == 0 #The distinction must be made for 0 penetration

        for MC_sample = 1:MC_size
            # println("   MC_sample = $MC_sample")
            for (i, load) in data["load"]
                bus = load["load_bus"]
                if bus == load_bus_gauss_1
    
                    data["load"][i]["pd"] = load_sample_gauss1[i][MC_sample]
                
                elseif bus == load_bus_gauss_2

                    data["load"][i]["pd"] = load_sample_gauss2[i][MC_sample]

                end
    
            end
    
            if solve_case
                global result_acdc["$pen_level"]["MC Results"]["$MC_sample"] = _PMACDC.run_acdcopf(data, _PM.ACPPowerModel, ipopt_solver; setting = s)
            end
        end

    else
        #Generate PV samples for >0 penetration
        PV_dist = Beta(α_g, β_g) * PV_size
        PV_sample["$pen_level"] = rand(PV_dist, MC_size)

        for MC_sample = 1:MC_size
            # println("   MC_sample = $MC_sample")
            for (i, load) in data["load"]
                bus = load["load_bus"]
                if bus == load_bus_gauss_1

                    data["load"][i]["pd"] = load_sample_gauss1[i][MC_sample] - (PV_pf * PV_sample["$pen_level"][MC_sample])

                elseif bus == load_bus_gauss_2
                
                    data["load"][i]["pd"] = load_sample_gauss2[i][MC_sample] - (PV_pf * PV_sample["$pen_level"][MC_sample])

                else

                    data["load"][i]["pd"] = data_org["load"][i]["pd"] - (PV_pf * PV_sample["$pen_level"][MC_sample])

                end

                data["load"][i]["qd"] = data_org["load"][i]["qd"] - (1-PV_pf) * PV_sample["$pen_level"][MC_sample]

            end

            if solve_case
                global result_acdc["$pen_level"]["MC Results"]["$MC_sample"] = _PMACDC.run_acdcopf(data, _PM.ACPPowerModel, ipopt_solver; setting = s)
            end
        end

    end

    

    if solve_case
        global feas_ctr = 0

        for i =1:MC_size

            if string(result_acdc["$pen_level"]["MC Results"]["$i"]["primal_status"]) == "FEASIBLE_POINT"
                
                #Count feasible cases for chance-constraint calculation
                global feas_ctr += 1

                #Store the objective values of feasible results
                global feas_cost["$pen_level"]["$i"] = result_acdc["$pen_level"]["MC Results"]["$i"]["objective"]

            end
            
        end
                
        result_acdc["$pen_level"]["CC"] = feas_ctr / MC_size; #Calculate the chance_constraint level

        if feas_ctr > 0
            
            #When calculating average cost only include the results from the feasible ones
            result_acdc["$pen_level"]["Expected Cost"] = mean([feas_cost["$pen_level"]["$i"] for i in keys(feas_cost["$pen_level"])])

        else 
            result_acdc["$pen_level"]["Expected Cost"] = 0
        end

        if result_acdc["$pen_level"]["CC"] < CC_level
            global solve_case = false #If violation in the CC level starts, stop the simulations
        end

    else
        delete!(result_acdc, "$pen_level") #Delete the initialized empty dictionary
    end



end

# Reporting
if contains(case, "ACDC")
    case = "ACDC";
else
    case = "AC";
end

exp_cost_export = Dict()
[exp_cost_export["$pen_level"] = result_acdc["$pen_level"]["Expected Cost"] for pen_level in keys(result_acdc)]
CC_export = Dict()
[CC_export["$pen_level"] = result_acdc["$pen_level"]["CC"] for pen_level in keys(result_acdc)]

file_name = "Results\\MC Results $case.xlsx"
fid    = XLSX.openxlsx(file_name, mode="w")
header = ["Penetration Level"; "Expected Cost"; "CC";]
export_data   = [[collect(keys(result_acdc))]; [collect(values(exp_cost_export))]; [collect(values(CC_export))]]
XLSX.writetable(file_name, export_data, header)

