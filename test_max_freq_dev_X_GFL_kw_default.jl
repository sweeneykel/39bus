using PowerSystems
using PowerSimulationsDynamics
using PowerFlows
using Sundials
using Plots
using LinearAlgebra
using Statistics
using Dates
using Printf
using Logging
using CSV
using DataFrames
using SparseArrays

const PSY  = PowerSystems
const PSID = PowerSimulationsDynamics

# Fixed Model Definitions
include("case39_matpower.jl")   # no changes made

# Updated parameter definitions and calculations
include("base_simulation_parameters.jl")
include("update_simulation_parameters.jl")

include("simulation.jl")
# A Single Simulation Run
#      - importing variables into local scope
#      - include("simulation_helper_funct.jl")
#      - include("model_parameters.jl") 

function run_experiment()
	simulation_test_type = "test_max_freq"
	independent_variable = "GFL_KW_default"

	test_name = splitext(basename(@__FILE__))[1]
	
	base_params =  get_default_params() 

	#for k in 0:35:5
	k = 35
	updated_base_params = merge(base_params, (GFL_KW_DEFAULT_test = float(k),))
	all_params = calculate_dependent_params(updated_base_params)
		
	run_simulation(k, #changing param value
					   all_params; 				   #updated Tuple of values
					   show_system = true, 
					   show_plot = false, 
					   simulation_test_type,	#test type
					   independent_variable	
					   )
	#end

end


run_experiment()