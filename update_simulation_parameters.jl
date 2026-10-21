# responsibility: calculates parameters that depend on parameters in simulation_parameters.jl

function calculate_dependent_params(params)
    GFL_KW_map_ext_test = Dict{Int,Float64}(
							30 => params.GFL_KW_DEFAULT_test,              
							31 => params.GFL_KW_DEFAULT_test,                   
							) 

	LOAD_CHANGE_EVENTS_test = [
					Dict(
						:buses        => [35],                            
						:event_time_s => 0.1,                             
						:p_pu_map     => Dict(35 => 70.0 / params.SYS_BASE_MVA_test), 
						:q_pu_map     => Dict(35 => 10.0 / params.SYS_BASE_MVA_test), )
						]
						
    GFL_DEFAULT_PARAMS_test = (
						kp_p     = 0.8,#Active-power PI proportional gain
						ki_p     = 8.0,#Active-power PI integral gain
						ωz_p     = 80.0,#Active-power PI zero shaping frequency
						kp_q     = 0.8,#Reactive-power PI proportional gain
						ki_q     = 12.0,#Reactive-power PI integral gain
						ωf_q     = 120.0,#Reactive loop filter cutoff frequency
						kω_droop = min(params.GFL_KW_DEFAULT_test, params.GFL_KW_DROOP_CAP_test),#Frequency-droop gain
						kpc      = 1.2,#Current controller proportional gain
						kic      = 120.0,#Current controller integral gain
						kffv     = 0.0,#Voltage feedforward gain
						pll_ωlp  = 80.0,#PLL low-pass filter bandwidth
						pll_kp   = params.GFL_PLL_KP_test,#PLL proportional gain
						pll_ki   = params.GFL_PLL_KI_test,#PLL integral gain
						lf       = 0.08,#LCL inverter-side inductance
						rf       = 0.05,#LCL inverter-side resistance
						cf       = 0.074,#LCL shunt capacitance
						lg       = 0.2,#LCL grid-side inductance
						rg       = 0.04,#LCL grid-side resistance
	    				)
	
	return merge(params, (GFL_KW_map_ext_test = GFL_KW_map_ext_test,
                          LOAD_CHANGE_EVENTS_test = LOAD_CHANGE_EVENTS_test,
                          GFL_DEFAULT_PARAMS_test = GFL_DEFAULT_PARAMS_test ))  # NamedTuple
end





						