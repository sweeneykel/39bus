# responsibility: all parameters that could be modified to test in simulation
function get_default_params()
	return ( 
		SYS_BASE_MVA_test = 100.0,
		SG_BUSES_test  = [33, 35, 37, 38, 39],      
		GFM_BUSES_test = [31, 36],       
		GFL_BUSES_test = [30, 34, 32], 
		GFL_PLL_KP_test = 0.03,
		GFL_PLL_KI_test = 1.5,
		GFL_KW_DEFAULT_test = 0.0,         
		GFL_KW_DROOP_CAP_test = 35.0,
	)
	end