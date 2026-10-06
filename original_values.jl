global test_name = splitext(basename(@__FILE__))[1]  #"testname"

GFL_PLL_KP_test = 0.03          #PLL Kp
GFL_PLL_KI_test = 1.5           #PLL Ki

# DON'T MODIFY!
const SYS_BASE_MVA = 100.0

#For larger sudden load steps, tighten solver tolerances for numerical stability
LOAD_CHANGE_EVENTS_test = [
    Dict(
        :buses        => [35],                            #Buses with load-step events
        :event_time_s => 0.1,                             #Load-step time[s]
        :p_pu_map     => Dict(35 => 70.0 / SYS_BASE_MVA), #Active-power load step[p.u.]
        :q_pu_map     => Dict(35 => 10.0 / SYS_BASE_MVA), #Reactive-power load step[p.u.]
    ),
]

"""
Removes one bus from SG_BUSES_test and adds it to GFL_BUSES_test
"""

global SG_BUSES_test  = [39]
global GFM_BUSES_test = [31, 36]
global GFL_BUSES_test = [30, 34, 32, 33, 35, 37, 38]

# 30;34;32;33;35;37;38,

include("ieee39_main_run.jl")



