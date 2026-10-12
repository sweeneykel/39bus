include("ieee39_main_run.jl")

test_name = splitext(basename(@__FILE__))[1]

const SYS_BASE_MVA = 100.0     #System base[MVA]

const SG_BUSES  = [33, 35, 37, 38, 39]      #Synchronous-generator buses
const GFM_BUSES = [31, 36]       #Grid-forming inverter buses
const GFL_BUSES = [30, 34, 32]         #Grid-following inverter buses
const PLOT_BUSES_EXT = vcat(SG_BUSES, GFM_BUSES, GFL_BUSES)              #Buses shown in plots

const GFL_PLL_KP = 0.03    #PLL Kp
const GFL_PLL_KI = 1.5          #PLL Ki

#For larger sudden load steps, tighten solver tolerances for numerical stability
const LOAD_CHANGE_EVENTS = [
    Dict(
        :buses        => [35],                            #Buses with load-step events
        :event_time_s => 0.1,                             #Load-step time[s]
        :p_pu_map     => Dict(35 => 70.0 / SYS_BASE_MVA), #Active-power load step[p.u.]
        :q_pu_map     => Dict(35 => 10.0 / SYS_BASE_MVA), #Reactive-power load step[p.u.]
    ),
]

#kω_droop = min(GFL_KW_DEFAULT, GFL_KW_DROOP_CAP),#Frequency-droop gain

# start:step:stop
for kw in 0:1:35
    GFL_KW_DEFAULT = kw         #GFL frequency-droop gain
    GFL_KW_map_ext = Dict{Int,Float64}(
        30 => GFL_KW_DEFAULT,                    #Bus 30 GFL frequency-droop gain
        31 => GFL_KW_DEFAULT,                    #Bus 31 GFL frequency-droop gain
    )

    indep_var_kw = GFL_KW_DEFAULT

    main(test_name, indep_var_kw; show_system = true, show_plot = false, test_max_freq = true)

end