include("plot_enumerate_num_GFL.jl")

global test_name = splitext(basename(@__FILE__))[1]  #"testname"

GFL_PLL_KP_test = 0.03          #PLL Kp
GFL_PLL_KI_test = 1.5           #PLL Ki

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

SG_BUSES_orig  = [33, 35, 37, 38, 39]
global GFM_BUSES_test = [31, 36]
GFL_BUSES_orig = [30, 34, 32]

# Bus 39 must always remain synchronous
convertible_SG = filter(bus -> bus != 39, SG_BUSES_orig)

for n_convert in 1:length(convertible_SG)

    for buses_to_convert in combinations(convertible_SG, n_convert)

        global SG_BUSES_test  = copy(SG_BUSES_orig)
        global GFL_BUSES_test = copy(GFL_BUSES_orig)

        append!(GFL_BUSES_test, buses_to_convert)

        filter!(
            bus -> bus ∉ buses_to_convert,
            SG_BUSES_test
        )
        
        println("Current iteration:----------------------------")
        println("SG_BUSES_test: ", SG_BUSES_test)
        println("GFL_BUSES_test: ", GFL_BUSES_test)
        println("GFM_BUSES_test: ", GFM_BUSES_test)
        
        #=
        print("Press Enter to continue, to type q to quit: ")
        user_input = readline()

        if lowercase(strip(user_input)) == "q"
            println("Exiting loop.")
            break
        end
        =#
        
        include("ieee39_main_run.jl")
    end
end

this_test_dir = joinpath(joinpath(pwd(), "plots"), test_name)
path_to_freq_csv = joinpath(this_test_dir, "max_min_freq.csv")

p_apex = plot_apex(path_to_freq_csv)
p_nadir = plot_nadir(path_to_freq_csv)

savefig(p_apex, joinpath(this_test_dir, "apex.png"))
savefig(p_nadir, joinpath(this_test_dir, "nadir.png"))