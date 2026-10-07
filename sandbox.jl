include("plot_enumerate_num_GFL.jl")

plt_ap_sandbox = plot_apex("sandbox_max_min_freq.csv")  # won't run without load summary passed
plt_na_sandbox = plot_nadir("sandbox_max_min_freq.csv")

savefig(plt_ap_sandbox, joinpath(pwd(),"_my_sandbox_test_apex_plot.png"))
savefig(plt_na_sandbox, joinpath(pwd(),"_my_sandbox_test_nadir_plot.png"))