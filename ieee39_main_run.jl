#!/usr/bin/env julia

#Main script for balanced EMT transient simulation of the IEEE 39-bus system.
#also Includes additinal script for  small-signal stability analysis.
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

#Case data for tge 39 bus 
include("case39_matpower.jl")

# Model and export parameters
include("ieee39_model_parameters.jl")

# All functions that run simulation and process results
include("ieee39_main_helper_funct.jl")

# test_name = splitext(basename(@__FILE__))[1]
test_name = "GFL_KW_DEFAULT"
changing_param = GFL_KW_DEFAULT

# Main run
function main(test_name::String, changing_param; show_system::Bool = true, show_plot::Bool = false, test_max_freq::Bool = true)
    validate_partitions!(SG_BUSES, GFM_BUSES, GFL_BUSES)

    outdir   = pwd()
    plotsdir = joinpath(outdir, "plots")
    csvdir   = joinpath(plotsdir, "csv")
    pfdir    = joinpath(plotsdir, "pf")
    freqpdfdir   = joinpath(plotsdir, "freq_bus30_39_pdfs")
    rotorpdfdir  = joinpath(plotsdir, "rotor_angle_bus30_39_pdfs")
    tgpdfdir     = joinpath(plotsdir, "tgtypei_bus_pdfs")
    pgpdfdir     = joinpath(plotsdir, "pg_bus_pdfs")
    gfmpdfdir    = joinpath(plotsdir, "gfm_power_pdfs")
    gflpdfdir    = joinpath(plotsdir, "gfl_pref_eff_pdfs")
    voltpdfdir   = joinpath(plotsdir, "voltage_bus30_39_pdfs")
    for d in (plotsdir, csvdir, pfdir, freqpdfdir, rotorpdfdir, tgpdfdir, pgpdfdir, gfmpdfdir, gflpdfdir, voltpdfdir)
        mkpath(d)
    end

    realized_steps = build_realized_load_steps_from_events()
    load_event_times = unique(sort(vcat(event_times_from_realized_steps(realized_steps), fault_event_times())))
    load_event_summary = LOAD_CHANGE_ENABLE ? event_summary_string(realized_steps) : fault_event_summary_string()
    if LOAD_CHANGE_ENABLE
        print_load_change_schedule(realized_steps)
    end
    if BUS_FAULT_ENABLE
        print_fault_disturbance_schedule()
    end

    casefile = joinpath(outdir, "modified_ieee39.m")
    write_modified_ieee39_case(casefile)
    @info "Wrote MATPOWER case" file = casefile

    sys = System(casefile)
    set_units_base_system!(sys, "SYSTEM_BASE")

    enforce_load_bases!(sys)

    enforce_min_line_r!(sys)

    add_dynamic_lines!(sys)

    pf_data = PowerFlowData(ACPowerFlow(), sys)
    solve_powerflow!(pf_data)
    write_pf_results(pf_data, pfdir; sys_base_mva = SYS_BASE_MVA)

    update_system_voltages!(sys, pf_data)

    loads_to_constant_impedance!(sys)

# Generator machine bases and SG machine/AVR parameters
    gen_mbase_map = generator_mbase_map()
    mach = sg_machine_params_by_bus()
    stable_avr = stable_avr_params()

    sync_static_gens_to_pf!(sys, pf_data; buses = collect(30:39))

    print_inverter_reference_summary(sys; gfm_buses = GFM_BUSES, gfl_buses = GFL_BUSES)

    for b in 30:39
        g = gen_at_bus(sys, b)
        try
            PSY.set_base_power!(g, gen_mbase_map[b])
        catch e
            @warn "Could not set generator base power" bus = b err = string(e)
        end
    end

    print_dynamic_reference_summary(sys; buses = collect(30:39))

# Attach dynamic models
    for b in SG_BUSES
        g = gen_at_bus(sys, b)
        attach_synchronous_generator!(sys, g;
            H = getH(b), D = getDsg(b), droop_R = R_SG_PU,
            machine_params = mach[b], avr_params = stable_avr)
    end

    for b in GFM_BUSES
        g = gen_at_bus(sys, b)
        attach_grid_forming_inverter!(sys, g;
            droop_R = R_GFM_PU, Ta = getGfmTa(b), kd = getGfmKd(b), kω = getGfmKw(b),
            kpv = GFM_DEFAULT_PARAMS.kpv, kiv = GFM_DEFAULT_PARAMS.kiv,
            kpc = GFM_DEFAULT_PARAMS.kpc, kic = GFM_DEFAULT_PARAMS.kic,
            ωad = GFM_DEFAULT_PARAMS.ωad)
    end

    for b in GFL_BUSES
        g = gen_at_bus(sys, b)
        attach_grid_following_inverter!(sys, g;
            kω_droop = min(getGflKw(b), GFL_KW_DROOP_CAP),
            kp_q = GFL_MAIN_CALL_OVERRIDES.kp_q,
            ki_q = GFL_MAIN_CALL_OVERRIDES.ki_q,
            ωf_q = GFL_MAIN_CALL_OVERRIDES.ωf_q,
            kic = GFL_MAIN_CALL_OVERRIDES.kic)
    end

    if show_system
        @info "System summary after attaching dynamic models"
        try
            show(sys)
        catch
        end
    end

# Build disturbances
    perts = PowerSimulationsDynamics.Perturbation[]
    if LOAD_CHANGE_ENABLE
        ensure_perturbable_loads_all_buses!(sys)
        for b in unique([s.bus for s in realized_steps])
            print_loads_on_bus(sys, b)
        end
        append!(perts, make_load_changes_from_event_schedule(sys, realized_steps))
    end
    if BUS_FAULT_ENABLE
        append!(perts, make_fault_perturbations(sys))
    end

    tspan = (0.0, TF_SIM)

    sim = Simulation!(ResidualModel, sys, pwd(), tspan, perts)

    try
        report_small_signal(sim)
    catch e
        @warn "small-signal analysis failed (continuing)" err = string(e)
    end

    execute!(sim, IDA(linear_solver = :Dense, max_order = SOLVER_MAX_ORDER);
             abstol = SOLVER_ABSTOL, reltol = SOLVER_RELTOL,
             dtinit = SOLVER_INIT_DT, dtmax = SOLVER_DTMAX)

    results = read_results(sim)

    gen_buses = collect(30:39)

# CSV exports
    for b in gen_buses
        save_frequency_csv_for_bus(results, sys, b, csvdir; t_step = minimum(load_event_times; init = 1.0))
        save_rotor_angle_csv_for_bus(results, sys, b, csvdir; t_step = minimum(load_event_times; init = 1.0))
        save_voltage_csv_for_bus(results, sys, b, csvdir; t_step = minimum(load_event_times; init = 1.0))
        save_terminal_pq_csv_for_bus(results, sys, b, csvdir; event_times = load_event_times)
    end

    for b in SG_BUSES
        save_mechanical_power_csv_for_bus(results, sys, b, csvdir; t_step = minimum(load_event_times; init = 1.0))
    end
    for b in GFM_BUSES
        save_gfm_power_csv_for_bus(results, sys, b, csvdir; t_step = minimum(load_event_times; init = 1.0))
    end
    for b in GFL_BUSES
        save_gfl_pref_eff_csv_for_bus(results, sys, b, csvdir; t_step = minimum(load_event_times; init = 1.0))
    end

#checknetwork P/Q cache
    net = nothing
    if SAVE_DIAGNOSTIC_NETWORK_PQ
        net = network_bus_pq_dev(results, sys, pf_data, realized_steps; buses = NETWORK_PQ_ALL_BUSES)
        print_network_pq_reconstruction_summary(net)
        for b in NETWORK_PQ_EXPORT_BUSES
            save_network_pq_csv(results, sys, pf_data, b, csvdir, realized_steps; net_cache = net)
        end
    end

#load-bus CSVs
    export_extra_load_bus_emt_csvs(results, sys, pf_data, csvdir, realized_steps;
                                   t_step = minimum(load_event_times; init = 1.0), net_cache = net)

    print_final_steady_state_freq(results, sys; buses = gen_buses)
    print_terminal_pq_final_summary(results, sys; buses = gen_buses)

# Max Frequency Study Per Iteration)
    if test_max_freq == true
        
        # Get the max frequency deviation for each bus (positive or negative)
        maxfreqdir = joinpath(plotsdir, "test_max_freq_x_$(test_name)")
        fname_bus = "each_bus_max_abs_freq_dev_x_$(test_name).csv"
        for b in gen_buses                  # 30...39             
            save_max_freq_for_bus(results, sys, b, maxfreqdir; fname=fname_bus, t_step = minimum(load_event_times; init = 1.0))
        end

        # Get the max frequency deviation for the entire system
        df = CSV.read(joinpath(maxfreqdir, fname_bus), DataFrame)
        row_idx = argmax(df.max_freq_dev)

        # format what should be stored
        df_save = DataFrame(independent_variable=changing_param, bus=df.bus_number[row_idx], max_abs_sys_freq_dev=df.max_freq_dev[row_idx])

        fname_sys = "system_max_abs_freq_deviation_x_$(test_name).csv"
        fpath = joinpath(maxfreqdir, fname_sys)

        # check if file exists and has content (>0 bytes)
        file_has_data = isfile(fpath) && filesize(fpath) > 0

        # append without re-writing the header if the file already contains data
        CSV.write(fpath, df_save; append = file_has_data, writeheader = !file_has_data)

        @info "Saved system max frequency deviation" file = fpath 

    end

#Plots
    plot_terminal_pq_panels(results, sys, PLOT_BUSES_EXT, plotsdir, load_event_times, load_event_summary)
    if net !== nothing
        plot_network_pq_panels(net, NETWORK_PQ_EXPORT_BUSES, plotsdir, realized_steps, load_event_times, load_event_summary)
    end
    plot_voltage_magnitude_panel(results, sys, gen_buses, load_event_times, load_event_summary, plotsdir)
    plot_voltage_deviation_panel(results, sys, gen_buses, load_event_times, load_event_summary, plotsdir)
    save_individual_voltage_plots(results, sys, gen_buses, voltpdfdir, load_event_times, load_event_summary;
                                  t_step = minimum(load_event_times; init = 0.25))

#Per-bus frequency plots
    for b in gen_buses
        series = frequency_deviation_hz_at_bus(results, sys, b)
        series === nothing && continue
        t, df, lab, sym = series
        if RESAMPLE_FOR_PLOTS
            tg = uniform_grid_with_event(TF_SIM, EXPORT_DT, minimum(load_event_times; init = 1.0))
            df = interp1_linear(t, df, tg)
            t = tg
        end
        p = plot(t, df; title = build_panel_title_julia(@sprintf("Frequency deviation at bus %d", b), load_event_summary),
            xlabel = "Time [s]", ylabel = "Δf [Hz]", lw = 2, color = bus_plot_color(b),
            label = lab, size = (1500, 700), titlefont = font(8), legend = :best, xlims = (0.0, TF_SIM), widen = false)
        hline!(p, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
        for (k, tev) in enumerate(load_event_times)
            vline!(p, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = (k == 1 ? "event time(s)" : ""))
        end
        savefig(p, joinpath(freqpdfdir, @sprintf("freq_bus%02d.pdf", b)))
        savefig(p, joinpath(freqpdfdir, @sprintf("freq_bus%02d.png", b)))
    end

# Per-bus rotor-angle plots
    for b in gen_buses
        series = rotor_angle_dev_rad_at_bus(results, sys, b; t_step = minimum(load_event_times; init = 1.0), for_plots = true)
        series === nothing && continue
        t, dδ, lab, sym = series
        p = plot(t, dδ; title = build_panel_title_julia(@sprintf("Rotor/angle deviation at bus %d", b), load_event_summary),
            xlabel = "Time [s]", ylabel = "Δδ [rad]", lw = 2, color = bus_plot_color(b),
            label = lab, size = (1500, 700), titlefont = font(8), legend = :best, xlims = (0.0, TF_SIM), widen = false)
        hline!(p, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
        for (k, tev) in enumerate(load_event_times)
            vline!(p, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = (k == 1 ? "event time(s)" : ""))
        end
        savefig(p, joinpath(rotorpdfdir, @sprintf("rotor_angle_bus%02d.pdf", b)))
        savefig(p, joinpath(rotorpdfdir, @sprintf("rotor_angle_bus%02d.png", b)))
    end

#Per-bus power plots
    for b in SG_BUSES
        plot_tgtypei_response!(results, sys, b, tgpdfdir, load_event_times, load_event_summary)
        out = get_mechanical_power_dev_pu_at_bus(results, sys, b; t_step = minimum(load_event_times; init = 1.0), for_plots = true)
        if out !== nothing
            t, ΔPm, Pm0, mBase, src_sym, label = out
            p = plot(t, ΔPm; title = build_panel_title_julia(@sprintf("SG mechanical power deviation at bus %d", b), load_event_summary),
                xlabel = "Time [s]", ylabel = "ΔP_m [pu, system base]", lw = 2, color = bus_plot_color(b),
                label = label, size = (1500, 700), titlefont = font(8), legend = :best, xlims = (0.0, TF_SIM), widen = false)
            hline!(p, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
            for (k, tev) in enumerate(load_event_times)
                vline!(p, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = (k == 1 ? "event time(s)" : ""))
            end
            savefig(p, joinpath(pgpdfdir, @sprintf("pg_bus%02d.pdf", b)))
            savefig(p, joinpath(pgpdfdir, @sprintf("pg_bus%02d.png", b)))
        end
    end

    for b in GFM_BUSES
        out = get_gfm_power_dev_pu_at_bus(results, sys, b; t_step = minimum(load_event_times; init = 1.0), for_plots = true)
        out === nothing && continue
        t, ΔP, P0, dev_name, label = out
        p = plot(t, ΔP; title = build_panel_title_julia(@sprintf("GFM active-power deviation at bus %d", b), load_event_summary),
            xlabel = "Time [s]", ylabel = "ΔP [pu, system base]", lw = 2, color = bus_plot_color(b),
            label = label, size = (1500, 700), titlefont = font(8), legend = :best, xlims = (0.0, TF_SIM), widen = false)
        hline!(p, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
        for (k, tev) in enumerate(load_event_times)
            vline!(p, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = (k == 1 ? "event time(s)" : ""))
        end
        savefig(p, joinpath(gfmpdfdir, @sprintf("gfm_power_bus%02d.pdf", b)))
        savefig(p, joinpath(gfmpdfdir, @sprintf("gfm_power_bus%02d.png", b)))
    end

    for b in GFL_BUSES
        out = get_gfl_pref_eff_series(results, sys, b; t_step = minimum(load_event_times; init = 1.0), for_plots = true)
        out === nothing && continue
        t, p_ref_eff, Δp_ref_eff, p_ref_const, ω_pll, Kω, ωsym, label = out
        p = plot(t, Δp_ref_eff; title = build_panel_title_julia(@sprintf("GFL effective P-ref deviation at bus %d", b), load_event_summary),
            xlabel = "Time [s]", ylabel = "Δp_ref_eff [pu]", lw = 2, color = bus_plot_color(b),
            label = label, size = (1500, 700), titlefont = font(8), legend = :best, xlims = (0.0, TF_SIM), widen = false)
        hline!(p, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
        for (k, tev) in enumerate(load_event_times)
            vline!(p, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = (k == 1 ? "event time(s)" : ""))
        end
        savefig(p, joinpath(gflpdfdir, @sprintf("gfl_pref_eff_bus%02d.pdf", b)))
        savefig(p, joinpath(gflpdfdir, @sprintf("gfl_pref_eff_bus%02d.png", b)))
    end

# Combined preview plot
    ordered = order_buses_for_plot_julia(PLOT_BUSES_EXT)

    p_freq = plot(title = build_panel_title_julia("Frequency deviation", load_event_summary),
        xlabel = "Time [s]", ylabel = "Δf [Hz]", legend = :outertopright)
    for b in ordered
        series = frequency_deviation_hz_at_bus(results, sys, b)
        series === nothing && continue
        t, df, lab, sym = series
        if RESAMPLE_FOR_PLOTS
            tg = uniform_grid_with_event(TF_SIM, EXPORT_DT, minimum(load_event_times; init = 1.0))
            df = interp1_linear(t, df, tg)
            t = tg
        end
        plot!(p_freq, t, df; lw = 2, color = bus_plot_color(b), label = @sprintf("Bus %d (%s)", b, model_label(b)))
    end
    hline!(p_freq, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    style_panel_plot!(p_freq)

    p_rotor = plot(title = build_panel_title_julia("Rotor angle deviation (Δδ from initial)", load_event_summary),
        xlabel = "Time [s]", ylabel = "Δδ [rad]", legend = :outertopright)
    for b in ordered
        series = rotor_angle_dev_rad_at_bus(results, sys, b; t_step = minimum(load_event_times; init = 1.0), for_plots = true)
        series === nothing && continue
        t, dδ, lab, sym = series
        plot!(p_rotor, t, dδ; lw = 2, color = bus_plot_color(b), label = @sprintf("Bus %d (%s)", b, model_label(b)))
    end
    hline!(p_rotor, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    style_panel_plot!(p_rotor)

    p_pm = plot(title = build_panel_title_julia("SG mechanical power deviation (ΔP_m from pre-event setpoint)", load_event_summary),
        xlabel = "Time [s]", ylabel = "ΔP_m [pu, system base]", legend = :outertopright)
    for b in order_buses_for_plot_julia(SG_BUSES)
        out = get_mechanical_power_dev_pu_at_bus(results, sys, b; t_step = minimum(load_event_times; init = 1.0), for_plots = true)
        out === nothing && continue
        t, ΔPm, Pm0, mBase, src_sym, label = out
        plot!(p_pm, t, ΔPm; lw = 2, color = bus_plot_color(b), label = @sprintf("Bus %d (SG)", b))
    end
    hline!(p_pm, [0.0]; lw = 0.8, ls = :dot, c = :gray, label = "")
    style_panel_plot!(p_pm)

    for (k, tev) in enumerate(load_event_times)
        vline!(p_freq,  [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = "")
        vline!(p_rotor, [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = "")
        vline!(p_pm,    [tev]; lw = 1, ls = :dash, c = :black, alpha = 0.35, label = "")
    end

    combined = plot(p_freq, p_rotor, p_pm; layout = (3, 1), size = (1500, 1500), dpi = 140,
        left_margin = 10Plots.mm, right_margin = 4Plots.mm, bottom_margin = 6Plots.mm)
    savefig(combined, joinpath(plotsdir, "ieee39_combined_freq_rotor_pm.png"))
    savefig(combined, joinpath(plotsdir, "ieee39_combined_freq_rotor_pm.pdf"))
    @info "Saved combined 3-panel preview" folder = plotsdir

    if show_plot
        display(combined)
    end

    @info "DONE. All CSVs and plots written." outdir = outdir plots = plotsdir csv = csvdir
    return sys, results
end



main(test_name, changing_param; show_system = true, show_plot = false, test_max_freq = true)

