using Combinatorics

using CSV, DataFrames
using Plots

function plot_apex(csv_file::String)

    df = CSV.read(csv_file, DataFrame)

    # DataFrame(list_GFL=list_GFL_busses,  bus_number=bus,  apex=max_val, nadir=min_val)


    # ---------------------------------------------------------
    # Prepare GFL configuration labels
    # ---------------------------------------------------------

    # Example CSV entry:
    # "[30, 34, 32, 33]"
    #
    # Convert to:
    # "30, 34, 32, 33"

    df.GFL_label = replace.(df.list_GFL, "[" => "", "]" => "")

    # Count how many buses are in each GFL configuration
    df.GFL_count = [
        length(split(label, ","))
        for label in df.GFL_label
    ]

    # Find unique GFL configurations
    groups = unique(
        select(df, [:GFL_label, :GFL_count])
    )

    # Sort first by number of GFL buses,
    # then by the actual configuration
    sort!(groups, [:GFL_count, :GFL_label])

    # ---------------------------------------------------------
    # Give each GFL configuration its own x-axis position
    # ---------------------------------------------------------

    group_to_x = Dict(
        row.GFL_label => i
        for (i, row) in enumerate(eachrow(groups))
    )

    df.x_position = [
        group_to_x[label]
        for label in df.GFL_label
    ]

    tick_positions = 1:nrow(groups)

    tick_labels = groups.GFL_label

    # ---------------------------------------------------------
    # Create plot
    # ---------------------------------------------------------

    p = plot(
        xlabel = "GFL Bus Configuration",
        ylabel = "Δf [Hz]",

        xticks = (
            tick_positions,
            tick_labels
        ),

        xrotation = 45,

        legend = :outerright,

        bottom_margin = 10Plots.mm,
        top_margin = 10Plots.mm,

        size = (1200, 700)
    )

    # ---------------------------------------------------------
    # Fixed colors based ONLY on generator type
    # ---------------------------------------------------------

    GFL_COLOR = :purple
    SG_COLOR  = :green
    GFM_COLOR = :orange

    # ---------------------------------------------------------
    # Plot data
    # ---------------------------------------------------------

    for bus in sort(unique(df.bus_number))

        bus_df = filter(
            row -> row.bus_number == bus,
            df
        )

        for row in eachrow(bus_df)

            # Convert the GFL configuration string back
            # into bus numbers
            #
            # "30, 34, 32, 33"
            # ->
            # [30, 34, 32, 33]

            gfl_buses = parse.(
                Int,
                strip.(split(row.GFL_label, ","))
            )

            # Determine generator type for THIS configuration

            if bus in gfl_buses

                point_color = GFL_COLOR

            elseif bus in GFM_BUSES_orig

                point_color = GFM_COLOR

            else

                point_color = SG_COLOR

            end

            scatter!(
                p,
                [row.x_position],
                [row.apex],

                marker = :circle,
                markersize = 6,

                color = point_color,

                label = false
            )

        end
    end


    # ---------------------------------------------------------
    # Add legend entries once
    # ---------------------------------------------------------

    scatter!(
        p,
        [NaN],
        [NaN],

        marker = :circle,
        markersize = 6,

        color = GFL_COLOR,

        label = "GFL"
    )

    scatter!(
        p,
        [NaN],
        [NaN],

        marker = :circle,
        markersize = 6,

        color = SG_COLOR,

        label = "SG"
    )

    scatter!(
        p,
        [NaN],
        [NaN],

        marker = :circle,
        markersize = 6,

        color = GFM_COLOR,

        label = "GFM"
    )


    # ---------------------------------------------------------
    # Add grouping based on number of GFL buses
    # ---------------------------------------------------------

    counts = sort(
        unique(groups.GFL_count)
    )

    # Current y-axis limits
    ylims_current = ylims(p)

    y_range =
        ylims_current[2] -
        ylims_current[1]

    # Leave some room above data for labels
    new_ymax =
        ylims_current[2] +
        0.10 * y_range

    ylims!(
        p,
        ylims_current[1],
        new_ymax
    )

    # Position of "4 GFL", "5 GFL", etc.
    label_y =
        ylims_current[2] +
        0.05 * y_range


    for count in counts

        # Find configurations with this number of GFL buses
        positions = findall(
            groups.GFL_count .== count
        )

        first_x = first(positions)
        last_x  = last(positions)

        midpoint =
            (first_x + last_x) / 2


        # Group title
        annotate!(
            p,
            midpoint,
            label_y,

            text(
                "$count GFL",
                10,
                :center
            )
        )


        # Vertical separator between groups
        if last_x < nrow(groups)

            vline!(
                p,
                [last_x + 0.5],

                color = :gray,
                linestyle = :dash,
                linewidth = 1,

                label = false
            )

        end

    end


    return p

end

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
        
        include("ieee39_main_run.jl")
    end
end

#path_to_freq_apex = joinpath(joinpath(joinpath(pwd(), "plots"), test_name), "max_min_freq.csv")

#plot_apex(path_to_freq_apex)


