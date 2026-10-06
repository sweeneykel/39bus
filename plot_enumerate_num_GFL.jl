using Combinatorics

using CSV, DataFrames
using Plots

function plot_apex(csv_file::String)

    df = CSV.read(csv_file, DataFrame)

    # Expected CSV format:
    #
    # list_GFL          bus_number    apex    nadir
    # "30;34;32;33"     30            ...     ...
    #
    # list_GFL is stored as a String in one CSV cell.


    # ---------------------------------------------------------
    # Prepare GFL configuration labels
    # ---------------------------------------------------------

    # Convert:
    # "30;34;32;33"
    #
    # into:
    # "30, 34, 32, 33"

    df.GFL_label = [
        join(strip.(split(config, ";")), ", ")
        for config in df.list_GFL
    ]

    # Count how many buses are in each GFL configuration
    df.GFL_count = [
        length(split(config, ";"))
        for config in df.list_GFL
    ]


    # ---------------------------------------------------------
    # Find unique GFL configurations
    # ---------------------------------------------------------

    groups = unique(
        select(df, [:list_GFL, :GFL_label, :GFL_count])
    )

    # Sort first by number of GFL buses,
    # then by configuration
    sort!(groups, [:GFL_count, :GFL_label])


    # ---------------------------------------------------------
    # Give each GFL configuration its own x-axis position
    # ---------------------------------------------------------

    group_to_x = Dict(
        row.list_GFL => i
        for (i, row) in enumerate(eachrow(groups))
    )

    df.x_position = [
        group_to_x[config]
        for config in df.list_GFL
    ]

    tick_positions = 1:nrow(groups)
    tick_labels = groups.GFL_label


    # ---------------------------------------------------------
    # Create plot
    # ---------------------------------------------------------

    p = plot(
        xlabel = "GFL Bus Configuration",
        ylabel = "Δf [Hz]",
        title = "Frequency apex per bus for different combinations of GFL and SG",

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
    # Fixed color for each bus
    # ---------------------------------------------------------

    BUS_COLORS = Dict(
        30 => :red,
        31 => :blue,
        32 => :green,
        33 => :orange,
        34 => :purple,
        35 => :brown,
        36 => :pink,
        37 => :cyan,
        38 => :magenta,
        39 => :black
    )


    # ---------------------------------------------------------
    # Plot data
    # ---------------------------------------------------------

    for bus in sort(unique(df.bus_number))

        bus_df = filter(
            row -> row.bus_number == bus,
            df
        )

        for row in eachrow(bus_df)

            scatter!(
                p,
                [row.x_position],
                [row.apex],

                marker = :circle,
                markersize = 6,

                color = BUS_COLORS[bus],

                label = false
            )

        end
    end


    # ---------------------------------------------------------
    # Add legend entries once
    # ---------------------------------------------------------

    for bus in 30:39

        scatter!(
            p,
            [NaN],
            [NaN],

            marker = :circle,
            markersize = 6,

            color = BUS_COLORS[bus],

            label = "Bus $bus"
        )

    end


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


function plot_nadir(csv_file::String)

    df = CSV.read(csv_file, DataFrame)

    # Expected CSV format:
    #
    # list_GFL          bus_number    apex    nadir
    # "30;34;32;33"     30            ...     ...
    #
    # list_GFL is stored as a String in one CSV cell.


    # ---------------------------------------------------------
    # Prepare GFL configuration labels
    # ---------------------------------------------------------

    # Convert:
    # "30;34;32;33"
    #
    # into:
    # "30, 34, 32, 33"

    df.GFL_label = [
        join(strip.(split(config, ";")), ", ")
        for config in df.list_GFL
    ]

    # Count how many buses are in each GFL configuration
    df.GFL_count = [
        length(split(config, ";"))
        for config in df.list_GFL
    ]


    # ---------------------------------------------------------
    # Find unique GFL configurations
    # ---------------------------------------------------------

    groups = unique(
        select(df, [:list_GFL, :GFL_label, :GFL_count])
    )

    # Sort first by number of GFL buses,
    # then by configuration
    sort!(groups, [:GFL_count, :GFL_label])


    # ---------------------------------------------------------
    # Give each GFL configuration its own x-axis position
    # ---------------------------------------------------------

    group_to_x = Dict(
        row.list_GFL => i
        for (i, row) in enumerate(eachrow(groups))
    )

    df.x_position = [
        group_to_x[config]
        for config in df.list_GFL
    ]

    tick_positions = 1:nrow(groups)
    tick_labels = groups.GFL_label


    # ---------------------------------------------------------
    # Create plot
    # ---------------------------------------------------------

    p = plot(
        title = "Frequency nadir per bus for different combinations of GFL and SG",
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
    # Fixed color for each bus
    # ---------------------------------------------------------

    BUS_COLORS = Dict(
        30 => :red,
        31 => :blue,
        32 => :green,
        33 => :orange,
        34 => :purple,
        35 => :brown,
        36 => :pink,
        37 => :cyan,
        38 => :magenta,
        39 => :black
    )


    # ---------------------------------------------------------
    # Plot data
    # ---------------------------------------------------------

    for bus in sort(unique(df.bus_number))

        bus_df = filter(
            row -> row.bus_number == bus,
            df
        )

        for row in eachrow(bus_df)

            scatter!(
                p,
                [row.x_position],
                [row.nadir],

                marker = :circle,
                markersize = 6,

                color = BUS_COLORS[bus],

                label = false
            )

        end
    end


    # ---------------------------------------------------------
    # Add legend entries once
    # ---------------------------------------------------------

    for bus in 30:39

        scatter!(
            p,
            [NaN],
            [NaN],

            marker = :circle,
            markersize = 6,

            color = BUS_COLORS[bus],

            label = "Bus $bus"
        )

    end


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